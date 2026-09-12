

# ==========================================================================
# cycles + energy budget
# ==========================================================================

def cycle_state_points(case):
    """ORC and heat-pump state points, and the two cycle efficiencies.

    Every state point is fixed by a prescribed approach temperature. Note the
    expansions and compressions carry no isentropic efficiency, so eta_ORC and
    COP are idealised -- only the electrical and mechanical efficiencies appear.
    """
    T_2d_C = case.T_m_C - case.DT_m_2D
    T = dict(
        T_1e=case.T_sink_C + case.DT_E_sink + 273.15,
        T_3e=T_2d_C - case.DT_2D_3E + 273.15,
        T_2d=T_2d_C + 273.15,
        T_3d=T_2d_C - case.DT_3C_2C + 273.15,     # DT_2D_3D is tied to the glide
        T_4c=case.T_m_C + case.DT_4C_M + 273.15,
        T_4a=case.T_source_C - case.DT_3A_4A + 273.15,
        T_13h=case.T_source_C - case.DT_3A_13H + 273.15,
    )
    T["T_3c"] = T["T_4c"]
    T["T_2c"] = T["T_3c"] - case.DT_3C_2C
    T["T_2h"] = T["T_3c"] + case.DT_2H_3C

    rank = double_stage_rankine(case.fluid, T["T_1e"], T["T_3e"])
    hp = two_stage_htheatpump_2regs(case.refrig, T["T_13h"], T["T_2h"],
                                    case.DT_sub)
    for name, v in (("rank_eff", rank.get("rank_eff")),
                    ("hp_cop", hp.get("hp_cop"))):
        if v is None or not np.isfinite(v):
            raise RuntimeError(f"{name} is not finite ({v})")
    return rank, hp, T


def energy_budget(case, rank_eff, hp_cop, T):
    """Work backwards from the specified electrical output to the charging duty."""
    W_dot_T = case.W_dot_el_out / case.Turb_eff / case.ElG_eff
    Q_dot_in_ORC = W_dot_T / rank_eff
    D_E_in_ORC = Q_dot_in_ORC * case.t_dc * 3600.0                 # kJ
    D_E_out_HP = D_E_in_ORC * (1.0 + case.loss_surplus)
    Q_dot_out_HP = D_E_out_HP / case.t_ch / 3600.0
    W_dot_C_HP = Q_dot_out_HP / hp_cop
    W_dot_el_in = W_dot_C_HP / case.Comp_eff / case.ElH_eff

    cp_w = CP.PropsSI("C", "T", 0.5 * (T["T_3c"] + T["T_2c"]), "P",
                      case.P, case.fluid2) / 1000.0
    m_dot_w_ch = Q_dot_out_HP / cp_w / case.DT_3C_2C
    return dict(Q_dot_in_ORC=Q_dot_in_ORC, D_E_in_ORC=D_E_in_ORC,
                D_E_out_HP=D_E_out_HP, Q_dot_out_HP=Q_dot_out_HP,
                W_dot_el_in=W_dot_el_in, m_dot_w_ch=m_dot_w_ch)


def melting_temperatures(case):
    """Layer melting temperatures for charging and for discharging."""
    dTm = case.DT_3C_2C / case.N_lay
    i = np.arange(case.N_lay, dtype=float)
    return ((case.T_m - i * dTm).tolist(),
            (case.T_m - case.DT_3C_2C + (i + 1.0) * dTm).tolist())


# ==========================================================================
# sizing
# ==========================================================================

class SizingResult(dict):
    """The sizing dict, with retired key names that fail loudly.

    THERE ARE ONLY TWO INDEPENDENT NUMBERS. The model once exposed four names
    for them, and two of those names were traps:

      N_inventory  was an ALIAS for N_capacity, so `N_inventory == N_capacity`
                   returned True and nobody could tell them apart -- but in
                   v0.3 the same name meant something ELSE entirely (the well
                   count at which eps_PCM(N) = 1, a marched quantity needing a
                   root solve). Same name, two meanings, no warning. Anyone
                   comparing a v0.3 printout with a v0.4 one was comparing
                   different quantities.

      N_heat       was the dict key, while the report said "charge-rate
                   criterion" and the sweep column said N_rate. Three labels,
                   one quantity, and `d['N_rate']` raised KeyError.

    Both are now retired. Reading either raises with an explanation rather
    than returning a number that may mean something other than you think.
    """

    _RETIRED = {
        'N_inventory':
            "'N_inventory' is retired. In v0.4 it was a silent alias for "
            "'N_capacity'; in v0.3 it meant the well count at which "
            "eps_PCM(N) = 1, which is NOT the same quantity. Use "
            "'N_capacity' if you want the closed-form capacity bound, and "
            "do not compare v0.3 printouts of N_inventory with it.",
        'N_heat':
            "'N_heat' is retired; the rate criterion is now 'N_rate' "
            "everywhere -- dict key, report text and sweep column.",
    }

    def __missing__(self, key):
        if key in self._RETIRED:
            raise KeyError(self._RETIRED[key])
        raise KeyError(key)

    def get(self, key, default=None):
        # .get() must not silently return None for a retired name: that is the
        # exact failure this class exists to prevent.
        if key in self._RETIRED and key not in self:
            raise KeyError(self._RETIRED[key])
        return super().get(key, default)


def _bisect(f, lo, hi, tol, max_iter, what):
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise RuntimeError(f"{what}: root outside bracket [{lo}, {hi}]")
    for it in range(1, max_iter + 1):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol or (hi - lo) < 1e-3:
            return mid, it
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    raise RuntimeError(f"{what}: no convergence in {max_iter} iterations")


def size_well_field(case, T_inlet, T_m_seg, m_dot_total, D_E_required_kJ,
                    times, k_wall, n_segments=60, T_discharge_in=None):
    """N_wells = max(N_charge_rate, N_capacity).

    TWO CRITERIA, and they answer different questions.

    N_charge_rate -- can the field ABSORB D_E_required within t_ch? A rate
    question. This is the loop v0.1 actually solved
    (`find_N_wells_for_Q_ratio_ch_fast`), and it is the only one it used.

    N_capacity -- can the field CONTAIN D_E_required at all? A capacity
    question, and the one the v0.1 notebook computed as `N_wells_ideal` and then
    spent only on costing:

        N_capacity = D_E_required / E_well,
        E_well = rho V_well [ h_m + cp_l (T_charge_in - T_m) ]

    It is a hard lower bound: no thermal resistance, no losses, so no design can
    do better. Because it is a closed form it costs nothing to evaluate, unlike
    the version this replaces.

    WHY NOT SIZE ON THE DISCHARGE. Sizing must act on whichever variable is
    free. During charging the flow is pinned by the specified glide
    (m_dot = Q_out_HP / cp / DT_3C_2C), so N is the only freedom and the charge
    can size it. During discharging N is already fixed, so the flow is the
    freedom and the discharge sizes that instead. Pinning the discharge flow to
    its glide as well leaves the delivered energy saturating at 0.996 of the
    requirement for ANY well count -- there is no root -- because the fluid
    leaves slightly short of the nominal outlet temperature. The v0.1
    arrangement is correct, and this is why.
    """
    def evaluate(N):
        m1 = m_dot_total / (N * case.num_tubes)
        r = march_h(case, T_inlet, T_m_seg, m1, k_wall, times,
                    n_segments=n_segments)
        Q_total = r["Q_cum_J"] * case.num_tubes * N / 1000.0        # kJ
        return r, Q_total

    lo, hi = case.N_wells_bracket

    # ---- capacity: closed form, no march needed ---------------------------
    E_well = well_capacity_kJ(case, T_inlet, T_discharge_in)
    N_cap = D_E_required_kJ / E_well

    # ---- charge rate ------------------------------------------------------
    try:
        N_rate, _ = _bisect(lambda N: D_E_required_kJ / evaluate(N)[1] - 1.0,
                            lo, hi, case.tol_Q_ratio, case.max_iterations,
                            "charge-rate sizing")
    except RuntimeError:
        if D_E_required_kJ / evaluate(lo)[1] <= 1.0:
            N_rate = lo
        else:
            raise

    N = max(N_rate, N_cap)
    r, Q_total = evaluate(N)
    A_avail, E_lat, _, _ = pcm_capacities(case)
    dz = case.L_tube / n_segments
    eps = float(np.sum(r["A_melt"] * dz)) * case.num_tubes / case.V_well

    # Which criterion actually set N, and by how much did it win? A margin of a
    # few percent means the design sits on the crossover: a parameter sweep will
    # switch criterion partway through, and reading either criterion ALONE
    # across such a sweep gives a trend that is not the trend in N_wells.
    binding = "capacity" if N_cap >= N_rate else "rate"
    margin = abs(N_cap - N_rate) / N

    return SizingResult(N_wells=N, N_rate=N_rate, N_capacity=N_cap,
                binding=binding, binding_margin=margin,
                near_crossover=bool(margin < 0.05),
                eps_pcm=eps, E_well_kJ=E_well,
                eps_local_max=float(r["eps_local"].max()),
                eps_local_min=float(r["eps_local"].min()),
                merge_proximity_max=r["merge_proximity_max"],
                f_near_merge=r["f_near_merge"],
                delta_merge=r["delta_merge"])


def sizing_report(case, detail, label=""):
    """Print the well count with the criterion that SET it named first.

    THERE ARE ONLY TWO INDEPENDENT NUMBERS, and one derived from them:

        N_rate      can the field ABSORB the energy within t_ch?   RATE
        N_capacity  can the field CONTAIN the energy at all?       CAPACITY
        N_wells     = max(N_rate, N_capacity)                      THE ANSWER

    N_capacity is the seductive one: it is a closed form, it is smooth, it
    responds to every geometric parameter, and it is wrong to read on its own.
    It is a LOWER BOUND -- the best any design could do with no thermal
    resistance and infinite time. Because N_wells is a MAXIMUM, whenever the
    rate criterion binds N_capacity moves while the answer does not, and it
    can move the opposite way.

    The failure this guards against is real and was made by a careful reader.
    Sweeping fin radial length, N_capacity fell monotonically (11.573, 11.239,
    10.985 for 7.5, 3.75, 0.75 mm) which reads as "fins make things worse".
    N_capacity = D_E / (rho V_well [h_m + cp dT]) contains no heat transfer at
    all: shrinking a fin returns its metal volume to the PCM, so the bound
    falls by pure volume bookkeeping. Over the same sweep the real answer rose
    11.573 -> 11.569 -> 14.654, because the rate criterion took over. Removing
    the fins entirely costs 42 % more wells.

    So: lead with N_wells, name the binding criterion, and label N_capacity as
    a bound rather than a result.
    """
    d = detail
    N = d["N_wells"]
    name = {"capacity": "CAPACITY (can the field CONTAIN the energy at all?)",
            "rate":     "RATE     (can the field ABSORB it within t_ch?)"}
    slack = "capacity" if d["binding"] == "rate" else "rate"

    head = f"WELL FIELD SIZING{('  --  ' + label) if label else ''}"
    print(head)
    print("=" * max(len(head), 66))
    print(f"  N_wells = {N:8.3f}        set by {name[d['binding']]}")
    print()
    print(f"    rate criterion          N_rate     = {d['N_rate']:8.3f}"
          f"   {'<-- BINDING' if d['binding'] == 'rate' else ''}")
    print(f"    capacity criterion      N_capacity = {d['N_capacity']:8.3f}"
          f"   {'<-- BINDING' if d['binding'] == 'capacity' else ''}")
    print(f"    N_wells = max(the two);  {slack} criterion is slack by "
          f"{100*d['binding_margin']:.1f} %")
    print()
    print("    N_capacity is a LOWER BOUND, not a result. It is the closed form")
    print("    D_E / (rho V_well [h_m + cp_l dT]) and contains no heat transfer.")
    print("    Do not read it across a parameter sweep on its own.")

    if d.get("near_crossover"):
        print()
        print(f"  ** the two criteria are within {100*d['binding_margin']:.1f} % "
              f"-- this design sits ON the crossover.")
        print("     A sweep of any geometric parameter will switch the binding")
        print("     criterion partway through. Plot N_wells, not either criterion.")

    # ---- H3: how close is the melt front to its neighbour? ---------------
    if "merge_proximity_max" in d:
        p, f = d["merge_proximity_max"], d["f_near_merge"]
        print()
        print(f"  H3 (concentric annulus, fronts do not interact):")
        print(f"    fronts merge at delta = {1000*d['delta_merge']:.2f} mm"
              f"   (the CELL radius, not the borehole wall)")
        print(f"    max delta / delta_merge          = {p:.3f}")
        print(f"    fraction of well above 0.90      = {100*f:.0f} %")
        if p >= 0.999:
            print("    -> the front SITS ON the merge line. Under Formulation C")
            print("       delta cannot EXCEED it: the enthalpy cap enforces")
            print("       A_melt <= A_avail, and A(delta_merge) = A_avail")
            print("       identically. So this is not a violation. It is the")
            print("       solution running at the exact limit of H3's validity.")
        elif f > 0.5:
            print("    -> most of the well runs near the merge line.")
        else:
            print("    -> comfortably inside H3.")
        if p > 0.90:
            print()
            print("       Two consequences, both biasing the model AGAINST fins:")
            print("       (a) the annular resistance ln(1+delta/r_e)/(2 pi k)")
            print("           assumes the full azimuth 2 pi (r_e+delta) conducts.")
            print("           Past merge it does not -- the neighbouring front")
            print("           has taken part of it, and the last solid sits in")
            print("           the corners between legs, on a longer and narrower")
            print("           path. Late-stage resistance is UNDERSTATED, which")
            print("           flatters whichever design is struggling to finish.")
            print("       (b) fins enter U_i only as added SURFACE at the tube")
            print("           (eta_f, P_T). They are not represented as radial")
            print("           conduction paths reaching into that corner PCM --")
            print("           which is the single thing fins do best in exactly")
            print("           this regime. The model cannot credit it.")


# ==========================================================================
# run   (charge, size, discharge, KPIs)
# ==========================================================================

def run_cycle(case):
    """The whole calculation: cycles, budget, sizing, charge, discharge, KPIs."""
    rank, hp, T = cycle_state_points(case)
    rank_eff, hp_cop = rank["rank_eff"], hp["hp_cop"]
    E = energy_budget(case, rank_eff, hp_cop, T)
    T_m_lay, T_m_lay_dc = melting_temperatures(case)

    times_ch = np.logspace(0.0, np.log10(case.t_ch * 3600.0), case.n_times)
    times_dc = np.logspace(0.0, np.log10(case.t_dc * 3600.0), case.n_times)
    gv = case.geom_vector()
    k_charge = case.k_wall

    detail = SizingResult()

    n_seg = case.n_segments
    T_m_seg = layer_map(T_m_lay, n_seg, case.N_lay)
    # THE DISCHARGE MARCHES FROM THE OTHER END. The PCM at a given physical
    # position has ONE melting temperature; it is the march INDEX that mirrors,
    # because the flow reverses between half-cycles. Building the discharge
    # cascade as the exact reverse of the charge cascade also removes a
    # one-layer (6.111 K) mismatch that layer_map introduces whenever
    # n_segments is not a multiple of N_lay.
    T_m_seg_dc = T_m_seg[::-1].copy()

    sz = size_well_field(case, T["T_4c"], T_m_seg, E["m_dot_w_ch"],
                         E["D_E_out_HP"], times_ch, k_charge, n_seg,
                         T_discharge_in=T["T_3d"])
    N_wells = sz["N_wells"]
    eps_pcm = sz["eps_pcm"]
    detail.update(sz)

    m1_ch = E["m_dot_w_ch"] / (N_wells * case.num_tubes)
    r_ch = march_h(case, T["T_4c"], T_m_seg, m1_ch, k_charge, times_ch,
                   n_segments=n_seg)
    E_stored = r_ch["Q_cum_J"] * case.num_tubes * N_wells / 1000.0

    def discharged(ratio):
        # E' is measured from a datum of SOLID AT THE LOCAL T_m, so the state
        # array mirrors along with the march index. Passing it unreversed shifts
        # the datum by up to 48.9 K and reports PCM hotter than any fluid that
        # ever touched it: 176.6 C against a 160.0 C charge inlet.
        r = march_h(case, T["T_3d"], T_m_seg_dc, ratio * m1_ch,
                    case.k_wall, times_dc, n_segments=n_seg,
                    E0=r_ch["E"][::-1])
        return abs(r["Q_cum_J"]) * case.num_tubes * N_wells / 1000.0, r

    ratio, _ = _bisect(lambda x: E["D_E_in_ORC"] / discharged(x)[0] - 1.0,
                       *case.ratio_bracket, case.tol_Q_ratio,
                       case.max_iterations, "discharge flow")
    Q_dc, r_dc = discharged(ratio)
    m_dot_d_well1 = ratio * m1_ch
    dz = case.L_tube / n_seg
    detail.update(
        closure_charge=r_ch["closure"], closure_discharge=r_dc["closure"],
        E_stored_kJ=E_stored, E_discharged_kJ=Q_dc,
        eta_storage=Q_dc / E_stored,
        eps_end_of_discharge=float(np.sum(r_dc["A_melt"] * dz))
        * case.num_tubes / case.V_well,
        E_sensible_charge_J=r_ch["E_sensible_J"],
        E_sensible_discharge_J=r_dc["E_sensible_J"])

    # ---- parasitics ----
    _, pp_dc = calculate_pressure_drop(gv, m_dot_d_well1, case.fluid2,
                                       T["T_2d"], T["T_3d"], case.P)
    pumping_dc = pp_dc * case.num_tubes * N_wells / 1000.0
    m1_ch = E["m_dot_w_ch"] / (N_wells * case.num_tubes)
    _, pp_ch = calculate_pressure_drop(gv, m1_ch, case.fluid2,
                                       T["T_4c"], T["T_4a"], case.P)
    pumping_ch = pp_ch * case.num_tubes * N_wells / 1000.0

    # ---- KPIs ----
    D_E_in_ORC_kWh = E["D_E_in_ORC"] / 3600.0
    E_well = D_E_in_ORC_kWh / N_wells / 1000.0                     # MWh
    kpis = {
        "eta_rte": ((case.W_dot_el_out - pumping_dc) * case.t_dc
                    / ((E["W_dot_el_in"] + pumping_ch) * case.t_ch)),
        "eta_rte_nopump": (case.W_dot_el_out * case.t_dc
                           / (E["W_dot_el_in"] * case.t_ch)),
        "cop_hp": hp_cop,
        "eta_orc": rank_eff,
        "eps_pcm": eps_pcm,
        "E_well": E_well,
        "rho_E": E_well * 1000.0 / case.V_well,
        "N_wells": N_wells,
        "wells_per_MW": N_wells / (case.W_dot_el_out / 1000.0),
        "f_pump": (pumping_ch + pumping_dc) / case.W_dot_el_out,
        "c_pcm": case.cost_per_kWh,
    }

    # ---- DoE-study indicators -------------------------------------------
    # Restored from the earlier factorial study. Three are new; eps_m, RTE and
    # dE_therm were already present as eps_pcm, eta_rte and E_well.
    #
    #   eta_T      thermal -> net electric conversion of the discharge side,
    #              W_el_out / Q_in_ORC. Identically rank_eff * Turb_eff * ElG_eff.
    #   eta_T_eff  the same after the discharge pumping parasitic is charged
    #              against it: eta_T - N W_f,dc / Q_in_ORC.
    #   eps_rte    RTE weighted by the melted fraction, RTE * eps_m.
    #   dE_therm   thermal energy stored per well  [kWh]
    #   dE_elec    net electric energy stored per well, dE_therm * eta_T_eff
    #
    # CAVEAT on eps_rte, which matters for reading a DoE table. eps_m saturates:
    # under the enthalpy formulation, once a segment has melted all its PCM the
    # further energy goes into SUPERHEAT, which eps_m cannot see. At this design
    # point eps_m = 1.0000 exactly, so eps_rte == RTE and carries no information.
    # It only discriminates across designs that do NOT saturate. `util_enthalpy`
    # below is the non-saturating companion: energy actually banked per well as a
    # fraction of the capacity ceiling of well_capacity_kJ.
    eta_T = case.W_dot_el_out / E["Q_dot_in_ORC"]
    eta_T_eff = eta_T - pumping_dc / E["Q_dot_in_ORC"]
    dE_therm = E["D_E_in_ORC"] / 3600.0 / N_wells                   # kWh/well
    kpis.update({
        "eta_T": eta_T,
        "eta_T_eff": eta_T_eff,
        "eps_rte": kpis["eta_rte"] * eps_pcm,
        "dE_therm_kWh": dE_therm,
        "dE_elec_kWh": dE_therm * eta_T_eff,
        "util_enthalpy": (detail["E_stored_kJ"] / N_wells) / detail["E_well_kJ"],
    })
    detail.update(flow_ratio_dc_ch=ratio, m_dot_d_well1=m_dot_d_well1,
                  pumping_ch_kW=pumping_ch, pumping_dc_kW=pumping_dc)
    return {"kpis": kpis, "detail": detail, "T": T}




# ==========================================================================
# cyclic steady state   (v0.5)
# ==========================================================================
# run_cycle() above answers "how many wells does one charge from a COLD START
# need?". That is a first-cycle question, and the cycle does not return to the
# state it started from, so the answer is not what a plant repeating a 24 h
# schedule would see.
#
# WHY SIZING CANNOT SIMPLY BE MOVED TO CSS. At cyclic steady state the state
# returns to itself, so the enthalpy change over a cycle is zero. This model has
# no loss path -- H2 no axial conduction, H6 adiabatic wall, H9 no azimuthal
# exchange -- so charge energy EQUALS discharge energy identically:
#
#       eta_storage(CSS) == 1     exactly, for every N and every flow.
#
# Two consequences follow, and both are structural rather than numerical.
#
#   1. The assumed loss surplus lambda is UNATTAINABLE at CSS. The budget
#      requires D_E_out_HP = (1+lambda) D_E_in_ORC, so the ratio of CSS charge
#      to requirement pins at 1/(1+lambda) for EVERY N. There is no root.
#      On the first cycle the surplus is not lost but PARKED in the store as
#      residual melt; at CSS there is nowhere left to park it.
#
#   2. Even with lambda = 0 the rate criterion goes DEGENERATE. Charge and
#      discharge are then the same number, so one equation is left for two
#      unknowns: it fixes the FLOW RATIO and says nothing about N.
#
# THE REFORMULATION. Stop asking the energy balance to determine N. Specify the
# hardware, march to CSS, and report what comes out:
#
#       N        from latent heat alone, D_E_out_HP / (rho V_well h_m)
#       m_ch     pinned by the charging glide
#       m_dc     pinned by the discharging glide
#
# Nothing is solved -- there is no root-find anywhere in this routine -- so the
# question of convergence does not arise. What was an unsatisfiable constraint
# becomes a reported output: the deviation of the delivered energy, and hence of
# the net electrical output, from its target.
#
# This also dissolves an older objection. size_well_field's docstring argues
# that the discharge flow cannot be pinned to its glide because the delivered
# energy then saturates below the requirement for any well count, leaving no
# root. Under this framing that saturation is not a failure; it is the answer.


def simulate_css(case, N=None, n_cycles=80, tol=1e-9, record=False):
    """March a specified well field to cyclic steady state.

    No root-finding. `N` defaults to the latent-heat-only well count and both
    flows are pinned by their glides, so the routine is a pure simulation.

    `lambda` is forced to zero regardless of `case.loss_surplus`, because a
    non-zero surplus cannot be absorbed at CSS in a model with no loss path
    (see the note above). The override is reported in the returned dict as
    `lambda_overridden` so it can never happen silently.

    Returns the converged half-cycle states, the convergence history, and the
    deviation of the delivered energy from the target.
    """
    lambda_overridden = case.loss_surplus != 0.0
    case = case.with_(loss_surplus=0.0) if lambda_overridden else case

    rank, hp, T = cycle_state_points(case)
    Eb = energy_budget(case, rank["rank_eff"], hp["hp_cop"], T)
    n = case.n_segments
    T_m_lay, T_m_lay_dc = melting_temperatures(case)
    T_m_seg = layer_map(T_m_lay, n, case.N_lay)
    T_m_seg_dc = T_m_seg[::-1].copy()     # flow reverses; see run_cycle
    times_ch = np.logspace(0.0, np.log10(case.t_ch * 3600.0), case.n_times)
    times_dc = np.logspace(0.0, np.log10(case.t_dc * 3600.0), case.n_times)

    # ---- hardware and flows: all SPECIFIED, none solved -------------------
    if N is None:
        N = Eb["D_E_out_HP"] * 1000.0 / (case.rho_latent * case.V_well * case.h_m)
    cp_ch = CP.PropsSI("C", "T", 0.5 * (T["T_3c"] + T["T_2c"]), "P",
                       case.P, case.fluid2) / 1000.0
    cp_dc = CP.PropsSI("C", "T", 0.5 * (T["T_2d"] + T["T_3d"]), "P",
                       case.P, case.fluid2) / 1000.0
    m1_ch = Eb["Q_dot_out_HP"] / cp_ch / case.DT_3C_2C / (N * case.num_tubes)
    m1_dc = Eb["Q_dot_in_ORC"] / cp_dc / case.DT_3C_2C / (N * case.num_tubes)

    # ---- march until the state repeats -----------------------------------
    E = np.zeros(n)
    history = []
    for k in range(n_cycles):
        ch = march_h(case, T["T_4c"], T_m_seg, m1_ch, case.k_wall, times_ch,
                     n_segments=n, E0=E, record=record)
        dc = march_h(case, T["T_3d"], T_m_seg_dc, m1_dc, case.k_wall, times_dc,
                     n_segments=n, E0=ch["E"][::-1], record=record)
        Q_ch = ch["Q_cum_J"] * case.num_tubes * N / 1000.0          # kJ, field
        Q_dc = abs(dc["Q_cum_J"]) * case.num_tubes * N / 1000.0
        drift = float(np.max(np.abs(dc["E"][::-1] - E)))
        history.append(dict(cycle=k + 1, Q_charge_kJ=Q_ch, Q_discharge_kJ=Q_dc,
                            drift=drift,
                            eps_end_charge=float(ch["eps_local"].mean()),
                            eps_end_discharge=float(dc["eps_local"].mean())))
        E = dc["E"][::-1]                 # back to charge-march indexing
        if drift < tol:
            break

    req = Eb["D_E_in_ORC"]
    deviation = Q_dc / req - 1.0
    return dict(
        N_wells=N, m1_ch=m1_ch, m1_dc=m1_dc, flow_ratio=m1_dc / m1_ch,
        cycles=len(history), converged=history[-1]["drift"] < tol,
        history=history, charge=ch, discharge=dc, E_css=E,
        Q_charge_kJ=Q_ch, Q_discharge_kJ=Q_dc, required_kJ=req,
        deviation=deviation,
        W_el_out_implied=case.W_dot_el_out * (1.0 + deviation),
        eta_storage=Q_dc / Q_ch, lambda_overridden=lambda_overridden,
        T=T, budget=Eb, T_m_seg=T_m_seg, T_m_seg_dc=T_m_seg_dc)


def css_report(case, r):
    """Print the cyclic-steady-state result."""
    print("CYCLIC STEADY STATE")
    print("=" * 66)
    if r["lambda_overridden"]:
        print("  NOTE  lambda forced to 0. A non-zero loss surplus cannot be")
        print("        absorbed at CSS: the state returns to itself and this")
        print("        model has no loss path, so charge == discharge exactly.")
        print()
    print(f"  N_wells          {r['N_wells']:10.4f}   from latent heat alone, not solved")
    print(f"  m_dot charge     {r['m1_ch']:10.4f} kg/s per leg-pair, pinned by the glide")
    print(f"  m_dot discharge  {r['m1_dc']:10.4f} kg/s per leg-pair, pinned by the glide")
    print(f"  flow ratio       {r['flow_ratio']:10.4f}   a consequence, not a solve")
    print()
    print(f"  converged in {r['cycles']} cycles"
          f"   (drift {r['history'][-1]['drift']:.1e})")
    print(f"  eta_storage      {r['eta_storage']:10.6f}   must be 1: adiabatic, closed cycle")
    print()
    print(f"  required   D_E_in_ORC  {r['required_kJ']:12.5e} kJ")
    print(f"  delivered  at CSS      {r['Q_discharge_kJ']:12.5e} kJ")
    print(f"  DEVIATION              {100*r['deviation']:+12.3f} %")
    print()
    print(f"  implied net output     {r['W_el_out_implied']/1000:10.4f} MWe"
          f"   (target {case.W_dot_el_out/1000:.4f})")
    ch, dc = r["charge"], r["discharge"]
    print()
    print("  melted fraction at CSS")
    print(f"    end of charge     mean {ch['eps_local'].mean():.4f}"
          f"   min {ch['eps_local'].min():.4f}   max {ch['eps_local'].max():.4f}")
    print(f"    end of discharge  mean {dc['eps_local'].mean():.4f}"
          f"   min {dc['eps_local'].min():.4f}   max {dc['eps_local'].max():.4f}")
    if dc["eps_local"].mean() > 0.05 and ch["eps_local"].mean() > 0.99:
        print()
        print("    The store melts completely but does not refreeze completely,")
        print("    so the shortfall is a DISCHARGE RATE limit, not an inventory")
        print("    limit. More PCM per well is not the lever.")
