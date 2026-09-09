

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
                    times, k_wall, n_segments=60):
    """N_wells = max(N_heat_transfer, N_inventory).

    v0.1 solved only the first: can the wells ABSORB the energy in time? That is
    a rate question. Nothing in the sizing required the wells to CONTAIN the
    energy, so with the k_w error -- which made the rate criterion the larger of
    the two -- the omission never showed.

    N_inventory here is NOT the same quantity as the 'Ideal number of wells'
    printed by the v0.1 notebook. That figure (12.215) comes from a separate,
    simpler calculation: the thermal load divided by the sensible plus latent
    capacity of the PCM, giving the PCM volume required and hence a well count.
    It carries no thermal resistance and no losses, so it is a genuine LOWER
    BOUND on N -- the best any design could do. N_inventory is computed from the
    marched model and must sit at or above it.

    The two should differ by roughly the sensible-heat contribution the marched
    model omits (H5). See the check in section 12.
    """
    def evaluate(N):
        m1 = m_dot_total / (N * case.num_tubes)
        r = march(case, T_inlet, T_m_seg, m1, k_wall, times,
                  n_segments=n_segments, mode="charge")
        dz = case.L_tube / n_segments
        eps = float(np.sum(r["A_melt"] * dz)) * case.num_tubes / case.V_well
        Q_total = r["Q_cum_J"] * case.num_tubes * N / 1000.0        # kJ
        return eps, Q_total

    lo, hi = case.N_wells_bracket

    try:                                    # inventory: smallest N with eps <= 1
        N_inv, _ = _bisect(lambda N: evaluate(N)[0] - 1.0, lo, hi, 1e-3,
                           case.max_iterations, "inventory constraint")
    except RuntimeError:
        if evaluate(lo)[0] <= 1.0:
            N_inv = lo                      # inventory never binds
        else:
            raise

    try:                                    # heat transfer: Q delivered >= required
        N_heat, _ = _bisect(lambda N: D_E_required_kJ / evaluate(N)[1] - 1.0,
                            lo, hi, case.tol_Q_ratio, case.max_iterations,
                            "heat-transfer sizing")
    except RuntimeError:
        if D_E_required_kJ / evaluate(lo)[1] <= 1.0:
            N_heat = lo
        else:
            raise

    N = max(N_heat, N_inv)
    eps, Q_total = evaluate(N)
    return dict(N_wells=N, N_heat=N_heat, N_inventory=N_inv,
                binding="inventory" if N_inv >= N_heat else "heat_transfer",
                eps_pcm=eps)


# ==========================================================================
# run   (charge, size, discharge, KPIs)
# ==========================================================================

def run_cycle(case):
    """The whole calculation. `case.front` selects the melt-front formulation."""
    rank, hp, T = cycle_state_points(case)
    rank_eff, hp_cop = rank["rank_eff"], hp["hp_cop"]
    E = energy_budget(case, rank_eff, hp_cop, T)
    T_m_lay, T_m_lay_dc = melting_temperatures(case)

    times_ch = np.logspace(0.0, np.log10(case.t_ch * 3600.0), case.n_times)
    times_dc = np.logspace(0.0, np.log10(case.t_dc * 3600.0), case.n_times)
    gv = case.geom_vector()

    # The k_w slot of the charging chain. False reproduces the v0.1 defect.
    k_charge = case.k_wall if case.charge_uses_wall_conductivity else case.k_l

    detail = {}

    if case.front == "energy_balance":
        n_seg = case.n_segments
        T_m_seg = layer_map(T_m_lay, n_seg, case.N_lay)
        T_m_seg_dc = layer_map(T_m_lay_dc, n_seg, case.N_lay)

        sz = size_well_field(case, T["T_4c"], T_m_seg, E["m_dot_w_ch"],
                             E["D_E_out_HP"], times_ch, k_charge, n_seg)
        N_wells = sz["N_wells"]
        eps_pcm = sz["eps_pcm"]
        detail.update(sz)

        m1_ch = E["m_dot_w_ch"] / (N_wells * case.num_tubes)
        r_ch = march(case, T["T_4c"], T_m_seg, m1_ch, k_charge, times_ch,
                     n_segments=n_seg, mode="charge")
        E_stored = r_ch["Q_cum_J"] * case.num_tubes * N_wells / 1000.0

        def discharged(ratio):
            r = march(case, T["T_3d"], T_m_seg_dc, ratio * m1_ch, case.k_wall,
                      times_dc, n_segments=n_seg, A_melt0=r_ch["A_melt"],
                      mode="discharge")
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
            * case.num_tubes / case.V_well)

    elif case.front == "closed_form":
        lo, hi = case.N_wells_bracket
        N_wells = _find_N_wells_closed_form(case, gv, T, T_m_lay, k_charge,
                                            times_ch, E)
        m1_ch = E["m_dot_w_ch"] / (N_wells * case.num_tubes)
        _, _, Vm, _ = time_profiles_melt(
            [case.t_ch * 3600.0], gv, T["T_4c"], case.N_lay, T_m_lay, k_charge,
            case.Rf_i, m1_ch, case.P, case.fluid2, case.k_l, case.cp_l,
            case.rho_l, case.h_m, case.n_segments, case.delta_max)
        V_tube = Vm[case.t_ch * 3600.0]
        eps_pcm = V_tube * case.num_tubes / case.V_well
        m_dot_d_well1, ratio = _find_discharge_closed_form(
            case, gv, T, T_m_lay_dc, times_dc, m1_ch, N_wells, E)
    else:
        raise ValueError(f"unknown front: {case.front!r}")

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
    detail.update(flow_ratio_dc_ch=ratio, m_dot_d_well1=m_dot_d_well1,
                  pumping_ch_kW=pumping_ch, pumping_dc_kW=pumping_dc,
                  front=case.front)
    return {"kpis": kpis, "detail": detail, "T": T}


def _find_N_wells_closed_form(case, gv, T, T_m_lay, k_charge, times_ch, E):
    """v0.1 sizing, using the ORIGINAL Illinois regula-falsi solver.

    Note the parameter this solver names `k_m_l` -- PCM liquid conductivity --
    is the slot the wall conductivity is passed into. The defect is visible in
    the signature itself.
    """
    lo, hi = case.N_wells_bracket
    N = find_N_wells_for_Q_ratio_ch_fast(
        1.0, case.tol_Q_ratio, lo, hi, E["m_dot_w_ch"], case.num_tubes, gv,
        times_ch, T["T_4c"], case.N_lay, T_m_lay, k_charge, case.Rf_i, case.P,
        case.fluid2, case.cp_l, case.rho_l, case.h_m, case.n_segments,
        case.delta_max, E["D_E_out_HP"], N_hint=None,
        max_iterations=case.max_iterations, verbose=False)
    if not np.isfinite(N) or N >= hi - 1.0:
        raise RuntimeError(f"well-field sizing did not converge (N={N})")
    return N


def _find_discharge_closed_form(case, gv, T, T_m_lay_dc, times_dc, m1_ch,
                                N_wells, E):
    """v0.1 discharge, original solver. The front restarts from solid PCM."""
    rlo, rhi = case.ratio_bracket
    m_dot_d = find_m_dot_d_well1_for_Q_ratio_dc_fast(
        1.0, case.tol_Q_ratio, rlo, rhi, m1_ch, N_wells, case.num_tubes, gv,
        times_dc, T["T_3d"], case.N_lay, T_m_lay_dc, case.k_wall, case.Rf_i,
        case.P, case.fluid2, case.k_s, case.cp_s, case.rho_s, case.h_m,
        case.n_segments, case.delta_max, E["D_E_in_ORC"], ratio_hint=None,
        max_iterations=case.max_iterations, verbose=False)
    if not np.isfinite(m_dot_d) or m_dot_d <= 0:
        raise RuntimeError(f"discharge flow did not converge ({m_dot_d})")
    return m_dot_d, m_dot_d / m1_ch
