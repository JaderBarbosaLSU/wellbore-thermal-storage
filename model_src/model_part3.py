

# ==========================================================================
# cycles + energy budget
# ==========================================================================

def T_m_bottom(case):
    """Melting temperature of the COLDEST cascade layer [K].

    The cascade is graded so that each layer sits DT_4C_M below the charging
    fluid at its own leading face, under a linear glide. The layers are spaced
    DT_3C_2C / N_lay apart, so the span of the melting temperatures is one
    layer width SHORT of the glide:

        T_m,bottom = T_m,top - DT_3C_2C * (N_lay - 1) / N_lay

    For the default case that is 150 - 55*8/9 = 101.111 C against a 55 K glide.
    The discharge enters the exchanger at this end, which is why this is the
    reference for DT_M_1D and not T_m,top.
    """
    return case.T_m - case.cascade_span


def _orc_cold_composite(rank, fluid, n=300):
    """The ORC-side composite curve of the evaporator: (Q, T), cold end first.

    Two cold streams share the water:

      * the MAIN stream, 1 kg per kg of evaporator flow, from state 10e to
        state 3e at p_evape -- liquid preheat, then boiling at the constant
        temperature T_3e;
      * the REHEAT stream, y_frac of the flow, from 5e to 6e at p_inte,
        superheated vapour throughout.

    Both are built on real enthalpy, so the latent plateau appears as a
    vertical segment rather than being smeared into a straight line in T.
    That matters: the plateau is where this evaporator pinches.
    """
    y = rank["y_frac"]
    p_e, p_i = rank["p_evape"], rank["p_5e"]
    h10, h3 = rank["h_10e"], rank["h_3e"]
    h5, h6 = rank["h_5e"], rank["h_6e"]
    T10, T3 = rank["T_10e"], rank["T_3e"]
    T5, T6 = rank["T_5e"], rank["T_6e"]
    h_bub = CP.PropsSI("H", "P", p_e, "Q", 0, fluid) / 1000.0

    # sample in ENTHALPY, not temperature: a (T, p) call exactly on the
    # saturation line raises, and the preheat ends exactly there
    h_pre = np.linspace(h10, h_bub, n)
    T_pre = np.array([CP.PropsSI("T", "H", x * 1000.0, "P", p_e, fluid)
                      for x in h_pre])
    h_rh = np.linspace(h5, h6, n)
    T_rh = np.array([CP.PropsSI("T", "H", x * 1000.0, "P", p_i, fluid)
                     for x in h_rh])

    def Q_main(t):
        if t <= T10:
            return 0.0
        if t >= T3:
            return h3 - h10                      # preheat AND latent
        return float(np.interp(t, T_pre, h_pre)) - h10

    def Q_reheat(t):
        if t <= T5:
            return 0.0
        if t >= T6:
            return y * (h6 - h5)
        return y * (float(np.interp(t, T_rh, h_rh)) - h5)

    lo, hi = min(T10, T5), max(T3, T6)
    grid = set(np.linspace(lo, hi, 4 * n))
    grid |= {T10, T5, T6, T3, np.nextafter(T3, lo)}   # break points
    Ts = np.array(sorted(t for t in grid if lo <= t <= hi))
    Qs = np.array([Q_main(t) + Q_reheat(t) for t in Ts])
    keep = np.concatenate(([True], np.diff(Qs) > 1e-12))
    return Qs[keep], Ts[keep]


def orc_pinch(rank, T_w_hot, T_w_cold, fluid, P_water, n=600):
    """Smallest water-minus-ORC temperature difference in the evaporator [K].

    Negative when the two composite curves cross, which is a second-law
    violation: the reported ORC efficiency is then unreachable, whatever the
    component efficiencies.

    The water is the hot stream and follows REAL enthalpy at `P_water`, not a
    straight line in temperature, and the search runs over an even grid in Q
    that includes every break point. The pinch is found wherever it lies; it
    is emphatically not assumed to sit at an end. For the case as shipped in
    v0.8 it sat at 23 % of the duty, where boiling begins.
    """
    Qc, Tc = _orc_cold_composite(rank, fluid)
    Q_total = Qc[-1]
    Tw = np.linspace(T_w_cold, T_w_hot, 400)
    hw = np.array([CP.PropsSI("H", "T", t, "P", P_water, "Water") / 1000.0
                   for t in Tw])
    Qh = (hw - hw[0]) * Q_total / (hw[-1] - hw[0])    # water flow scales out
    Qq = np.unique(np.concatenate([np.linspace(0.0, Q_total, n), Qc, Qh]))
    Qq = Qq[(Qq >= 0.0) & (Qq <= Q_total)]
    return float(np.min(np.interp(Qq, Qh, Tw) - np.interp(Qq, Qc, Tc)))


def feasible_rankine(case, T_1e, T_w_hot, T_w_cold):
    """The ORC evaporating temperature, as the MINIMUM of two constraints.

        T_3e = min( T_2d - DT_2D_3E ,  T_3e^pinch )

    The first is a minimum approach at the HOT END, which is what anybody
    writes first and what this model used alone until v0.9. The second is the
    largest boiling temperature for which the composite curves keep
    `case.DT_pinch_ORC` everywhere, found by bisection. Neither subsumes the
    other, so the binding one is whichever is lower, and the returned dict
    says which it was in `T_3e_binding`.

    At the default DT_2D_3E = 12 K the PINCH binds across the whole useful
    range of glides -- the hot-end rule never becomes active -- so that field
    is inert in practice while remaining a real constraint if raised. Writing
    it as a minimum rather than as "start from the old rule and correct it"
    is the same arithmetic and a more honest statement of intent.
    """
    T_3e_hot_end = T_w_hot - case.DT_2D_3E

    def pinch_at(t):
        return orc_pinch(double_stage_rankine(case.fluid, T_1e, t,
                                      case.eta_turb_s, case.eta_pump_s),
                         T_w_hot, T_w_cold, case.fluid, case.P)

    # is the hot-end value already pinch-feasible?
    if pinch_at(T_3e_hot_end) >= case.DT_pinch_ORC:
        rank = double_stage_rankine(case.fluid, T_1e, T_3e_hot_end,
                                    case.eta_turb_s, case.eta_pump_s)
        rank["T_3e_uncorrected"] = T_3e_hot_end
        rank["T_3e_hot_end"] = T_3e_hot_end
        rank["T_3e_binding"] = "hot-end approach DT_2D_3E"
        rank["pinch"] = orc_pinch(rank, T_w_hot, T_w_cold, case.fluid, case.P)
        return rank

    lo = T_1e + 10.0                     # floor: no cycle worth the name below
    hi = T_3e_hot_end
    if pinch_at(lo) < case.DT_pinch_ORC:
        raise RuntimeError(
            f"no ORC boiling temperature above {lo - 273.15:.1f} C clears the "
            f"required pinch of {case.DT_pinch_ORC:.1f} K against water at "
            f"{T_w_hot - 273.15:.1f} -> {T_w_cold - 273.15:.1f} C "
            f"(best is {pinch_at(lo):+.2f} K). The water glide is too large "
            f"for a cycle that boils at one temperature: reduce the "
            f"discharging glide, or raise the water temperatures, or relax "
            f"DT_pinch_ORC.")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if pinch_at(mid) >= case.DT_pinch_ORC:
            lo = mid
        else:
            hi = mid
    rank = double_stage_rankine(case.fluid, T_1e, lo,
                                case.eta_turb_s, case.eta_pump_s)
    rank["T_3e_uncorrected"] = T_3e_hot_end      # kept: reported in css_report
    rank["T_3e_hot_end"] = T_3e_hot_end
    rank["T_3e_binding"] = "evaporator pinch DT_pinch_ORC"
    rank["pinch"] = orc_pinch(rank, T_w_hot, T_w_cold, case.fluid, case.P)
    return rank


def _stream_curve(h_cold_end, h_hot_end, p, fluid, n=400):
    """(Q, T) for a single-pressure stream, from its COLD end.

    Q is duty per unit mass of that stream. Sampling is in enthalpy, so a
    latent plateau comes out flat and a (T, p) call never lands exactly on
    the saturation line.
    """
    hs = np.linspace(h_cold_end, h_hot_end, n)
    Ts = np.array([CP.PropsSI("T", "H", x * 1000.0, "P", p, fluid)
                   for x in hs])
    return hs - h_cold_end, Ts


def _water_curve(T_cold, T_hot, P, Q_total, n=400):
    """(Q, T) for the water, scaled so its total duty matches Q_total."""
    Ts = np.linspace(T_cold, T_hot, n)
    hs = np.array([CP.PropsSI("H", "T", t, "P", P, "Water") / 1000.0
                   for t in Ts])
    if abs(hs[-1] - hs[0]) < 1e-12:                 # isothermal reservoir
        return np.array([0.0, Q_total]), np.array([T_cold, T_hot])
    return (hs - hs[0]) * Q_total / (hs[-1] - hs[0]), Ts


def _UA_from_curves(Qh, Th, Qc, Tc, Q_duty_kW, n=2000):
    """UA and the effective mean approach, from two composite curves.

    Both curves are (Q, T) measured from the COLD END of the exchanger, on
    whatever duty basis; only their SHAPE is used. The conductance follows
    from the local driving difference:

        UA = int dQ / dT      ->     UA = Q_duty * <1/dT>

    and the effective mean approach is the number that reproduces it,
    dT_eff = Q_duty / UA. For a counterflow exchanger with both streams
    sensible and constant cp, dT_eff is the log-mean difference exactly;
    with a phase change on one side it is the correct generalisation, which
    the LMTD is not.

    Returns (UA [kW/K], dT_eff [K], pinch [K]). UA is infinite and dT_eff
    zero if the curves touch -- which is the honest answer, and is what the
    heat-pump evaporator returned before v0.9a gave it an approach.
    """
    Q_top = min(Qh[-1], Qc[-1])
    q = np.linspace(0.0, Q_top, n)
    dT = np.interp(q, Qh, Th) - np.interp(q, Qc, Tc)
    pinch = float(np.min(dT))
    if pinch <= 0.0:
        return float("inf"), 0.0, pinch
    # trapezoid on 1/dT, in units of the curves' own duty, then rescaled
    inv = np.trapezoid(1.0 / dT, q) if hasattr(np, "trapezoid") \
        else np.trapz(1.0 / dT, q)
    UA = Q_duty_kW * inv / Q_top
    return float(UA), float(Q_duty_kW / UA), pinch


def exchanger_UA(case, rank, hp, T, Eb):
    """Conductance UA [kW/K] required by each of the four exchangers.

    The model reports efficiencies but has never reported what they cost in
    hardware. That gap makes every approach temperature a free lunch: tighten
    `DT_pinch_ORC` to 1 K and eta_ORC rises with nothing to pay for it. UA is
    the missing axis, and it is the one that makes the approach temperatures
    genuine design variables rather than free parameters.

    Each exchanger is integrated over its real composite curves -- the same
    construction as `orc_pinch`, so a latent plateau is a flat segment and
    the water follows real enthalpy. Duties are plant-level, from the energy
    budget, so the UA values are for the whole field and add up.

    Returns a dict of dicts, one per exchanger, each with:
        duty_kW, UA_kW_per_K, dT_eff_K, pinch_K
    """
    F, R = case.fluid, case.refrig
    out = {}

    # ---- ORC evaporator: water (hot) against the two ORC streams ----------
    Q_orc = Eb["Q_dot_in_ORC"]
    Qc, Tc = _orc_cold_composite(rank, F)
    Qh, Th = _water_curve(T["T_3d"], T["T_2d"], case.P, Qc[-1])
    out["ORC evaporator"] = dict(zip(
        ("UA_kW_per_K", "dT_eff_K", "pinch_K"),
        _UA_from_curves(Qh, Th, Qc, Tc, Q_orc))) | {"duty_kW": Q_orc}

    # ---- ORC condenser: ORC fluid (hot) against the sink ------------------
    Q_cond_orc = Q_orc * (1.0 - rank["rank_eff"])
    Qh, Th = _stream_curve(rank["h_1e"], rank["h_4e"], rank["p_conde"], F)
    Qc, Tc = _water_curve(case.T_sink_C + 273.15,
                          case.T_sink_C + case.DT_sink_glide + 273.15,
                          case.P, Qh[-1])
    out["ORC condenser"] = dict(zip(
        ("UA_kW_per_K", "dT_eff_K", "pinch_K"),
        _UA_from_curves(Qh, Th, Qc, Tc, Q_cond_orc))) | {"duty_kW": Q_cond_orc}

    # ---- HTHP condenser: refrigerant (hot) against the charging water -----
    Q_cond_hp = Eb["Q_dot_out_HP"]
    Qh, Th = _stream_curve(hp["h_3h"], hp["h_2h"], hp["p_condh"], R)
    Qc, Tc = _water_curve(T["T_2c"], T["T_3c"], case.P, Qh[-1])
    out["HTHP condenser"] = dict(zip(
        ("UA_kW_per_K", "dT_eff_K", "pinch_K"),
        _UA_from_curves(Qh, Th, Qc, Tc, Q_cond_hp))) | {"duty_kW": Q_cond_hp}

    # ---- HTHP evaporator: source water (hot) against the refrigerant ------
    # Duty is the condenser heat less the compressor work: what the source
    # must actually give up.
    Q_evap_hp = Q_cond_hp * (1.0 - 1.0 / hp["hp_cop"])
    Qc, Tc = _stream_curve(hp["h_4h"], hp["h_13h"], hp["p_evaph"], R)
    Qh, Th = _water_curve(T["T_4a"], case.T_source_C + 273.15,
                          case.P, Qc[-1])
    out["HTHP evaporator"] = dict(zip(
        ("UA_kW_per_K", "dT_eff_K", "pinch_K"),
        _UA_from_curves(Qh, Th, Qc, Tc, Q_evap_hp))) | {"duty_kW": Q_evap_hp}

    finite = [v["UA_kW_per_K"] for v in out.values()
              if np.isfinite(v["UA_kW_per_K"])]
    out["TOTAL"] = dict(duty_kW=sum(v["duty_kW"] for v in out.values()),
                        UA_kW_per_K=sum(finite) if len(finite) == 4
                        else float("inf"),
                        dT_eff_K=float("nan"), pinch_K=min(
                            v["pinch_K"] for v in out.values()))
    return out


def ua_report(case, r):
    """Print the conductance each exchanger needs, and what it buys."""
    ua = exchanger_UA(case, r["rank"], r["hp"], r["T"], r["budget"])
    print("=" * 74)
    print("EXCHANGER CONDUCTANCE  (plant totals, all wells)")
    print("=" * 74)
    print(f"  {'exchanger':<18}{'duty':>10}{'UA':>12}"
          f"{'dT_eff':>10}{'pinch':>9}")
    print(f"  {'':<18}{'kW':>10}{'kW/K':>12}{'K':>10}{'K':>9}")
    print("  " + "-" * 60)
    for name in ("HTHP evaporator", "HTHP condenser",
                 "ORC evaporator", "ORC condenser"):
        v = ua[name]
        print(f"  {name:<18}{v['duty_kW']:>10.1f}{v['UA_kW_per_K']:>12.2f}"
              f"{v['dT_eff_K']:>10.2f}{v['pinch_K']:>9.2f}")
    print("  " + "-" * 60)
    t = ua["TOTAL"]
    print(f"  {'TOTAL':<18}{t['duty_kW']:>10.1f}{t['UA_kW_per_K']:>12.2f}"
          f"{'':>10}{t['pinch_K']:>9.2f}")
    print()
    print("  dT_eff is the single approach that would need the same UA. For a")
    print("  counterflow exchanger with both streams sensible it IS the")
    print("  log-mean difference; with a phase change on one side it is the")
    print("  correct generalisation and the LMTD is not.")
    print()
    print("  Read this beside eta_RTE. An approach temperature bought cheaply")
    print("  in efficiency is paid for here, and until v0.10 nothing said so.")
    return ua


def cycle_state_points(case, T_2d=None):
    """ORC and heat-pump state points, and the two cycle efficiencies.

    Note the expansions and compressions carry no isentropic efficiency, so
    eta_ORC and COP are idealised -- only the electrical and mechanical
    efficiencies appear.

    THE DISCHARGE CLOSURE (changed in v0.6). Only the two exchanger INLETS are
    boundary conditions on the march; both outlets are results of the heat
    transfer and cannot be prescribed. Until v0.5a the discharge was closed the
    other way round -- the outlet was pinned at T_2d = T_m,top - DT_m_2D with
    DT_m_2D = 0, and the inlet derived by subtracting the glide. That is not
    merely arbitrary, it is contradicted by the model's own solution: at the
    start of a discharge at CSS the water leaves at 157.93 C, nearly 8 K above
    the top layer's melting point, because the PCM is
    superheated. It is also silently tied to N_lay -- the implied inlet
    approach was DT_3C_2C / N_lay, one layer width, chosen by nobody.

    The inlets are now prescribed symmetrically:

        charge:    T_4c = T_m,top    + DT_4C_M     (superheat, drives melting)
        discharge: T_3d = T_m,bottom - DT_M_1D     (subcooling, drives freezing)

    and T_2d is DEMOTED to a provisional estimate, the value it would take if
    the water achieved the full glide. It is provisional because the ORC needs
    an evaporating temperature before the store has been marched: T_3e is tied
    to T_2d. Pass the realised outlet back in as `T_2d` for the correction pass
    (see `simulate_css`); the estimate is good enough that one pass suffices.

    Setting DT_M_1D = DT_3C_2C / N_lay recovers the retired closure exactly,
    which is how the change is verified to be a no-op.
    """
    T_3d_C = T_m_bottom(case) - 273.15 - case.DT_M_1D
    # provisional unless the caller supplies the realised value
    T_2d_C = (T_3d_C + case.glide_dc_spec) if T_2d is None else (T_2d - 273.15)
    T = dict(
        T_1e=case.T_sink_C + case.DT_sink_glide + case.DT_E_sink + 273.15,
        T_3e=T_2d_C - case.DT_2D_3E + 273.15,
        T_2d=T_2d_C + 273.15,
        T_3d=T_3d_C + 273.15,                     # the borehole inlet, state 1d
        T_4c=case.T_m_C + case.DT_4C_M + 273.15,
        # The source stream, and the evaporating temperature DERIVED from its
        # outlet. Until v0.9 these were two independent fields and the
        # approach between them was their difference, unchecked: the default
        # pair made it exactly zero, and reversing them made it negative with
        # no symptom anywhere, because T_4a was written here and never read
        # again. See DN-17.
        T_4a=case.T_source_C - case.DT_3A_4A + 273.15,
        T_13h=case.T_source_C - case.DT_3A_4A - case.DT_pinch_HPE + 273.15,
    )
    T["T_3c"] = T["T_4c"]
    T["T_2c"] = T["T_3c"] - case.DT_3C_2C
    T["T_2h"] = T["T_3c"] + case.DT_2H_3C

    # T["T_3e"] above is only the OLD hot-end rule, kept so the uncorrected
    # value stays visible. What the water can actually deliver is decided by
    # the whole evaporator, so the boiling temperature comes back from the
    # pinch check and is written into T -- otherwise the reported state points
    # and the efficiency would describe different cycles. This runs on the
    # correction pass too, with the realised T_2d.
    rank = feasible_rankine(case, T["T_1e"], T["T_2d"], T["T_3d"])
    T["T_3e_hot_end_rule"] = T["T_3e"]
    T["T_3e"] = rank["T_3e"]
    hp = two_stage_htheatpump_2regs(case.refrig, T["T_13h"], T["T_2h"],
                                    case.DT_sub, case.eta_comp_s)
    for name, v in (("rank_eff", rank.get("rank_eff")),
                    ("hp_cop", hp.get("hp_cop"))):
        if v is None or not np.isfinite(v):
            raise RuntimeError(f"{name} is not finite ({v})")
    return rank, hp, T


def water_h(case, T_K):
    """Specific enthalpy of the heat-transfer water, kJ/kg."""
    return CP.PropsSI("H", "T", T_K, "P", case.P, case.fluid2) / 1000.0


def water_s(case, T_K):
    """Specific entropy of the heat-transfer water, kJ/kg-K."""
    return CP.PropsSI("S", "T", T_K, "P", case.P, case.fluid2) / 1000.0


def stream_exergy_rate(case, m_dot, T_hot, T_cold):
    """Exergy given up by a water stream cooling from T_hot to T_cold, kW.

        m_dot [ (h_hot - h_cold) - T_0 (s_hot - s_cold) ],   T_0 = T_sink

    Positive when T_hot > T_cold >= T_0. Used for the process water and for
    the geothermal stream, which are the same fluid at the same pressure.
    Passing T_cold = T_0 gives the stream's exergy relative to the dead
    state, which is what the 'resource' convention charges for.
    """
    T_0 = case.T_sink_C + 273.15
    return m_dot * ((water_h(case, T_hot) - water_h(case, T_cold))
                    - T_0 * (water_s(case, T_hot) - water_s(case, T_cold)))


def energy_budget(case, rank_eff, hp_cop, T):
    """Work backwards from the specified electrical output to the charging duty."""
    # rank_eff and hp_cop now ALREADY carry the machine losses, so only
    # the electrical efficiencies remain here. Dividing by Turb_eff as
    # well would double-count them.
    W_dot_T = case.W_dot_el_out / case.ElG_eff
    Q_dot_in_ORC = W_dot_T / rank_eff
    D_E_in_ORC = Q_dot_in_ORC * case.t_dc * 3600.0                 # kJ
    D_E_out_HP = D_E_in_ORC * (1.0 + case.loss_surplus)
    Q_dot_out_HP = D_E_out_HP / case.t_ch / 3600.0
    W_dot_C_HP = Q_dot_out_HP / hp_cop
    W_dot_el_in = W_dot_C_HP / case.ElH_eff

    # v0.13: the ENTHALPY difference, not c_p at the mean temperature times
    # the glide. Over the 105 -> 160 C charging rise c_p varies about 3 %, so
    # the linearised form gives a flow rate 0.096 % inconsistent with the
    # model's own CoolProp data. That is invisible in the energy balance --
    # both sides use the same wrong flow -- but it fabricates 7.17 kW of
    # exergy and stopped `exergy_audit` closing. The identical formula was
    # written a SECOND time in `run_cycle`, where it is live; see the note
    # there. DN-24.
    m_dot_w_ch = Q_dot_out_HP / (water_h(case, T["T_3c"])
                                 - water_h(case, T["T_2c"]))
    return dict(Q_dot_in_ORC=Q_dot_in_ORC, D_E_in_ORC=D_E_in_ORC,
                D_E_out_HP=D_E_out_HP, Q_dot_out_HP=Q_dot_out_HP,
                W_dot_el_in=W_dot_el_in, m_dot_w_ch=m_dot_w_ch)


def melting_temperatures(case):
    """Layer melting temperatures for charging and for discharging."""
    dTm = case.cascade_param / case.N_lay
    i = np.arange(case.N_lay, dtype=float)
    return ((case.T_m - i * dTm).tolist(),
            (case.T_m - case.cascade_param + (i + 1.0) * dTm).tolist())


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


def simulate_css(case, N=None, n_cycles=80, tol=1e-9, record=False, T_2d=None):
    """March a specified well field to cyclic steady state.

    One pass. `T_2d` overrides the provisional exchanger outlet used to set the
    ORC evaporating temperature; `simulate_css_corrected` is the wrapper that
    supplies the realised value. Called directly, this is the uncorrected
    result.

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

    rank, hp, T = cycle_state_points(case, T_2d=T_2d)
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
    # v0.13: enthalpy differences, not c_p at the mean temperature times the
    # glide. This is the LIVE copy of the formula fixed in `energy_budget`;
    # the two were written independently and both were wrong the same way,
    # which is why patching only one of them produced an encouraging no-op
    # and left the defect in place. See DN-24.
    m1_ch = (Eb["Q_dot_out_HP"] / (water_h(case, T["T_3c"])
                                   - water_h(case, T["T_2c"]))
             / (N * case.num_tubes))
    # NOTE the discharge flow stays pinned by the SPECIFIED glide, as it has
    # been since v0.6, and not by the realised T_2d that the correction pass
    # knows about. Using the realised outlet here is tempting -- it cuts the
    # reported deviation from +4.46 % to +1.45 % and lifts eta_RTE by 0.28 %
    # -- but that is partly circular: feeding last pass's outlet into this
    # pass's flow definition drives the delivered energy toward the target by
    # construction, without iterating to a fixed point. It would also change
    # the model's SPECIFICATION ("the flow is pinned by the glide"), which is
    # a different decision from fixing a linearisation. Measured and left
    # alone. See DN-25.
    m1_dc = (Eb["Q_dot_in_ORC"]
             / (water_h(case, T["T_3d"] + case.glide_dc_spec)
                - water_h(case, T["T_3d"]))
             / (N * case.num_tubes))

    # ---- march until the state repeats -----------------------------------
    def cycle(E0, rec=False):
        ch = march_h(case, T["T_4c"], T_m_seg, m1_ch, case.k_wall, times_ch,
                     n_segments=n, E0=E0, record=rec)
        dc = march_h(case, T["T_3d"], T_m_seg_dc, m1_dc, case.k_wall, times_dc,
                     n_segments=n, E0=ch["E"][::-1], record=rec)
        return ch, dc

    E = np.zeros(n)
    history = []
    for k in range(n_cycles):
        ch, dc = cycle(E)
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

    if record:
        # One extra cycle from the converged state, this time recording. The
        # state repeats to `tol`, so this reproduces the last cycle; doing it
        # here rather than inside the loop means the profile sweep that closes
        # each recorded march is paid once instead of `cycles` times.
        ch, dc = cycle(E, rec=True)

    # The discharge marches from the far end, so its segment index runs
    # backwards. Everything returned to the caller is in CHARGE (depth)
    # indexing, which is why there is one cascade, `T_m_seg`, and not two.
    dc = unmirror_march(dc)

    req = Eb["D_E_in_ORC"]
    deviation = Q_dc / req - 1.0

    # ---- the outlet temperatures, which are results and not inputs ---------
    T_2d_realised = mixed_mean_outlet(dc)      # what the ORC evaporator sees
    T_2c_realised = mixed_mean_outlet(ch)      # what returns to the HP condenser
    glide_dc = T_2d_realised - T["T_3d"]
    glide_ch = T["T_4c"] - T_2c_realised
    return dict(
        N_wells=N, m1_ch=m1_ch, m1_dc=m1_dc, flow_ratio=m1_dc / m1_ch,
        cycles=len(history), converged=history[-1]["drift"] < tol, tol=tol,
        history=history, charge=ch, discharge=dc, E_css=E,
        Q_charge_kJ=Q_ch, Q_discharge_kJ=Q_dc, required_kJ=req,
        deviation=deviation,
        W_el_out_implied=case.W_dot_el_out * (1.0 + deviation),
        eta_storage=Q_dc / Q_ch, lambda_overridden=lambda_overridden,
        T_2d_assumed=T["T_2d"], T_2d_realised=T_2d_realised,
        T_2c_assumed=T["T_2c"], T_2c_realised=T_2c_realised,
        glide_dc=glide_dc, glide_ch=glide_ch,
        glide_ratio_dc=glide_dc / case.glide_dc_spec,
        T_2d_was_provisional=T_2d is None,
        T=T, budget=Eb, T_m_seg=T_m_seg, rank=rank, hp=hp)
        # no `T_m_seg_dc`: the discharge is returned in depth indexing, where
        # the cascade is the same physical object as on charge. The mirrored
        # copy is an internal detail of the march and does not leave here.


def geothermal_resource(case, r):
    """The 60 C source as a FINITE, priced resource.

    Until v0.13 the source was effectively unlimited: `DT_3A_4A` said how far
    the stream was cooled, the flow followed from the duty, and no output ever
    recorded how much flow that was. Nothing in the model objected to a design
    that quietly demanded five times the geothermal field.

    The flow is DERIVED, not prescribed -- the same choice already made for
    `N_wells`. Prescribing a plant-level flow would put the uncertainty in a
    quantity nobody has measured; deriving it puts the uncertainty in the
    per-well yield, which is field data. See DN-22.

    Note that the exergy of the source is homogeneous of degree one in the
    flow, and m_geo = Q_src / (c_p dT), so BOTH conventions below reduce to
    the duty, the glide and the two temperatures. Deriving the flow rather
    than prescribing it costs the second-law analysis nothing.
    """
    hp, T = r["hp"], r["T"]
    T_0 = case.T_sink_C + 273.15
    T_src = case.T_source_C + 273.15
    T_rei = T["T_4a"]                          # T_source - DT_3A_4A

    x = hp["x_8h"]
    m_ref = r["budget"]["Q_dot_out_HP"] / (hp["h_2h"] - hp["h_3h"])
    Q_src = m_ref * (1.0 - x) * (hp["h_13h"] - hp["h_4h"])          # kW
    m_geo = Q_src / (water_h(case, T_src) - water_h(case, T_rei))   # kg/s

    ex_full = stream_exergy_rate(case, m_geo, T_src, T_0)
    ex_rei = stream_exergy_rate(case, m_geo, T_rei, T_0)
    ex_stripped = stream_exergy_rate(case, m_geo, T_src, T_rei)

    N_geo = m_geo / case.m_dot_geo_well
    return dict(
        Q_src_kW=Q_src, m_geo=m_geo, m_ref=m_ref,
        T_reinject_C=T_rei - 273.15,
        Ex_full_kW=ex_full, Ex_reinject_kW=ex_rei, Ex_stripped_kW=ex_stripped,
        utilisation=ex_stripped / ex_full if ex_full else float("nan"),
        N_geo=N_geo, N_geo_per_MWe=N_geo / (case.W_dot_el_out / 1000.0),
        # which of the two the audit will charge for
        Ex_charged_kW=(ex_full if case.exergy_convention == "resource"
                       else ex_stripped),
        Ex_loss_kW=(ex_rei if case.exergy_convention == "resource" else 0.0),
        convention=case.exergy_convention)


def borehole_pumping(case, r):
    """Field pumping power for both half-cycles, kW, at the realised flows.

    ONE copy, called by `performance_indices` and by `exergy_audit`. It was
    briefly written twice -- which is the DN-24 trap exactly, a formula
    duplicated so that a later fix can land in one copy and miss the other.
    """
    gv = case.geom_vector()
    N, T = r["N_wells"], r["T"]
    _, pp_dc = calculate_pressure_drop(gv, r["m1_dc"], case.fluid2,
                                       r["T_2d_realised"], T["T_3d"], case.P)
    _, pp_ch = calculate_pressure_drop(gv, r["m1_ch"], case.fluid2,
                                       T["T_4c"], r["T_2c_realised"], case.P)
    f = case.num_tubes * N / 1000.0
    return pp_ch * f, pp_dc * f


def exergy_audit(case, r):
    """Exergy destruction by component over one cycle at CSS, in kWh.

    THE DEAD STATE is T_0 = T_sink. It is the lowest temperature in the
    system, so every exergy is positive, and the heat dumped at the ORC
    condenser is genuinely unrecoverable. It is NOT the heat-pump evaporator
    source: that is the 60 C geothermal stream, a second and quite separate
    reservoir, which is why `geothermal_resource` has to price it.

    THE CLOSURE, and why this function is a gate rather than a report:

        W_el_in + Ex_geo  =  W_el_out + Ex_reinject + sum(I_j)

    The left and right sides are computed from DIFFERENT things. The boundary
    terms come from four stream states. Every I_j comes from the entropy
    generated inside one component, from its own state points and flow
    fractions. If the cycle state points, the flow splits and the duties are
    mutually consistent, the two agree to machine precision; if any one of
    them is wrong, they do not. That is what caught the c_p linearisation in
    the water flow rates (DN-24), which was invisible to every energy balance
    in the model because both sides of those balances used the same wrong
    flow.

    WHAT IS A DESTRUCTION AND WHAT IS A LOSS. The reinjected geothermal
    stream leaves the boundary intact -- its exergy is not destroyed, it is
    discarded. It is reported on its own line and attributed to NO component,
    because no component is responsible for it; the decision to size the
    field this way is. Under the 'stripped' convention it is not charged at
    all and the line is zero.

    THE BOREHOLE is one bucket here, obtained from the difference of the two
    water-stream exergies over a full cycle. That is exact at CSS -- the PCM
    returns to its own initial state, so its net exergy change is zero, and no
    PCM entropy function is needed. It does not say WHERE in the well or WHEN
    the destruction happens; that needs s(E') and is the next step.
    """
    T_0 = case.T_sink_C + 273.15
    hp, rank, T, Eb = r["hp"], r["rank"], r["T"], r["budget"]
    t_ch, t_dc = case.t_ch, case.t_dc                       # hours
    geo = geothermal_resource(case, r)
    m_ref = geo["m_ref"]
    x = hp["x_8h"]

    # ---- heat pump, per kg of condenser flow ------------------------------
    s = lambda k: hp["s_" + k] / 1000.0                     # kJ/kg-K
    h = lambda k: hp["h_" + k]
    gen = {}
    gen["HTHP compressor HP"] = 1.0 * (s("2h") - s("6h"))
    gen["HTHP compressor LP"] = (1.0 - x) * (s("5h") - s("1h"))
    gen["HTHP throttle 7h-8h"] = 1.0 * (s("8h") - s("7h"))
    gen["HTHP throttle 9h-4h"] = (1.0 - x) * (s("4h") - s("9h"))
    gen["HTHP flash separator"] = (x * s("10h") + (1.0 - x) * s("9h")
                                   - s("8h"))
    gen["HTHP IHX-1"] = 1.0 * (s("12h") - s("3h")) + x * (s("11h") - s("10h"))
    gen["HTHP IHX-2"] = (1.0 * (s("7h") - s("12h"))
                         + (1.0 - x) * (s("1h") - s("13h")))
    m_w_ch = Eb["m_dot_w_ch"]
    gen["HTHP condenser"] = (1.0 * (s("3h") - s("2h"))
                             + m_w_ch * (water_s(case, T["T_3c"])
                                         - water_s(case, T["T_2c"])) / m_ref)
    gen["HTHP evaporator"] = ((1.0 - x) * (s("13h") - s("4h"))
                              + geo["m_geo"]
                              * (water_s(case, T["T_4a"])
                                 - water_s(case, case.T_source_C + 273.15))
                              / m_ref)
    I = {k: T_0 * m_ref * v * t_ch for k, v in gen.items()}     # kWh/cycle

    W_comp = m_ref * ((h("2h") - h("6h")) + (1.0 - x) * (h("5h") - h("1h")))
    W_el_in = Eb["W_dot_el_in"]
    I["HTHP motor"] = (W_el_in - W_comp) * t_ch

    # ---- ORC, per kg through the first turbine stage ----------------------
    y = rank["y_frac"]
    se = lambda k: rank["s_" + k] / 1000.0
    he = lambda k: rank["h_" + k]
    m_orc = Eb["Q_dot_in_ORC"] / ((he("3e") - he("10e"))
                                  + y * (he("6e") - he("5e")))
    ge = {}
    ge["ORC turbine 1"] = 1.0 * (se("5e") - se("3e"))
    ge["ORC turbine 2"] = y * (se("4e") - se("6e"))
    ge["ORC pump 1"] = y * (se("2e") - se("1e"))
    ge["ORC pump 2"] = (1.0 - y) * (se("8e") - se("7e"))
    ge["ORC regenerator"] = ((1.0 - y) * (se("7e") - se("5e"))
                             + y * (se("9e") - se("2e")))
    ge["ORC mixer"] = se("10e") - (y * se("9e") + (1.0 - y) * se("8e"))
    m_w_dc = Eb["Q_dot_in_ORC"] / (water_h(case, T["T_2d"])
                                   - water_h(case, T["T_3d"]))
    ge["ORC evaporator"] = ((1.0 * (se("3e") - se("10e"))
                             + y * (se("6e") - se("5e")))
                            + m_w_dc * (water_s(case, T["T_3d"])
                                        - water_s(case, T["T_2d"])) / m_orc)
    # the sink is AT the dead state, so the rejected heat carries no exergy
    # and the whole of it is destroyed in the condenser
    Q_rej = m_orc * y * (he("4e") - he("1e"))
    ge["ORC condenser"] = y * (se("1e") - se("4e")) + Q_rej / T_0 / m_orc
    for k, v in ge.items():
        I[k] = T_0 * m_orc * v * t_dc

    W_turb = m_orc * ((he("3e") - he("5e")) + y * (he("6e") - he("4e"))
                      - y * (he("2e") - he("1e"))
                      - (1.0 - y) * (he("8e") - he("7e")))
    I["ORC generator"] = (W_turb - case.W_dot_el_out) * t_dc

    # ---- THE TWO GATES THAT CAN ACTUALLY FAIL -----------------------------
    # Each machine bucket is summed from entropy generated INSIDE its
    # components and checked against a boundary balance built from stream
    # states only. These are independent computations and they are what
    # caught DN-24.
    I_hp = sum(v for k, v in I.items() if k.startswith("HTHP"))
    I_orc = sum(v for k, v in I.items() if k.startswith("ORC"))
    dEx_ch = stream_exergy_rate(case, m_w_ch, T["T_3c"], T["T_2c"]) * t_ch
    dEx_dc = stream_exergy_rate(case, m_w_dc, T["T_2d"], T["T_3d"]) * t_dc
    bnd_hp = (W_comp * t_ch + geo["Ex_stripped_kW"] * t_ch - dEx_ch
              + (W_el_in - W_comp) * t_ch)
    bnd_orc = dEx_dc - W_turb * t_dc + (W_turb - case.W_dot_el_out) * t_dc
    gates = {"HTHP": (I_hp, bnd_hp), "ORC": (I_orc, bnd_orc)}

    # ---- borehole: BY DIFFERENCE, and therefore not a gate ----------------
    # T_3c = T_4c and T_2d/T_3d are shared with the surface exchangers, so
    # this term is algebraically the residual of the water loop: it cancels
    # the two exchanger terms exactly and the GLOBAL sum below is an
    # identity, not a test. An earlier draft printed that identity as a
    # closure check reading 0.00e+00 and called it a pass. It cannot fail,
    # so it cannot pass. The independent value comes from the resolved
    # integral over z and t once s(E') exists; until then this number is
    # correct but unverified.
    I["borehole"] = dEx_ch - dEx_dc

    # The water loop does NOT in fact close at CSS: the well returns water at
    # T_2c_realised, not at the T_2c the condenser assumes, and likewise on
    # discharge. Reported rather than absorbed silently.
    mism = float("nan")
    if r.get("T_2c_realised") is not None:
        mism = (stream_exergy_rate(case, m_w_ch, T["T_3c"],
                                   r["T_2c_realised"]) * t_ch - dEx_ch
                - (stream_exergy_rate(case, m_w_dc, r["T_2d_realised"],
                                      T["T_3d"]) * t_dc - dEx_dc))

    # ---- REALISED BASIS, and the same boundary as eta_RTE ----------------
    # Until v0.15 this function evaluated the whole audit on the TARGET
    # duties out of `energy_budget`, and charged no pumping. `eta_RTE` in
    # `performance_indices` does neither: it uses the energy the field
    # actually moved, and it charges the borehole parasitics. So psi and
    # eta_RTE were computed on different bases and different boundaries, and
    # quoting "0.317 becomes 0.246" compared two things that were not
    # like for like. The manuscript caught this before the code did.
    #
    # At CSS the storage efficiency is identically 1 and lambda is forced to
    # zero, so Q_charge = Q_discharge and D_E_out_HP = D_E_in_ORC: both half
    # cycles scale by the SAME factor, and the rebase is one multiplication.
    scale = r["Q_discharge_kJ"] / Eb["D_E_in_ORC"]
    for key in I:
        I[key] *= scale
    W_in = W_el_in * t_ch * scale
    W_out = case.W_dot_el_out * t_dc * scale
    Ex_geo = geo["Ex_charged_kW"] * t_ch * scale
    Ex_loss = geo["Ex_loss_kW"] * t_ch * scale

    # The borehole parasitics are work in, dissipated as friction in the
    # well, so they are an input AND a destruction. Charged with eta_RTE's
    # sign convention: the charge pump adds to the input, the discharge pump
    # subtracts from the output. They are computed at the realised flows
    # already, so they are not scaled.
    pump_ch, pump_dc = borehole_pumping(case, r)
    P_ch, P_dc = pump_ch * t_ch, pump_dc * t_dc
    I["borehole pumping"] = P_ch + P_dc
    W_in += P_ch
    W_out -= P_dc

    boundary = W_in + Ex_geo - W_out - Ex_loss
    total = sum(I.values())

    # A THIRD GATE, and the only one that spans two functions. eta_RTE is
    # computed here from the scaled budget plus the parasitics, and in
    # `performance_indices` from the energy the field actually moved. Two
    # routes, one answer. It is what would have caught the basis mismatch
    # that the manuscript found before the code did, and it can only hold
    # while psi and eta_RTE share a boundary -- which is the point.
    # Recomputed the way `performance_indices` does it -- from the energy the
    # field moved and the same parasitics -- rather than from the budget.
    Q_in_ORC = r["Q_discharge_kJ"] / (t_dc * 3600.0)
    Q_out_HP = r["Q_charge_kJ"] / (t_ch * 3600.0)
    rk_eff = case.W_dot_el_out / case.ElG_eff / Eb["Q_dot_in_ORC"]
    cop_r = Eb["Q_dot_out_HP"] / (Eb["W_dot_el_in"] * case.ElH_eff)
    W_out_pi = (Q_in_ORC * rk_eff * case.ElG_eff - pump_dc) * t_dc
    W_in_pi = (Q_out_HP / cop_r / case.ElH_eff + pump_ch) * t_ch
    gates["eta_RTE"] = (W_out / W_in * W_in_pi, W_out_pi)

    return dict(
        T_0=T_0, components_kWh=I, total_kWh=total, boundary_kWh=boundary,
        scale=scale, pumping_ch_kWh=P_ch, pumping_dc_kWh=P_dc,
        gates={k: dict(sum_kWh=a, boundary_kWh=b, residual_kWh=a - b,
                       residual_rel=abs(a - b) / abs(b) if b else float("nan"))
               for k, (a, b) in gates.items()},
        gates_pass=all(abs(a - b) <= 1e-9 * abs(b) for a, b in gates.values()),
        global_is_identity=True,
        loop_mismatch_kWh=mism,
        W_el_in_kWh=W_in, W_el_out_kWh=W_out,
        Ex_geo_kWh=Ex_geo, Ex_reinject_kWh=Ex_loss,
        eta_RTE=W_out / W_in,
        psi=W_out / (W_in + Ex_geo),
        psi_stripped=W_out / (W_in + geo["Ex_stripped_kW"] * t_ch * scale),
        shares={k: v / total for k, v in I.items()},
        geo=geo)


def exergy_report(case, r, a=None):
    """Print the destruction table, largest first, and the closure check."""
    a = exergy_audit(case, r) if a is None else a
    g = a["geo"]
    print(f"EXERGY AUDIT   dead state T_0 = {a['T_0'] - 273.15:.1f} C "
          f"(the sink), convention '{g['convention']}'")
    print(f"  IN   electrical           {a['W_el_in_kWh']:10.1f} kWh")
    print(f"  IN   geothermal           {a['Ex_geo_kWh']:10.1f} kWh   "
          f"({g['m_geo']:.1f} kg/s, {g['N_geo']:.2f} producers)")
    print(f"  OUT  electrical           {a['W_el_out_kWh']:10.1f} kWh")
    print(f"  OUT  reinjected (LOSS)    {a['Ex_reinject_kWh']:10.1f} kWh   "
          f"at {g['T_reinject_C']:.1f} C, "
          f"{100 * (1 - g['utilisation']):.1f} % of the resource discarded")
    print(f"  {'-' * 56}")
    for k, v in sorted(a["components_kWh"].items(), key=lambda kv: -kv[1]):
        print(f"    {k:<26}{v:10.1f} kWh{100 * a['shares'][k]:7.1f} %")
    print(f"  {'-' * 56}")
    print(f"    {'SUM of components':<26}{a['total_kWh']:10.1f} kWh")
    print(f"    {'boundary balance':<26}{a['boundary_kWh']:10.1f} kWh"
          "   (identity -- see below)")
    print("\n  GATES. Only these can fail. Each sums the entropy generated "
          "inside\n  one machine's components and checks it against a "
          "balance built from\n  stream states alone.")
    for k, g in a["gates"].items():
        ok = g["residual_rel"] <= 1e-9
        print(f"    {k:<8}{g['sum_kWh']:11.2f} vs {g['boundary_kWh']:11.2f} "
              f"kWh   rel {g['residual_rel']:.2e}   "
              f"{'OK' if ok else '*** FAIL ***'}")
    print("    borehole  by difference -- the water-loop residual. It cancels"
          "\n              the two exchanger terms exactly, so the global sum"
          "\n              above is an IDENTITY and cannot fail. Verified only"
          "\n              when the resolved s(E') integral lands.")
    if a["loop_mismatch_kWh"] == a["loop_mismatch_kWh"]:
        print(f"    loop      {a['loop_mismatch_kWh']:+.1f} kWh between the "
              "outlets the surface assumes\n              and the outlets the "
              "well actually delivers at CSS")
    print(f"\n  eta_RTE {a['eta_RTE']:.4f}    psi (resource) {a['psi']:.4f}"
          f"    psi (stripped) {a['psi_stripped']:.4f}")
    return a


def performance_indices(case, r):
    """The performance indicators of the factorial study, evaluated at CSS.

    These were defined for the first-cycle framing, where the field delivered
    the target energy by construction. At cyclic steady state it delivers
    whatever it delivers, so every one of them is computed from the ENERGY THE
    FIELD ACTUALLY MOVES rather than from the target. That matters: charging
    with the target input while crediting the realised output would flatter the
    round-trip efficiency by exactly the deviation.

        eta_T       thermal -> net electric conversion on the discharge side,
                    W_el_out / Q_in_ORC. Identically eta_ORC * eta_turb * eta_gen,
                    so it is a property of the power block alone.
        eta_T_eff   the same with the discharge pumping parasitic charged
                    against it. This is the "effective discharge efficiency".
        eta_RTE     round-trip: net electricity out over electricity in, both
                    including their pumping parasitics.
        eps_RTE     eta_RTE weighted by the fraction of the store that actually
                    CYCLES, eps(end of charge) - eps(end of discharge).
        dE_therm    thermal energy delivered per well per cycle   [kWh]
        dE_elec     net electric energy delivered per well,
                    dE_therm * eta_T_eff                          [kWh]

    A note on eps_RTE, because it changed meaning and the change matters when
    reading a factorial table. The first-cycle version weighted by eps_m, the
    melted fraction at end of charge. Under the enthalpy formulation that
    SATURATES -- once a segment has melted all its PCM the further energy goes
    into superheat, which eps_m cannot see -- and at this design point it is
    1.0000 exactly, so eps_RTE was identical to eta_RTE and carried no
    information at all. The cycled fraction does not saturate: a store that
    fills completely but only half empties returns 0.5, which is the useful
    statement.

    Note also that eta_RTE is INDEPENDENT of the deviation, exactly. A field
    that delivers 1 % more than target also drew 1 % more in, because the charge
    flow is pinned by its glide in the same way. The deviation is a statement
    about plant size, not about efficiency.
    """
    T, Eb = r["T"], r["budget"]
    N = r["N_wells"]
    t_ch_s, t_dc_s = case.t_ch * 3600.0, case.t_dc * 3600.0

    # ---- what the field actually moved, per cycle ------------------------
    Q_dot_in_ORC = r["Q_discharge_kJ"] / t_dc_s          # kW, field
    Q_dot_out_HP = r["Q_charge_kJ"] / t_ch_s             # kW, field
    rank_eff = Eb["Q_dot_in_ORC"] and (case.W_dot_el_out / case.ElG_eff
                                       / Eb["Q_dot_in_ORC"])
    W_el_out = Q_dot_in_ORC * rank_eff * case.ElG_eff
    cop = Eb["Q_dot_out_HP"] / (Eb["W_dot_el_in"] * case.ElH_eff)
    W_el_in = Q_dot_out_HP / cop / case.ElH_eff

    # ---- parasitics, at the realised temperatures ------------------------
    pumping_ch, pumping_dc = borehole_pumping(case, r)

    # ---- the indicators --------------------------------------------------
    eta_T = W_el_out / Q_dot_in_ORC
    eta_T_eff = eta_T - pumping_dc / Q_dot_in_ORC
    eta_RTE = ((W_el_out - pumping_dc) * case.t_dc
               / ((W_el_in + pumping_ch) * case.t_ch))
    eps_cycled = float(r["charge"]["eps_local"].mean()
                       - r["discharge"]["eps_local"].mean())
    dE_therm = r["Q_discharge_kJ"] / 3600.0 / N          # kWh per well
    out = dict(
        eta_T=eta_T, eta_T_eff=eta_T_eff,
        eta_RTE=eta_RTE, eps_RTE=eta_RTE * eps_cycled,
        eps_cycled=eps_cycled,
        dE_therm_kWh=dE_therm, dE_elec_kWh=dE_therm * eta_T_eff,
        pumping_ch_kW=pumping_ch, pumping_dc_kW=pumping_dc,
        f_pump=(pumping_ch + pumping_dc) / W_el_out,
        W_el_out_kW=W_el_out, W_el_in_kW=W_el_in,
        rho_E_kWh_m3=dE_therm / case.V_well,
        # ---- POWER density, which is not energy density -----------------
        # rho_E says how much energy a cubic metre of store holds; it says
        # nothing about how fast it can be delivered, and the two are traded
        # against each other by the discharge window. A store with twice the
        # energy density and half the discharge rate is a different machine.
        # Per unit DEPTH as well, because drilling is priced by the metre.
        rho_P_kW_m3=(W_el_out / N) / case.V_well,
        rho_P_kW_per_m=(W_el_out / N) / case.L_well,
        t_discharge_h=case.t_dc)

    # ---- what the efficiencies above cost in hardware --------------------
    # Without these, every approach temperature is a free lunch: tighten
    # DT_pinch_ORC from 10 K to 1 K and eta_RTE rises 11 % with nothing to
    # pay, while the evaporator conductance more than doubles. Reported per
    # unit of delivered electricity as well, which is the form that lets two
    # designs of different size be compared.
    if "rank" in r and "hp" in r:
        ua = exchanger_UA(case, r["rank"], r["hp"], r["T"], r["budget"])
        out["UA_total_kW_K"] = ua["TOTAL"]["UA_kW_per_K"]
        out["UA_per_MWe"] = ua["TOTAL"]["UA_kW_per_K"] / (W_el_out / 1000.0)
        out["pinch_min_K"] = ua["TOTAL"]["pinch_K"]
        for key, tag in (("HTHP evaporator", "hpe"), ("HTHP condenser", "hpc"),
                         ("ORC evaporator", "orce"), ("ORC condenser", "orcc")):
            out[f"UA_{tag}_kW_K"] = ua[key]["UA_kW_per_K"]
            out[f"dTeff_{tag}_K"] = ua[key]["dT_eff_K"]
            # The four pinches individually, not just their minimum. Which
            # exchanger is tightest is design information; `pinch_min_K`
            # throws it away, and a design can move the binding one without
            # the minimum changing at all.
            out[f"pinch_{tag}_K"] = ua[key]["pinch_K"]

        # ---- second law, and the resource the first law books as free ----
        # eta_RTE treats the 60 C geothermal stream as costless. It is not:
        # a dedicated producer is drilled and paid for whether or not its
        # exergy is used. psi prices it, and comes out well BELOW eta_RTE.
        # The well count is the number the operator actually cares about,
        # and until v0.13 the model reported only half of it.
        ex = exergy_audit(case, r)
        out["psi"] = ex["psi"]
        out["psi_stripped"] = ex["psi_stripped"]
        out["I_total_kWh"] = ex["total_kWh"]
        out["Ex_reinject_kWh"] = ex["Ex_reinject_kWh"]
        out["resource_utilisation"] = ex["geo"]["utilisation"]
        out["m_geo_kg_s"] = ex["geo"]["m_geo"]
        out["N_geo"] = ex["geo"]["N_geo"]
        out["T_reinject_C"] = ex["geo"]["T_reinject_C"]
        N_st = r.get("N_wells")
        if N_st is not None:
            out["N_wells_total"] = N_st + ex["geo"]["N_geo"]
            out["wells_per_MWe"] = ((N_st + ex["geo"]["N_geo"])
                                    / (W_el_out / 1000.0))
        for name, v in ex["components_kWh"].items():
            out["I_" + name.replace(" ", "_").replace("-", "")] = v
        out["exergy_gates_pass"] = ex["gates_pass"]

        # ---- per-well productivity, and one scalar for the whole audit ----
        # rho_P says power per cubic metre of store; this says power per WELL,
        # which is the number an operator compares against a drilling cost.
        if N_st is not None:
            out["kW_per_well"] = W_el_out / N_st
        out["I_per_MWh"] = ex["total_kWh"] / (ex["W_el_out_kWh"] / 1000.0)

    # ---- REGIME LABELS, not performance -------------------------------
    # These do not say whether a design is good. They say whether two runs
    # are the same plant, and a screen that averages across them will report
    # main effects taken over a regime change. A run whose ORC switched from
    # pinch-bound to hot-end-bound is a different question, not a worse
    # answer. Carried in the KPI dict so a factorial or a Morris screen
    # cannot lose them.
    out["deviation"] = r["deviation"]
    out["merge_proximity_max"] = float(max(r["charge"]["merge_proximity_max"],
                                           r["discharge"]["merge_proximity_max"]))
    tb = r["rank"].get("T_3e_binding", "")
    out["T_3e_binding"] = tb
    out["regime_pinch_bound"] = 1.0 if "pinch" in tb else 0.0
    msgs = validate_case(case, verbose=False)
    out["feasible"] = 0.0 if any(k == "error" for k, _ in msgs) else 1.0
    out["n_warnings"] = float(sum(1 for k, _ in msgs if k == "warn"))
    return out


def kpi_report(case, r, k=None):
    """Print the factorial-study indicators."""
    k = k or performance_indices(case, r)
    print("PERFORMANCE INDICATORS   (at cyclic steady state, per cycle)")
    print("=" * 66)
    print(f"  eta_T       discharge thermal -> electric   {k['eta_T']:9.4f}")
    print(f"  eta_T,eff   the same, less discharge pumping{k['eta_T_eff']:9.4f}"
          f"   <- effective")
    print(f"  eta_RTE     round trip, both parasitics     {k['eta_RTE']:9.4f}")
    print(f"  eps_RTE     eta_RTE x cycled fraction       {k['eps_RTE']:9.4f}")
    print()
    print(f"  dE_therm    thermal energy per well      {k['dE_therm_kWh']:12.1f} kWh")
    print(f"  dE_elec     net electric energy per well {k['dE_elec_kWh']:12.1f} kWh")
    print(f"  rho_E       thermal energy density       {k['rho_E_kWh_m3']:12.2f} kWh/m3")
    print(f"  rho_P       electric POWER density      {k['rho_P_kW_m3']:12.4f} kW/m3"
          f"   ({k['t_discharge_h']:.0f} h window)")
    print(f"              per metre of well           {k['rho_P_kW_per_m']:12.4f} kW/m")
    print("      rho_E and rho_P are traded against each other by the")
    print("      discharge window; a store with twice the energy density and")
    print("      half the rate is a different machine, not a better one.")
    print()
    print(f"  cycled fraction of the store   {k['eps_cycled']:8.4f}"
          f"   (melted {r['charge']['eps_local'].mean():.4f},"
          f" residual {r['discharge']['eps_local'].mean():.4f})")
    print(f"  pumping, charge / discharge    {k['pumping_ch_kW']:8.2f}"
          f" / {k['pumping_dc_kW']:.2f} kW"
          f"   = {100*k['f_pump']:.2f} % of gross output")
    if "UA_total_kW_K" in k:
        print()
        print(f"  UA total, all four exchangers  "
              f"{k['UA_total_kW_K']:8.1f} kW/K"
              f"   = {k['UA_per_MWe']:.1f} per MWe")
        print(f"      HTHP evap {k['UA_hpe_kW_K']:7.1f}  "
              f"HTHP cond {k['UA_hpc_kW_K']:7.1f}  "
              f"ORC evap {k['UA_orce_kW_K']:7.1f}  "
              f"ORC cond {k['UA_orcc_kW_K']:7.1f}")
        print(f"  the four pinches, individually (K)")
        print(f"      HTHP evap {k['pinch_hpe_K']:7.2f}  "
              f"HTHP cond {k['pinch_hpc_K']:7.2f}  "
              f"ORC evap {k['pinch_orce_K']:7.2f}  "
              f"ORC cond {k['pinch_orcc_K']:7.2f}")
        print(f"      minimum   {k['pinch_min_K']:7.2f}  <- which exchanger is"
              f" tightest is design information")
        print( "                         that the minimum alone throws away.")
        print("      Read UA beside eta_RTE. Every approach temperature is")
        print("      bought here; `ua_report` breaks it down.")
    print()
    print("  eta_RTE is independent of the deviation: a field that delivers")
    print("  1 % above target also drew 1 % more in. The deviation is a")
    print("  statement about plant SIZE, not about efficiency.")


def simulate_css_corrected(case, N=None, record=False, **kw):
    """CSS with a single ORC correction pass. The v0.6 default path.

    The ORC needs an evaporating temperature before the store has been marched,
    because T_3e is tied to T_2d. Since v0.6 T_2d is a result, so the first pass
    uses the provisional estimate -- the outlet the water would reach if it
    achieved the full glide -- and the second re-evaluates the ORC, the energy
    budget, the well count and both flows at the outlet the first pass actually
    produced.

    ONE pass is enough, and that is a property of the problem rather than a
    hopeful choice: the outlet is pinned by the store, not by the plant. Over a
    sweep of the discharge subcooling from 6 to 15 K -- which moves the inlet by
    9 K and the deviation by twelve percentage points -- the realised outlet
    moves by less than 0.3 K. `dT_2d_pass` is returned so the assumption is
    visible rather than assumed; if it is ever large, iterate.

    This keeps the plant level decoupled from the store in the sense that
    matters for Figure "procedure": there is no outer loop, only a second
    evaluation.
    """
    first = simulate_css(case, N=N, record=False, **kw)
    second = simulate_css(case, N=N, record=record,
                          T_2d=first["T_2d_realised"], **kw)
    second["pass1"] = first
    second["dT_2d_pass"] = second["T_2d_realised"] - first["T_2d_realised"]
    second["corrected"] = True
    return second


def css_report(case, r):
    """Print the cyclic-steady-state result."""
    print("CYCLIC STEADY STATE")
    print("=" * 66)
    if r["lambda_overridden"]:
        print(f"  NOTE  lambda forced to 0, from case.loss_surplus = "
              f"{case.loss_surplus:.3f}. A non-zero loss surplus cannot be")
        print("        absorbed at CSS: the state returns to itself and this")
        print("        model has no loss path, so charge == discharge exactly.")
        print("        The ZERO is what drove this run; the 0.050 in the Case")
        print("        is inert here and affects only energy_budget called on")
        print("        its own.")
        print()
    if r["T"].get("T_3e_hot_end_rule") is not None:
        old = r["T"]["T_3e_hot_end_rule"] - 273.15
        new = r["T"]["T_3e"] - 273.15
        if new < old - 1e-6:
            print(f"  NOTE  ORC boiling temperature set by the EVAPORATOR "
                  f"PINCH, not by DT_2D_3E:")
            print(f"        hot-end rule would give {old:7.2f} C; the pinch "
                  f"allows {new:7.2f} C.")
            print()
    print(f"  N_wells          {r['N_wells']:10.4f}   from latent heat alone, not solved")
    print(f"  m_dot charge     {r['m1_ch']:10.4f} kg/s per leg-pair, pinned by the glide")
    print(f"  m_dot discharge  {r['m1_dc']:10.4f} kg/s per leg-pair, pinned by the glide")
    print(f"  flow ratio       {r['flow_ratio']:10.4f}   a consequence, not a solve")
    print()
    if r["converged"]:
        print(f"  converged in {r['cycles']} cycles"
              f"   (drift {r['history'][-1]['drift']:.1e})")
    else:
        # Loud, because the result is still printed and looks perfectly normal.
        # Slow convergence is not rare: it scales with the well count, and a
        # case with low latent heat or many wells can need three times the
        # default. Re-run with a larger n_cycles before believing anything
        # below this line.
        print(f"  *** NOT CONVERGED *** stopped at the n_cycles limit of"
              f" {r['cycles']} with drift {r['history'][-1]['drift']:.2e},"
              f" tolerance {r['tol']:.0e}")
        print("      Re-run with a larger n_cycles. Everything below assumes")
        print("      the state has returned to itself, and it has not.")
    print(f"  eta_storage      {r['eta_storage']:10.6f}   must be 1: adiabatic, closed cycle")
    print()
    print(f"  required   D_E_in_ORC  {r['required_kJ']:12.5e} kJ")
    print(f"  delivered  at CSS      {r['Q_discharge_kJ']:12.5e} kJ")
    print(f"  DEVIATION              {100*r['deviation']:+12.3f} %")
    print()
    print(f"  implied net output     {r['W_el_out_implied']/1000:10.4f} MWe"
          f"   (target {case.W_dot_el_out/1000:.4f})")
    print()
    print("  the deviation IS the glide shortfall")
    print(f"    realised discharge glide / assumed   {r['glide_ratio_dc']:10.6f}")
    print(f"    delivered energy / required          "
          f"{r['Q_discharge_kJ']/r['required_kJ']:10.6f}")
    print("    The flow is pinned by the ASSUMED glide, so the energy")
    print("    delivered is just the ratio of realised to assumed glide.")
    print()
    print("  exchanger outlets are RESULTS, not inputs (v0.6)")
    print(f"    discharge, to the ORC       assumed {r['T_2d_assumed']-273.15:8.3f} C"
          f"   realised {r['T_2d_realised']-273.15:8.3f} C")
    print(f"    charge, to the HP condenser assumed {r['T_2c_assumed']-273.15:8.3f} C"
          f"   realised {r['T_2c_realised']-273.15:8.3f} C")
    if r.get("corrected"):
        print(f"    ORC corrected at the realised outlet; second pass moved it"
              f" {r['dT_2d_pass']:+.4f} K")
    elif r["T_2d_was_provisional"]:
        print("    UNCORRECTED: the ORC still runs on the provisional outlet.")
        print("    Use simulate_css_corrected for the reported result.")
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
