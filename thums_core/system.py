"""`run_cycle` — the strangler boundary.

This module changes NO physics. It translates:

    Case object            ->  the legacy positional argument lists
    legacy return tuples   ->  named KPIs from thums_core.kpi.REGISTRY

Everything inside still calls `thums_core._legacy_physics`, i.e. the code that
produced the conference-paper numbers, unmodified. That is deliberate: while the
port is in progress this function must reproduce the frozen fixture exactly, so
that when a ported module later changes a number, the change is attributable to
that module and nothing else.

Once `network.py`, `stefan.py`, `well.py` and `cycle.py` exist, their calls are
substituted here one at a time, checking the fixture at each swap. When the last
legacy call is gone, `_legacy_physics` stops being imported.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from . import _legacy_physics as L
from . import kpi as kpi_mod
from . import sizing
from . import well_marched as W
from .config import Case
from .errors import ThumsConvergenceError, ThumsPropertyError

CP = L.CP


@dataclass
class CycleResult:
    kpis: dict = field(default_factory=dict)
    detail: dict = field(default_factory=dict)
    converged: bool = False
    warnings: list = field(default_factory=list)
    wall_time_s: float = float("nan")


def _cycles(c: Case):
    """HTHP and ORC state points. Wraps the legacy cycle functions."""
    cy, pcm = c.cycle, c.pcm
    T_2d_C = pcm.T_m_C - cy.DT_m_2D
    T = dict(
        T_1e=cy.T_sink_C + cy.DT_E_sink + 273.15,
        T_3e=T_2d_C - cy.DT_2D_3E + 273.15,
        T_2d=T_2d_C + 273.15,
        T_3d=T_2d_C - cy.DT_2D_3D + 273.15,
        T_4c=pcm.T_m_C + cy.DT_4C_M + 273.15,
        T_4a=cy.T_source_C - cy.DT_3A_4A + 273.15,
        T_13h=cy.T_source_C - cy.DT_3A_13H + 273.15,
    )
    T["T_3c"] = T["T_4c"]
    T["T_2c"] = T["T_3c"] - cy.DT_3C_2C
    T["T_2h"] = T["T_3c"] + cy.DT_2H_3C

    try:
        if cy.orc_stages == 2:
            rank = L.double_stage_rankine(cy.fluid, T["T_1e"], T["T_3e"])
        else:
            rank = L.single_stage_rankine(cy.fluid, T["T_1e"], T["T_3e"])
        if cy.hthp_stages == 2:
            hp = L.two_stage_htheatpump_2regs(cy.refrig, T["T_13h"], T["T_2h"], cy.DT_sub)
        else:
            hp = L.single_stage_htheatpump(cy.refrig, T["T_13h"], T["T_2h"], cy.DT_sub)
    except ValueError as e:
        raise ThumsPropertyError(f"cycle state points failed: {e}") from e

    for name, v, d in (("rank_eff", rank.get("rank_eff"), rank),
                       ("hp_cop", hp.get("hp_cop"), hp)):
        if v is None or not np.isfinite(v):
            # legacy returns float('inf') here; inf passes every isnan-only guard
            raise ThumsConvergenceError(f"{name} is not finite ({v})")
    return rank, hp, T



def _legacy_well_field(case, *, T, T_m_lay, T_m_lay_dc, m_dot_w_ch, D_E_out_HP,
                       times_ch, times_dc, gv, k_charge_wall, res):
    """The v0.1 path, verbatim. Retained so the two formulations can be run
    side by side rather than argued about; `Case(front="legacy")` selects it."""
    geom, pcm, cy, num = case.geometry, case.pcm, case.cycle, case.numerics
    N_lay = case.N_lay
    lo, hi = num.N_wells_bracket
    N_wells = L.find_N_wells_for_Q_ratio_ch_fast(
        1.0, num.tol_Q_ratio, lo, hi, m_dot_w_ch, geom.num_tubes, gv, times_ch,
        T["T_4c"], N_lay, T_m_lay, k_charge_wall, geom.Rf_i, cy.P, cy.fluid2,
        pcm.cp_l, pcm.rho_l, pcm.h_m, num.n_segments, num.delta_max, D_E_out_HP,
        N_hint=None, max_iterations=num.max_iterations, verbose=False)

    if not np.isfinite(N_wells) or N_wells >= hi - 1.0:
        # The legacy bisection returns the bracket bound and reports it as an
        # answer; cell 22 silently `continue`s. Here it is a recorded failure.
        raise ThumsConvergenceError("well-field sizing (charging)",
                                    iterations=num.max_iterations, last=N_wells)

    m_dot_c_well1 = m_dot_w_ch / N_wells / geom.num_tubes
    rlo, rhi = num.ratio_bracket
    D_E_in_ORC = D_E_out_HP / (1.0 + cy.loss_surplus)
    m_dot_d_well1 = L.find_m_dot_d_well1_for_Q_ratio_dc_fast(
        1.0, num.tol_Q_ratio, rlo, rhi, m_dot_c_well1, N_wells, geom.num_tubes,
        gv, times_dc, T["T_3d"], N_lay, T_m_lay_dc, geom.k_wall, geom.Rf_i,
        cy.P, cy.fluid2, pcm.k_s, pcm.cp_s, pcm.rho_s, pcm.h_m,
        num.n_segments, num.delta_max, D_E_in_ORC,
        ratio_hint=None, max_iterations=num.max_iterations, verbose=False)

    if not np.isfinite(m_dot_d_well1) or m_dot_d_well1 <= 0:
        raise ThumsConvergenceError("discharge flow rate", last=m_dot_d_well1)
    ratio = m_dot_d_well1 / m_dot_c_well1
    if not (rlo * 1.001 < ratio < rhi * 0.999):
        raise ThumsConvergenceError("discharge flow rate hit its bracket",
                                    last=ratio)

    # Melted volume at the end of charging, from the independent Stefan front.
    t_end = [cy.t_ch * 3600.0]
    _, _, V_melt_map, _ = L.time_profiles_melt(
        t_end, gv, T["T_4c"], N_lay, T_m_lay, k_charge_wall, geom.Rf_i,
        m_dot_c_well1, cy.P, cy.fluid2, pcm.k_l, pcm.cp_l, pcm.rho_l, pcm.h_m,
        num.n_segments, num.delta_max)
    V_melt_tube = V_melt_map.get(t_end[0], np.nan) if isinstance(V_melt_map, dict) \
        else float(np.asarray(V_melt_map).ravel()[-1])
    eps_pcm = V_melt_tube * geom.num_tubes / geom.V_well
    if np.isfinite(eps_pcm) and eps_pcm > 1.0:
        res.warnings.append(f"eps_pcm = {eps_pcm:.3f} > 1: melt fronts are not "
                            "bounded by the borehole (audit 2.5)")
    return N_wells, m_dot_d_well1, ratio, eps_pcm, {}


def _layer_map(T_m_lay, n_segments: int, N_lay: int) -> np.ndarray:
    """Per-segment melting temperature, matching the legacy node convention.

    `_legacy_physics.temperature_profile_melt` builds T_m on the n+1 nodes of
    `linspace(0, L, n+1)`, assigns node j the layer containing it, and takes the
    melt area at the segment *ends* (nodes 1..n). This reproduces that, so the
    marched and legacy fronts see an identical PCM layer profile and the
    difference between them is attributable to the front formulation alone.
    """
    z = np.linspace(0.0, 1.0, n_segments + 1)
    z_lay = np.linspace(0.0, 1.0, N_lay + 1)
    T_m = np.empty(n_segments + 1)
    for j in range(N_lay):
        mask = (z >= z_lay[j]) & (z < z_lay[j + 1])
        T_m[mask] = T_m_lay[j]
    T_m[-1] = T_m_lay[-1]
    return T_m[1:]


def _marched_well_field(case, *, T, T_m_lay, T_m_lay_dc, m_dot_w_ch,
                        D_E_out_HP, D_E_in_ORC, times_ch, times_dc,
                        k_charge_wall, res):
    """Size the field and solve the discharge with the marched front.

    Replaces two legacy calls:

        find_N_wells_for_Q_ratio_ch_fast   -> sizing.size_well_field
        find_m_dot_d_well1_for_Q_ratio_dc  -> a bisection over a marched discharge
                                              that STARTS FROM the charge state

    The second substitution is the one that changes the character of the model.
    v0.1 discharges from a bare tube in fully solid PCM regardless of what
    charging produced, so the energy recovered is not bounded by the energy
    stored. Here `state=r_ch.state` carries the melt distribution across, and
    `stefan.advance` will not freeze PCM that is already solid.
    """
    geom, pcm, cy, num = case.geometry, case.pcm, case.cycle, case.numerics
    n_seg = num.n_segments
    N_lay = case.N_lay

    T_m_seg_ch = _layer_map(T_m_lay, n_seg, N_lay)
    T_m_seg_dc = _layer_map(T_m_lay_dc, n_seg, N_lay)

    # ---- charging: N = max(heat transfer, inventory) ------------------------
    sz = sizing.size_well_field(
        geom=geom, pcm=pcm, cy=cy, num=num, T_inlet=T["T_4c"],
        T_m_seg=T_m_seg_ch, m_dot_total=m_dot_w_ch,
        D_E_required_kJ=D_E_out_HP, times=times_ch, k_wall=k_charge_wall,
        n_segments=n_seg)
    N_wells = sz.N_wells
    res.warnings.extend(sz.warnings)
    res.warnings.append(
        f"sizing: N_heat={sz.N_heat:.3f}, N_inventory={sz.N_inventory:.3f}, "
        f"binding={sz.binding}")

    # ---- the charge end-state at the chosen N -------------------------------
    m_dot_ch_well1 = m_dot_w_ch / (N_wells * geom.num_tubes)
    r_ch = W.march(geom, pcm, T_inlet=T["T_4c"], T_m_seg=T_m_seg_ch,
                   m_dot=m_dot_ch_well1, P=cy.P, fluid=cy.fluid2,
                   k_wall=k_charge_wall, Rf_i=geom.Rf_i, times=times_ch,
                   n_segments=n_seg, mode="charge",
                   strict_bounds=case.strict)
    dz = geom.L_tube / n_seg
    eps_pcm = r_ch.state.volume(dz, geom.num_tubes) / geom.V_well
    E_stored_kJ = r_ch.Q_cum_J * geom.num_tubes * N_wells / 1000.0

    # ---- discharge: bisect the flow ratio, starting from the charge state ----
    rlo, rhi = num.ratio_bracket

    def discharged_kJ(ratio):
        r = W.march(geom, pcm, T_inlet=T["T_3d"], T_m_seg=T_m_seg_dc,
                    m_dot=ratio * m_dot_ch_well1, P=cy.P, fluid=cy.fluid2,
                    k_wall=geom.k_wall, Rf_i=geom.Rf_i, times=times_dc,
                    n_segments=n_seg, state=r_ch.state, mode="discharge")
        return abs(r.Q_cum_J) * geom.num_tubes * N_wells / 1000.0, r

    def f(ratio):
        Q, _ = discharged_kJ(ratio)
        return D_E_in_ORC / Q - 1.0 if Q > 0 else 1e6

    flo, fhi = f(rlo), f(rhi)
    if flo * fhi > 0:
        # The bracket does not contain the root. With the inventory now bounding
        # the discharge this is a physical statement rather than a numerical
        # accident: if even the largest admissible flow cannot deliver
        # D_E_in_ORC, the store does not hold it. v0.1 returned the bracket
        # bound and reported it as an answer.
        Q_hi, _ = discharged_kJ(rhi)
        raise ThumsConvergenceError(
            f"discharge cannot deliver {D_E_in_ORC:.0f} kJ within the flow "
            f"bracket; at ratio={rhi} it delivers {Q_hi:.0f} kJ "
            f"(stored {E_stored_kJ:.0f} kJ)", last=Q_hi)

    lo, hi = rlo, rhi
    for _ in range(num.max_iterations):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < num.tol_Q_ratio or (hi - lo) < 1e-6:
            lo = hi = mid
            break
        if flo * fm < 0:
            hi = mid
        else:
            lo, flo = mid, fm
    ratio = 0.5 * (lo + hi)
    Q_dc_kJ, r_dc = discharged_kJ(ratio)
    m_dot_d_well1 = ratio * m_dot_ch_well1

    # ---- what the store actually returned -----------------------------------
    eps_after = r_dc.state.volume(dz, geom.num_tubes) / geom.V_well
    eta_storage = Q_dc_kJ / E_stored_kJ if E_stored_kJ else float("nan")
    for r_, tag in ((r_ch, "charge"), (r_dc, "discharge")):
        if r_.warnings_energy:
            res.warnings.append(f"{tag}: {r_.warnings_energy}")
        res.warnings.extend(f"{tag}: {w}" for w in r_.warnings)

    detail = {
        "N_heat_transfer": sz.N_heat,
        "N_inventory": sz.N_inventory,
        "binding_constraint": sz.binding,
        "closure_charge": r_ch.closure,
        "closure_discharge": r_dc.closure,
        "E_stored_kJ": E_stored_kJ,
        "E_discharged_kJ": Q_dc_kJ,
        "eta_storage": eta_storage,
        "eps_pcm_end_of_discharge": eps_after,
        "Q_rejected_charge_J": r_ch.Q_rejected_J,
        "Q_rejected_discharge_J": r_dc.Q_rejected_J,
    }
    return N_wells, m_dot_d_well1, ratio, eps_pcm, detail


def run_cycle(case: Case) -> CycleResult:
    """Charge, size the well field, discharge, and return named KPIs."""
    t_start = time.time()
    res_warn = case.validate()
    cy, pcm, geom, num = case.cycle, case.pcm, case.geometry, case.numerics
    res = CycleResult()
    res.warnings.extend(res_warn)

    rank, hp, T = _cycles(case)
    rank_eff, hp_cop = rank["rank_eff"], hp["hp_cop"]

    # ---- energy bookkeeping (legacy cell 21) --------------------------------
    W_dot_T = cy.W_dot_el_out / cy.Turb_eff / cy.ElG_eff
    Q_dot_in_ORC = W_dot_T / rank_eff
    D_E_in_ORC = Q_dot_in_ORC * cy.t_dc * 3600.0                 # kJ
    D_E_out_HP = D_E_in_ORC * (1.0 + cy.loss_surplus)
    Q_dot_out_HP = D_E_out_HP / cy.t_ch / 3600.0
    W_dot_C_HP = Q_dot_out_HP / hp_cop
    W_dot_el_in = W_dot_C_HP / cy.Comp_eff / cy.ElH_eff

    cp_w = CP.PropsSI("C", "T", 0.5 * (T["T_3c"] + T["T_2c"]), "P", cy.P, cy.fluid2) / 1000.0
    m_dot_w_ch = Q_dot_out_HP / cp_w / cy.DT_3C_2C

    # ---- layer melting temperatures ----------------------------------------
    N_lay = case.N_lay
    dTm = cy.DT_3C_2C / N_lay
    i = np.arange(N_lay, dtype=float)
    T_m_lay = (pcm.T_m - i * dTm).tolist()
    T_m_lay_dc = (pcm.T_m - cy.DT_3C_2C + (i + 1.0) * dTm).tolist()

    times_ch = np.logspace(np.log10(1.0), np.log10(cy.t_ch * 3600.0), num.n_times)
    times_dc = np.logspace(np.log10(1.0), np.log10(cy.t_dc * 3600.0), num.n_times)
    gv = geom.legacy_vector()

    # The k_w slot of the charging chain. False = reproduce the legacy defect.
    k_charge_wall = geom.k_wall if case.charge_uses_wall_conductivity else pcm.k_l
    if not case.charge_uses_wall_conductivity:
        res.warnings.append("charge path uses PCM liquid conductivity as k_wall (legacy defect)")

    # ---- well field and discharge ------------------------------------------
    if case.front == "marched":
        N_wells, m_dot_d_well1, ratio, eps_pcm, march_detail = \
            _marched_well_field(
                case, T=T, T_m_lay=T_m_lay, T_m_lay_dc=T_m_lay_dc,
                m_dot_w_ch=m_dot_w_ch, D_E_out_HP=D_E_out_HP,
                D_E_in_ORC=D_E_in_ORC, times_ch=times_ch, times_dc=times_dc,
                k_charge_wall=k_charge_wall, res=res)
    elif case.front == "legacy":
        N_wells, m_dot_d_well1, ratio, eps_pcm, march_detail = \
            _legacy_well_field(
                case, T=T, T_m_lay=T_m_lay, T_m_lay_dc=T_m_lay_dc,
                m_dot_w_ch=m_dot_w_ch, D_E_out_HP=D_E_out_HP,
                times_ch=times_ch, times_dc=times_dc, gv=gv,
                k_charge_wall=k_charge_wall, res=res)
    else:
        raise ValueError(f"unknown front formulation: {case.front!r}")

    m_dot_c_well1 = m_dot_w_ch / N_wells / geom.num_tubes

    # ---- parasitics --------------------------------------------------------
    _, pp_dc_tube = L.calculate_pressure_drop(gv, m_dot_d_well1, cy.fluid2,
                                              T["T_2d"], T["T_3d"], cy.P)
    pumping_dc = pp_dc_tube * geom.num_tubes * N_wells / 1000.0

    m_dot_ch_well1 = m_dot_w_ch / (N_wells * geom.num_tubes)
    _, pp_ch_tube = L.calculate_pressure_drop(gv, m_dot_ch_well1, cy.fluid2,
                                              T["T_4c"], T["T_4a"], cy.P)
    pumping_ch = pp_ch_tube * geom.num_tubes * N_wells / 1000.0
    res.warnings.append("charging pressure drop evaluated between T_4c and T_4a "
                        "(source-loop pair); T_2c is the borehole return")

    # ---- KPIs --------------------------------------------------------------
    eta_rte = ((cy.W_dot_el_out - pumping_dc) * cy.t_dc
               / ((W_dot_el_in + pumping_ch) * cy.t_ch))
    eta_rte0 = cy.W_dot_el_out * cy.t_dc / (W_dot_el_in * cy.t_ch)
    D_E_in_ORC_kWh = D_E_in_ORC / 3600.0
    E_well = D_E_in_ORC_kWh / N_wells / 1000.0                    # MWh
    W_gross = cy.W_dot_el_out
    res.kpis = {
        "eta_rte": eta_rte,
        "eta_rte_nopump": eta_rte0,
        "cop_hp": hp_cop,
        "eta_orc": rank_eff,
        "eps_pcm": eps_pcm,
        "E_well": E_well,
        "rho_E": E_well * 1000.0 / geom.V_well,
        "N_wells": N_wells,
        "wells_per_MW": N_wells / (cy.W_dot_el_out / 1000.0),
        "f_pump": (pumping_ch + pumping_dc) / W_gross,
        "c_pcm": pcm.cost_per_kWh,
    }
    res.detail = {
        "Q_dot_out_HP_kW": Q_dot_out_HP, "Q_dot_in_ORC_kW": Q_dot_in_ORC,
        "W_dot_el_in_kW": W_dot_el_in, "m_dot_w_ch_kg_s": m_dot_w_ch,
        "m_dot_d_well1_kg_s": m_dot_d_well1, "flow_ratio_dc_ch": ratio,
        "pumping_ch_kW": pumping_ch, "pumping_dc_kW": pumping_dc,
        "D_E_in_ORC_kWh": D_E_in_ORC_kWh, "T_m_lay": T_m_lay,
        "stefan_charge": pcm.stefan(T["T_4c"] - pcm.T_m),
        "front": case.front,
    }
    res.detail.update(march_detail)
    res.converged = True
    res.wall_time_s = time.time() - t_start
    return res
