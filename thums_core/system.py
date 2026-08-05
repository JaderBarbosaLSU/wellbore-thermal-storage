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

    # ---- size the well field (charging) ------------------------------------
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

    # ---- discharge flow ----------------------------------------------------
    m_dot_c_well1 = m_dot_w_ch / N_wells / geom.num_tubes
    rlo, rhi = num.ratio_bracket
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

    # ---- storage utilisation ----------------------------------------------
    # Melted volume at the end of charging. The legacy helper
    # `mass_effectiveness_end_of_charge_cached` lives in a driver cell, not in the
    # physics; its arithmetic is reproduced here so the package does not depend on
    # a notebook cell. eps = V_melt(per tube) * num_tubes / V_well, i.e. the same
    # ratio the legacy code forms from masses (rho_l cancels).
    t_end = [cy.t_ch * 3600.0]
    _, _, V_melt_map, _ = L.time_profiles_melt(
        t_end, gv, T["T_4c"], N_lay, T_m_lay, k_charge_wall, geom.Rf_i,
        m_dot_ch_well1, cy.P, cy.fluid2, pcm.k_l, pcm.cp_l, pcm.rho_l, pcm.h_m,
        num.n_segments, num.delta_max)
    V_melt_tube = V_melt_map.get(t_end[0], np.nan) if isinstance(V_melt_map, dict) \
        else float(np.asarray(V_melt_map).ravel()[-1])
    eps_pcm = V_melt_tube * geom.num_tubes / geom.V_well
    if np.isfinite(eps_pcm) and eps_pcm > 1.0:
        res.warnings.append(f"eps_pcm = {eps_pcm:.3f} > 1: melt fronts are not "
                            "bounded by the borehole (audit 2.5)")

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
    }
    res.converged = True
    res.wall_time_s = time.time() - t_start
    return res
