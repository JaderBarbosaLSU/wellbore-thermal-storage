"""Measure the energy-closure error of the v0.1 wellbore model.

Compares, over one charging period:

    Q_in    = integral of the heat given up by the fluid   (from the enthalpy drop)
    Q_latent= rho * h_m * V_melt                            (from the melt-front position)

These are two routes to the same quantity and must agree. They do not, because
the fluid-side conductance uses the FINNED perimeter while the Stefan front
solves a BARE concentric cylinder. This script quantifies the gap and attributes it.

Usage:  PYTHONPATH=. python tools/diagnose_closure.py
"""

import numpy as np

from thums_core import _legacy_physics as L
from thums_core.config import BASELINE
from thums_core.system import _cycles

PI = np.pi


def main(N_wells=12.0, k_wall_charge=None):
    case = BASELINE
    cy, pcm, geom, num = case.cycle, case.pcm, case.geometry, case.numerics
    rank, hp, T = _cycles(case)

    W_dot_T = cy.W_dot_el_out / cy.Turb_eff / cy.ElG_eff
    D_E_in_ORC = (W_dot_T / rank["rank_eff"]) * cy.t_dc * 3600.0
    D_E_out_HP = D_E_in_ORC * 1.05
    Q_out_HP = D_E_out_HP / cy.t_ch / 3600.0
    cp_w = L.CP.PropsSI("C", "T", 0.5 * (T["T_3c"] + T["T_2c"]), "P", cy.P, cy.fluid2) / 1000.0
    m_dot_w_ch = Q_out_HP / cp_w / cy.DT_3C_2C
    m_dot_1 = m_dot_w_ch / (N_wells * geom.num_tubes)

    k_w = pcm.k_l if k_wall_charge is None else k_wall_charge
    N_lay = case.N_lay
    dTm = cy.DT_3C_2C / N_lay
    T_m_lay = (pcm.T_m - np.arange(N_lay, dtype=float) * dTm).tolist()
    times = np.logspace(np.log10(1.0), np.log10(cy.t_ch * 3600.0), num.n_times)
    gv = geom.legacy_vector()

    profs, Qs, Vmelts, Q_cum_kJ = L.time_profiles_melt(
        times, gv, T["T_4c"], N_lay, T_m_lay, k_w, geom.Rf_i, m_dot_1,
        cy.P, cy.fluid2, pcm.k_l, pcm.cp_l, pcm.rho_l, pcm.h_m,
        num.n_segments, num.delta_max)

    t_end = times[-1]
    V_melt = Vmelts[t_end] if isinstance(Vmelts, dict) else float(np.asarray(Vmelts)[-1])
    Q_latent_kJ = pcm.rho_l * V_melt * pcm.h_m / 1000.0

    print(f"charging, per tube, N_wells = {N_wells:g}, k_wall(charge) = {k_w} W/m.K")
    print(f"  m_dot per tube          {m_dot_1:12.5f} kg/s")
    print(f"  Q_in  (integral of Q)   {Q_cum_kJ:12.4e} kJ")
    print(f"  Q_latent (rho h_m V)    {Q_latent_kJ:12.4e} kJ")
    print(f"  CLOSURE RATIO Q_in/Q_latent = {Q_cum_kJ / Q_latent_kJ:6.3f}"
          f"   (should be 1.00 + a few % of sensible)")

    # ---- attribution -------------------------------------------------------
    perim_bare = 2 * PI * geom.r_e
    perim_fin = 2 * PI * geom.r_e + 2 * geom.num_fins * geom.fin_L
    print("\nattribution")
    print(f"  fluid side sees perim_T = 2*pi*r_e + 2*n_fin*fin_L = {perim_fin:.4f} m")
    print(f"  Stefan front sees        2*pi*r_e                   = {perim_bare:.4f} m")
    print(f"  ratio                                              = {perim_fin/perim_bare:6.3f}")

    sens = pcm.cp_l * (T["T_4c"] - pcm.T_m) / pcm.h_m
    print(f"  neglected sensible heat of the melt (cp*dT/h_m)    = {sens*100:6.2f} %")

    df = profs[t_end] if isinstance(profs, dict) else profs
    d = df["delta [m]"].to_numpy()
    print(f"\n  delta at t_end: min {d.min()*1000:7.2f} mm  "
          f"max {d.max()*1000:7.2f} mm  (borehole radius {geom.D_well/2*1000:.1f} mm)")
    print(f"  segments with delta > borehole radius: "
          f"{int((d > geom.D_well/2).sum())} of {len(d)}")


if __name__ == "__main__":
    for N in (12.0, 24.0):
        for kw in (None, 45.0):
            main(N, kw)
            print("-" * 72)
