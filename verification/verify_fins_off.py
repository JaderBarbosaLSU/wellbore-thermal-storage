"""Does removing the fins collapse the closure discrepancy?

FIX_energy_closure.md attributes the energy-balance failure to a surface
mismatch: the fluid side uses a finned perimeter (0.4924 m), the melt front a
bare cylinder (0.1324 m), a factor of 3.72.

If that diagnosis is complete, then setting num_fins = 0 makes the two surfaces
identical, and the legacy front and the marched front must then agree to within
the neglected sensible heat of the melt layer (~5%).

If it is incomplete -- if the quasi-steady closed-form solve is ALSO wrong for
reasons unrelated to the fins -- the discrepancy survives with the fins off.

This is the check that distinguishes them. The 1e-16 closure reported in the FIX
doc is true by construction (`advance` integrates the same q' that `closure_error`
differences) and therefore tests nothing.
"""
import numpy as np

from thums_core import _legacy_physics as L
from thums_core import well_marched as W
from thums_core.config import BASELINE
from thums_core.system import _cycles


def run(num_fins, k_wall_charge, N_wells, label):
    case = BASELINE.with_(num_fins=num_fins)
    geom, pcm, cy, num = case.geometry, case.pcm, case.cycle, case.numerics
    rank, hp, T = _cycles(case)

    # charging mass flow, reproducing run_cycle's bookkeeping
    W_dot_T = cy.W_dot_el_out / cy.Turb_eff / cy.ElG_eff
    Q_dot_in_ORC = W_dot_T / rank["rank_eff"]
    D_E_in_ORC = Q_dot_in_ORC * cy.t_dc * 3600.0
    D_E_out_HP = D_E_in_ORC * (1.0 + cy.loss_surplus)
    Q_dot_out_HP = D_E_out_HP / cy.t_ch / 3600.0
    cp_w = L.CP.PropsSI("C", "T", 0.5 * (T["T_3c"] + T["T_2c"]), "P", cy.P,
                        cy.fluid2) / 1000.0
    m_dot_w_ch = Q_dot_out_HP / cp_w / cy.DT_3C_2C
    m1 = m_dot_w_ch / (N_wells * geom.num_tubes)

    N_lay = case.N_lay
    dTm = cy.DT_3C_2C / N_lay
    i = np.arange(N_lay, dtype=float)
    T_m_lay = (pcm.T_m - i * dTm).tolist()

    n_seg = num.n_segments
    times = np.logspace(0.0, np.log10(cy.t_ch * 3600.0), num.n_times)

    # ---- legacy: independent closed-form Stefan front ----------------------
    gv = geom.legacy_vector()
    _, Qs, Vm, Q_cum_kJ = L.time_profiles_melt(
        times, gv, T["T_4c"], N_lay, T_m_lay, k_wall_charge, geom.Rf_i, m1,
        cy.P, cy.fluid2, pcm.k_l, pcm.cp_l, pcm.rho_l, pcm.h_m,
        n_seg, num.delta_max)
    V_leg = Vm[times[-1]]
    Q_leg_J = Q_cum_kJ * 1000.0
    lat_leg = pcm.rho_l * pcm.h_m * V_leg
    closure_leg = Q_leg_J / lat_leg

    # ---- marched: front from the delivered heat ----------------------------
    # per-segment melting temperature, same layer map the legacy march uses
    seg_per_lay = n_seg / N_lay
    T_m_seg = np.array([T_m_lay[min(int(s / seg_per_lay), N_lay - 1)]
                        for s in range(n_seg)])
    r = W.march(geom, pcm, T_inlet=T["T_4c"], T_m_seg=T_m_seg, m_dot=m1,
                P=cy.P, fluid=cy.fluid2, k_wall=k_wall_charge, Rf_i=geom.Rf_i,
                times=times, n_segments=n_seg, mode="charge")
    dz = geom.L_tube / n_seg
    V_mar = float(np.sum(r.state.A_melt * dz))
    Q_mar_J = r.Q_cum_J

    print(f"{label:34s} fins={num_fins:2d} k_w={k_wall_charge:5.2f}")
    print(f"   legacy   Q_in={Q_leg_J/1e9:8.3f} GJ  V_melt={V_leg:7.4f} m3  "
          f"Q_in/Q_latent={closure_leg:6.3f}")
    print(f"   marched  Q_in={Q_mar_J/1e9:8.3f} GJ  V_melt={V_mar:7.4f} m3  "
          f"Q_in/Q_latent={Q_mar_J/(pcm.rho_l*pcm.h_m*V_mar):6.3f}")
    print(f"   V_marched / V_legacy = {V_mar/V_leg:6.3f}    "
          f"Q_marched / Q_legacy = {Q_mar_J/Q_leg_J:6.3f}")
    print()
    return closure_leg, V_mar / V_leg


if __name__ == "__main__":
    print("=" * 74)
    print("Closure of the LEGACY front, fins on vs fins off")
    print("If the fin/bare mismatch is the whole story, Q_in/Q_latent -> ~1.05")
    print("when num_fins = 0.")
    print("=" * 74)
    for fins in (24, 0):
        for kw in (0.45, 45.0):
            run(fins, kw, 24.123429921, f"fins={fins} k_w={kw}")
