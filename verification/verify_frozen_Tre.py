"""Is the closed-form Stefan solve algebraically wrong, or is its ASSUMPTION wrong?

Two separate questions:

(a) Given a genuinely constant T_re, does `compute_delta2_fast` agree with a
    direct integration of the same quasi-steady bare-cylinder ODE?
        -> if yes, the algebra is right.

(b) In the actual well, is T_re constant from t = 0?
        -> if no, applying the closed form with the CURRENT T_re over-predicts
           melting, because T_re rises monotonically from ~T_m as the melt layer
           thickens and it is being applied retroactively to the whole history.

Together these say whether the fin/bare-surface mismatch identified in
FIX_energy_closure.md is the whole story or only part of it.
"""
import numpy as np

from thums_core import _legacy_physics as L
from thums_core.config import BASELINE

geom, pcm = BASELINE.geometry, BASELINE.pcm
r_e = geom.r_e
k_m, cp_m, rho_m, h_m = pcm.k_l, pcm.cp_l, pcm.rho_l, pcm.h_m
T_m = pcm.T_m

print("(a) closed form vs direct ODE integration, T_re held constant")
print("    d/dt[rho*h_m*pi*(s^2-r_e^2)] = 2*pi*k_m*(T_re-T_m)/ln(s/r_e)")
print()
for dT in (2.0, 10.0, 30.0):
    T_re = T_m + dT
    t_end = 10.0 * 3600.0
    d_closed = L.compute_delta2_fast(r_e, k_m, cp_m, rho_m, h_m, T_m, T_re,
                                     t_end, delta_max=0.5)
    # The ODE integrates exactly to  r_e^2*[x^2/2*ln x - x^2/4 + 1/4] = k*dT*t/(rho*h).
    # Check the returned root satisfies it. (Direct Euler stepping fails here:
    # q' ~ 1/ln(s/r_e) is singular at delta = 0, so the first step is unbounded.
    # The singularity is integrable, which is exactly why the closed form exists.)
    x = 1.0 + d_closed / r_e
    lhs = 0.5 * x ** 2 * np.log(x) - 0.25 * x ** 2 + 0.25
    rhs = k_m * dT * t_end / (rho_m * h_m * r_e ** 2)
    print(f"    dT={dT:5.1f} K   delta={d_closed*1000:8.3f} mm   "
          f"term(x)={lhs:.6e}   k*dT*t/(rho*h*r_e^2)={rhs:.6e}   "
          f"residual={abs(lhs-rhs)/rhs:.2e}")

print()
print("(b) is T_re actually constant? profile of T_re over the charging window,")
print("    fins OFF (so the fluid-side and front-side surfaces are identical)")
print()

case = BASELINE.with_(num_fins=0)
g = case.geometry
gv = g.legacy_vector()
cy = case.cycle
N_lay = case.N_lay
dTm = cy.DT_3C_2C / N_lay
T_m_lay = (pcm.T_m - np.arange(N_lay, dtype=float) * dTm).tolist()
T_4c = pcm.T_m_C + cy.DT_4C_M + 273.15
m1 = 23.40677583278621 / (24.123429921 * g.num_tubes)

times = np.logspace(0.0, np.log10(cy.t_ch * 3600.0), 8)
print(f"    {'t [h]':>8s} {'T_re-T_m at mid-well [K]':>26s} {'delta [mm]':>12s}")
for t in times:
    df, Q, V = L.temperature_profile_melt(
        gv, T_4c, N_lay, T_m_lay, g.k_wall, g.Rf_i, m1, cy.P, cy.fluid2, t,
        pcm.k_l, pcm.cp_l, pcm.rho_l, pcm.h_m, 100, 0.5)
    col_Tre = [c for c in df.columns if "T_re" in c or "re" in c.lower()]
    mid = len(df) // 2
    tre = df[col_Tre[0]].iloc[mid] if col_Tre else float("nan")
    dcol = [c for c in df.columns if "delta" in c.lower()]
    dv = df[dcol[0]].iloc[mid] if dcol else float("nan")
    print(f"    {t/3600:8.3f} {tre - T_m_lay[min(mid//(100//N_lay), N_lay-1)]:26.3f} "
          f"{dv*1000:12.3f}")
