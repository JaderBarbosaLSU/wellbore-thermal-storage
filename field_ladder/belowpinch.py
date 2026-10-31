"""Handoff 2(a): how much of the discharge arrives below the ORC's reach.

The plant model evaluates the ORC at the flow-averaged T_2d with a fixed
evaporating temperature. If the outlet falls below T_3e + DT_pinch late
in the discharge, an evaporator at that pressure cannot take that heat,
and eta_RTE is overstated by whatever share of dE_dc arrives there.

Reported two ways on purpose. The ENERGY share is what bounds the
overstatement; the TIME share is what a reader sees on the profile
figure, and the two are very different here.
"""
import json, os, pickle
import numpy as np
import boreas as B, boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
L_MD = 1250.0
THUMS = dict(D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
             k_cement=0.7)
CASE = B.Case(n_segments=99, L_ft=L_MD/0.3048, t_operation_yr=10.0, **THUMS)
trapz = getattr(np, 'trapezoid', None) or np.trapz


def march(key):
    if key == 'adiabatic':
        return B.simulate_css_corrected(CASE.with_(k_rock=1e-9), record=True)
    traj = 'vertical' if key == 'R3' else 's_curve'
    md, Rg, tvd, _ = G.cached_reference_field(
        B, CASE, CSV, trajectory=traj, kop_rule='staggered', theta_max=55.0,
        kop_ft=300.0, L_md=L_MD, ds=G.DS_STATION, verbose=False)
    with G.field_ground(B, md, Rg, tvd):
        return B.simulate_css_corrected(CASE, record=True)


out = {}
for key in ('adiabatic', 'R4', 'R3'):
    rec = march(key)
    T_lim = rec['rank']['T_3e'] - 273.15 + CASE.DT_pinch_ORC
    H = rec['discharge']['history']
    z = rec['discharge']['z']
    th = np.asarray(rec['discharge']['t_hist'], float)
    q = np.abs(np.asarray(H['q_prime'], float).sum(axis=1) * (z[1] - z[0]))
    T = np.asarray(rec['discharge']['T_out'], float) - 273.15
    if len(T) == len(th) - 1:
        T = np.r_[T[0], T]
    below = T < T_lim
    e_share = trapz(np.where(below, q, 0.0), th) / trapz(q, th)
    t_share = trapz(below.astype(float), th) / (th[-1] - th[0])
    out[key] = dict(T_lim_C=float(T_lim), T_out_end_C=float(T[-1]),
                    T2d_mean_C=float(rec['T_2d_realized'] - 273.15),
                    energy_share_pct=100 * float(e_share),
                    time_share_pct=100 * float(t_share),
                    rate_at_crossing_kW=float(q[below][0] / 1e3)
                    if below.any() else 0.0)
    print(f"{key:10s} T_lim {T_lim:6.1f} C   outlet ends at {T[-1]:6.1f} C"
          f"   energy below {100*e_share:5.2f} %   time below "
          f"{100*t_share:5.1f} %", flush=True)
json.dump(out, open('belowpinch.json', 'w'), indent=1)
print('\nThe energy share is what bounds the overstatement of eta_RTE.')
