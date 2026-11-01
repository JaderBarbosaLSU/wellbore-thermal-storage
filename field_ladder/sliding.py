"""Handoff v2.6 item 2: the ORC evaluated along the discharge.

The v2.5 answer used the wrong threshold. T_3e + DT_pinch = 95.0 C is
the outlet below which the water cannot heat the ORC AT ALL; it is not
the condition for a FIXED-PRESSURE ORC to meet its pinch. The pinch sits
at the onset of boiling, 15.3 % of the evaporator duty, where the water
is at

    T_3d + 0.153 (T_2d - T_3d)

so at the design pressure the pinch holds only while the instantaneous
outlet is at or above the AVERAGED outlet it was set from. That is
roughly half the discharge, not 0.01 % of it.

What this computes, for the adiabatic limit, R4 and R3:

  * eta_ORC re-solved at the instantaneous outlet of every recorded time
    level (T_3e by the model's own bisection on the pinch), which is the
    sliding-pressure plant;
  * the energy-weighted mean of it, and the eta_RTE that follows;
  * the share of dE_dc delivered while the outlet is below the mean.

Infeasible levels -- where no evaporating temperature clears the pinch
at all -- are counted as eta_ORC = 0, which is what a plant would get.
"""
import json
import os
import pickle

import numpy as np

import boreas as B
import boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
L_MD = 1250.0
THUMS = dict(D_in=201.2 / 25.4, r_casing_o=0.21910 / 2,
             r_bore=0.31120 / 2, k_cement=0.7)
CASE = B.Case(n_segments=99, L_ft=L_MD / 0.3048, t_operation_yr=10.0,
              **THUMS)
STORE = 'sliding.json'
trapz = getattr(np, 'trapezoid', None) or np.trapz


def march(key):
    if key == 'R4' and os.path.exists('rec_field.pkl'):
        return pickle.load(open('rec_field.pkl', 'rb'))['rec']
    if key == 'adiabatic':
        return B.simulate_css_corrected(CASE.with_(k_rock=1e-9), record=True)
    traj = 'vertical' if key == 'R3' else 's_curve'
    md, Rg, tvd, _ = G.cached_reference_field(
        B, CASE, CSV, trajectory=traj, kop_rule='staggered', theta_max=55.0,
        kop_ft=300.0, L_md=L_MD, ds=G.DS_STATION, verbose=False)
    with G.field_ground(B, md, Rg, tvd):
        return B.simulate_css_corrected(CASE, record=True)


def one(key):
    rec = march(key)
    kk = B.performance_indices(CASE, rec)
    T_bar = rec['T_2d_realized']
    T_3d = rec['T']['T_3d']
    eta_fixed = rec['rank']['rank_eff']

    H = rec['discharge']['history']
    z = rec['discharge']['z']
    th = np.asarray(rec['discharge']['t_hist'], float)
    q = np.abs(np.asarray(H['q_prime'], float).sum(axis=1) * (z[1] - z[0]))
    T = np.asarray(rec['discharge']['T_out'], float)
    if len(T) == len(th) - 1:
        T = np.r_[T[0], T]

    # eta_ORC at each instantaneous outlet, from the curve built once by
    # etacurve.py. eta_ORC depends on T_2d alone -- everything else in
    # the cycle is fixed -- so one curve serves all three cases, and the
    # 9 s pinch bisection is paid 13 times instead of 120.
    cur = json.load(open('eta_curve.json'))
    xs = np.array(sorted(float(k) for k in cur))
    ys = np.array([cur[f'{x:.3f}']['eta_ORC'] for x in xs])
    eta_t = np.interp(T - 273.15, xs, ys)
    eta_t = np.where(T <= T_3d + 1e-6, 0.0, eta_t)
    feasible = int((eta_t > 0).sum())

    E_tot = trapz(q, th)
    eta_slide = trapz(eta_t * q, th) / E_tot          # energy-weighted
    below = T < T_bar
    share_below = trapz(np.where(below, q, 0.0), th) / E_tot

    # eta_RTE follows: the only thing that changes is the electrical
    # output, so scale the fixed-pressure numerator by the ratio of the
    # two conversion efficiencies and subtract the same pumping.
    W_out_fixed = kk['W_el_out_kW']
    W_out_slide = W_out_fixed * eta_slide / eta_fixed
    pump_dc, pump_ch = kk['pumping_dc_kW'], kk['pumping_ch_kW']
    W_in = kk['W_el_in_kW']
    rte_fixed = kk['eta_RTE']
    rte_slide = ((W_out_slide - pump_dc) * CASE.t_dc
                 / ((W_in + pump_ch) * CASE.t_ch))

    out = dict(
        T2d_mean_C=float(T_bar - 273.15),
        T2d_start_C=float(T[0] - 273.15), T2d_end_C=float(T[-1] - 273.15),
        eta_ORC_fixed=float(eta_fixed), eta_ORC_sliding=float(eta_slide),
        eta_ORC_ratio=float(eta_slide / eta_fixed),
        eta_RTE_fixed=float(rte_fixed), eta_RTE_sliding=float(rte_slide),
        energy_below_mean_pct=100.0 * float(share_below),
        time_below_mean_pct=100.0 * float(
            trapz(below.astype(float), th) / (th[-1] - th[0])),
        levels=int(len(T)), levels_feasible=int(feasible))
    print(f"{key:10s} eta_ORC {eta_fixed:.4f} fixed -> {eta_slide:.4f} "
          f"sliding ({100*(eta_slide/eta_fixed-1):+.1f} %)   "
          f"eta_RTE {rte_fixed:.4f} -> {rte_slide:.4f}   "
          f"energy below the mean outlet {100*share_below:.1f} %",
          flush=True)
    return out


def main():
    res = json.load(open(STORE)) if os.path.exists(STORE) else {}
    for key in ('R4', 'R3', 'adiabatic'):
        if key in res:
            continue
        res[key] = one(key)
        json.dump(res, open(STORE, 'w'), indent=1)
    if len(res) == 3:
        print('\nThe fixed-pressure plant is evaluated at the AVERAGED '
              'outlet, so it\nover-rates the cold half of the discharge '
              'and under-rates the hot half.')


if __name__ == '__main__':
    main()
