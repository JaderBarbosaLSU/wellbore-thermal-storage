"""eta_ORC as a function of the store outlet, built once.

Everything else in the cycle is fixed, so eta_ORC depends on T_2d alone
and one curve serves the adiabatic limit, R3 and R4. Each point costs a
pinch bisection, about 9 s, so 40 time levels x 3 cases would be 18
minutes of re-solving the same function; 16 points and an interpolation
is the same answer in two.
"""
import json, os, sys, time
import numpy as np
import boreas as B

OUT = 'eta_curve.json'
CASE = B.Case(n_segments=99, L_ft=1250 / 0.3048, t_operation_yr=10.0,
              D_in=201.2 / 25.4, r_casing_o=0.21910 / 2,
              r_bore=0.31120 / 2, k_cement=0.7)
GRID = [90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 116.544, 120.0, 125.0,
        130.0, 135.0, 139.544, 141.6, 145.0]

d = json.load(open(OUT)) if os.path.exists(OUT) else {}
budget = float(sys.argv[1]) if len(sys.argv) > 1 else 150.0
t0 = time.time()
for T in GRID:
    k = f'{T:.3f}'
    if k in d:
        continue
    if time.time() - t0 > budget:
        print('BUDGET', flush=True)
        break
    try:
        r, _, _ = B.cycle_state_points(CASE, T_2d=T + 273.15)
        d[k] = dict(eta_ORC=float(r['rank_eff']),
                    T_3e_C=float(r['T_3e'] - 273.15))
    except Exception as e:
        d[k] = dict(eta_ORC=0.0, T_3e_C=None,
                    error=f'{type(e).__name__}: {e}'[:120])
    json.dump(d, open(OUT, 'w'), indent=1)
    print(f"  T_2d {T:7.3f} C -> eta_ORC {d[k]['eta_ORC']:.4f}", flush=True)
print(f'{len(d)} of {len(GRID)} points', flush=True)
