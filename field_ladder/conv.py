"""Station-spacing convergence, WITHOUT the symmetrization (handoff 5.1).

The criterion is now convergence of the result, not the absence of
negative strengths: the negatives were an artifact of symmetrizing a
matrix that is not symmetric, and with that removed there are none at any
spacing. Choose the largest ds whose loss is within 0.1 points, and e'
within 0.2 %, of the next finer one.
"""
import json, time
import numpy as np
import boreas as B, boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
THUMS = dict(D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
             k_cement=0.7)
P = G.grissom_positions(CSV)
KOP, _ = G.staggered_kop(P)
N = len(P)
DISC = G.disc_positions(N, B0=G.B0_MEASURED)
ROWS = G.rows_positions(N, B0=G.B0_MEASURED, w=G.W_DEF)
GEOM = {'R1': (DISC, 'vertical', None), 'R2': (ROWS, 'vertical', None),
        'R3': (P, 'vertical', None), 'R4': (P, 's_curve', KOP)}
DS_GRID = (100.0, 50.0, 25.0, 12.5, 6.25)

import os
out = json.load(open('station_convergence.json')) if os.path.exists('station_convergence.json') else []
DONE = {(x['step'], x['L_md'], x['nst']) for x in out}
for L in (1250.0, 400.0):
    c = B.Case(n_segments=99, L_ft=L/0.3048, t_operation_yr=10.0, **THUMS)
    Rn = float(np.atleast_1d(B.formation_resistance(c, np.array([0.0])))[0]
               - B.ground_resistance(c))
    for step in ('R4', 'R3', 'R1'):
        Pp, traj, kk = GEOM[step]
        md, X = G.trajectories(Pp, L, kind=traj, kop=kk, theta_max=55.0)
        for ds in DS_GRID:
            nst = int(max(3, round(L / ds)))
            if (step, L, nst) in DONE:
                continue
            gr = G.Ground(k=c.k_rock, rhocp=c.rho_cp_rock, t_yr=10.0,
                          rb=c.r_bore)
            t0 = time.time()
            try:
                pr = gr.profiles(md, X, Rn, nst=nst)
            except Exception as e:
                print(f'{step} L={L:.0f} ds={ds} FAILED {e}', flush=True)
                continue
            r = G.run_boreas(B, L, pr['md'], pr['temp'], tvd_prof=pr['tvd'],
                             n_segments=99, t_operation_yr=10.0, **THUMS)
            out.append(dict(step=step, L_md=L, ds_m=L/nst, ds_target=ds,
                            nst=nst, ds_over_B0=(L/nst)/G.B0_MEASURED,
                            n_negative=int(gr.n_negative),
                            off_station_weight=gr.off_station_weight,
                            Rg_mean=float(pr['temp'].mean()),
                            loss=r['loss'], e_prime=r['e_prime'],
                            RTE=r['RTE'], secs=time.time()-t0))
            print(f"{step} L={L:6.0f} ds={L/nst:6.2f} nst={nst:4d} "
                  f"offW={100*gr.off_station_weight:5.1f}% "
                  f"neg={gr.n_negative:5d} Rg={pr['temp'].mean():8.4f} "
                  f"loss={100*r['loss']:7.3f}% e'={r['e_prime']:.4f} "
                  f"({time.time()-t0:.0f} s)", flush=True)
            json.dump(out, open('station_convergence.json', 'w'), indent=1)
print('ALLDONE')
