"""Station-count convergence (2c.3), with the structural result that makes
most of the grid irrelevant.

A target on well i at station s couples to well k at station s_nb, where
s_nb is the station of the nearest point on k's trajectory. For PARALLEL
VERTICAL wells that nearest point is at the same depth, so s_nb == s for
every entry: the matrix is exactly block diagonal in station, every block
is the same, and R_g,eff cannot depend on nst at all. Measured below as
the off-station fraction.

For the S-curve the wells converge and diverge, so s_nb differs by at most
one station: the matrix is block TRIDIAGONAL, 5 % of the kernel weight sits
off the diagonal block, and the count does matter. That is the case the
grid below actually tests.
"""
import json, time
import numpy as np
import boreas as B, boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
THUMS = dict(D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
             k_cement=0.7)
P_SLOT = G.grissom_positions(CSV)
KOP, _ = G.staggered_kop(P_SLOT)
N = len(P_SLOT)
DISC = G.disc_positions(N, B0=G.B0_MEASURED)

PLAN = [('R4', P_SLOT, 's_curve', KOP, (25, 50, 100, 200)),
        ('R3', P_SLOT, 'vertical', None, (25, 50)),
        ('R1', DISC,   'vertical', None, (25, 50))]

out = []
for L in (1250.0, 400.0):
    c = B.Case(n_segments=99, L_ft=L/0.3048, t_operation_yr=10.0, **THUMS)
    Rn = float(np.atleast_1d(B.formation_resistance(c, np.array([0.0])))[0]
               - B.ground_resistance(c))
    for step, Pp, traj, kk, grid in PLAN:
        md, X = G.trajectories(Pp, L, kind=traj, kop=kk, theta_max=55.0)
        for nst in grid:
            gr = G.Ground(k=c.k_rock, rhocp=c.rho_cp_rock, t_yr=10.0,
                          rb=c.r_bore)
            t0 = time.time()
            pr = gr.profiles(md, X, Rn, nst=nst)
            nneg, off_w = gr.n_negative, gr.off_station_weight
            r = G.run_boreas(B, L, pr['md'], pr['temp'], tvd_prof=pr['tvd'],
                             n_segments=99, t_operation_yr=10.0, **THUMS)
            out.append(dict(step=step, L_md=L, nst=nst, ds_m=L/nst,
                            n_negative=int(nneg),
                            off_station_weight=off_w,
                            Rg_mean=float(np.mean(pr['temp'])),
                            loss=r['loss'], e_prime=r['e_prime'],
                            RTE=r['RTE'], secs=time.time()-t0))
            print(f"{step} L={L:6.0f} nst={nst:4d} ds={L/nst:6.2f} m "
                  f"offW={100*off_w:5.2f}% neg={nneg:6d} "
                  f"Rg={out[-1]['Rg_mean']:8.3f} loss={100*r['loss']:6.2f}% "
                  f"e'={r['e_prime']:.4f} ({out[-1]['secs']:.0f} s)",
                  flush=True)
            json.dump(out, open('station_convergence.json', 'w'), indent=1)
print('ALLDONE')
