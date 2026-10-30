import numpy as np, boreas as B, boreas_ground as G
CSV='THUMS_wells_WellSTAR_clean.csv'
CASE = B.Case(T_m_C=150.0, h_m=380_000.0, k_s=0.60, k_l=0.45, N_lay=9,
              DT_3C_2C=55.0, L_ft=5000.0, n_ft=3, num_fins=24, fin_L=0.0035,
              fin_t=0.0015, t_ins=0.008, k_ins=0.007,
              D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
              k_cement=0.7, T_surface_C=18.0, grad_K_per_km=45.0, k_rock=2.5,
              rho_cp_rock=2.34e6, t_operation_yr=10.0, W_dot_el_out=1000.0,
              t_ch=10.0, t_dc=10.0, T_source_C=60.0, T_sink_C=20.0, P=10.0e5,
              DT_pinch_ORC=5.0, DT_2H_3C=10.0, n_segments=99)
c = CASE.with_(L_ft=1250.0/0.3048)
P = G.grissom_positions(CSV)
kop, _ = G.staggered_kop(P)
Rn = float(np.atleast_1d(B.formation_resistance(c, np.array([0.0])))[0]
           - B.ground_resistance(c))
for traj, nst in (('s_curve', 25), ('vertical', 50)):
    md, X = G.trajectories(P, 1250.0, kind=traj,
                           kop=(kop if traj == 's_curve' else None),
                           theta_max=55.0)
    gr = G.Ground(k=c.k_rock, rhocp=c.rho_cp_rock, t_yr=c.t_operation_yr,
                  rb=c.r_bore)
    pr = gr.profiles(md, X, Rn, nst=nst)
    assert getattr(gr, 'n_negative', 0) == 0, (traj, gr.n_negative)
    info = dict(N=len(P), L_md=1250.0, n_stations=nst, R_near=Rn,
                Rg_mean=float(np.mean(pr['temp'])),
                Rg_top=float(pr['temp'][0]), Rg_bottom=float(pr['temp'][-1]),
                tvd_bottom=float(pr['tvd'][-1]), island='Grissom',
                trajectory=traj, kop_rule='staggered', theta_max=55.0)
    key = G._ref_key(c, 'Grissom', traj, 'staggered', 55.0, 300.0, 1250.0)
    np.savez(key + '.npz', md=pr['md'], Rg=pr['temp'], tvd=pr['tvd'],
             info=np.array(info, dtype=object))
    print(traj, nst, 'Rg_mean', info['Rg_mean'], '->', key + '.npz', flush=True)
print('DONE')
