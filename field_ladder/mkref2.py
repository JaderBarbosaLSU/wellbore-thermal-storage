import numpy as np, boreas as B, boreas_ground as G
CSV='THUMS_wells_WellSTAR_clean.csv'
CASE = B.Case(T_m_C=150.0, h_m=380_000.0, k_s=0.60, k_l=0.45, N_lay=9,
              DT_3C_2C=55.0, n_ft=3, num_fins=24, fin_L=0.0035,
              fin_t=0.0015, t_ins=0.008, k_ins=0.007,
              D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
              k_cement=0.7, T_surface_C=18.0, grad_K_per_km=45.0, k_rock=2.5,
              rho_cp_rock=2.34e6, t_operation_yr=10.0, W_dot_el_out=1000.0,
              t_ch=10.0, t_dc=10.0, T_source_C=60.0, T_sink_C=20.0, P=10.0e5,
              DT_pinch_ORC=5.0, DT_2H_3C=10.0, n_segments=99)
REF = CASE.with_(L_ft=1250.0/0.3048)
for traj in ('s_curve','vertical'):
    G.cached_reference_field(B, REF, CSV, fetch=False, trajectory=traj,
                             kop_rule='staggered', theta_max=55.0,
                             kop_ft=300.0, L_md=1250.0, ds=50.0)
print('DONE')
