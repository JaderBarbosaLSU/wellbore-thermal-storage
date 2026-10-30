"""Does ONE fixed station SPACING keep V11 (q > 0) over the whole sweep?"""
import sys, json
import numpy as np, boreas as B, boreas_ground as G
DS = float(sys.argv[1]); step = sys.argv[2]
CSV='THUMS_wells_WellSTAR_clean.csv'
THUMS=dict(D_in=201.2/25.4,r_casing_o=0.21910/2,r_bore=0.31120/2,k_cement=0.7)
P=G.grissom_positions(CSV); K,_=G.staggered_kop(P); N=len(P)
GEOM={'R1':(G.disc_positions(N,B0=G.B0_MEASURED),'vertical',None),
      'R2':(G.rows_positions(N,B0=G.B0_MEASURED,w=G.W_DEF),'vertical',None),
      'R3':(P,'vertical',None),'R4':(P,'s_curve',K)}
Pp,traj,kk=GEOM[step]
for L in (100.,200.,300.,400.,500.,600.,800.,1000.,1250.,1600.,2000.):
    nst=max(3,int(round(L/DS)))
    c=B.Case(n_segments=99,L_ft=L/0.3048,t_operation_yr=10.0,**THUMS)
    Rn=float(np.atleast_1d(B.formation_resistance(c,np.array([0.0])))[0]-B.ground_resistance(c))
    md,X=G.trajectories(Pp,L,kind=traj,kop=kk,theta_max=55.0)
    gr=G.Ground(k=c.k_rock,rhocp=c.rho_cp_rock,t_yr=10.0,rb=c.r_bore)
    pr=gr.profiles(md,X,Rn,nst=nst)
    print(f"{step} ds={DS:4.0f} L={L:6.0f} nst={nst:4d} actual_ds={L/nst:6.1f} "
          f"neg={gr.n_negative:6d} offW={100*gr.off_station_weight:5.1f}% "
          f"Rg={pr['temp'].mean():8.3f}", flush=True)
