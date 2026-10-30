import json, sys, time
import numpy as np, boreas as B, boreas_ground as G
step, L, nst = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
CSV='THUMS_wells_WellSTAR_clean.csv'
THUMS=dict(D_in=201.2/25.4,r_casing_o=0.21910/2,r_bore=0.31120/2,k_cement=0.7)
P=G.grissom_positions(CSV); K,_=G.staggered_kop(P); N=len(P)
GEOM={'R1':(G.disc_positions(N,B0=G.B0_MEASURED),'vertical',None),
      'R2':(G.rows_positions(N,B0=G.B0_MEASURED,w=G.W_DEF),'vertical',None),
      'R3':(P,'vertical',None),'R4':(P,'s_curve',K)}
Pp,traj,kk=GEOM[step]
c=B.Case(n_segments=99,L_ft=L/0.3048,t_operation_yr=10.0,**THUMS)
Rn=float(np.atleast_1d(B.formation_resistance(c,np.array([0.0])))[0]-B.ground_resistance(c))
md,X=G.trajectories(Pp,L,kind=traj,kop=kk,theta_max=55.0)
gr=G.Ground(k=c.k_rock,rhocp=c.rho_cp_rock,t_yr=10.0,rb=c.r_bore)
t0=time.time(); pr=gr.profiles(md,X,Rn,nst=nst)
r=G.run_boreas(B,L,pr['md'],pr['temp'],tvd_prof=pr['tvd'],n_segments=99,
               t_operation_yr=10.0,**THUMS)
rec=dict(step=step,L_md=L,nst=nst,ds_m=L/nst,ds_over_B0=(L/nst)/G.B0_MEASURED,
         n_negative=int(gr.n_negative),off_station_weight=gr.off_station_weight,
         Rg_mean=float(pr['temp'].mean()),loss=r['loss'],e_prime=r['e_prime'],
         RTE=r['RTE'],secs=time.time()-t0)
try: all_=json.load(open('station_convergence.json'))
except Exception: all_=[]
all_=[x for x in all_ if not (x['step']==step and x['L_md']==L and x['nst']==nst)]
all_.append(rec); json.dump(all_,open('station_convergence.json','w'),indent=1)
print(f"{step} L={L:6.0f} nst={nst:4d} ds={L/nst:6.2f} ds/B0={rec['ds_over_B0']:6.1f} "
      f"offW={100*rec['off_station_weight']:5.2f}% neg={rec['n_negative']:5d} "
      f"Rg={rec['Rg_mean']:8.3f} loss={100*r['loss']:6.2f}% e'={r['e_prime']:.4f} "
      f"({rec['secs']:.0f} s)", flush=True)
