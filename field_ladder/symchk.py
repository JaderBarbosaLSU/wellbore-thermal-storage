"""Check 5.1 independently: is G asymmetric, and does the symmetrization
cause the negative strengths?"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import boreas as B, boreas_ground as G
CSV='THUMS_wells_WellSTAR_clean.csv'
THUMS=dict(D_in=201.2/25.4,r_casing_o=0.21910/2,r_bore=0.31120/2,k_cement=0.7)
P=G.grissom_positions(CSV); K,_=G.staggered_kop(P); N=len(P)
print(f"{'L':>6}{'nst':>5}{'ds':>7}{'|G-Gt|/|G|':>12}"
      f"{'sym neg':>9}{'sym Rg':>9}{'raw neg':>9}{'raw Rg':>9}"
      f"{'diag':>8}{'max offrow':>11}")
for L, nsts in ((1250.0,(25,50,100)), (400.0,(8,25,50))):
    c=B.Case(n_segments=99,L_ft=L/0.3048,t_operation_yr=10.0,**THUMS)
    Rn=float(np.atleast_1d(B.formation_resistance(c,np.array([0.0])))[0]
             -B.ground_resistance(c))
    md,X=G.trajectories(P,L,kind='s_curve',kop=K,theta_max=55.0)
    for nst in nsts:
        gr=G.Ground(k=c.k_rock,rhocp=c.rho_cp_rock,t_yr=10.0,rb=c.r_bore)
        st,r,cc,d,Nn=gr.pairs(md,X,nst)
        g=gr.kern(d); M=N*nst; self_t=gr.kern(np.array([gr.rb]))[0]
        Graw=sp.csr_matrix((g,(r,cc)),shape=(M,M))
        asym=spla.norm(Graw-Graw.T)/spla.norm(Graw)
        b=np.ones(M)
        Asym=(Graw+Graw.T)*0.5+sp.identity(M)*(Rn+self_t)
        Araw=Graw+sp.identity(M)*(Rn+self_t)
        res={}
        for tag,A in (('sym',Asym),('raw',Araw)):
            q=spla.spsolve(A.tocsc(),b).reshape(N,nst)
            res[tag]=(int((q<0).sum()), float(1.0/q.mean(0).mean()-Rn))
        offrow=float(np.asarray(Graw.sum(1)).ravel().max())
        print(f"{L:6.0f}{nst:5d}{L/nst:7.1f}{asym:12.2f}"
              f"{res['sym'][0]:9d}{res['sym'][1]:9.4f}"
              f"{res['raw'][0]:9d}{res['raw'][1]:9.4f}"
              f"{Rn+self_t:8.3f}{offrow:11.2f}", flush=True)
