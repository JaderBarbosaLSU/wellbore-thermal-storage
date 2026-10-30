"""Section 2b: the model is exactly linear in plant size.

Run the reference case at 1 MW and at 10 MW and compare everything the
paper reports per well, per metre or per MWe.
"""
import numpy as np, boreas as B, boreas_ground as G
CSV='THUMS_wells_WellSTAR_clean.csv'
THUMS=dict(D_in=201.2/25.4,r_casing_o=0.21910/2,r_bore=0.31120/2,k_cement=0.7)
CASE=B.Case(n_segments=99,L_ft=1250/0.3048,t_operation_yr=10.0,**THUMS)
md,Rg,tvd,_=G.cached_reference_field(B,CASE,CSV,trajectory='s_curve',
    kop_rule='staggered',theta_max=55.0,kop_ft=300.0,L_md=1250.0,ds=50.0,
    verbose=False)
out={}
for W in (1000.0, 10000.0):
    c=CASE.with_(W_dot_el_out=W)
    with G.field_ground(B,md,Rg,tvd):
        r=B.simulate_css_corrected(c); k=B.performance_indices(c,r)
    out[W]=dict(m1_ch=r['m1_ch'], m1_dc=r['m1_dc'],
                e_prime=k['kWh_per_m_well'], RTE=k['eta_RTE'],
                loss=k['f_formation_loss'], eta_storage=r['eta_storage'],
                T2d=r['T_2d_realized'], N_wells=r['N_wells'],
                wells_per_MWe=k['wells_per_MWe'],
                N_geo_per_MWe=k['N_geo_per_MWe'],
                UA_per_MWe=k['UA_per_MWe'],
                m_w_dc_per_MWe=k['m_w_dc_per_MWe'],
                W_el_in_per_MWe=k['W_el_in_per_MWe'],
                Q_orce_over_Wnet=k['Q_orce_over_Wnet'])
a,b=out[1000.0],out[10000.0]
print(f"{'quantity':22s}{'1 MW':>16s}{'10 MW':>16s}{'rel diff':>12s}")
worst=0.0
for q in a:
    rel=abs(b[q]/a[q]-1) if a[q] else 0.0
    if q!='N_wells': worst=max(worst,rel)
    print(f'{q:22s}{a[q]:16.6g}{b[q]:16.6g}{rel:12.2e}')
print(f"\nN_wells scales 10x as it must: {b['N_wells']/a['N_wells']:.6f}")
print(f"worst relative difference, everything else: {worst:.2e}")
assert worst < 1e-8, worst
print('PASS -- the 1 MW runs in this pack are valid at any plant size.')
