"""V1-V13. Every one an assertion, and every one printed with its margin."""
import json, math, os
import numpy as np
import boreas as B, boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
THUMS = dict(D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
             k_cement=0.7)
LEGACY = dict(D_in=7.0, r_casing_o=0.098, r_bore=0.110, k_cement=1.0)
RES, FAIL = [], 0

def check(tag, what, ok, detail):
    global FAIL
    RES.append(dict(id=tag, check=what, pass_=bool(ok), detail=detail))
    FAIL += (not ok)
    print(f"  {tag:5s} {'PASS' if ok else 'FAIL'}  {what:<52s} {detail}",
          flush=True)

def case(L, casing=THUMS, **kw):
    return B.Case(n_segments=99, L_ft=L/0.3048, t_operation_yr=10.0,
                  **casing, **kw)

P = G.grissom_positions(CSV)
KOP, _ = G.staggered_kop(P)
N = len(P)

print('BOREAS verification suite, ds =', G.DS_STATION, 'm')

# ---- V1: energy closure and the loss-sign gate -------------------------
# The v1.3 "7.93 % / 0.2983" reference is RETIRED, not failed: it was the
# shielded case under the packing/perimeter field model, and that model no
# longer exists. Its replacement as a regression is V7, which reproduces
# the prototype's legacy-casing rungs. What V1 checks now is the thing no
# geometry change can excuse: the march closes, and the loss carries the
# sign it claims to.
c = case(1250.0)
r = B.simulate_css_corrected(c)
cl = abs(r['charge']['closure'])
cz = r['charge']['loss_clausius']
check('V1', 'charge closes and the Clausius gate is positive',
      cl < 1e-6 and cz > 0,
      f"closure {cl:.2e}, clausius {cz:+.3e}")

# ---- V2/V3/V4: the ground model, as standalone assertions --------------
c = case(1250.0)
gr = G.Ground(k=c.k_rock, rhocp=c.rho_cp_rock, t_yr=c.t_operation_yr,
              rb=c.r_bore)
check('V2', 'G_ii degenerates to the shipped isolated R_g',
      abs(gr.kern(c.r_bore) - B.ground_resistance(c)) < 1e-6,
      f"{gr.kern(c.r_bore):.9f} vs {B.ground_resistance(c):.9f}")

Rn = float(np.atleast_1d(B.formation_resistance(c, np.array([0.0])))[0]
           - B.ground_resistance(c))
md, X = G.trajectories(P[:1], 1250.0, kind='vertical')
pr1 = gr.profiles(md, X, Rn, nst=25)
# tolerance 1e-5, not 1e-6: V2 measures the gap between the E1 kernel and
# the closed-form log approximation at r_b as 5.7e-7 m K/W absolute, which
# is 1.7e-6 relative. V3 inherits exactly that, so a 1e-6 relative bound
# would be testing the approximation rather than the solver.
check('V3', 'one well in the field solver = isolated well',
      abs(pr1['temp'].mean() - B.ground_resistance(c)) / B.ground_resistance(c)
      < 1e-5,
      f"{pr1['temp'].mean():.6f} vs {B.ground_resistance(c):.6f} m K/W, "
      f"rel {abs(pr1['temp'].mean()/B.ground_resistance(c)-1):.1e}")

md, X = G.trajectories(P, 1250.0, kind='vertical')
prT = gr.profiles(md, X, Rn, nst=25)
# The direction: equal flux is the OPTIMISTIC model, and optimistic means
# a HIGHER effective resistance, hence a lower loss. It averages a wall
# temperature that is not uniform, under-weighting the edge wells, which
# are the ones doing the losing. The cache agrees: at R3, flux gives
# 12.97 % loss against 13.37 % for equal temperature.
check('V4', 'equal flux is the optimistic model (R_flux > R_temp)',
      prT['flux'].mean() > prT['temp'].mean(),
      f"flux {prT['flux'].mean():.4f} > temp {prT['temp'].mean():.4f} m K/W")

# ---- V5: adiabatic limit ----------------------------------------------
r0 = B.simulate_css_corrected(case(1250.0).with_(k_rock=1e-9))
check('V5', 'k_rock -> 0 gives zero loss and eta_storage = 1',
      r0['Q_loss_kJ'] == 0.0 and abs(r0['eta_storage'] - 1) < 1e-9,
      f"Q_loss {r0['Q_loss_kJ']:.1e} kJ, eta_storage {r0['eta_storage']:.9f}")

# ---- V6: the signed discharge ------------------------------------------
r6 = B.simulate_css_corrected(case(300.0, LEGACY))
k6 = B.performance_indices(case(300.0, LEGACY), r6)
check('V6', 'short isolated legacy well reports non-delivery',
      (not r6['delivers']) and k6['kWh_per_m_well'] == 0.0,
      f"delivers {r6['delivers']}, e' {k6['kWh_per_m_well']:.4f}, "
      f"T_2d {r6['T_2d_realized']-273.15:.1f} C")

# ---- V7/V8: the ladder against the cache --------------------------------
runs = [v for v in json.load(open('pack_runs.json')).values()
        if 'error' not in v]
def pick(step, L=1250.0, ground='superposition'):
    c_ = [v for v in runs if v['step'] == step and v['ground'] == ground
          and abs(v['L_md'] - L) < 1e-6]
    return c_[0] if c_ else None

check('V7', 'legacy-casing rungs regress to the prototype',
      abs(100*pick('R3_legacy')['loss'] - 17.4) < 0.15
      and abs(100*pick('R4_legacy')['loss'] - 70.7) < 0.3,
      f"R3 {100*pick('R3_legacy')['loss']:.2f}% (17.4), "
      f"R4 {100*pick('R4_legacy')['loss']:.2f}% (70.7)")

check('V8', 'R4 short-length agreement with the prototype wrapper',
      abs(100*pick('R4', 100.)['loss'] - 19.3) < 1.0
      and abs(100*pick('R4', 300.)['loss'] - 26.8) < 1.0,
      f"L=100 {100*pick('R4', 100.)['loss']:.1f}% (19.3), "
      f"L=300 {100*pick('R4', 300.)['loss']:.1f}% (26.8)")

# ---- V9: the ladder is monotonic in shielding ---------------------------
seq = [100*pick(s)['loss'] for s in ('R1', 'R2', 'R3')]
check('V9', 'disc < double row < real slots, all vertical',
      seq[0] < seq[1] < seq[2],
      ' < '.join(f'{v:.2f}%' for v in seq))

# ---- V10: R_g,eff at a given MD does not depend on L --------------------
# With a FIXED station spacing the stations sit at the same measured
# depths whatever L is, so this compares like with like. It is the check
# that caught the short-L defect in v2.0.
vals = {}
for L in (200.0, 2000.0):
    md, X = G.trajectories(P, L, kind='s_curve', kop=KOP, theta_max=55.0)
    c_ = case(L)
    Rn_ = float(np.atleast_1d(B.formation_resistance(c_, np.array([0.0])))[0]
                - B.ground_resistance(c_))
    g_ = G.Ground(k=c_.k_rock, rhocp=c_.rho_cp_rock, t_yr=10.0, rb=c_.r_bore)
    pr_ = g_.profiles(md, X, Rn_, nst=G.Ground.stations_for(L))
    vals[L] = float(np.interp(50.0, pr_['md'], pr_['temp']))
# Compared at L = 200 m rather than 100 m, and that is not cherry-picking:
# with ds = 50 m the stations sit at 25, 75, 125 ... m for EVERY length, so
# L = 200 and L = 2000 sample R_g,eff at exactly the same measured depths
# and the comparison is like for like. L = 100 m is the one length where
# the floor of three stations bites (ds = 33.3 m there), so its stations
# are offset from every other run's and the 1.3 % it shows is that offset,
# not a depth dependence. Reported below rather than hidden.
rel = abs(vals[200.0] / vals[2000.0] - 1)
check('V10', 'R_g,eff at MD = 50 m is the same for L = 200 and 2000 m',
      rel < 0.01,
      f"{vals[200.0]:.4f} vs {vals[2000.0]:.4f} m K/W, {100*rel:.2f} %")

md, X = G.trajectories(P, 100.0, kind='s_curve', kop=KOP, theta_max=55.0)
c_ = case(100.0)
Rn_ = float(np.atleast_1d(B.formation_resistance(c_, np.array([0.0])))[0]
            - B.ground_resistance(c_))
g_ = G.Ground(k=c_.k_rock, rhocp=c_.rho_cp_rock, t_yr=10.0, rb=c_.r_bore)
pr_ = g_.profiles(md, X, Rn_, nst=G.Ground.stations_for(100.0))
v100 = float(np.interp(50.0, pr_['md'], pr_['temp']))
check('V10b', 'the same at L = 100 m, where the 3-station floor bites',
      abs(v100 / vals[2000.0] - 1) < 0.02,
      f"{v100:.4f} vs {vals[2000.0]:.4f} m K/W, "
      f"{100*abs(v100/vals[2000.0]-1):.2f} % (floor: 3 stations, ds 33.3 m)")

# ---- V11: every solved strength positive, every cached run -------------
neg = [v for v in runs if v.get('n_negative', 0)]
check('V11', 'no negative well strengths in any cached run',
      not neg, f"{len(runs)} runs, {len(neg)} with negatives")

# ---- V12/V13: the normalization basis ----------------------------------
c = case(1250.0)
r = B.simulate_css_corrected(c)
k = B.performance_indices(c, r)
check('V12', "wells per MWe = t_dc / (e' L)",
      k['V12_rel_error'] < 1e-9,
      f"{k['wells_per_MWe']:.3f} vs {k['wells_per_MWe_from_eprime']:.3f}, "
      f"rel {k['V12_rel_error']:.2e}")

k10 = B.performance_indices(c.with_(W_dot_el_out=10000.0),
                            B.simulate_css_corrected(
                                c.with_(W_dot_el_out=10000.0)))
keys = ['wells_per_MWe', 'N_geo_per_MWe', 'UA_per_MWe', 'kWh_per_m_well',
        'eta_RTE', 'f_formation_loss', 'm_w_dc_per_MWe']
worst = max((abs(k10[q]/k[q] - 1) if k[q] else 0.0) for q in keys)
check('V13', 'per-MWe outputs identical at 1 MW and 10 MW',
      worst < 1e-9, f"worst relative difference {worst:.2e} over {len(keys)} keys")

# ---- V14: the system is solved as assembled, not symmetrized ----------
# This is a regression guard, not a physics check. Symmetrizing G was the
# v2.1 defect: it is not a symmetric operator -- row (i,s) is about the
# rise at a different point from row (j,s') -- and averaging produced
# negative strengths that were then misread as a grid problem.
import scipy.sparse as sp
import scipy.sparse.linalg as spla
c = case(1250.0)
gr = G.Ground(k=c.k_rock, rhocp=c.rho_cp_rock, t_yr=10.0, rb=c.r_bore)
Rn = float(np.atleast_1d(B.formation_resistance(c, np.array([0.0])))[0]
           - B.ground_resistance(c))
md, X = G.trajectories(P, 1250.0, kind='s_curve', kop=KOP, theta_max=55.0)
nst = G.Ground.stations_for(1250.0)
st, rr, cc, dd, Nn = gr.pairs(md, X, nst)
g = gr.kern(dd)
M = Nn * nst
self_t = gr.kern(np.array([gr.rb]))[0]
Graw = sp.csr_matrix((g, (rr, cc)), shape=(M, M))
asym = spla.norm(Graw - Graw.T) / spla.norm(Graw)
qs = spla.spsolve(((Graw + Graw.T) * 0.5
                   + sp.identity(M) * (Rn + self_t)).tocsc(), np.ones(M))
qr = spla.spsolve((Graw + sp.identity(M) * (Rn + self_t)).tocsc(), np.ones(M))
check('V14', 'G is asymmetric and is NOT symmetrized before solving',
      asym > 0.1 and (qr < 0).sum() == 0,
      f'||G-G^T||/||G|| = {asym:.2f}; negatives {int((qr<0).sum())} as '
      f'assembled against {int((qs<0).sum())} symmetrized')

check('V15', 'the diagonal is NOT dominant, so CG would be invalid',
      float(np.asarray(Graw.sum(1)).ravel().max()) > (Rn + self_t),
      f'diagonal {Rn + self_t:.3f} against a largest off-diagonal row sum '
      f'of {float(np.asarray(Graw.sum(1)).ravel().max()):.2f} m K/W')

# ---- V16: every generated file was built from THIS cache -------------
# v2.2 shipped a tab_ladder.tex that was byte-for-byte the v2.1 file,
# while the return note claimed every derived table had been recomputed.
# Stale LaTeX is still valid LaTeX, so nothing caught it until the paper
# printed 97.7 MWh beside a CSV that said 97.592. Each generator now
# stamps the SHA-256 of the cache it read, and this compares them.
# ---- V17: the plotted rates integrate to the budget ------------------
# v2.3 divided the per-tube rate by n_ft on top of q_prime, which is
# ALREADY per fin tube, and the figures came out a factor of three low.
# Nothing caught it because a rate curve has no independent reference.
# This gives it one: n_ft times the exact cumulative must be the energy
# the budget reports per well.
import pickle as _pk
if os.path.exists('rec_field.pkl'):
    _rec = _pk.load(open('rec_field.pkl', 'rb'))['rec']
    _c = case(1250.0)
    _worst = 0.0
    _det = []
    for _half, _key in (('charge', 'Q_charge_kJ'),
                        ('discharge', 'Q_discharge_kJ')):
        _tube = abs(_rec[_half]['Q_cum_J']) / 3.6e9          # MWh per tube
        _well = abs(_rec[_key]) / 3.6e6 / _rec['N_wells']    # MWh per well
        _rel = abs(_c.n_ft * _tube / _well - 1.0)
        _worst = max(_worst, _rel)
        _det.append(f'{_half} {_c.n_ft}x{_tube:.4f} vs {_well:.4f} MWh '
                    f'({_rel:.1e})')
    check('V17', 'n_ft x the per-tube integral equals dE per well',
          _worst < 1e-9, '; '.join(_det))

import packstamp
stale, unstamped, fresh = packstamp.check('results_pack_v2')
detail = f'{len(fresh)} current'
if stale:
    detail += '; STALE: ' + ', '.join(f'{f} ({h})' for f, h in stale)
if unstamped:
    detail += '; unstamped: ' + ', '.join(unstamped)
check('V16', 'every generated .tex was built from the current cache',
      not stale and not unstamped,
      detail + f'  [cache {packstamp.pack_hash()}]')

json.dump(RES, open('results_pack_v2/verification.json', 'w'), indent=1)
print(f"\n{len(RES) - FAIL}/{len(RES)} pass")
raise SystemExit(1 if FAIL else 0)
