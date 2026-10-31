"""v2.3 extra runs: the cutoff robustness sweep (6.5) and the three
casing sensitivities (6.2c). Resumable, own cache."""
import json, os, sys, time
import numpy as np
import boreas as B, boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
CACHE = 'extra_runs.json'
THUMS = dict(D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
             k_cement=0.7)
P = G.grissom_positions(CSV)
KOP, _ = G.staggered_kop(P)

# ---- 6.5: does the cutoff survive a finer station grid? ---------------
# The margin is 5 %: the marginal is 0.124 over 400-500 m against a
# threshold of 0.131, and at L = 400 m the discretization drift is about
# 0.4 points of loss. So the question is fair and the answer has to be
# measured, not argued.
CUTOFF_L = (100., 200., 300., 400., 500., 600., 1250.)
CUTOFF_DS = (50.0, 25.0, 12.5)

# ---- 6.2c: the casing unknowns ----------------------------------------
# k_cement is the ANNULUS fill. Top of cement at the surface means cement
# all the way (1.0); at the shoe it means mud all the way, and the mud
# range 0.6-1.5 brackets that. Casing weight enters as the inside
# diameter, which is what bounds the PCM: 8-5/8 in 32# is 201.2 mm,
# 36# is 198.8 mm.
CASING = [('ref', dict()),
          ('toc_surface', dict(k_cement=1.0)),
          ('mud_low', dict(k_cement=0.6)),
          ('mud_high', dict(k_cement=1.5)),
          ('weight_36', dict(D_in=198.8/25.4))]


def field(L, traj, ds):
    c = B.Case(n_segments=99, L_ft=L/0.3048, t_operation_yr=10.0, **THUMS)
    return G.cached_reference_field(
        B, c, CSV, trajectory=traj, kop_rule='staggered', theta_max=55.0,
        kop_ft=300.0, L_md=L, ds=ds, verbose=False)[:3]


def run(spec):
    L, ds, traj = spec['L_md'], spec['ds'], spec['traj']
    md, Rg, tvd = field(L, traj, ds)
    kw = dict(THUMS); kw.update(spec.get('case', {}))
    r = G.run_boreas(B, L, md, Rg, tvd_prof=tvd, n_segments=99,
                     t_operation_yr=10.0, **kw)
    return dict(spec, loss=r['loss'], e_prime=r['e_prime'], RTE=r['RTE'],
                E_out_kWh=r['E_out_kWh'], E_in_kWh=r['E_in_kWh'],
                eta_storage=r['eta_storage'], T2d_C=r['T2d_C'],
                Rg_mean=float(np.mean(Rg)))


def specs():
    S = []
    for ds in CUTOFF_DS:
        for L in CUTOFF_L:
            S.append(dict(kind='cutoff', L_md=L, ds=ds, traj='s_curve',
                          case={}))
    for tag, kw in CASING:
        for step, traj in (('R3', 'vertical'), ('R4', 's_curve')):
            S.append(dict(kind='casing', tag=tag, step=step, L_md=1250.0,
                          ds=50.0, traj=traj, case=kw))
    return S


def key(s):
    return json.dumps(s, sort_keys=True, default=str)


def main(budget=150.0):
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    todo = [s for s in specs() if key(s) not in cache]
    print(f'{len(cache)} cached, {len(todo)} to run', flush=True)
    t0 = time.time()
    for s in todo:
        if time.time() - t0 > budget:
            print('BUDGET', flush=True); break
        try:
            r = run(s)
            lab = s.get('tag', '') or f"ds={s['ds']}"
            print(f"  {s['kind']:7s} {lab:12s} {s.get('step',''):3s} "
                  f"L={s['L_md']:6.0f} loss={100*r['loss']:6.3f}% "
                  f"e'={r['e_prime']:.4f}", flush=True)
        except Exception as e:
            r = dict(s, error=f'{type(e).__name__}: {e}'[:200])
            print(f"  ERROR {r['error'][:80]}", flush=True)
        cache[key(s)] = r
        json.dump(cache, open(CACHE, 'w'), indent=0, default=float)
    left = sum(1 for s in specs() if key(s) not in cache)
    print(f'{len(cache)} cached, {left} remaining', flush=True)
    return left


if __name__ == '__main__':
    sys.exit(0 if main(float(sys.argv[1]) if len(sys.argv) > 1 else 150.) == 0
             else 1)
