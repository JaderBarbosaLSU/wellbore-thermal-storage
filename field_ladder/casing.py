"""R3 (Grissom slots, vertical) at L = 1250 m under casing variants.

The THUMS upper-zone program maps onto three existing Case fields, so the
uniform case needs no new code:
    D_in        = 201.2/25.4 in   (8-5/8 in 32# INSIDE diameter -> PCM bound)
    r_casing_o  = 0.10955 m       (219.1 mm OD)
    r_bore      = 0.1556  m       (311.2 mm, 12-1/4 in hole)
    k_cement    = k_mud or k_cem  (TOC at the shoe = mud everywhere)
A depth-varying TOC would need a profile; the spread below shows why that
is not worth building yet.
"""
import json, os, sys, time, numpy as np, boreas as B, boreas_ground as G

CSV, L, CACHE = 'THUMS_wells_WellSTAR_clean.csv', 1250.0, 'casing_v2.json'
THUMS = dict(D_in=201.2 / 25.4, r_casing_o=0.21910 / 2, r_bore=0.31120 / 2)
VARIANTS = {
    'legacy':             {},
    'thums, mud 0.7':     {**THUMS, 'k_cement': 0.7},
    'thums, cement 1.0':  {**THUMS, 'k_cement': 1.0},
    'thums, mud 0.6':     {**THUMS, 'k_cement': 0.6},
    'thums, mud 1.5':     {**THUMS, 'k_cement': 1.5},
    'thums 36#, mud 0.7': {**THUMS, 'D_in': 198.8 / 25.4, 'k_cement': 0.7},
}


def run(kw):
    c0 = B.Case(n_segments=99, L_ft=L / 0.3048, **kw)
    R_near = float(np.atleast_1d(B.formation_resistance(c0, np.array([0.0])))[0]
                   - B.ground_resistance(c0))
    P = G.grissom_positions(CSV)
    md, X = G.trajectories(P, L, kind='vertical')
    gr = G.Ground(k=c0.k_rock, rhocp=c0.rho_cp_rock, t_yr=c0.t_operation_yr,
                  rb=c0.r_bore)                      # r_b follows the hole
    pr = gr.profiles(md, X, R_near, nst=50)
    r = G.run_boreas(B, L, pr['md'], pr['temp'], tvd_prof=pr['tvd'],
                     n_segments=99, **kw)
    A = B.pcm_capacities(c0)[0]
    return dict(loss=r['loss'], RTE=r['RTE'], e_prime=r['e_prime'],
                T2d=r['T2d_C'], R_near=R_near, pcm_cm2=A * 1e4,
                Rg_mean=float(np.mean(pr['temp'])))


if __name__ == '__main__':
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 110.0
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    t0 = time.time()
    for name, kw in VARIANTS.items():
        if name in cache:
            continue
        if time.time() - t0 > budget:
            print('budget spent; re-run', flush=True); break
        cache[name] = run(kw)
        json.dump(cache, open(CACHE, 'w'), indent=0)
        print(f'  done {name} ({time.time()-t0:.0f} s)', flush=True)
    print(f"\n{'variant':22}{'R_near':>9}{'PCM cm2':>9}{'loss%':>8}{'RTE':>8}{'e_prime':>9}")
    for n in VARIANTS:
        if n in cache:
            v = cache[n]
            print(f"{n:22}{v['R_near']:>9.4f}{v['pcm_cm2']:>9.1f}"
                  f"{100*v['loss']:>8.1f}{v['RTE']:>8.3f}{v['e_prime']:>9.3f}")
