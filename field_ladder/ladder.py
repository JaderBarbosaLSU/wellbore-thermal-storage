"""Ladder rungs at L = 1250 m, equal-temperature superposition.

R1 and R2 are built at the MEASURED Grissom spacing (1.65 m), not the
nominal 6 ft, so the step R2 -> R3 isolates the LAYOUT rather than mixing
it with a 10 % spacing change. Pass B0_NOMINAL to reproduce the prototype.
"""
import json, os, time, numpy as np, boreas as B, boreas_ground as G

CACHE = 'ladder_v2.json'


def load():
    return json.load(open(CACHE)) if os.path.exists(CACHE) else {}


def save(d):
    json.dump(d, open(CACHE, 'w'), indent=0)

CSV, L, N = 'THUMS_wells_WellSTAR_clean.csv', 1250.0, 363


def rung(name, P, kind, kop=None, tag=''):
    c0 = B.Case(n_segments=99, L_ft=L / 0.3048)
    R_near = float(np.atleast_1d(B.formation_resistance(c0, np.array([0.0])))[0]
                   - B.ground_resistance(c0))
    md, X = G.trajectories(P, L, kind=kind, kop=kop)
    gr = G.Ground(k=c0.k_rock, rhocp=c0.rho_cp_rock,
                  t_yr=c0.t_operation_yr, rb=c0.r_bore)
    pr = gr.profiles(md, X, R_near, nst=50)
    r = G.run_boreas(B, L, pr['md'], pr['temp'], tvd_prof=pr['tvd'], n_segments=99)
    return dict(rung=name, tag=tag, loss=r['loss'], RTE=r['RTE'],
                e_prime=r['e_prime'], T2d=r['T2d_C'], delivers=r['delivers'],
                Rg_mean=float(np.mean(pr['temp'])), R_near=R_near)


JOBS = []
for _b0, _tag in ((G.B0_MEASURED, '1.65 measured'), (G.B0_NOMINAL, '1.83 nominal')):
    JOBS.append(('R1_disc', lambda b=_b0: G.disc_positions(N, B0=b), 'vertical', _tag))
    JOBS.append(('R2_rows', lambda b=_b0: G.rows_positions(N, B0=b, w=G.W_DEF),
                 'vertical', _tag))
JOBS.append(('R3_grissom_vert', lambda: G.grissom_positions(CSV), 'vertical',
             'real slots'))

if __name__ == '__main__':
    import sys
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
    cache, t0 = load(), time.time()
    for name, mk, kind, tag in JOBS:
        key = f'{name}|{tag}'
        if key in cache:
            continue
        if time.time() - t0 > budget:
            print('budget spent; re-run to continue', flush=True)
            break
        cache[key] = rung(name, mk(), kind, tag=tag)
        save(cache)
        print(f'  done {key} ({time.time()-t0:.0f} s)', flush=True)
    rows = [cache[f'{n}|{t}'] for n, _, _, t in JOBS if f'{n}|{t}' in cache]
    print(f'\n{len(rows)}/{len(JOBS)} rungs cached\n')
    print(f"{'rung':18}{'spacing':>16}{'loss%':>8}{'RTE':>8}{'e_prime':>9}"
          f"{'Rg_mean':>10}{'T2d':>8}")
    for o in rows:
        print(f"{o['rung']:18}{o['tag']:>16}{100*o['loss']:>8.1f}{o['RTE']:>8.3f}"
              f"{o['e_prime']:>9.3f}{o['Rg_mean']:>10.2f}{o['T2d']:>8.1f}")
