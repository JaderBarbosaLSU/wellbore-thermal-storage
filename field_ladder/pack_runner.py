"""Results pack v2.0 -- the field ladder, run once and consistently.

ONE definition of the geometry, used by every rung:

  * idealised rungs (disc, double row) are built at the MEASURED Grissom
    spacing, 1.65 m, not the nominal 6 ft. Mixing the two is what made the
    prototype's R2 -> R3 step understate the layout penalty by 40 %;
  * R3 onwards use the real WellSTAR slot positions, 363 of them;
  * every rung carries the THUMS upper-zone casing, which is the §1.5
    default and is worth +36 % on e' against the legacy string;
  * the ground model is equal-temperature superposition unless stated.

Cached per run, so compute accumulates across sessions and nothing is
recomputed. Re-run until `remaining` is 0.
"""
import json
import os
import sys
import time

import numpy as np

import boreas as B
import boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
CACHE = 'pack_runs.json'
N_FIELD = 363
T_OP_DEFAULT = 10.0
# (geometry) -> a station count already known to solve cleanly. Persisted,
# because the search costs an extra sparse solve per rejected count and the
# sensitivities repeat the same geometries with a different t_op or theta.
_NST_MEMO_FILE = 'nst_memo.json'
try:
    _NST_MEMO = {tuple(k.split('|')): v
                 for k, v in json.load(open(_NST_MEMO_FILE)).items()}
except Exception:
    _NST_MEMO = {}

# The THUMS upper-zone string (SPE-2562) expressed in existing Case fields.
# D_in is the casing INSIDE diameter: it is what bounds the PCM, and the
# 201.2 mm ID is why this casing carries 46 % more material per tube.
THUMS_CASING = dict(D_in=201.2 / 25.4, r_casing_o=0.21910 / 2,
                    r_bore=0.31120 / 2, k_cement=0.7)
LEGACY_CASING = {}

L_SHOE = 1250.0          # §1.6 reference bound; sweeps may exceed it
LENGTHS = (100., 200., 300., 400., 500., 600., 800., 1000., 1250.)
LONG = (1600., 2000.)    # only for the rungs §5 asks to extend


def _positions(kind):
    if kind == 'disc':
        return G.disc_positions(N_FIELD, B0=G.B0_MEASURED)
    if kind == 'rows':
        return G.rows_positions(N_FIELD, B0=G.B0_MEASURED, w=G.W_DEF)
    if kind == 'slots':
        return G.grissom_positions(CSV)
    raise ValueError(kind)


def _case(L, casing, t_op):
    return B.Case(n_segments=99, L_ft=L / 0.3048,
                  t_operation_yr=t_op, **casing)


def _R_near(c):
    return float(np.atleast_1d(B.formation_resistance(c, np.array([0.0])))[0]
                 - B.ground_resistance(c))


def run_one(spec):
    """One ladder point. `spec` is the cache key's content."""
    L = spec['L_md']
    casing = THUMS_CASING if spec['casing'] == 'thums' else LEGACY_CASING
    t_op = spec.get('t_op_yr', T_OP_DEFAULT)
    c0 = _case(L, casing, t_op)
    R_near = _R_near(c0)

    if spec['positions'] == 'isolated':
        # No field: R_g alone, constant in MD. The two R0 rungs differ ONLY
        # in the TVD mapping -- an inclined well reaches cooler rock at the
        # same measured depth, which is the whole point of R0i.
        md = np.arange(0.0, L + 1e-9, 5.0)
        Rg = np.full(md.shape, B.ground_resistance(c0))
        if spec['trajectory'] == 'vertical':
            tvd = md.copy()
        else:
            inc = np.radians(G.s_curve(md, spec['kop_ft'] * 0.3048, 4.0,
                                       spec['theta_max'], 3000 * 0.3048))
            incm = 0.5 * (inc[1:] + inc[:-1])
            tvd = np.r_[0.0, np.cumsum(np.cos(incm) * 5.0)]
        q_spread = 0.0
        nst_used = len(md)
    else:
        P = _positions(spec['positions'])
        kop = None
        if spec['trajectory'] == 's_curve':
            if spec['kop_rule'] == 'staggered':
                kop, _ = G.staggered_kop(P)
            else:
                kop = np.full(len(P), spec['kop_ft'] * 0.3048)
        md, X = G.trajectories(P, L, kind=spec['trajectory'], kop=kop,
                               theta_max=spec['theta_max'])
        gr = G.Ground(k=c0.k_rock, rhocp=c0.rho_cp_rock, t_yr=t_op,
                      rb=c0.r_bore)
        # The station count is set by the GEOMETRY, not by the trajectory
        # type: the superposition needs the station spacing to be large
        # against the well spacing or the solve returns negative strengths.
        # See Ground.stations_for. A vertical well's profile is constant in
        # MD so it is insensitive to the count anyway; the S-curve is not,
        # which is why the rule is applied to both rather than special-cased.
        # STATION COUNT, chosen by the geometry and then verified.
        #
        # The superposition needs the station spacing to be large against
        # the well spacing (Ground.stations_for). That rule removes the
        # catastrophic failure -- at L = 100 m with 50 stations, 47 % of
        # the strengths go negative and R_g,eff reads -0.18 m K/W -- but it
        # does not by itself guarantee a clean solve: the count at which a
        # few negatives appear is erratic, because it depends on which
        # station of a neighbour happens to be nearest. Measured at
        # L = 100 m, the negative count runs 0, 28, 0, 38, 59, 24, 115, 0
        # for nst = 3..10 while R_g,eff stays within 3.11-3.18, so those
        # few negatives are discretisation noise and not a physical error.
        #
        # So: take the rule's count, then step DOWN to the first count that
        # solves with no negatives at all. Any of them gives the same
        # answer; this one is also clean, and V11 then holds by
        # construction rather than by tolerance.
        # Memoised by geometry AND FIELD AGE. The first version of this
        # key left t_op out, on the argument that the age enters only
        # through R_near, which shifts the diagonal uniformly. That was
        # wrong: t_op is inside the kernel, E1(d^2/4 a t), so it changes
        # which neighbours are coupled at all and therefore which station
        # count solves cleanly. The t_op = 30 yr sweep failed outright at
        # two lengths because it was handed the 10 yr count as its only
        # candidate -- which is also why the memo is now a STARTING POINT
        # and never the whole candidate list.
        memo_key = tuple(str(x) for x in
                         (spec['positions'], spec['trajectory'], L,
                          spec['theta_max'], spec['kop_rule'], t_op))
        nst0 = _NST_MEMO.get(memo_key) or G.Ground.stations_for(L)
        # GEOMETRIC stride, not one at a time. Stepping down singly from 50
        # costs up to 47 sparse solves -- over an hour for one ladder point,
        # and the clean counts are scattered anyway, so the fine search buys
        # nothing. The same stride reference_field uses, so the notebook and
        # the pack can land on the same count.
        cands = sorted({max(3, int(round(nst0 * f)))
                        for f in (1.0, .85, .7, .6, .5, .4, .3, .2, .1)},
                       reverse=True)
        pr, nst, n_neg = None, None, None
        for cand in cands:
            try:
                pr = gr.profiles(md, X, R_near, nst=cand)
            except MemoryError:
                continue          # a resource failure, not a result
            n_neg = getattr(gr, 'n_negative', 0)
            if n_neg == 0:
                nst = cand
                break
        if nst is not None:
            _NST_MEMO[memo_key] = nst
        if nst is None:
            raise RuntimeError(
                f'no clean station count among {cands} at L={L}')
        key = {'superposition': 'temp', 'flux': 'flux'}[spec['ground']]
        md, Rg, tvd = pr['md'], pr[key], pr['tvd']
        q_spread = float(np.mean(pr['q_spread']))
        nst_used = nst

    r = G.run_boreas(B, L, md, Rg, tvd_prof=tvd, n_segments=99,
                     t_operation_yr=t_op, **casing)
    out = dict(spec)
    out.update(loss=r['loss'], eta_storage=r['eta_storage'], RTE=r['RTE'],
               e_prime=r['e_prime'], E_out_kWh=r['E_out_kWh'],
               E_in_kWh=r['E_in_kWh'], T2d_C=r['T2d_C'],
               delivers=bool(r['delivers']), N_field=N_FIELD,
               tvd_bottom=float(tvd[-1]), R_near=R_near,
               Rg_mean=float(np.mean(Rg)), q_spread=q_spread,
               n_stations=int(nst_used), n_negative=0,
               field_MWh_per_cycle=N_FIELD * r['e_prime'] * L / 1000.0,
               beyond_shoe=bool(L > L_SHOE))
    return out


# ------------------------------------------------------------- the matrix
def build_specs():
    S = []

    def add(step, positions, trajectory, lengths, ground='superposition',
            kop_rule='staggered', kop_ft=300.0, theta_max=55.0,
            t_op=T_OP_DEFAULT, casing='thums'):
        for L in lengths:
            S.append(dict(step=step, positions=positions,
                          trajectory=trajectory, L_md=float(L),
                          ground=ground, kop_rule=kop_rule, kop_ft=kop_ft,
                          theta_max=theta_max, t_op_yr=t_op, casing=casing))

    # --- the ladder, §5
    add('R0v', 'isolated', 'vertical', LENGTHS)
    add('R0i', 'isolated', 's_curve', LENGTHS + LONG, kop_rule='fixed')
    add('R1', 'disc', 'vertical', LENGTHS)
    add('R2', 'rows', 'vertical', LENGTHS)
    add('R3', 'slots', 'vertical', LENGTHS)
    add('R4', 'slots', 's_curve', LENGTHS + LONG)

    # --- comparison column: other ground models, L = 1250 only
    for step, pos, traj in (('R1', 'disc', 'vertical'), ('R2', 'rows', 'vertical'),
                            ('R3', 'slots', 'vertical'), ('R4', 'slots', 's_curve')):
        add(step, pos, traj, (1250.,), ground='flux')

    # --- sensitivities on R4, §5
    add('R4_kop300', 'slots', 's_curve', LENGTHS + LONG, kop_rule='fixed',
        kop_ft=300.0)
    add('R4_th45', 'slots', 's_curve', LENGTHS, theta_max=45.0)
    add('R4_th70', 'slots', 's_curve', LENGTHS, theta_max=70.0)
    add('R4_top1', 'slots', 's_curve', LENGTHS, t_op=1.0)
    add('R4_top30', 'slots', 's_curve', LENGTHS, t_op=30.0)

    # --- V1/V7 regression needs the legacy casing at the reference
    add('R3_legacy', 'slots', 'vertical', (1250.,), casing='legacy')
    add('R4_legacy', 'slots', 's_curve', (1250.,), casing='legacy')
    return S


def key_of(spec):
    return json.dumps({k: spec[k] for k in sorted(spec)}, sort_keys=True)


def _worker(spec):
    try:
        return spec, run_one(spec)
    except Exception as e:
        return spec, dict(spec, error=f'{type(e).__name__}: {e}'[:200])


def main(budget=110.0, workers=4):
    """Parallel across runs. A field point costs ~85 s -- the superposition
    solve over 363 wells x 50 stations dominates and barely depends on L --
    so serial execution would be two hours. Set OMP_NUM_THREADS=1: the
    workers otherwise contend over BLAS threads and the speedup collapses
    from 4x to barely 1.5x."""
    import multiprocessing as mp
    specs = build_specs()
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    todo = [s for s in specs if key_of(s) not in cache]
    if not todo:
        print(f'nothing to do; {len(specs)} runs cached')
        return 0
    t0, done = time.time(), 0
    pool = mp.Pool(processes=workers)
    try:
        it = pool.imap_unordered(_worker, todo)
        while True:
            left_s = budget - (time.time() - t0)
            if left_s <= 0:
                break
            try:
                spec, res = it.next(timeout=left_s)
            except (mp.TimeoutError, StopIteration):
                break
            cache[key_of(spec)] = res
            done += 1
            json.dump(cache, open(CACHE, 'w'), indent=0, default=float)
            print(f"  {spec['step']:10s} L={spec['L_md']:6.0f} "
                  f"{spec['ground'][:4]:5s} ({time.time()-t0:.0f} s)", flush=True)
    finally:
        pool.terminate(); pool.join()
    left = len(todo) - done
    print(f'\ndone {done}, remaining {left}, total {len(specs)}')
    return left


if __name__ == '__main__':
    sys.exit(0 if main(float(sys.argv[1]) if len(sys.argv) > 1 else 110.0) == 0
             else 1)
