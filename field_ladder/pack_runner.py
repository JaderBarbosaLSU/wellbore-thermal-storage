"""Results pack v2.0 -- the field ladder, run once and consistently.

ONE definition of the geometry, used by every rung:

  * idealized rungs (disc, double row) are built at the MEASURED Grissom
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
# The station spacing is part of the SPECIFICATION, so it is in the cache
# key: a run at a different discretization is a different run, not a
# replacement for one already there.
DS_STATION = G.DS_STATION
N_FIELD = 363
T_OP_DEFAULT = 10.0

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
        off_w = 0.0
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
        # ONE FIXED STATION SPACING for every run (Ground.stations_for).
        # v2.0 chose the count per run, as the first count that solved
        # with no negative strengths, and that made runs of the same sweep
        # incomparable: at L = 500 m the reference case landed on 25
        # stations and the 45 deg variant on 6, which moved e' by 2.8 %
        # against the 0.4 % the inclination was actually worth.
        nst = G.Ground.stations_for(L, spec.get('ds_station'))
        pr = gr.profiles(md, X, R_near, nst=nst)
        n_neg = getattr(gr, 'n_negative', 0)
        if n_neg:
            # V11 is an assertion now. No fallback.
            raise RuntimeError(
                f'V11 FAILED: {n_neg} of {len(P) * nst} strengths negative '
                f'at L = {L:.0f} m, ds = {spec.get("ds_station") or G.DS_STATION:.0f} m')
        key = {'superposition': 'temp', 'flux': 'flux'}[spec['ground']]
        md, Rg, tvd = pr['md'], pr[key], pr['tvd']
        q_spread = float(np.mean(pr['q_spread']))
        off_w = gr.off_station_weight
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
               ds_station=float(spec.get('ds_station') or G.DS_STATION),
               off_station_weight=float(off_w),
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
                          theta_max=theta_max, t_op_yr=t_op, casing=casing,
                          ds_station=float(DS_STATION)))

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
