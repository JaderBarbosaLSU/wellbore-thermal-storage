"""tvd_bottom meant the last STATION, not the bottom of the store.

The station grid at ds = 50 m puts its last centre at MD = 1225 m, so
pr['tvd'][-1] is the TVD there -- 890.9 m -- while the store ends at
MD = 1250 m, where the trajectory has reached 908.6 m. Every 'the
S-curve reaches 891 m TVD at 1250 m MD' is therefore about 18 m, 2 %,
too shallow.

This is pure geometry: it needs no march, so the cache is corrected in
place rather than recomputed, and tvd_last_station is kept beside it so
the old number is still traceable.
"""
import json
import numpy as np
import boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
P = G.grissom_positions(CSV)
KOP, _ = G.staggered_kop(P)
N = len(P)
POS = {'slots': P, 'disc': G.disc_positions(N, B0=G.B0_MEASURED),
       'rows': G.rows_positions(N, B0=G.B0_MEASURED, w=G.W_DEF)}
_cache = {}


def tvd_at_bottom(spec):
    """Mean TVD of the field at MD = L, from the trajectory itself."""
    L = float(spec['L_md'])
    traj = spec['trajectory']
    key = (spec['positions'], traj, L, spec['theta_max'], spec['kop_rule'],
           spec.get('kop_ft'))
    if key in _cache:
        return _cache[key]
    if spec['positions'] == 'isolated':
        if traj == 'vertical':
            v = L
        else:
            inc = np.radians(G.s_curve(np.arange(0.0, L + 1e-9, 5.0),
                                       spec['kop_ft'] * G.FT, 4.0,
                                       spec['theta_max'], 3000 * G.FT))
            im = 0.5 * (inc[1:] + inc[:-1])
            v = float(np.sum(np.cos(im) * 5.0))
    else:
        Pp = POS[spec['positions']]
        kop = None
        if traj == 's_curve':
            kop = (G.staggered_kop(Pp)[0] if spec['kop_rule'] == 'staggered'
                   else np.full(len(Pp), spec['kop_ft'] * G.FT))
        md, X = G.trajectories(Pp, L, kind=traj, kop=kop,
                               theta_max=spec['theta_max'])
        v = float(X[:, -1, 2].mean())
    _cache[key] = v
    return v


c = json.load(open('pack_runs.json'))
n = 0
for k, v in c.items():
    if 'error' in v or 'tvd_bottom' not in v:
        continue
    v.setdefault('tvd_last_station', v['tvd_bottom'])
    new = tvd_at_bottom(v)
    if abs(new - v['tvd_bottom']) > 1e-9:
        n += 1
    v['tvd_bottom'] = new
json.dump(c, open('pack_runs.json', 'w'), indent=0, default=float)
print(f'corrected tvd_bottom on {n} of {len(c)} rows')
for step in ('R4', 'R3', 'R0i'):
    r = [v for v in c.values() if v.get('step') == step
         and abs(v.get('L_md', 0) - 1250.0) < 1e-6
         and v.get('ground') == 'superposition']
    if r:
        print(f"  {step:4s} last station {r[0]['tvd_last_station']:7.1f} m"
              f"   bottom of store {r[0]['tvd_bottom']:7.1f} m")
