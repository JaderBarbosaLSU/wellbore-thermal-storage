"""Run the model from dbhe_master.ipynb without opening the notebook.

Loads every code cell except the demo/dashboard/chart cells into a namespace,
so sweeps and regression checks call the notebook's own functions rather than a
re-implementation. If the notebook and this harness ever disagree, the harness
is wrong by definition.

    python harness.py baseline out.json          # default case, annual samples
    python harness.py sweep 4500,4677 3,4,5 o.json
    python harness.py check reference/v10.2_baseline.json
"""
import contextlib
import io
import json
import os
import re
import sys
import warnings

import matplotlib
matplotlib.use('Agg')

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.join(HERE, '..', 'geothermal_cooling', 'dbhe_master.ipynb')

# cells containing any of these run the demo or draw figures -> skip on import
SKIP = ('run_facility(Lbh_default', 'psychro_background')

# what a fixture records, per year
TRACK = ('T_ret', 'T_wr', 'T_wi', 'T_acc', 'T_ri', 'T_ro', 'T_rd',
         'T_pd', 'T_pi', 'w_ri', 'w_rd', 'w_pd', 'w_pi_req',
         'UA_cool', 'd_rA_ctrl', 'dW_process')


def load(nb_path=None):
    # resolved at call time, not at import time: binding NB as a default argument
    # would silently ignore any later override and make a regression check that
    # cannot fail.
    nb = json.load(open(nb_path or NB))
    g = {'__name__': 'nbmodel'}
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        s = ''.join(c['source'])
        if any(k in s for k in SKIP):
            continue
        s = '\n'.join(l for l in s.splitlines() if not re.match(r'\s*[!%]', l))
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(s, '<nb>', 'exec'), g)
    return g


def notebook_version(nb_path=None):
    nb = json.load(open(nb_path or NB))
    for c in nb['cells']:
        m = re.search(r'NB_VERSION\s*=\s*"([^"]+)"', ''.join(c['source']))
        if m:
            return m.group(1)
    return 'unknown'


def run(Lbh, flow_m3_h, kt_key='kt_insulated', g=None):
    g = g or load()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with contextlib.redirect_stdout(io.StringIO()):
            out, dbhe, acc = g['run_facility'](Lbh, flow_m3_h, g[kt_key], verbose=False)
    n = len(out['T_ret'])
    spy = 365.25 / g['n_days_per_step']
    years = list(range(1, int(g['n_years']) + 1))
    idx = [min(n - 1, int(round(y * spy)) - 1) for y in years]
    return dict(
        case=dict(Lbh=Lbh, flow_m3_h=flow_m3_h, kt=kt_key,
                  BHT=g['T_surf'] + g['Gg'] * Lbh,
                  mdot_loop=float(dbhe['mdot_loop'])),
        years=years,
        series={k: [round(float(out[k][i]), 6) for i in idx] for k in TRACK if k in out},
        wheel=dict(calls=int(g['DW_REGIME']['calls']),
                   reversed=int(g['DW_REGIME']['reversed']),
                   worst_dw=float(g['DW_REGIME']['worst_dw'])),
    )


def compare(ref, new, rtol=1e-9, atol=1e-9):
    """Return a list of (key, year, ref, new, absdiff) for every value that moved."""
    diffs = []
    for k, ref_v in ref['series'].items():
        new_v = new['series'].get(k)
        if new_v is None:
            diffs.append((k, None, 'present', 'MISSING', float('nan')))
            continue
        for y, a, bb in zip(ref['years'], ref_v, new_v):
            if abs(a - bb) > atol + rtol * abs(a):
                diffs.append((k, y, a, bb, abs(a - bb)))
    for k in set(new['series']) - set(ref['series']):
        diffs.append((k, None, 'ABSENT', 'new series', float('nan')))
    return diffs


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'baseline':
        g = load()
        r = run(g['Lbh_default'], g['flow_rate_m3_h_demo'], 'kt_demo_key'
                if 'kt_demo_key' in g else 'kt_insulated', g)
        r['notebook_version'] = notebook_version()
        json.dump(r, open(sys.argv[2], 'w'), indent=1)
        print(f"baseline written for notebook v{r['notebook_version']}: "
              f"T_ret[20]={r['series']['T_ret'][-1]:.3f} C")
    elif cmd == 'sweep':
        from multiprocessing import Pool
        depths = [float(x) for x in sys.argv[2].split(',')]
        flows = [float(x) for x in sys.argv[3].split(',')]
        grid = [(L, q) for L in depths for q in flows]
        G = {}

        def _init():
            G['g'] = load()

        def _run(a):
            try:
                return run(a[0], a[1], 'kt_insulated', G['g'])
            except Exception as e:
                return dict(case=dict(Lbh=a[0], flow_m3_h=a[1]),
                            error=f'{type(e).__name__}: {e}')
        with Pool(2, initializer=_init) as p:
            res = list(p.imap_unordered(_run, grid))
        json.dump(res, open(sys.argv[4], 'w'), indent=1)
        print(f'{len(res)} cases written')
    elif cmd == 'check':
        ref = json.load(open(sys.argv[2]))
        g = load()
        new = run(ref['case']['Lbh'], ref['case']['flow_m3_h'], ref['case']['kt'], g)
        d = compare(ref, new)
        cur = notebook_version()
        if not d:
            print(f'PASS  notebook v{cur} reproduces {sys.argv[2]} '
                  f"(recorded under v{ref.get('notebook_version', '?')}) exactly.")
        else:
            print(f'{len(d)} value(s) moved vs {sys.argv[2]}  '
                  f"(ref v{ref.get('notebook_version','?')} -> now v{cur})")
            print(f"{'series':<12}{'yr':>4}{'reference':>14}{'now':>14}{'|diff|':>12}")
            for k, y, a, bb, dd in d[:40]:
                ys = '' if y is None else str(y)
                fa = f'{a:14.6f}' if isinstance(a, float) else f'{a:>14}'
                fb = f'{bb:14.6f}' if isinstance(bb, float) else f'{bb:>14}'
                fd = '' if dd != dd else f'{dd:12.3e}'
                print(f'{k:<12}{ys:>4}{fa}{fb}{fd}')
            if len(d) > 40:
                print(f'... and {len(d)-40} more')
            sys.exit(1)
    else:
        print(__doc__)
        sys.exit(2)
