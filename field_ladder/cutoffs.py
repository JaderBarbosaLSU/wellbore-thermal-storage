"""The cutoff criterion, applied as DEFINED.

'marginal RTE below half that of the first 100 m'. The first 100 m is the
interval 0 -> 100, and since E_out(0) = E_in(0) = 0 that marginal IS
RTE(100). v2.0 used the 100 -> 200 interval instead, which is a different
and larger number, and the error was compounded because B1 had inflated it.
"""
import csv, json, sys
import packstamp  # noqa: F401 -- CSV, but keep the import contract
import numpy as np

c = json.load(open('pack_runs.json'))
rows = [v for v in c.values() if 'loss' in v and v['ground'] == 'superposition']
steps = sorted({r['step'] for r in rows if r['step'].startswith('R4')
                or r['step'] in ('R0v', 'R0i', 'R1', 'R2', 'R3')})
out = []
print(f"{'step':12}{'RTE(100)':>10}{'threshold':>11}{'cutoff':>14}")
for s in steps:
    rr = sorted([r for r in rows if r['step'] == s], key=lambda r: r['L_md'])
    if len(rr) < 3:
        continue
    L = np.array([r['L_md'] for r in rr])
    Eo = np.array([r['E_out_kWh'] for r in rr])
    Ei = np.array([r['E_in_kWh'] for r in rr])
    if abs(L[0] - 100.0) > 1e-6:
        continue
    first = Eo[0] / Ei[0]                 # the 0 -> 100 m interval
    thr = 0.5 * first
    marg = np.diff(Eo) / np.diff(Ei)
    # FIRST crossing only. A marginal that comes back above the threshold
    # later does not undo the cutoff -- and at R4 it does come back, over
    # 1250-1600 m, which is beyond the shoe and therefore an
    # extrapolation of the geometry rather than a design. Reported in its
    # own column so it cannot be read as a second cutoff.
    lo = hi = None
    recross = ''
    below = False
    for a, b, m in zip(L[:-1], L[1:], marg):
        if m < thr and not below:
            lo, hi, below = a, b, True
        elif m >= thr and below and not recross:
            recross = f'{a:.0f}-{b:.0f}'
    out.append(dict(step=s, rte_first_100m=first, threshold=thr,
                    cutoff_lo=lo, cutoff_hi=hi,
                    recrosses_at=recross,
                    beyond_shoe=bool(recross and float(recross.split('-')[0])
                                     >= 1250.0),
                    marginals=';'.join(f'{a:.0f}-{b:.0f}:{m:.3f}'
                                       for a, b, m in zip(L[:-1], L[1:], marg))))
    print(f"{s:12}{first:>10.3f}{thr:>11.3f}"
          f"{(f'{lo:.0f}-{hi:.0f} m' if lo else 'none'):>14}"
          f"{('  re-crosses ' + recross if recross else ''):>22}")
with open('results_pack_v2/cutoffs.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)
print(f'\nwrote results_pack_v2/cutoffs.csv ({len(out)} steps)')
