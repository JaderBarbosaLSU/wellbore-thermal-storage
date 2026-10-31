"""Handoff 6.5: does the 400-500 m cutoff survive a finer station grid?

The margin is 5 % -- the marginal is 0.124 over 400-500 m against a
threshold of 0.131 -- and the discretization drift at L = 400 m is about
0.4 points of loss. So this is measured, not argued.
"""
import csv, json
import numpy as np

rows = [v for v in json.load(open('extra_runs.json')).values()
        if v.get('kind') == 'cutoff' and 'error' not in v]
out, verdict = [], True
print(f"{'ds':>7}{'first 100 m':>13}{'threshold':>11}  marginals by interval")
for ds in sorted({r['ds'] for r in rows}, reverse=True):
    rr = sorted([r for r in rows if r['ds'] == ds], key=lambda r: r['L_md'])
    L = np.array([r['L_md'] for r in rr])
    Eo = np.array([r['E_out_kWh'] for r in rr])
    Ei = np.array([r['E_in_kWh'] for r in rr])
    first = Eo[0] / Ei[0]            # E(0) = 0, so this IS RTE(100)
    thr = 0.5 * first
    marg = np.diff(Eo) / np.diff(Ei)
    lo = hi = None
    for a, b, m in zip(L[:-1], L[1:], marg):
        if m < thr:
            lo, hi = a, b
            break
    txt = '  '.join(f'{a:.0f}-{b:.0f}:{m:.3f}'
                    for a, b, m in zip(L[:-1], L[1:], marg))
    print(f'{ds:7.2f}{first:13.4f}{thr:11.4f}  {txt}')
    print(f'{"":7s}{"cutoff":>13s}{"":11s}  '
          f'{(f"{lo:.0f}-{hi:.0f} m MD" if lo else "none in range")}')
    out.append(dict(ds_m=ds, rte_first_100m=first, threshold=thr,
                    cutoff_lo=lo, cutoff_hi=hi, marginals=txt))
    if (lo, hi) != (400.0, 500.0):
        verdict = False
with open('results_pack_v2/cutoff_robustness.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)
print(f'\nCUTOFF ROBUST: {"yes" if verdict else "no"} '
      f'-- 400-500 m at every station spacing tested'
      if verdict else '\nCUTOFF ROBUST: no')
json.dump(dict(robust='yes' if verdict else 'no', rows=out),
          open('cutoff_robust.json', 'w'), indent=1)
