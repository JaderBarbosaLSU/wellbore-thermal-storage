"""Handoff 6.2c: the three casing unknowns, quantified at R3 and R4.

The paper has been carrying 'less than 0.2 points or 2.2 %' as a
placeholder with no source in the pack. These are the numbers, and one of
them is larger than that placeholder: at R4 the unknowns are worth up to
2.7 points of loss and 5.8 % of e'. The v2.0 estimate was made at R3,
where the field is well shielded and the near-well resistance is a small
share of the total; at R4 the shielding is weaker and the same annulus
change matters more.
"""
import csv, json
rows = [v for v in json.load(open('extra_runs.json')).values()
        if v.get('kind') == 'casing' and 'error' not in v]
ref = {r['step']: r for r in rows if r['tag'] == 'ref'}
GROUP = {'toc': ['toc_surface'],        # cement to surface vs mud to shoe
         'mud': ['mud_low', 'mud_high'],
         'weight': ['weight_36']}
out, worst_l, worst_e = [], 0.0, 0.0
print(f"{'unknown':9s}{'step':5s}{'variant':13s}{'loss %':>9}"
      f"{'d loss [pt]':>13}{'d e-prime [%]':>15}")
for g, tags in GROUP.items():
    for step in ('R3', 'R4'):
        b = ref[step]
        for t in tags:
            r = [x for x in rows if x['tag'] == t and x['step'] == step][0]
            dl = 100 * (r['loss'] - b['loss'])
            de = 100 * (r['e_prime'] / b['e_prime'] - 1)
            out.append(dict(unknown=g, step=step, variant=t,
                            loss=r['loss'], e_prime=r['e_prime'],
                            d_loss_pt=dl, d_eprime_pct=de))
            print(f'{g:9s}{step:5s}{t:13s}{100*r["loss"]:9.3f}'
                  f'{dl:+13.3f}{de:+15.2f}')
        mx_l = max(abs(x['d_loss_pt']) for x in out if x['unknown'] == g)
        mx_e = max(abs(x['d_eprime_pct']) for x in out if x['unknown'] == g)
        worst_l, worst_e = max(worst_l, mx_l), max(worst_e, mx_e)
with open('results_pack_v2/casing_sensitivity.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)
summ = {}
for g in GROUP:
    sub = [x for x in out if x['unknown'] == g]
    summ[g] = dict(max_d_loss_pt=max(abs(x['d_loss_pt']) for x in sub),
                   max_d_eprime_pct=max(abs(x['d_eprime_pct']) for x in sub))
summ['max'] = dict(max_d_loss_pt=worst_l, max_d_eprime_pct=worst_e)
json.dump(summ, open('casing_sens.json', 'w'), indent=1)
print(f'\nlargest over all three: {worst_l:.2f} points of loss, '
      f'{worst_e:.2f} % of e-prime')
for g, v in summ.items():
    if g != 'max':
        print(f"  {g:7s} {v['max_d_loss_pt']:5.2f} pt   "
              f"{v['max_d_eprime_pct']:5.2f} %")
