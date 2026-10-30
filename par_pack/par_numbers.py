"""par_numbers.tex -- macros for the parametric prose. Same rule as the
results pack: nothing is retyped. Raises if any macro would be unusable."""
import json
import numpy as np
import par_grids as A

OUT = 'par_pack'
R = {}
for k, v in json.load(open('par_runs.json')).items():
    R.setdefault(json.loads(k)[0], []).append(v)
ok = lambda rs: [x for x in rs if 'error' not in x]

M = {}
def add(n, v, fmt='{:.3f}', note=''):
    M[n] = (fmt.format(v) if not isinstance(v, str) else v, note)

# ---- h_m, the dominant lever
row = sorted([x for x in ok(R['hm_Tm']) if x['T_m_C'] == 150.],
             key=lambda x: x['h_m'])
for x in row:
    t = f'Hm{int(x["h_m"]/1000)}'
    add(f'{t}Eprime', x['e_prime'], '{:.3f}', f"e' at h_m={x['h_m']/1000:.0f}")
    add(f'{t}Wells', x['wells_per_MWe'], '{:.1f}', '')
    add(f'{t}Loss', x['loss'], '{:.3f}', '')
add('HmRangeFactor', row[-1]['h_m'] / row[0]['h_m'], '{:.1f}',
    'h_m range as a factor')
add('HmEprimeFactor', row[-1]['e_prime'] / row[0]['e_prime'], '{:.1f}',
    "e' response to it")
col = sorted([x for x in ok(R['hm_Tm']) if abs(x['h_m'] - 380e3) < 1],
             key=lambda x: x['T_m_C'])
add('TmEprimeFactor', max(x['e_prime'] for x in col)
    / min(x['e_prime'] for x in col), '{:.1f}', "e' response to T_m alone")

# ---- cascade
c = ok(R['cascade'])
ref = [x for x in c if x['T_m_C'] == 150. and x['DT_3C_2C'] == 55.][0]
be = max(c, key=lambda x: x['e_prime'])
br = max(c, key=lambda x: x['RTE'])
for tag, x in (('CascRef', ref), ('CascBestEp', be), ('CascBestRTE', br)):
    add(f'{tag}Tm', x['T_m_C'], '{:.0f}', '')
    add(f'{tag}Grading', x['DT_3C_2C'], '{:.0f}', '')
    add(f'{tag}Eprime', x['e_prime'], '{:.3f}', '')
    add(f'{tag}RTE', x['RTE'], '{:.3f}', '')
add('CascInfeasible', float(sum(1 for x in R['cascade'] if 'error' in x)),
    '{:.0f}', 'hatched cells')
add('TmCap', max(A.TM), '{:.0f}', 'T_m cap, C')

# ---- geometry
nft1 = sorted([x for x in ok(R['L_nft']) if x['n_ft'] == 1],
              key=lambda x: x['L_md'])
add('NftOnePumpAtShoe', nft1[-1]['f_pump'] * 100, '{:.0f}',
    'pumping at n_ft=1, L=shoe, % of gross')
add('NftOneRTEAtShoe', nft1[-1]['RTE'], '{:+.3f}', '')
add('LShoeFt', A.L_SHOE_FT, '{:.0f}', 'map L cap, ft')
add('LShoeM', max(A.LW_M), '{:.0f}', 'the same in m MD')
add('LRefM', A.L_REF, '{:.0f}', '')

# ---- N_lay confound
ef = np.array([x['e_prime'] for x in sorted(ok(R['Nlay_free']),
                                            key=lambda r: r['N_lay'])])
ep = np.array([x['e_prime'] for x in sorted(ok(R['Nlay_pinned']),
                                            key=lambda r: r['N_lay'])])
sf = 100 * (ef.max() - ef.min()) / ef.mean()
sp = 100 * (ep.max() - ep.min()) / ep.mean()
add('NlaySpreadFree', sf, '{:.2f}', "e' spread, span free, %")
add('NlaySpreadPinned', sp, '{:.2f}', "e' spread, span pinned, %")
add('NlayConfoundShare', 100 * (1 - sp / sf), '{:.0f}',
    '% of the apparent effect that is the span')
add('NlaySpanPinned', A.SPAN0, '{:.1f}', 'pinned span, K')
add('ParNPoints', float(sum(len(v) for v in R.values())), '{:.0f}', '')

BAD = [n for n, (v, _) in M.items()
       if (not str(v).strip()) or str(v).strip().lower()
       in ('nan', 'inf', '-inf', 'none', '--')]
if BAD:
    raise SystemExit('unusable macros: ' + ', '.join(sorted(BAD)))
with open(f'{OUT}/par_numbers.tex', 'w') as f:
    f.write('% par_numbers.tex -- GENERATED from par_runs.json. Do not edit\n'
            '% and do not retype: \\input this file and use the macros.\n'
            f'% {len(M)} macros\n\n')
    for n, (v, note) in sorted(M.items()):
        pad = ' ' * max(1, 32 - len(n))
        f.write(f'\\newcommand{{\\Par{n}}}{{{v}}}{pad}'
                + (f'% {note}\n' if note else '\n'))
print(f'wrote {OUT}/par_numbers.tex ({len(M)} macros)')
