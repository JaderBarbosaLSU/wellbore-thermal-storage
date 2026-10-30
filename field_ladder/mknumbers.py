"""boreas_numbers.tex -- every number the prose quotes, as a macro.

Section 0 of the handoff says numbers are generated and never retyped.
Tables already obey that. Prose cannot, unless the numbers it quotes are
macros, which is what this emits.
"""
import csv, json
import numpy as np

OUT = 'results_pack_v2'
runs = [v for v in json.load(open('pack_runs.json')).values()
        if 'error' not in v]
cut = {r['step']: r for r in csv.DictReader(open(f'{OUT}/cutoffs.csv'))}
ver = {v['id']: v for v in json.load(open(f'{OUT}/verification.json'))}
conv = list(csv.DictReader(open(f'{OUT}/station_convergence.csv')))

def pick(step, L=1250.0, ground='superposition'):
    c = [v for v in runs if v['step'] == step and v['ground'] == ground
         and abs(v['L_md'] - L) < 1e-6]
    return c[0] if c else None

M = {}
def add(name, value, fmt='{:.1f}', note=''):
    M[name] = (fmt.format(value) if not isinstance(value, str) else value, note)

REF, L = pick('R4'), 1250.0
# --- the ladder at the reference length
for step, tag in (('R0v', 'RZeroV'), ('R0i', 'RZeroI'), ('R1', 'ROne'),
                  ('R2', 'RTwo'), ('R3', 'RThree'), ('R4', 'RFour')):
    r = pick(step)
    add(f'LadderLoss{tag}', 100 * r['loss'], '{:.1f}', f'{step} loss, %')
    add(f'LadderEp{tag}', r['e_prime'], '{:.3f}', f"{step} e', kWh/m")
    add(f'LadderRTE{tag}', r['RTE'], '{:.3f}', f'{step} eta_RTE')
    add(f'LadderTvd{tag}', r['tvd_bottom'], '{:.0f}', f'{step} TVD at the shoe, m')
# --- the steps between rungs
for a, b, tag in (('R0v', 'R1', 'FieldStep'), ('R1', 'R2', 'LayoutStep'),
                  ('R2', 'R3', 'SlotStep'), ('R3', 'R4', 'TrajectoryStep')):
    add(f'Ladder{tag}', 100 * (pick(b)['loss'] - pick(a)['loss']), '{:+.1f}',
        f'{a} -> {b}, points of loss')
# --- the reference case
add('RefLmd', L, '{:.0f}', 'reference storage length, m MD')
add('RefTvd', REF['tvd_bottom'], '{:.0f}', 'reference TVD, m')
add('RefNwells', REF['N_field'], '{:.0f}', 'Grissom slots')
add('RefTtwod', REF['T2d_C'], '{:.1f}', 'realized T_2d, C')
add('RefFieldMWh', REF['field_MWh_per_cycle'], '{:.1f}', 'field MWh per cycle')
add('RefRgMean', REF['Rg_mean'], '{:.3f}', 'mean R_g,eff, m K/W')
add('RefStations', REF['n_stations'], '{:.0f}', 'stations at the reference')
add('DsStation', REF['ds_station'], '{:.0f}', 'station spacing, m')
add('OffStationWeight', 100 * REF['off_station_weight'], '{:.1f}',
    'off-station kernel weight at the reference, %')
# --- cutoffs
c4 = cut['R4']
add('CutoffRfourLo', float(c4['cutoff_lo']), '{:.0f}', 'cutoff lower bound, m MD')
add('CutoffRfourHi', float(c4['cutoff_hi']), '{:.0f}', 'cutoff upper bound, m MD')
add('CutoffFirstHundred', float(c4['rte_first_100m']), '{:.3f}',
    'marginal RTE over the first 100 m')
add('CutoffThreshold', float(c4['threshold']), '{:.3f}', 'half of it')
for s_, tag in (('R4_th45', 'ThFortyFive'), ('R4_th70', 'ThSeventy'),
                ('R4_top1', 'TopOne')):
    if cut[s_]['cutoff_lo']:
        add(f'Cutoff{tag}Lo', float(cut[s_]['cutoff_lo']), '{:.0f}', '')
        add(f'Cutoff{tag}Hi', float(cut[s_]['cutoff_hi']), '{:.0f}', '')
add('CutoffKopThreeHundred', 'none within the sweep', None,
    'R4_kop300 has no cutoff')
# --- sensitivities at the reference length
for s_, tag in (('R4_kop300', 'KopThreeHundred'), ('R4_th45', 'ThFortyFive'),
                ('R4_th70', 'ThSeventy'), ('R4_top1', 'TopOne'),
                ('R4_top30', 'TopThirty')):
    r = pick(s_)
    add(f'Sens{tag}Loss', 100 * r['loss'], '{:.1f}', f'{s_} loss, %')
    add(f'Sens{tag}Ep', r['e_prime'], '{:.3f}', f"{s_} e'")
add('TopSpreadLoss', 100 * (pick('R4_top1')['loss'] - pick('R4_top30')['loss']),
    '{:.1f}', 'field-age spread, points of loss')
# --- casing
add('CasingGainEp',
    100 * (pick('R4')['e_prime'] / pick('R4_legacy')['e_prime'] - 1), '{:.0f}',
    "THUMS casing gain on e', %")
add('CasingGainEpRthree',
    100 * (pick('R3')['e_prime'] / pick('R3_legacy')['e_prime'] - 1), '{:.0f}',
    "same at R3, %")
add('LegacyLossRthree', 100 * pick('R3_legacy')['loss'], '{:.1f}', '')
add('LegacyLossRfour', 100 * pick('R4_legacy')['loss'], '{:.1f}', '')
# --- ground model
for s_, tag in (('R1', 'ROne'), ('R3', 'RThree'), ('R4', 'RFour')):
    add(f'FluxLoss{tag}', 100 * pick(s_, ground='flux')['loss'], '{:.1f}',
        f'{s_}, equal-flux model, %')
# --- geometry
add('SpacingMeasured', 1.65, '{:.2f}', 'measured median spacing, m')
add('SpacingNominal', 1.83, '{:.2f}', 'nominal 6 ft, m')
# --- convergence
r25 = [c for c in conv if c['step'] == 'R4' and c['L_md'] == '1250.0'
       and c['nst'] == '25'][0]
r200 = [c for c in conv if c['step'] == 'R4' and c['L_md'] == '1250.0'
        and c['nst'] == '200'][0]
add('ConvLossSpread',
    abs(100 * (float(r200['loss']) - float(r25['loss']))), '{:.2f}',
    'loss spread, 25 -> 200 stations, points')
add('ConvEpSpread',
    abs(100 * (float(r200['e_prime']) / float(r25['e_prime']) - 1)), '{:.2f}',
    "e' spread over the same, %")
add('ConvNegAtTwoHundred', float(r200['n_negative']), '{:.0f}',
    'negative strengths at 200 stations')

with open(f'{OUT}/boreas_numbers.tex', 'w') as f:
    f.write('% boreas_numbers.tex -- GENERATED. Do not edit, and do not\n'
            '% retype any of these into the text: \\input this file and use\n'
            '% the macros. Regenerate with mknumbers.py.\n')
    f.write(f'% {len(M)} macros, from pack_runs.json\n\n')
    for name, (val, note) in sorted(M.items()):
        pad = ' ' * max(1, 34 - len(name))
        f.write(f'\\newcommand{{\\{name}}}{{{val}}}{pad}'
                + (f'% {note}\n' if note else '\n'))
print(f'wrote {OUT}/boreas_numbers.tex ({len(M)} macros)')
