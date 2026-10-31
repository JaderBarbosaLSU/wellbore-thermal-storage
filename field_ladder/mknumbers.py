"""boreas_numbers.tex -- every number the prose quotes, as a macro.

Section 0 of the handoff says numbers are generated and never retyped.
Tables already obey that. Prose cannot, unless the numbers it quotes are
macros, which is what this emits.
"""
import csv, json
import numpy as np
import packstamp

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

# --- section 2b: the intensive tables -----------------------------------
TAG = {'adiabatic': 'Adiabatic', 'R4': 'RefFour', 'R3': 'RefThree',
       'R0v': 'RefZeroV'}
try:
    I = json.load(open(f'{OUT}/intensive.json'))
except FileNotFoundError:
    I = {}
for key, tag in TAG.items():
    if key not in I:
        continue
    d = I[key]
    add(f'{tag}TTwod', d['T2d_C'], '{:.1f}',
        f'{key}: mean T_2d, C')
    add(f'{tag}TThreee', d['T_3e_C'], '{:.1f}', f'{key}: T_3e, C')
    add(f'{tag}EtaORC', d['eta_ORC'], '{:.4f}', f'{key}: eta_ORC')
    add(f'{tag}EtaRTE', d['eta_RTE'], '{:.4f}', f'{key}: eta_RTE')
    add(f'{tag}EtaStorage', d['eta_storage'], '{:.4f}', '')
    add(f'{tag}CyclesToCSS', d['cycles'], '{:.0f}', '')
    add(f'{tag}EnergyConv', d['energy_conv_last5'], '{:.1e}',
        'relative spread of dE_dc over the last 5 cycles')
    add(f'{tag}WellsPerMWe', d['wells_per_MWe'], '{:.1f}',
        f'{key}: STORAGE wells per MWe, realized net basis')
    add(f'{tag}GeoPerMWe', d['N_geo_per_MWe'], '{:.1f}', '')
    add(f'{tag}AllWellsPerMWe', d['wells_total_per_MWe'], '{:.1f}', '')
    add(f'{tag}UAPerMWe', d['UA_per_MWe'], '{:.0f}',
        f'{key}: total UA, kW/K per MWe')
    add(f'{tag}Eprime', d['e_prime'], '{:.4f}', '')
    add(f'{tag}DEch', d['dE_ch_MWh'], '{:.3f}', '')
    add(f'{tag}DEdc', d['dE_dc_MWh'], '{:.3f}', '')
    add(f'{tag}Pumping', d['pumping_pct'], '{:.2f}',
        f'{key}: pumping, % of gross')
    add(f'{tag}Leakage',
        d['leakage_pct'], '{:.2f}', '')
    add(f'{tag}MdotTubeCh', d['m1_ch'], '{:.4f}', '')
    add(f'{tag}MdotTubeDc', d['m1_dc'], '{:.4f}', '')
if 'R4' in I:
    d = I['R4']
    add('RefHTHPInputPerMWe', d['HTHP_in_MW_per_MWe'], '{:.3f}',
        'HTHP electrical input at R4, MW per MWe')
    add('AdiabaticHTHPInputPerMWe', I['adiabatic']['HTHP_in_MW_per_MWe'],
        '{:.3f}', 'the same with no formation at all')
    add('RefChargeElecRatio',
        d['HTHP_in_MW_per_MWe'] / I['adiabatic']['HTHP_in_MW_per_MWe'],
        '{:.1f}', 'charged electricity per MWh delivered, R4 / adiabatic')
    for nm, tg in (('HTHP evaporator', 'hpe'), ('HTHP condenser', 'hpc'),
                   ('ORC evaporator', 'orce'), ('ORC condenser', 'orcc')):
        add(f'Duty{tg.upper()}', d[f'Q_{tg}_over_Wnet'], '{:.2f}',
            f'{nm}: Q/W_net at R4')
        add(f'UA{tg.upper()}PerMWe', d[f'UA_{tg}_per_MWe'], '{:.0f}', '')
    add('RefFieldMWePlant', 363.0 / d['wells_per_MWe'], '{:.1f}',
        'MWe from the 363 Grissom slots at the reference case')
try:
    C = json.load(open(f'{OUT}/composites.json'))
    for key, tag in (('adiabatic', 'Adiabatic'), ('R4', 'RefFour'),
                     ('R3', 'RefThree')):
        add(f'{tag}Pinch', C[key]['pinch_K'], '{:.2f}', f'{key}: ORC pinch, K')
        add(f'{tag}PinchAt', C[key]['pinch_at_pct'], '{:.1f}',
            f'{key}: pinch location, % of evaporator duty')
except FileNotFoundError:
    pass

# No macro may expand to nan, inf or nothing. \RefHTHPInputPerMWe shipped
# as 'nan' in v2.1 because the value it read was never written to
# intensive.json, and nothing caught it: LaTeX sets 'nan' happily.
# ---- 6.2: cycle indicators, deviation, casing sensitivities ----------
for key, tag in (('adiabatic', 'Adiabatic'), ('R4', 'RefFour'),
                 ('R3', 'RefThree'), ('R0v', 'RefZeroV')):
    if key in I and I[key].get('deviation') is not None:
        add(f'{tag}Deviation', I[key]['deviation'], '{:+.2f}',
            f'{key}: deviation D, %')
try:
    CY = json.load(open(f'{OUT}/cycle_states.json'))
    add('RefCOP', CY['hp_cop'], '{:.3f}', 'HTHP COP at R4')
    add('RefBleedY', CY['y_frac'], '{:.4f}', 'ORC reheat bleed fraction')
    add('RefFlashX', CY['x_8h'], '{:.4f}', 'flash vapor fraction x_8h')
    add('RefEpsRone', CY['epsilon_IHX_1'], '{:.4f}', '')
    add('RefEpsRtwo', CY['epsilon_IHX_2'], '{:.4f}', '')
except FileNotFoundError:
    pass
try:
    CS = json.load(open('casing_sens.json'))
    for g, tag in (('toc', 'Toc'), ('mud', 'Mud'), ('weight', 'Weight')):
        add(f'Sens{tag}Loss', CS[g]['max_d_loss_pt'], '{:.2f}',
            f'{g}: largest |d loss| over R3 and R4, points')
        add(f'Sens{tag}Ep', CS[g]['max_d_eprime_pct'], '{:.2f}',
            f'{g}: largest |d e-prime|, %')
    add('SensCasingMaxLoss', CS['max']['max_d_loss_pt'], '{:.2f}',
        'largest over all three unknowns, points of loss')
    add('SensCasingMaxEp', CS['max']['max_d_eprime_pct'], '{:.2f}',
        'the same for e-prime, %')
except FileNotFoundError:
    pass
try:
    CR = json.load(open('cutoff_robust.json'))
    add('CutoffRobust', CR['robust'], None,
        '400-500 m at ds = 50, 25 and 12.5 m')
except FileNotFoundError:
    pass
_rc = (cut['R4'].get('recrosses_at') or '')
if _rc:
    add('CutoffRecrossLo', float(_rc.split('-')[0]), '{:.0f}',
        'R4 marginal returns above the threshold here, beyond the shoe')
    add('CutoffRecrossHi', float(_rc.split('-')[1]), '{:.0f}', '')

BAD = [n for n, (v, _) in M.items()
       if (not str(v).strip()) or str(v).strip().lower()
       in ('nan', 'inf', '-inf', 'none', '--')]
if BAD:
    raise SystemExit('macros with no usable value: ' + ', '.join(sorted(BAD)))
# A LaTeX control sequence is letters only. A digit in the name does not
# fail here -- it fails in the paper's preamble, which is the worst place
# for it to fail. (Handoff 6.8, found in par_numbers.)
import re as _re
NONALPHA = [n for n in M if not _re.fullmatch(r'[A-Za-z]+', n)]
if NONALPHA:
    raise SystemExit('macro names must be letters only: '
                     + ', '.join(sorted(NONALPHA)))

with open(f'{OUT}/boreas_numbers.tex', 'w') as f:
    f.write(packstamp.stamp())
    f.write('% boreas_numbers.tex -- GENERATED. Do not edit, and do not\n'
            '% retype any of these into the text: \\input this file and use\n'
            '% the macros. Regenerate with mknumbers.py.\n')
    f.write(f'% {len(M)} macros, from pack_runs.json\n\n')
    for name, (val, note) in sorted(M.items()):
        pad = ' ' * max(1, 34 - len(name))
        f.write(f'\\newcommand{{\\{name}}}{{{val}}}{pad}'
                + (f'% {note}\n' if note else '\n'))
print(f'wrote {OUT}/boreas_numbers.tex ({len(M)} macros)')
