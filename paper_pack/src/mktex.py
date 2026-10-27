"""Generate the LaTeX tables for the enthalpy paper from papercalc.json."""
import json
import boreas as B

d = json.load(open('papercalc.json'))
ref = d['array_rows']
rank, hp = ref['rank'], ref['hpd']
OUT = '/sessions/festive-determined-cray/mnt/outputs/paper_pack/tex/'
import os; os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------- state points
lines = [r'% Generated from boreas.py v1.3, reference case = 1000-well vertical',
         r'% array in double rows on four islands. Do not hand-edit: regenerate',
         r'% with mktex.py so the table cannot drift from the model.',
         r'\begin{table}[htbp]', r'\centering',
         r'\caption{\emph{Supplementary material.} State points of the two cycles at the reference case. '
         r'The ORC is a double-stage reheat cycle with a regenerator '
         r'(7e--10e); the HTHP is a two-stage flash-economizer cycle, in '
         r'which 7h throttles into a drum at $p_\mathrm{int}$, the saturated '
         r'liquid 9h throttles on to 4h and the flash vapour 10h rejoins the '
         r'second stage. States 5h, 6h and 11h coincide.}',
         r'\label{tab:states}', r'\small',
         r'\begin{tabular}{@{}l r r r r | l r r r r@{}}', r'\toprule',
         r'ORC & $T$ [\si{\celsius}] & $p$ [\si{\bar}] & $h$ '
         r'[\si{\kilo\joule\per\kilogram}] & $s$ '
         r'[\si{\kilo\joule\per\kilogram\per\kelvin}] & HTHP & '
         r'$T$ [\si{\celsius}] & $p$ [\si{\bar}] & $h$ & $s$ \\', r'\midrule']
for i in range(1, 14):
    e = f'{i}e'; h = f'{i}h'
    if i <= 10:
        L = (f'{e} & {rank["T_"+e]-273.15:.2f} & {rank["p_"+e]/1e5:.3f} & '
             f'{rank["h_"+e]:.2f} & {rank["s_"+e]/1e3:.4f}')
    else:
        L = ' & & & & '
    R = (f'{h} & {hp["T_"+h]-273.15:.2f} & {hp["p_"+h]/1e5:.3f} & '
         f'{hp["h_"+h]:.2f} & {hp["s_"+h]/1e3:.4f}')
    lines.append(L + ' & ' + R + r' \\')
lines += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']
open(OUT + 'tab_states.tex', 'w').write('\n'.join(lines))

# ---------------------------------------------------- cycle summary
rows = [
    (r'$\etaORC$', f'{rank["rank_eff"]:.4f}', '--'),
    (r'$\COP$', '--', f'{hp["hp_cop"]:.4f}'),
    (r'$y$, reheat mass fraction', f'{rank["y_frac"]:.4f}', '--'),
    (r'$x_\mathrm{8h}$, flash vapour fraction', '--', f'{hp["x_8h"]:.4f}'),
    (r'$\varepsilon_\mathrm{R1}$', '--', f'{hp["epsilon_IHX_1"]:.4f}'),
    (r'$\varepsilon_\mathrm{R2}$', '--', f'{hp["epsilon_IHX_2"]:.4f}'),
    (r'$p_\mathrm{evap}$ [\si{\bar}]', f'{rank["p_evape"]/1e5:.3f}',
     f'{hp["p_evaph"]/1e5:.3f}'),
    (r'$p_\mathrm{int}$ [\si{\bar}]', f'{rank["p_inte"]/1e5:.3f}',
     f'{hp["p_inth"]/1e5:.3f}'),
    (r'$p_\mathrm{cond}$ [\si{\bar}]', f'{rank["p_conde"]/1e5:.3f}',
     f'{hp["p_condh"]/1e5:.3f}'),
    (r'$T_\mathrm{3e}$ [\si{\celsius}] (pinch-limited)',
     f'{rank["T_3e"]-273.15:.2f}', '--'),
    (r'evaporator pinch [\si{\kelvin}]', f'{rank["pinch"]:.2f}', '--'),
    (r'duty [\si{\kilo\watt}], field',
     f'{ref["Q_dot_in_ORC_kW"]:.1f}', f'{ref["Q_dot_out_HP_kW"]:.1f}'),
]
t = [r'\begin{table}[htbp]', r'\centering',
     r'\caption{Cycle indicators at the reference case. The plant models are '
     r'those of the hairpin study; only $T_\mathrm{3e}$, which follows the '
     r'realized outlet $\overline{T}_\mathrm{2d}$, depends on the store. '
     r'Electrical input to the heat pump is '
     rf'\SI{{{ref["W_dot_el_in_kW"]:.1f}}}{{\kilo\watt}} and the charging '
     rf'water flow \SI{{{ref["m_dot_w_ch"]:.2f}}}{{\kilogram\per\second}}; '
     rf'per fin tube the store sees \SI{{{ref["m1_ch"]:.4f}}} and '
     rf'\SI{{{ref["m1_dc"]:.4f}}}{{\kilogram\per\second}} on charge and '
     r'discharge.}', r'\label{tab:cycles}', r'\small',
     r'\begin{tabular}{@{}l r r@{}}', r'\toprule',
     r'quantity & ORC & HTHP \\', r'\midrule']
t += [f'{a} & {b} & {c} ' + r'\\' for a, b, c in rows]
t += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']
open(OUT + 'tab_cycles.tex', 'w').write('\n'.join(t))

# ---------------------------------------------------- exchanger UA
ua = ref['UA']
order = ['ORC evaporator', 'ORC condenser', 'HTHP condenser',
         'HTHP evaporator', 'TOTAL']
t = [r'\begin{table}[htbp]', r'\centering',
     r'\caption{Conductance required by the four plant exchangers at the '
     r'reference case, integrated over the real composite curves. '
     r'$\Delta T_\mathrm{eff}=\dot Q/UA$ is the approach that reproduces the '
     r'conductance; with a phase change on one side it is the correct '
     r'generalization of the log-mean difference, which the LMTD is not. '
     r'The recovered pinches reproduce the specified approaches exactly, '
     r'which is the check on the construction.}',
     r'\label{tab:UA}', r'\small',
     r'\begin{tabular}{@{}l r r r r@{}}', r'\toprule',
     r'exchanger & $\dot Q$ [\si{\kilo\watt}] & $UA$ '
     r'[\si{\kilo\watt\per\kelvin}] & $\Delta T_\mathrm{eff}$ '
     r'[\si{\kelvin}] & pinch [\si{\kelvin}] \\', r'\midrule']
for k in order:
    v = ua[k]
    dt = '--' if v['dT_eff_K'] != v['dT_eff_K'] else f'{v["dT_eff_K"]:.2f}'
    if k == 'TOTAL':
        t.append(r'\midrule')
    t.append(f'{k.lower()} & {v["duty_kW"]:.1f} & {v["UA_kW_per_K"]:.1f} & '
             f'{dt} & {v["pinch_K"]:.2f} ' + r'\\')
t += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']
open(OUT + 'tab_UA.tex', 'w').write('\n'.join(t))

# ---------------------------------------------------- T_3e by store
t = [r'\begin{table}[htbp]', r'\centering',
     r'\caption{The store reaches the ORC through one number. Loss to the '
     r'formation lowers the realized water outlet '
     r'$\overline{T}_\mathrm{2d}$, the pinch-limited evaporating temperature '
     r'follows and $\etaORC$ with it.}',
     r'\label{tab:T3e}', r'\small',
     r'\begin{tabular}{@{}l r r r@{}}', r'\toprule',
     r'store & $\overline{T}_\mathrm{2d}$ [\si{\celsius}] & '
     r'$T_\mathrm{3e}$ [\si{\celsius}] & $\etaORC$ \\', r'\midrule']
for key, lab in (('adiabatic', r'adiabatic, $k_r\to0$'),
                 ('array_rows', 'vertical array, 1000 wells, rows'),
                 ('fan1000_rows', 'fan, 1000 wells, rows'),
                 ('isolated', 'isolated well')):
    v = d[key]
    t.append(f'{lab} & {v["T_2d_realised_C"]:.2f} & {v["T_3e_C"]:.2f} & '
             f'{v["eta_ORC"]:.4f} ' + r'\\')
t += [r'\bottomrule', r'\end{tabular}', r'\end{table}', '']
open(OUT + 'tab_T3e.tex', 'w').write('\n'.join(t))
print('wrote 4 tables to', OUT)
