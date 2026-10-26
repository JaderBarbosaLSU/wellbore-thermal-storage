"""ORC evaporator composites for three store configurations.

The store reaches the ORC through ONE number, the realised water outlet
T_2d. Loss to the formation lowers it, the pinch-limited evaporating
temperature follows, and eta_ORC with it. Three panels on a common axis make
that chain visible in a way the table cannot.
"""
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
import boreas as B, boreas_plots as BP
BP.use_paper_style()
d = json.load(open('papercalc.json'))
O = '/sessions/festive-determined-cray/mnt/outputs/'

show = [('adiabatic', r'adiabatic, $k_r\to0$'),
        ('array_rows', 'vertical array, rows'),
        ('isolated', 'isolated well')]
fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.4), sharey=True)
case = B.Case(n_segments=99)
for ax, (key, lab) in zip(axes, show):
    v = d[key]
    rank, T = v['rank'], v['Tdict']
    Qc, Tc = B._orc_cold_composite(rank, case.fluid)
    Qh, Th = B._water_curve(T['T_3d'], T['T_2d'], case.P, Qc[-1])
    pc = 100 * np.asarray(Qc) / Qc[-1]
    ph = 100 * np.asarray(Qh) / Qc[-1]
    q = np.linspace(0, 100, 3000)
    Tw = np.interp(q, ph, np.asarray(Th) - 273.15)
    Tf = np.interp(q, pc, np.asarray(Tc) - 273.15)
    i = int(np.argmin(Tw - Tf))
    ax.plot(ph, np.asarray(Th) - 273.15, color=BP.RED, lw=1.9,
            label=r'water, $\overline{T}_\mathrm{2d}\to T_\mathrm{3d}$')
    ax.plot(pc, np.asarray(Tc) - 273.15, color=BP.NAVY, lw=1.9,
            label='ORC working fluid')
    ax.annotate('', xy=(q[i], Tw[i]), xytext=(q[i], Tf[i]),
                arrowprops=dict(arrowstyle='<->', color='#1a1a1a', lw=1.1))
    ax.annotate(rf'${Tw[i] - Tf[i]:.2f}$ K', (q[i], 0.5 * (Tw[i] + Tf[i])),
                textcoords='offset points', xytext=(7, -3), fontsize=8.6)
    ax.axhline(rank['T_3e'] - 273.15, color='0.6', lw=0.8, ls=(0, (1, 2.5)))
    ax.set_title(rf'{lab}' '\n'
                 rf'$\overline{{T}}_\mathrm{{2d}}={T["T_2d"] - 273.15:.1f}$, '
                 rf'$T_\mathrm{{3e}}={rank["T_3e"] - 273.15:.2f}\,^\circ$C, '
                 rf'$\eta_\mathrm{{ORC}}={rank["rank_eff"]:.4f}$', pad=8,
                 fontsize=9.8)
    ax.set_xlabel('cumulative duty  /  % of evaporator load')
axes[0].set_ylabel(r'temperature, $T$  /  $^{\circ}$C')
axes[0].legend(loc='lower right', handlelength=1.8)
fig.subplots_adjust(left=0.055, right=0.985, top=0.80, bottom=0.15, wspace=0.08)
fig.savefig(O + 'BOREAS_ORC_composite_stores.pdf')
fig.savefig(O + 'BOREAS_ORC_composite_stores.png', dpi=170)
print('written')
for key, _ in show:
    v = d[key]
    print(f"{key:14s} T_2d {v['T_2d_realised_C']:7.2f}  T_3e {v['T_3e_C']:7.2f}  "
          f"etaORC {v['eta_ORC']:.4f}  binding {v['T_3e_binding']}")
