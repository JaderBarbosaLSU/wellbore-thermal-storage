"""ORC evaporator composites, one INDEPENDENT panel per store.

Three panels for a subfigure environment, sized for 0.48\\linewidth in a
single-column layout. The temperature range is computed across all three
and applied to each, so the panels stay directly comparable even though
they are now separate figures -- that comparability is the whole point of
the figure and a per-panel autoscale would destroy it.
"""
import json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
import boreas as B, boreas_plots as BP
BP.use_paper_style(scale=1.0)
d = json.load(open('/sessions/festive-determined-cray/mnt/outputs/paper_pack/data/papercalc.json'))
O = '/sessions/festive-determined-cray/mnt/outputs/panels/'
case = B.Case(n_segments=99)
SHOW = [('adiabatic', 'comp_a', r'adiabatic, $k_r\to0$'),
        ('array_rows', 'comp_b', 'vertical array, rows'),
        ('isolated', 'comp_c', 'isolated well')]

def curves(key):
    v = d[key]; rank, T = v['rank'], v['Tdict']
    Qc, Tc = B._orc_cold_composite(rank, case.fluid)
    Qh, Th = B._water_curve(T['T_3d'], T['T_2d'], case.P, Qc[-1])
    return (100 * np.asarray(Qc) / Qc[-1], np.asarray(Tc) - 273.15,
            100 * np.asarray(Qh) / Qc[-1], np.asarray(Th) - 273.15, rank, T)

data = {k: curves(k) for k, _, _ in SHOW}
lo = min(min(v[1].min(), v[3].min()) for v in data.values())
hi = max(max(v[1].max(), v[3].max()) for v in data.values())
pad = 0.06 * (hi - lo)

for key, name, lab in SHOW:
    pc, Tc, ph, Th, rank, T = data[key]
    fig, ax = plt.subplots(figsize=(3.45, 2.95))
    q = np.linspace(0, 100, 3000)
    Tw = np.interp(q, ph, Th); Tf = np.interp(q, pc, Tc)
    i = int(np.argmin(Tw - Tf))
    ax.plot(ph, Th, color=BP.RED, lw=1.8,
            label=r'water, $\overline{T}_\mathrm{2d}\to T_\mathrm{3d}$')
    ax.plot(pc, Tc, color=BP.NAVY, lw=1.8, label='ORC working fluid')
    ax.annotate('', xy=(q[i], Tw[i]), xytext=(q[i], Tf[i]),
                arrowprops=dict(arrowstyle='<->', color='#1a1a1a', lw=1.1))
    ax.annotate(r'$\Delta T_\mathrm{pinch,ORC}$',
                (q[i], 0.5 * (Tw[i] + Tf[i])), textcoords='offset points',
                xytext=(7, -4), fontsize=7.8, bbox=BP.MASK, zorder=9)
    ax.axhline(rank['T_3e'] - 273.15, color='0.6', lw=0.8, ls=(0, (1, 2.5)))
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel('cumulative duty / % of load')
    ax.set_ylabel(r'temperature / $^\circ$C')
    if key == 'adiabatic':
        ax.legend(loc='lower right', handlelength=1.6, fontsize=7.4)
    fig.subplots_adjust(left=0.175, right=0.98, top=0.97, bottom=0.185)
    fig.savefig(O + name + '.pdf'); fig.savefig(O + name + '.png', dpi=190)
    plt.close(fig)
    print(f"{lab:26s} T_2d {T['T_2d']-273.15:7.2f}  T_3e {rank['T_3e']-273.15:7.2f}"
          f"  etaORC {rank['rank_eff']:.4f}")
print('3 composite panels written, common T range %.1f..%.1f C' % (lo, hi))
