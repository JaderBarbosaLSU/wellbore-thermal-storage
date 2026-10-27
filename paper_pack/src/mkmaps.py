"""Heat-map panels as INDEPENDENT figures, for a subfigure environment.

Each panel is its own PDF sized for a single-column subfigure at about
0.48\\linewidth, so it is reproduced near 1:1 and the cell values stay
legible. Panels are no longer tiled into one wide figure.
"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
import boreas_plots as BP, boreas_sweep as S, arr
BP.use_paper_style(scale=1.0)
O = '/sessions/festive-determined-cray/mnt/outputs/panels/'
import os; os.makedirs(O, exist_ok=True)
PW, PH = 3.45, 3.05            # inches, for width=0.48\linewidth

def panel(name, rows, xk, yk, zk, xs, ys, xlabel, ylabel, **kw):
    fig, ax = plt.subplots(figsize=(PW, PH))
    _, Z = BP.heatmap(rows, xk, yk, zk, xs, ys, ax=ax, fontsize=6.4, **kw)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    fig.subplots_adjust(left=0.165, right=0.99, top=0.995, bottom=0.155)
    fig.savefig(O + name + '.pdf'); fig.savefig(O + name + '.png', dpi=190)
    plt.close(fig)
    return Z

GLIDE = r'$\Delta T_\mathrm{3C,2C}$ / K'
TM    = r'$T_\mathrm{m}$ / $^\circ$C'
LFT   = r'$L$ / ft'
NLAY  = r'$N_\mathrm{lay}$'
NFT   = r'$n_\mathrm{ft}$'
NW    = r'$N$ wells'
FA    = r'$\theta$ / deg'
GR    = r'$G$ / K km$^{-1}$'

# ---- Fig 13: cascade design map
cas = S.collect('css', arr.CASCADE, fixed=arr.ARRAY)
panel('map_cascade_a', cas, 'DT_3C_2C', 'T_m_C', 'eta_RTE', arr.DG, arr.TM,
      GLIDE, TM, cmap='viridis', cbar_label=r'$\eta_\mathrm{RTE}$',
      mark=(55.0, 150.0), best='max')
panel('map_cascade_b', cas, 'DT_3C_2C', 'T_m_C', 'kWh_per_m_well', arr.DG, arr.TM,
      GLIDE, TM, cmap='magma', cbar_label=r"$e'$ / kWh m$^{-1}$",
      mark=(55.0, 150.0), best='max')
panel('map_cascade_c', cas, 'DT_3C_2C', 'T_m_C', 'loss_frac', arr.DG, arr.TM,
      GLIDE, TM, cmap='cividis', cbar_label='loss fraction',
      mark=(55.0, 150.0), best='min')

# ---- Fig 14: geometry maps
ln = S.collect('css', arr.G_LN, fixed=arr.ARRAY_NOSEG)
ft = S.collect('css', arr.G_FT, fixed=arr.ARRAY)
panel('map_geom_a', ln, 'L_ft', 'N_lay', 'kWh_per_m_well', arr.LW, arr.NL,
      LFT, NLAY, cmap='magma', cbar_label=r"$e'$ / kWh m$^{-1}$",
      mark=(5000.0, 9), best='max')
panel('map_geom_b', ln, 'L_ft', 'N_lay', 'eta_RTE', arr.LW, arr.NL,
      LFT, NLAY, cmap='viridis', cbar_label=r'$\eta_\mathrm{RTE}$',
      mark=(5000.0, 9), best='max')
panel('map_geom_c', ft, 'L_ft', 'n_ft', 'eta_RTE', arr.LW, arr.NFT,
      LFT, NFT, cmap='viridis', cbar_label=r'$\eta_\mathrm{RTE}$',
      mark=(5000.0, 3), best='max')

# ---- Fig 16: field maps
wf = S.collect('css', arr.G_WF, fixed=arr.ARRAY_NOFAN)
gl = S.collect('css', arr.G_GL, fixed=arr.ARRAY)
panel('map_field_a', wf, 'n_wells_field', 'fan_angle_deg', 'loss_frac',
      arr.NW, arr.FA, NW, FA, cmap='cividis', cbar_label='loss fraction',
      mark=(1000.0, 0.0), best='min')
panel('map_field_b', wf, 'n_wells_field', 'fan_angle_deg', 'eta_RTE',
      arr.NW, arr.FA, NW, FA, cmap='viridis',
      cbar_label=r'$\eta_\mathrm{RTE}$', mark=(1000.0, 0.0), best='max')
panel('map_field_c', gl, 'L_ft', 'grad_K_per_km', 'loss_frac', arr.LW, arr.GR,
      LFT, GR, cmap='cividis', cbar_label='loss fraction',
      mark=(5000.0, 56.5), best='min')
print('9 panels written to', O)
