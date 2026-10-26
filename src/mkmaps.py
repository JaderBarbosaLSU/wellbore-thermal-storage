import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
import boreas_plots as BP, boreas_sweep as S, arr
BP.use_paper_style()
O = '/sessions/festive-determined-cray/mnt/outputs/'

def panel(ax, rows, xk, yk, zk, xs, ys, **kw):
    _, Z = BP.heatmap(rows, xk, yk, zk, xs, ys, ax=ax, **kw)
    return Z

# ---------------- geometry: L_ft x N_lay, and n_ft x L_ft
ln = S.collect('css', arr.G_LN, fixed=arr.ARRAY_NOSEG)
ft = S.collect('css', arr.G_FT, fixed=arr.ARRAY)
fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.5))
Z1 = panel(ax[0], ln, 'L_ft', 'N_lay', 'kWh_per_m_well', arr.LW, arr.NL,
           cmap='magma', title="(a) $e'$: flat in well length",
           cbar_label=r"$e'$ / kWh per m of well", fmt='{:.3f}',
           mark=(5000.0, 9))
ax[0].set_ylabel(r'layers in the cascade, $N_\mathrm{lay}$')
Z2 = panel(ax[1], ln, 'L_ft', 'N_lay', 'eta_RTE', arr.LW, arr.NL,
           cmap='viridis', title=r'(b) $\eta_\mathrm{RTE}$',
           cbar_label=r'$\eta_\mathrm{RTE}$', fmt='{:.3f}', mark=(5000.0, 9))
for a in ax[:2]:
    a.set_xlabel(r'useful well length, $L$ / ft')
Z3 = panel(ax[2], ft, 'L_ft', 'n_ft', 'eta_RTE', arr.LW, arr.NFT,
           cmap='viridis', title=r'(c) $\eta_\mathrm{RTE}$ vs fin-tube count',
           cbar_label=r'$\eta_\mathrm{RTE}$', fmt='{:.3f}', mark=(5000.0, 3))
ax[2].set_ylabel(r'fin tubes per well, $n_\mathrm{ft}$')
ax[2].set_xlabel(r'useful well length, $L$ / ft')
fig.subplots_adjust(left=0.05, right=0.98, top=0.88, bottom=0.14, wspace=0.46)
fig.savefig(O + 'BOREAS_map_geometry.pdf'); fig.savefig(O + 'BOREAS_map_geometry.png', dpi=170)
print('geometry  e-prime spread %.1f %%   eta_RTE %.1f %%   n_ft map %.1f %%'
      % (BP.spread(Z1), BP.spread(Z2), BP.spread(Z3)))

# ---------------- field: n_wells x fan angle, gradient x L_ft
wf = S.collect('css', arr.G_WF, fixed=arr.ARRAY_NOFAN)
gl = S.collect('css', arr.G_GL, fixed=arr.ARRAY)
fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.5))
Z4 = panel(ax[0], wf, 'n_wells_field', 'fan_angle_deg', 'loss_frac',
           arr.NW, arr.FA, cmap='cividis',
           title='(a) shielding: count against divergence',
           cbar_label='formation loss fraction', fmt='{:.3f}',
           mark=(1000.0, 0.0))
ax[0].set_ylabel(r'fan half-angle, $\theta$ / deg')
ax[0].set_xlabel(r'wells in the field, $N$')
Z5 = panel(ax[1], wf, 'n_wells_field', 'fan_angle_deg', 'eta_RTE',
           arr.NW, arr.FA, cmap='viridis', title=r'(b) $\eta_\mathrm{RTE}$',
           cbar_label=r'$\eta_\mathrm{RTE}$', fmt='{:.3f}', mark=(1000.0, 0.0))
ax[1].set_xlabel(r'wells in the field, $N$')
Z6 = panel(ax[2], gl, 'L_ft', 'grad_K_per_km', 'loss_frac', arr.LW, arr.GR,
           cmap='cividis', title='(c) hotter rock against more boundary',
           cbar_label='formation loss fraction', fmt='{:.3f}',
           mark=(5000.0, 56.5))
ax[2].set_ylabel(r'geothermal gradient / K km$^{-1}$')
ax[2].set_xlabel(r'useful well length, $L$ / ft')
fig.subplots_adjust(left=0.05, right=0.98, top=0.88, bottom=0.14, wspace=0.46)
fig.savefig(O + 'BOREAS_map_field.pdf'); fig.savefig(O + 'BOREAS_map_field.png', dpi=170)
print('field     loss spread %.1f %%   eta_RTE %.1f %%   grad map %.1f %%'
      % (BP.spread(Z4), BP.spread(Z5), BP.spread(Z6)))

# ---------------- the confound figure
pin = sorted(S.collect('css', arr.G_NLPIN, fixed=arr.ARRAY_NOSEG), key=lambda r: r['N_lay'])
fm = {(r['L_ft'], r['N_lay']): r for r in ln}
free = [fm[(5000.0, n)] for n in arr.NL]
fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.2))
N = np.array(arr.NL)
for a, key, lab in ((ax[0], 'eta_RTE', r'$\eta_\mathrm{RTE}$'),
                    (ax[1], 'kWh_per_m_well', r"$e'$ / kWh per m")):
    a.plot(N, [r[key] for r in free], 'o-', color=BP.RED, lw=1.8, ms=5,
           label=r'$\Delta T_\mathrm{3C,2C}$ fixed (span varies)')
    a.plot(N, [r[key] for r in pin], 's--', color=BP.NAVY, lw=1.8, ms=5,
           label='span pinned at 48.9 K')
    a.set_xlabel(r'layers in the cascade, $N_\mathrm{lay}$')
    a.set_ylabel(lab)
    a.legend(handlelength=2.0)
ax[0].set_title(r'(a) the $N_\mathrm{lay}$ trend is the span, not the count', pad=8)
ax[1].set_title(r"(b) $e'$ converges by $N_\mathrm{lay}\approx 7$", pad=8)
fig.tight_layout()
fig.savefig(O + 'BOREAS_Nlay_confound.pdf'); fig.savefig(O + 'BOREAS_Nlay_confound.png', dpi=180)
print('confound figure written')
