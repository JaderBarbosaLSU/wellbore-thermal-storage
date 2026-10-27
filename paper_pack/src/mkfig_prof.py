import pickle, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
import boreas as B, boreas_plots as BP
BP.use_paper_style(scale=0.95)
O = '/sessions/festive-determined-cray/mnt/outputs/'
r = pickle.load(open('rec.pkl', 'rb'))
c = B.Case(n_segments=99, n_wells_field=1000., fan_angle_deg=0.,
           packing='rows', n_islands=4)
A = BP.pcm_capacities_A(c)
dseg, dnode = BP.depth_axis(B, c, 99)
FS = (6.5, 4.6)          # single column, width=\linewidth

def save(name, **kw):
    fig, ax = plt.subplots(1, 2, figsize=FS, sharex=True)
    BP.well_profiles(r, c, axes=ax, vertical=True,
                     depth_seg=dseg, depth_node=dnode, **kw)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.925, bottom=0.145, wspace=0.08)
    fig.savefig(O + name + '.pdf'); fig.savefig(O + name + '.png', dpi=220)
    plt.close(fig)

save('BOREAS_prof_Tfluid', field='T_fluid', staircase=True, offset=-273.15,
     n_times=6, overlay='T_pcm', overlay_label=r'$T_\mathrm{pcm}$',
     ylabel=r'temperature / $^\circ$C')
save('BOREAS_prof_melt', field='A_melt', scale=1.0 / A,
     ylabel=r'melted fraction, $A_\mathrm{melt}/A_\mathrm{av}$')
save('BOREAS_prof_Ui', field='U_i',
     ylabel=r'$U_i$ / W m$^{-2}$ K$^{-1}$')
fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.9))
BP.well_integrals(r, c, axes=ax)
fig.subplots_adjust(left=0.105, right=0.895, top=0.89, bottom=0.21, wspace=0.62)
fig.savefig(O + 'BOREAS_prof_integrals.pdf'); fig.savefig(O + 'BOREAS_prof_integrals.png', dpi=220)
print('profiles re-authored at single-column width')
