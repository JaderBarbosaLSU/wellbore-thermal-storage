import pickle, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
import boreas as B, boreas_plots as BP
BP.use_paper_style()
O = '/sessions/festive-determined-cray/mnt/outputs/'
r = pickle.load(open('rec.pkl', 'rb'))
c = B.Case(n_segments=99, n_wells_field=1000., fan_angle_deg=0.,
           packing='rows', n_islands=4)
A = BP.pcm_capacities_A(c)
dseg, dnode = BP.depth_axis(B, c, 99)
print(f'depth: seg[0]={dseg[0]:.1f} seg[-1]={dseg[-1]:.1f}  node[0]={dnode[0]:.1f} node[-1]={dnode[-1]:.1f}')

def save(name, **kw):
    fig, ax = plt.subplots(1, 2, figsize=(8.8, 5.8), sharex=True)
    BP.well_profiles(r, c, axes=ax, vertical=True,
                     depth_seg=dseg, depth_node=dnode, **kw)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.93, bottom=0.12, wspace=0.08)
    fig.savefig(O + name + '.pdf'); fig.savefig(O + name + '.png', dpi=170)
    plt.close(fig)

save('BOREAS_prof_Tfluid', field='T_fluid', staircase=True, offset=-273.15,
     n_times=6, overlay='T_pcm', overlay_label=r'$T_\mathrm{pcm}$',
     ylabel=r'temperature / $^\circ$C')
save('BOREAS_prof_melt', field='A_melt', scale=1.0 / A,
     ylabel=r'local melted fraction, $A_\mathrm{melt}/A_\mathrm{av}$')
save('BOREAS_prof_Ui', field='U_i',
     ylabel=r'overall coefficient, $U_i$ / W m$^{-2}$ K$^{-1}$')

fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.2))
BP.well_integrals(r, c, axes=ax)
fig.subplots_adjust(left=0.07, right=0.93, top=0.90, bottom=0.15, wspace=0.42)
fig.savefig(O + 'BOREAS_prof_integrals.pdf'); fig.savefig(O + 'BOREAS_prof_integrals.png', dpi=170)
print('4 figures written (true depth, T_pcm overlaid)')
