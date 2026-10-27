import pickle, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
import boreas as B, boreas_plots as BP
BP.use_paper_style()
O = '/sessions/festive-determined-cray/mnt/outputs/'
r = pickle.load(open('rec.pkl', 'rb'))
c = B.Case(n_segments=99, n_wells_field=1000., fan_angle_deg=0.,
           packing='rows', n_islands=4)
A = BP.pcm_capacities_A(c)

# 1. fluid temperature
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.4), sharey=True)
BP.well_profiles(r, c, 'T_fluid', axes=ax, staircase=True, offset=-273.15,
                 ylabel=r'secondary fluid temperature, $T_f$ / $^\circ$C')
fig.subplots_adjust(left=0.07, right=0.985, top=0.92, bottom=0.14, wspace=0.07)
fig.savefig(O + 'BOREAS_prof_Tfluid.pdf'); fig.savefig(O + 'BOREAS_prof_Tfluid.png', dpi=170)
plt.close(fig)

# 2. local melt fraction
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.4), sharey=True)
BP.well_profiles(r, c, 'A_melt', axes=ax, scale=1.0 / A,
                 ylabel=r'local melted fraction, $A_\mathrm{melt}/A_\mathrm{av}$')
fig.subplots_adjust(left=0.07, right=0.985, top=0.92, bottom=0.14, wspace=0.07)
fig.savefig(O + 'BOREAS_prof_melt.pdf'); fig.savefig(O + 'BOREAS_prof_melt.png', dpi=170)
plt.close(fig)

# 3. overall heat transfer coefficient
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.4), sharey=True)
BP.well_profiles(r, c, 'U_i', axes=ax,
                 ylabel=r'overall coefficient, $U_i$ / W m$^{-2}$ K$^{-1}$')
fig.subplots_adjust(left=0.07, right=0.985, top=0.92, bottom=0.14, wspace=0.07)
fig.savefig(O + 'BOREAS_prof_Ui.pdf'); fig.savefig(O + 'BOREAS_prof_Ui.png', dpi=170)
plt.close(fig)

# 4. integrals
fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.2))
BP.well_integrals(r, c, axes=ax)
fig.subplots_adjust(left=0.07, right=0.93, top=0.90, bottom=0.15, wspace=0.42)
fig.savefig(O + 'BOREAS_prof_integrals.pdf'); fig.savefig(O + 'BOREAS_prof_integrals.png', dpi=170)
print('4 figures written')
for half in ('charge', 'discharge'):
    H = r[half]['history']
    print(f"{half:10s} U_i range {np.nanmin(H['U_i']):7.1f} - {np.nanmax(H['U_i']):7.1f} W/m2K"
          f"   melt frac end {np.mean(H['A_melt'][-1]) / A:.4f}")
