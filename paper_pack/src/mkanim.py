import pickle, matplotlib; matplotlib.use('Agg')
import boreas as B, boreas_plots as BP
BP.use_paper_style()
r = pickle.load(open('rec.pkl', 'rb'))
c = B.Case(n_segments=99, n_wells_field=1000., fan_angle_deg=0.,
           packing='rows', n_islands=4)
O = '/sessions/festive-determined-cray/mnt/outputs/'
print('mp4 frames', BP.animate_well(r, c, O + 'BOREAS_cycle.mp4', fps=12,
                                    frames_per_half=60))
print('gif frames', BP.animate_well(r, c, O + 'BOREAS_cycle.gif', fps=10,
                                    frames_per_half=45, dpi=85))
