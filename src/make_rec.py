"""Produce rec.pkl: one recorded cycle at CSS for the reference case.

The profile figures and the animation are built from a RECORDED march,
which `simulate_css_corrected(..., record=True)` produces at the cost of
one extra cycle. It is cached to disk because it takes about a minute and
three separate scripts consume it.
"""
import pickle
import boreas as B

CASE = dict(n_segments=99, n_wells_field=1000., fan_angle_deg=0.,
            packing='rows', n_islands=4)

if __name__ == '__main__':
    c = B.Case(**CASE)
    r = B.simulate_css_corrected(c, record=True)
    with open('rec.pkl', 'wb') as f:
        pickle.dump(r, f)
    print('rec.pkl written:', r['charge']['history']['T_fluid'].shape,
          'time levels x nodes;  eta_storage %.4f' % r['eta_storage'])
