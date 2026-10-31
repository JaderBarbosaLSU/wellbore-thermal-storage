"""Section 4.3: the profile figures at R4, one PDF per panel, and the
numbers the text quotes (handoff 6.6 and 6.7).

The march runs INSIDE field_ground, so every profile is on the measured
Grissom field: the superposition resistance R_g,eff(MD) and the rock
temperature on TRUE VERTICAL DEPTH. The S-curve reaches 891 m TVD at
1250 m MD, so running this outside the context manager puts the bottom of
the store in rock about 15 K too hot.
"""
import json
import os
import pickle

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import boreas as B
import boreas_ground as G
import boreas_plots as BP
import packstamp

OUT = 'results_pack_v2'
CSV = 'THUMS_wells_WellSTAR_clean.csv'
L_MD = 1250.0
THUMS = dict(D_in=201.2 / 25.4, r_casing_o=0.21910 / 2,
             r_bore=0.31120 / 2, k_cement=0.7)
CASE = B.Case(n_segments=99, L_ft=L_MD / 0.3048, t_operation_yr=10.0,
              **THUMS)
PANEL = (4.9, 4.4)


def march():
    if os.path.exists('rec_field.pkl'):
        return pickle.load(open('rec_field.pkl', 'rb'))['rec']
    md, Rg, tvd, _ = G.cached_reference_field(
        B, CASE, CSV, trajectory='s_curve', kop_rule='staggered',
        theta_max=55.0, kop_ft=300.0, L_md=L_MD, ds=G.DS_STATION,
        verbose=False)
    with G.field_ground(B, md, Rg, tvd):
        rec = B.simulate_css_corrected(CASE, record=True)
    pickle.dump(dict(rec=rec), open('rec_field.pkl', 'wb'))
    return rec


def depth_axes():
    """Depth from the TRAJECTORY, not from the solver coordinate.

    The solver's z = 0 is the BOTTOM of the well on charge, and on an
    S-curve 1250 m of hole reaches only 891 m of depth, so neither the
    raw coordinate nor measured depth is the right vertical axis.
    """
    md, Rg, tvd, _ = G.cached_reference_field(
        B, CASE, CSV, trajectory='s_curve', kop_rule='staggered',
        theta_max=55.0, kop_ft=300.0, L_md=L_MD, ds=G.DS_STATION,
        verbose=False)
    n = CASE.n_segments
    xi = (np.arange(n) + 0.5) / n
    mdv = CASE.L_well * (1.0 - xi)
    dseg = np.interp(mdv, np.r_[0.0, md], np.r_[0.0, tvd])
    dnode = np.interp(CASE.L_well * (1.0 - np.linspace(0, 1, n + 1)),
                      np.r_[0.0, md], np.r_[0.0, tvd])
    return dseg, dnode


def two_panel(name, draw_pair, size=PANEL):
    """Render a charge/discharge pair as two separate files.

    The two axes live on different figures, so each panel is drawn
    natively rather than cropped out of a two-panel image afterwards.
    """
    fa, axa = plt.subplots(figsize=size)
    fb, axb = plt.subplots(figsize=size)
    draw_pair([axa, axb])
    for ax in (axa, axb):
        ax.set_title('')
    if axb.get_legend() is not None:
        axb.get_legend().remove()
    for f, tag in ((fa, 'a'), (fb, 'b')):
        f.tight_layout()
        f.savefig(f'{OUT}/{name}_{tag}.pdf')
        plt.close(f)
    print(f'  wrote {OUT}/{name}_a.pdf and {name}_b.pdf')


def main():
    os.makedirs(OUT, exist_ok=True)
    BP.use_paper_style()
    rec = march()
    dseg, dnode = depth_axes()
    print(f'depth runs {dseg[0]:.0f} m (bottom) -> {dseg[-1]:.0f} m (top)')

    two_panel('BOREAS_prof_Tfluid', lambda ax: BP.well_profiles(
        rec, CASE, field='T_fluid', axes=ax, vertical=True,
        depth_seg=dseg, depth_node=dnode, overlay='T_pcm',
        overlay_label=r'$T_\mathrm{pcm}$'))
    two_panel('BOREAS_prof_melt', lambda ax: BP.well_profiles(
        rec, CASE, field='A_melt', axes=ax, vertical=True,
        depth_seg=dseg, depth_node=dnode,
        scale=1.0 / BP.pcm_capacities_A(CASE),
        ylabel='melted fraction'))
    two_panel('BOREAS_prof_Ui', lambda ax: BP.well_profiles(
        rec, CASE, field='U_i', axes=ax, vertical=True,
        depth_seg=dseg, depth_node=dnode,
        ylabel=r'$U_i$ / W m$^{-2}$ K$^{-1}$'))
    two_panel('BOREAS_prof_integrals', lambda ax: BP.well_integrals(
        rec, CASE, axes=ax, per_tube=True, titles=False), size=(4.9, 3.9))

    # ---- the numbers section 4.3 quotes -------------------------------
    M = {}

    def add(n, v, fmt='{:.2f}', note=''):
        M[n] = (fmt.format(v), note)

    Tm = rec['T_m_seg'] - 273.15
    out_ch = np.asarray(rec['charge']['T_out'], float) - 273.15
    out_dc = np.asarray(rec['discharge']['T_out'], float) - 273.15
    add('TtwocLo', out_ch.min(), '{:.1f}', 'T_2c range over the charge, C')
    add('TtwocHi', out_ch.max(), '{:.1f}', '')
    add('TtwodLo', out_dc.min(), '{:.1f}', 'T_2d range over the discharge, C')
    add('TtwodHi', out_dc.max(), '{:.1f}', '')

    Tp = np.asarray(rec['discharge']['history']['T_pcm'], float)[-1] - 273.15
    sub = Tm - Tp                      # positive = PCM below its melting point
    add('SubcoolMax', float(sub.max()), '{:.1f}',
        'largest subcooling at the end of discharge, K')
    add('SubcoolMean', float(sub.mean()), '{:.1f}', '')

    add('EpsCh', float(rec['charge']['eps_local'].mean()), '{:.4f}',
        'mean cycled fraction, charge')
    add('EpsDc', float(rec['discharge']['eps_local'].mean()), '{:.4f}', '')

    U = {h: np.asarray(rec[h]['history']['U_i'], float) for h in
         ('charge', 'discharge')}
    allU = np.concatenate([u.ravel() for u in U.values()])
    add('UiLo', float(allU.min()), '{:.1f}', 'U_i range, W/m2/K')
    add('UiHi', float(allU.max()), '{:.1f}', '')
    # mean against the value at the surface end (last segment, index -1)
    add('UiMean', float(allU.mean()), '{:.1f}', '')
    add('UiSurface', float(U['charge'][-1, -1]), '{:.1f}',
        'U_i at the surface end, end of charge')
    add('UiMeanToSurface', float(allU.mean() / U['charge'][-1, -1]),
        '{:.2f}', 'mean over surface-end value')

    for half, tag in (('charge', 'Ch'), ('discharge', 'Dc')):
        H = rec[half]['history']
        z = rec[half]['z']
        dz = z[1] - z[0]
        q = np.asarray(H['q_prime'], float).sum(axis=1) * dz / 1e3 / CASE.n_ft
        add(f'Rate{tag}Start', float(q[0]), '{:.1f}',
            f'{half} rate per fin tube at the start, kW')
        add(f'Rate{tag}End', float(q[-1]), '{:.1f}', 'at the end')
        add(f'Rate{tag}Mean', float(q.mean()), '{:.1f}',
            'half-cycle average')

    bad = [n for n, (v, _) in M.items() if v.strip().lower()
           in ('nan', 'inf', '-inf', '')]
    if bad:
        raise SystemExit('unusable: ' + ', '.join(bad))
    import re
    nonalpha = [n for n in M if not re.fullmatch(r'[A-Za-z]+', n)]
    if nonalpha:
        raise SystemExit('names must be letters only: ' + ', '.join(nonalpha))

    with open(f'{OUT}/boreas_profile_numbers.tex', 'w') as f:
        f.write(packstamp.stamp())
        f.write('% GENERATED by mkprofiles.py -- section 4.3 numbers, at R4.\n'
                f'% {len(M)} macros\n\n')
        for n, (v, note) in sorted(M.items()):
            pad = ' ' * max(1, 30 - len(n))
            f.write(f'\\newcommand{{\\Prof{n}}}{{{v}}}{pad}'
                    + (f'% {note}\n' if note else '\n'))
    json.dump({n: v for n, (v, _) in M.items()},
              open(f'{OUT}/profile_numbers.json', 'w'), indent=1)
    print(f'  wrote {OUT}/boreas_profile_numbers.tex ({len(M)} macros)')
    for n, (v, note) in sorted(M.items()):
        print(f'    {n:18s}{v:>10s}  {note}')


if __name__ == '__main__':
    main()
