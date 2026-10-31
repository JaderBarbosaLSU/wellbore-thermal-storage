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

    # Temperatures in DEGREES CELSIUS, depth labelled TVD (it reaches
    # 891 m at 1250 m MD, so a reader must not take it for the storage
    # length), and the time key as a colorbar on panel (b) only.
    two_panel('BOREAS_prof_Tfluid', lambda ax: BP.well_profiles(
        rec, CASE, field='T_fluid', axes=ax, vertical=True,
        depth_seg=dseg, depth_node=dnode, overlay='T_pcm',
        overlay_label=r'$T_\mathrm{pcm}$', offset=-273.15,
        ylabel=r'temperature / $^\circ$C', colorbar=True))
    two_panel('BOREAS_prof_melt', lambda ax: BP.well_profiles(
        rec, CASE, field='A_melt', axes=ax, vertical=True,
        depth_seg=dseg, depth_node=dnode,
        scale=1.0 / BP.pcm_capacities_A(CASE),
        ylabel='melted fraction of the PCM', colorbar=True))
    two_panel('BOREAS_prof_Ui', lambda ax: BP.well_profiles(
        rec, CASE, field='U_i', axes=ax, vertical=True,
        depth_seg=dseg, depth_node=dnode,
        ylabel=r'$U_i$ / W m$^{-2}$ K$^{-1}$', colorbar=True))
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

    # NOT the cycled fraction. eps is the RESIDUAL melted fraction at the
    # end of each half-cycle; the fraction actually cycled is the
    # difference, about 0.58 here. A small eps_dc means the store refroze
    # almost completely, which is the intended behaviour, not a failure.
    e_ch = float(rec['charge']['eps_local'].mean())
    e_dc = float(rec['discharge']['eps_local'].mean())
    add('EpsCh', e_ch, '{:.4f}', 'melted fraction at the END OF CHARGE')
    add('EpsDc', e_dc, '{:.4f}',
        'residual melted fraction at the END OF DISCHARGE')
    add('EpsCycled', e_ch - e_dc, '{:.4f}',
        'the fraction actually cycled, eps_ch - eps_dc')

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

    # RATES PER FIN TUBE. q_prime is already per unit length of fin tube
    # (Eq. ode is written per tube, K_j is per tube length), so the
    # per-tube rate is the bare sum -- v2.3 divided by n_ft on top of
    # that and came out a factor of three low.
    #
    # The MEAN is step-weighted. The record grid is logarithmic, 41
    # points whose steps run from 0.3 s to 8491 s, so an arithmetic mean
    # over it weights the first minute as heavily as the last two hours.
    # For the discharge, which collapses to near zero in the last
    # quarter, that is the difference between -125 and -62 kW.
    trapz = getattr(np, 'trapezoid', None) or np.trapz
    for half, tag in (('charge', 'Ch'), ('discharge', 'Dc')):
        H = rec[half]['history']
        z = rec[half]['z']
        dz = z[1] - z[0]
        th = np.asarray(rec[half]['t_hist'], float)
        q = np.asarray(H['q_prime'], float).sum(axis=1) * dz / 1e3
        add(f'Rate{tag}Start', float(q[0]), '{:.1f}',
            f'{half} rate per fin tube at the start, kW')
        add(f'Rate{tag}End', float(q[-1]), '{:.1f}', 'at the end')
        add(f'Rate{tag}Mean', float(trapz(q, th) / (th[-1] - th[0])),
            '{:.1f}', 'half-cycle average, STEP-WEIGHTED')
        add(f'Rate{tag}MeanPerWell',
            float(CASE.n_ft * trapz(q, th) / (th[-1] - th[0])), '{:.0f}',
            'the same per well')

    # ---- how much of the discharge arrives below the ORC's reach ------
    # The plant model evaluates the ORC at the flow-averaged T_2d with a
    # fixed evaporating temperature. Late in the discharge the outlet
    # falls below T_3e + DT_pinch, and an evaporator at that pressure
    # cannot take the heat. This is the share of dE_dc delivered there.
    rank = rec['rank']
    T_lim = rank['T_3e'] - 273.15 + CASE.DT_pinch_ORC
    add('TlimORC', T_lim, '{:.1f}',
        'T_3e + DT_pinch_ORC, the lowest outlet the ORC can use, C')
    H = rec['discharge']['history']
    z = rec['discharge']['z']
    th = np.asarray(rec['discharge']['t_hist'], float)
    q = np.abs(np.asarray(H['q_prime'], float).sum(axis=1) * (z[1] - z[0]))
    # T_out is recorded once per step, so it is one shorter than t_hist;
    # pad with the initial outlet so the mask lines up with the rate.
    Tout = np.asarray(rec['discharge']['T_out'], float) - 273.15
    if len(Tout) == len(th) - 1:
        Tout = np.r_[Tout[0], Tout]
    below = Tout < T_lim
    tot = trapz(q, th)
    share = trapz(np.where(below, q, 0.0), th) / tot if tot else float('nan')
    add('BelowPinchShare', 100 * share, '{:.1f}',
        'share of dE_dc delivered with T_2d(t) below that limit, %')
    add('ToutDcEnd', float(Tout[-1]), '{:.1f}',
        'discharge outlet at the end, C')

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
