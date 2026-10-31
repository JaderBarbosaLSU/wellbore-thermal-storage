"""Results pack v2.0 -- CSV, profiles and the four §6 figures.

Reads pack_runs.json (produced by pack_runner.py) and emits exactly what
the handoff's §6 asks for. Nothing here computes physics: if a number is
wrong it is wrong in the runner, and this file will not hide it.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

import boreas as B
import boreas_ground as G

OUT = 'results_pack_v2'
CSV_IN = 'THUMS_wells_WellSTAR_clean.csv'
COLS = ['step', 'ground', 'L_md', 'tvd_bottom', 'loss', 'eta_storage', 'RTE',
        'e_prime', 'E_out_kWh', 'E_in_kWh', 'T2d_C', 'delivers', 'N_field',
        'field_MWh_per_cycle', 'kop_rule', 'theta_max', 't_op_yr',
        'casing', 'positions', 'trajectory', 'R_near', 'Rg_mean', 'q_spread',
        'n_stations', 'ds_station', 'off_station_weight', 'n_negative',
        'beyond_shoe']

# TrueType embedded, serif, per §6
rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42,
                 'font.family': 'serif',
                 'font.serif': ['Times New Roman', 'Nimbus Roman',
                                'Liberation Serif', 'DejaVu Serif'],
                 'mathtext.fontset': 'stix',
                 'axes.spines.top': False, 'axes.spines.right': False,
                 'savefig.bbox': 'tight', 'savefig.dpi': 300})

RED, NAVY, GREEN, GREY = '#b03028', '#23356b', '#2f6f4f', '#6a6a6a'
LADDER = ['R0v', 'R0i', 'R1', 'R2', 'R3', 'R4']
LABEL = {'R0v': 'isolated,\nvertical', 'R0i': 'isolated,\nS-curve',
         'R1': 'disc,\nvertical', 'R2': 'double row,\nvertical',
         'R3': 'Grissom slots,\nvertical',
         'R4': 'Grissom slots,\nS-curve (reference)'}
L_SHOE = 1250.0


def load():
    return list(json.load(open('pack_runs.json')).values())


def write_csv(rows):
    import csv
    os.makedirs(OUT, exist_ok=True)
    rows = sorted(rows, key=lambda r: (r['step'], r['ground'], r['L_md']))
    with open(f'{OUT}/ladder_results.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return len(rows)


def pick(rows, step, L=None, ground='superposition'):
    c = [r for r in rows if r['step'] == step and r['ground'] == ground
         and 'error' not in r]
    if L is not None:
        c = [r for r in c if abs(r['L_md'] - L) < 1e-6]
    return sorted(c, key=lambda r: r['L_md'])


# ---------------------------------------------------- one panel per file
# Handoff 6.7: the paper sets each panel as its own subfigure and supplies
# the (a), (b) labels itself. So nothing below draws a panel letter or a
# suptitle into the image, panels that sit side by side are the same size,
# and only the first of a row carries a legend.
PANEL = (4.9, 3.9)


def panel(name, draw, size=PANEL):
    fig, ax = plt.subplots(figsize=size)
    draw(ax)
    fig.tight_layout()
    fig.savefig(f'{OUT}/{name}')
    plt.close(fig)
    print(f'  wrote {OUT}/{name}')


# ------------------------------------------------------------- geometry
def fig_geometry():
    P = G.grissom_positions(CSV_IN)
    kop, is_out = G.staggered_kop(P)
    md, X = G.trajectories(P, 2000.0, kind='s_curve', kop=kop)

    def a(ax):
        for grp, col, lab in ((is_out, RED, 'outer row'),
                              (~is_out, NAVY, 'inner row')):
            ax.plot(P[grp, 0], P[grp, 1], 'o', ms=1.8, color=col, label=lab)
        ax.set_aspect('equal')
        ax.set_xlabel('east / m'); ax.set_ylabel('north / m')
        ax.legend(frameon=False, fontsize=7.5, markerscale=3)

    def b(ax):
        rad = np.hypot(X[:, :, 0], X[:, :, 1])
        for i in range(0, len(P), 7):
            ax.plot(rad[i], X[i, :, 2], lw=0.4, color=GREY, alpha=0.7)
        ax.invert_yaxis()
        ax.set_xlabel('radius from centroid / m'); ax.set_ylabel('TVD / m')

    def c(ax):
        from scipy.spatial import cKDTree
        for kind, col, lab in (('vertical', NAVY, 'vertical'),
                               ('s_curve', RED, 'S-curve')):
            m2, X2 = G.trajectories(P, 2000.0, kind=kind, kop=kop)
            nn = [np.median(cKDTree(X2[:, j, :2]).query(X2[:, j, :2], k=2)[0][:, 1])
                  for j in range(0, len(m2), 8)]
            ax.plot(m2[::8], nn, color=col, lw=1.6, label=lab)
        ax.set_yscale('log')
        ax.set_xlabel('MD / m'); ax.set_ylabel('median spacing / m')
        ax.legend(frameon=False, fontsize=8)

    def d(ax):
        r = np.hypot(P[:, 0] - P[:, 0].mean(), P[:, 1] - P[:, 1].mean())
        ax.hist(r / r.max(), bins=20, color=NAVY, alpha=0.85)
        ax.set_xlabel(r'radius / $R_\mathrm{max}$')
        ax.set_ylabel('slots')

    for tag, fn in (('a', a), ('b', b), ('c', c), ('d', d)):
        panel(f'BOREAS_grissom_geometry_{tag}.pdf', fn)


# ------------------------------------------------------------- shielding
def fig_Rg(profiles):
    def a(ax):
        for step, col in (('R3', NAVY), ('R4', RED)):
            if step in profiles:
                md, rg = profiles[step]
                ax.plot(rg, md, color=col, lw=1.8,
                        label=f'{step}: ' + LABEL[step].replace('\n', ' '))
        ax.invert_yaxis(); ax.set_xscale('log')
        ax.set_xlabel(r'$R_{g,\mathrm{eff}}$  /  m K W$^{-1}$')
        ax.set_ylabel('MD / m')
        ax.legend(frameon=False, fontsize=8)

    def b(ax):
        for step, col in (('R4_kop300', '#0a9396'), ('R4_th45', GREEN),
                          ('R4_th70', RED)):
            if step in profiles:
                md, rg = profiles[step]
                ax.plot(rg, md, color=col, lw=1.6, label=step)
        ax.invert_yaxis(); ax.set_xscale('log')
        ax.set_xlabel(r'$R_{g,\mathrm{eff}}$  /  m K W$^{-1}$')
        ax.set_ylabel('MD / m')
        ax.legend(frameon=False, fontsize=8)

    panel('BOREAS_grissom_Rg_a.pdf', a)
    panel('BOREAS_grissom_Rg_b.pdf', b)


# ---------------------------------------------------------------- ladder
def fig_ladder(rows, L=1250.0):
    steps = [s for s in LADDER if pick(rows, s, L)]
    loss = [100 * pick(rows, s, L)[0]['loss'] for s in steps]
    ep = [pick(rows, s, L)[0]['e_prime'] for s in steps]
    lab = [LABEL[s] for s in steps]
    x = np.arange(len(steps))

    def mk(vals, ylab, colors, fmt):
        def draw(ax):
            ax.bar(x, vals, color=colors)
            ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=7)
            ax.set_ylabel(ylab)
            for i, v in enumerate(vals):
                ax.text(i, v + 0.02 * max(vals), fmt.format(v), ha='center',
                        fontsize=7.5)
        return draw

    panel('BOREAS_grissom_ladder_a.pdf',
          mk(loss, 'formation loss / %',
             [RED if v > 50 else NAVY for v in loss], '{:.1f}'), (5.4, 4.0))
    panel('BOREAS_grissom_ladder_b.pdf',
          mk(ep, r"$e'$ / kWh per m of well", GREEN, '{:.3f}'), (5.4, 4.0))


# ----------------------------------------------------------------- sweep
def fig_sweep(rows):
    cols = dict(zip(LADDER, (GREY, '#9c27a0', NAVY, '#0a9396', GREEN, RED)))
    series = {}
    for s in LADDER:
        rr = pick(rows, s)
        if len(rr) < 2:
            continue
        L = np.array([r['L_md'] for r in rr])
        Eo = np.array([r['E_out_kWh'] for r in rr])
        Ei = np.array([r['E_in_kWh'] for r in rr])
        series[s] = (L, np.array([r['e_prime'] for r in rr]),
                     np.array([r['RTE'] for r in rr]),
                     0.5 * (L[1:] + L[:-1]), np.diff(Eo) / np.diff(Ei),
                     Eo[0] / Ei[0] if abs(L[0] - 100.) < 1e-6 else np.nan)

    def shade(ax):
        ax.axvspan(L_SHOE, 2000.0, color='#cccccc', alpha=.35, lw=0, zorder=0)
        ax.axvline(L_SHOE, color=GREY, lw=0.9, ls=(0, (2, 2)))
        ax.set_xlabel('storage length $L$ / m MD')

    def a(ax):
        for s, v in series.items():
            ax.plot(v[0], v[1], 'o-', ms=3, lw=1.4, color=cols[s],
                    label=f'{s}: ' + LABEL[s].replace('\n', ' '))
        shade(ax); ax.set_ylabel(r"$e'$ / kWh m$^{-1}$")
        ax.legend(frameon=False, fontsize=6.6, ncol=1)

    def b(ax):
        for s, v in series.items():
            ax.plot(v[0], v[2], 'o-', ms=3, lw=1.4, color=cols[s])
        shade(ax); ax.set_ylabel(r'$\eta_\mathrm{RTE}$')

    def c(ax):
        for s, v in series.items():
            ax.plot(v[3], v[4], 'o-', ms=3, lw=1.4, color=cols[s])
            if np.isfinite(v[5]):
                ax.plot([0.5 * v[0][0]], [v[5]], 'o', ms=3, color=cols[s])
        shade(ax); ax.set_ylim(bottom=0)
        ax.set_ylabel(r'marginal $\Delta E_\mathrm{out}/\Delta E_\mathrm{in}$')

    for tag, fn in (('a', a), ('b', b), ('c', c)):
        panel(f'BOREAS_grissom_sweep_{tag}.pdf', fn, (5.2, 4.0))


def write_profiles(L=1250.0):
    """R_g,eff(MD) and TVD(MD) for EVERY rung and every sensitivity.

    v2.0 shipped two of them. The sensitivities are the rows a reader is
    most likely to want to re-plot, and they cannot be reconstructed from
    the scalar CSV.
    """
    THUMS = dict(D_in=201.2 / 25.4, r_casing_o=0.21910 / 2,
                 r_bore=0.31120 / 2, k_cement=0.7)
    LEGACY = dict(D_in=7.0, r_casing_o=0.098, r_bore=0.110, k_cement=1.0)
    P = G.grissom_positions(CSV_IN)
    kop, _ = G.staggered_kop(P)
    N = len(P)
    DISC = G.disc_positions(N, B0=G.B0_MEASURED)
    ROWS = G.rows_positions(N, B0=G.B0_MEASURED, w=G.W_DEF)
    SPEC = [
        ('R1', DISC, 'vertical', None, 55.0, 10.0, THUMS),
        ('R2', ROWS, 'vertical', None, 55.0, 10.0, THUMS),
        ('R3', P, 'vertical', None, 55.0, 10.0, THUMS),
        ('R4', P, 's_curve', kop, 55.0, 10.0, THUMS),
        ('R4_kop300', P, 's_curve', np.full(N, 300 * 0.3048), 55.0, 10.0, THUMS),
        ('R4_th45', P, 's_curve', kop, 45.0, 10.0, THUMS),
        ('R4_th70', P, 's_curve', kop, 70.0, 10.0, THUMS),
        ('R4_top1', P, 's_curve', kop, 55.0, 1.0, THUMS),
        ('R4_top30', P, 's_curve', kop, 55.0, 30.0, THUMS),
        ('R3_legacy', P, 'vertical', None, 55.0, 10.0, LEGACY),
        ('R4_legacy', P, 's_curve', kop, 55.0, 10.0, LEGACY),
    ]
    store, profiles = {}, {}
    nst = G.Ground.stations_for(L)
    for step, Pp, traj, kk, th, top, casing in SPEC:
        c0 = B.Case(n_segments=99, L_ft=L / 0.3048, t_operation_yr=top,
                    **casing)
        Rn = float(np.atleast_1d(
            B.formation_resistance(c0, np.array([0.0])))[0]
            - B.ground_resistance(c0))
        md, X = G.trajectories(Pp, L, kind=traj, kop=kk, theta_max=th)
        gr = G.Ground(k=c0.k_rock, rhocp=c0.rho_cp_rock, t_yr=top,
                      rb=c0.r_bore)
        pr = gr.profiles(md, X, Rn, nst=nst)
        assert gr.n_negative == 0, (step, gr.n_negative)
        store[f'{step}_md'] = pr['md']
        store[f'{step}_Rg'] = pr['temp']
        store[f'{step}_Rg_flux'] = pr['flux']
        store[f'{step}_tvd'] = pr['tvd']
        store[f'{step}_q_spread'] = pr['q_spread']
        profiles[step] = (pr['md'], pr['temp'])
        print(f'  profile {step:10s} nst={nst} '
              f'Rg {pr["temp"].mean():7.3f} m K/W', flush=True)
    store['_ds_station'] = np.array([G.DS_STATION])
    store['_n_stations'] = np.array([nst])
    store['_L_md'] = np.array([L])
    np.savez(f'{OUT}/profiles.npz', **store)
    return profiles


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = [r for r in load() if 'error' not in r]
    n = write_csv(rows)
    profiles = write_profiles()
    fig_geometry(); fig_Rg(profiles)
    fig_ladder(rows); fig_sweep(rows)
    print(f'pack written to {OUT}/  ({n} rows, 11 single-panel '
      f'figures, profiles.npz)')


if __name__ == '__main__':
    main()
