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
         'R3': 'Grissom slots,\nvertical', 'R4': 'Grissom slots,\nS-curve'}


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


# ------------------------------------------------------------- figure 1
def fig_geometry():
    P = G.grissom_positions(CSV_IN)
    kop, is_out = G.staggered_kop(P)
    fig, ax = plt.subplots(1, 3, figsize=(10.6, 3.5))

    for grp, col, lab in ((is_out, RED, 'outer row'), (~is_out, NAVY, 'inner row')):
        ax[0].plot(P[grp, 0], P[grp, 1], 'o', ms=1.8, color=col, label=lab)
    ax[0].set_aspect('equal')
    ax[0].set_xlabel('east / m'); ax[0].set_ylabel('north / m')
    ax[0].set_title(f'(a) {len(P)} Grissom slots', pad=6)
    ax[0].legend(frameon=False, fontsize=7.5, markerscale=3)

    md, X = G.trajectories(P, 2000.0, kind='s_curve', kop=kop)
    rad = np.hypot(X[:, :, 0], X[:, :, 1])
    for i in range(0, len(P), 7):
        ax[1].plot(rad[i], X[i, :, 2], lw=0.4, color=GREY, alpha=0.7)
    ax[1].invert_yaxis()
    ax[1].set_xlabel('radius from centroid / m'); ax[1].set_ylabel('TVD / m')
    ax[1].set_title('(b) S-curve trajectories', pad=6)

    from scipy.spatial import cKDTree
    for kind, col, lab in (('vertical', NAVY, 'vertical'),
                           ('s_curve', RED, 'S-curve')):
        m2, X2 = G.trajectories(P, 2000.0, kind=kind, kop=kop)
        nn = []
        for j in range(0, len(m2), 8):
            d, _ = cKDTree(X2[:, j, :2]).query(X2[:, j, :2], k=2)
            nn.append(np.median(d[:, 1]))
        ax[2].plot(m2[::8], nn, color=col, lw=1.6, label=lab)
    ax[2].set_yscale('log')
    ax[2].set_xlabel('MD / m'); ax[2].set_ylabel('median spacing / m')
    ax[2].set_title('(c) nearest-neighbour spacing', pad=6)
    ax[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(f'{OUT}/grissom_geometry.pdf')
    plt.close(fig)


# ------------------------------------------------------------- figure 2
def fig_Rg(profiles):
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for step, col in (('R3', NAVY), ('R4', RED)):
        if step not in profiles:
            continue
        md, rg = profiles[step]
        ax.plot(rg, md, color=col, lw=1.8, label=f'{step}: ' +
                LABEL[step].replace('\n', ' '))
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.set_xlabel(r'$R_{g,\mathrm{eff}}$  /  m K W$^{-1}$')
    ax.set_ylabel('MD / m')
    ax.set_title('Shielding along the well, equal-temperature superposition',
                 pad=7, fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(f'{OUT}/grissom_Rg.pdf')
    plt.close(fig)


# ------------------------------------------------------------- figure 3
def fig_ladder(rows, L=1250.0):
    steps, loss, ep = [], [], []
    for s in LADDER:
        r = pick(rows, s, L)
        if not r:
            continue
        steps.append(s); loss.append(100 * r[0]['loss']); ep.append(r[0]['e_prime'])
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6))
    x = np.arange(len(steps))
    ax[0].bar(x, loss, color=[RED if v > 50 else NAVY for v in loss])
    ax[0].set_xticks(x); ax[0].set_xticklabels([LABEL[s] for s in steps], fontsize=7.5)
    ax[0].set_ylabel('formation loss / %')
    for i, v in enumerate(loss):
        ax[0].text(i, v + 1, f'{v:.1f}', ha='center', fontsize=7.5)
    ax[1].bar(x, ep, color=GREEN)
    ax[1].set_xticks(x); ax[1].set_xticklabels([LABEL[s] for s in steps], fontsize=7.5)
    ax[1].set_ylabel(r"$e'$ / kWh per m of well")
    for i, v in enumerate(ep):
        ax[1].text(i, v + 0.005, f'{v:.3f}', ha='center', fontsize=7.5)
    fig.suptitle(f'The field ladder at L = {L:.0f} m MD, '
                 f'{len(steps)} rungs, Grissom (N = 363)', y=1.02, fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{OUT}/grissom_ladder.pdf')
    plt.close(fig)


# ------------------------------------------------------------- figure 4
def fig_sweep(rows):
    fig, ax = plt.subplots(1, 3, figsize=(12.0, 3.6))
    cols = dict(zip(LADDER, (GREY, '#9c27a0', NAVY, '#0a9396', GREEN, RED)))
    for s in LADDER:
        rr = pick(rows, s)
        if len(rr) < 2:
            continue
        L = np.array([r['L_md'] for r in rr])
        ax[0].plot(L, [r['e_prime'] for r in rr], 'o-', ms=3, lw=1.4,
                   color=cols[s], label=s)
        ax[1].plot(L, [r['RTE'] for r in rr], 'o-', ms=3, lw=1.4, color=cols[s])
        Eo = np.array([r['E_out_kWh'] for r in rr])
        Ei = np.array([r['E_in_kWh'] for r in rr])
        dE = np.diff(Ei)
        marg = np.where(np.abs(dE) > 1e-12, np.diff(Eo) / np.where(dE == 0, 1, dE),
                        np.nan)
        ax[2].plot(0.5 * (L[1:] + L[:-1]), marg, 'o-', ms=3, lw=1.4, color=cols[s])
    for a, lab in zip(ax, (r"$e'$ / kWh m$^{-1}$", r'$\eta_\mathrm{RTE}$',
                           r'marginal $\Delta E_\mathrm{out}/\Delta E_\mathrm{in}$')):
        a.set_xlabel('storage length $L$ / m MD'); a.set_ylabel(lab)
        a.axvline(1250.0, color=GREY, lw=0.9, ls=(0, (2, 2)))
    ax[2].set_ylim(bottom=0)
    ax[0].legend(frameon=False, fontsize=7.5, ncol=2)
    ax[0].annotate('shoe (provisional)', (1250, ax[0].get_ylim()[1]),
                   fontsize=7, color=GREY, ha='right', va='top', rotation=90,
                   xytext=(-3, -2), textcoords='offset points')
    fig.tight_layout()
    fig.savefig(f'{OUT}/grissom_sweep.pdf')
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = [r for r in load() if 'error' not in r]
    n = write_csv(rows)
    profiles = {}
    THUMS = dict(D_in=201.2 / 25.4, r_casing_o=0.21910 / 2,
                 r_bore=0.31120 / 2, k_cement=0.7)
    P = G.grissom_positions(CSV_IN)
    kop, _ = G.staggered_kop(P)
    c0 = B.Case(n_segments=99, L_ft=1250 / 0.3048, **THUMS)
    Rn = float(np.atleast_1d(B.formation_resistance(c0, np.array([0.0])))[0]
               - B.ground_resistance(c0))
    for step, kind, kk in (('R3', 'vertical', None), ('R4', 's_curve', kop)):
        md, X = G.trajectories(P, 1250.0, kind=kind, kop=kk)
        gr = G.Ground(k=c0.k_rock, rhocp=c0.rho_cp_rock, t_yr=10.0, rb=c0.r_bore)
        pr = gr.profiles(md, X, Rn, nst=15 if kind == 'vertical' else 50)
        profiles[step] = (pr['md'], pr['temp'])
        np.savez(f'{OUT}/profiles.npz', **{f'{k}_md': v[0] for k, v in profiles.items()},
                 **{f'{k}_Rg': v[1] for k, v in profiles.items()})
    fig_geometry(); fig_Rg(profiles); fig_ladder(rows); fig_sweep(rows)
    print(f'pack written to {OUT}/  ({n} rows, 4 figures, profiles.npz)')


if __name__ == '__main__':
    main()
