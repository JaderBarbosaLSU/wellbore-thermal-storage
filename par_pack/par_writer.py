"""The parametric pack: maps, the confound figure and the tables.

Reads par_runs.json and emits only views of it. Nothing here computes
physics, so a wrong number here is wrong in par_runner.
"""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import boreas_plots as BP
import par_grids as A

OUT = 'par_pack'
RED, NAVY, GREEN, GREY = '#b03028', '#23356b', '#2f6f4f', '#6a6a6a'


def load():
    raw = json.load(open('par_runs.json'))
    out = {}
    for k, v in raw.items():
        name = json.loads(k)[0]
        out.setdefault(name, []).append(v)
    return out


def ok(rows):
    """Feasible points only, tagged so heatmap hatches the rest."""
    for r in rows:
        r['feasible'] = 'error' not in r
    return rows


def fig_cascade(G):
    rows = ok(G['cascade'])
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))
    for ax, z, lab, fmt in (
            (axes[0], 'e_prime', r"$e'$ / kWh per m of well", '{:.3f}'),
            (axes[1], 'RTE', r'$\eta_\mathrm{RTE}$', '{:.3f}'),
            (axes[2], 'wells_per_MWe', 'storage wells per MWe', '{:.0f}')):
        BP.heatmap(rows, 'DT_3C_2C', 'T_m_C', z, A.DG, A.TM, ax=ax,
                   cbar_label=lab, fmt=fmt,
                   best='min' if z == 'wells_per_MWe' else 'max')
        ax.set_xlabel(r'cascade grading $\Delta T_\mathrm{3C,2C}$ / K')
        ax.set_ylabel(r'$T_m$ of the hottest layer / $^\circ$C')
    fig.suptitle(f'Cascade, at the reference case (R4), '
                 f'L = {A.L_REF:.0f} m MD. Hatched cells are infeasible: '
                 f'the bottom layer cannot drive the ORC.', y=1.02,
                 fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/map_cascade.pdf'); plt.close(fig)


def fig_hm(G):
    rows = ok(G['hm_Tm'])
    for r in rows:
        r['h_m_kJ'] = r['h_m'] / 1000.0
    HK = tuple(h / 1000.0 for h in A.HM)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))
    for ax, z, lab, fmt, best in (
            (axes[0], 'e_prime', r"$e'$ / kWh per m of well", '{:.3f}', 'max'),
            (axes[1], 'wells_per_MWe', 'storage wells per MWe', '{:.0f}', 'min'),
            (axes[2], 'loss', 'formation loss', '{:.3f}', 'min')):
        BP.heatmap(rows, 'h_m_kJ', 'T_m_C', z, HK, A.TM, ax=ax,
                   cbar_label=lab, fmt=fmt, best=best)
        ax.set_xlabel(r'latent heat $h_m$ / kJ kg$^{-1}$')
        ax.set_ylabel(r'$T_m$ of the hottest layer / $^\circ$C')
    fig.suptitle('Choosing the PCM: latent heat against melting '
                 f'temperature, at R4, L = {A.L_REF:.0f} m MD', y=1.02,
                 fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/map_hm.pdf'); plt.close(fig)


def fig_geometry(G):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))
    BP.heatmap(ok(G['L_Nlay']), 'L_md', 'N_lay', 'e_prime', A.LW_M, A.NL,
               ax=axes[0], cbar_label=r"$e'$ / kWh m$^{-1}$", fmt='{:.3f}')
    axes[0].set_xlabel('storage length $L$ / m MD')
    axes[0].set_ylabel(r'$N_\mathrm{lay}$')
    BP.heatmap(ok(G['L_Nlay']), 'L_md', 'N_lay', 'wells_per_MWe', A.LW_M,
               A.NL, ax=axes[1], cbar_label='storage wells per MWe',
               fmt='{:.0f}', best='min')
    axes[1].set_xlabel('storage length $L$ / m MD')
    axes[1].set_ylabel(r'$N_\mathrm{lay}$')
    BP.heatmap(ok(G['L_nft']), 'L_md', 'n_ft', 'RTE', A.LW_M, A.NFT,
               ax=axes[2], cbar_label=r'$\eta_\mathrm{RTE}$', fmt='{:.3f}')
    axes[2].set_xlabel('storage length $L$ / m MD')
    axes[2].set_ylabel(r'fin tubes per well $n_\mathrm{ft}$')
    fig.suptitle('Well geometry at R4. The $L$ axis stops at the 5000 ft '
                 'casing shoe.', y=1.02, fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/map_geometry.pdf'); plt.close(fig)


def fig_grad(G):
    rows = ok(G['L_grad'])
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3))
    BP.heatmap(rows, 'L_md', 'grad_K_per_km', 'loss', A.LW_M, A.GR,
               ax=axes[0], cbar_label='formation loss', fmt='{:.3f}',
               best='min')
    BP.heatmap(rows, 'L_md', 'grad_K_per_km', 'e_prime', A.LW_M, A.GR,
               ax=axes[1], cbar_label=r"$e'$ / kWh m$^{-1}$", fmt='{:.3f}')
    for ax in axes:
        ax.set_xlabel('storage length $L$ / m MD')
        ax.set_ylabel(r'geothermal gradient / K km$^{-1}$')
    fig.suptitle('Hotter rock against more boundary, at R4. The gradient '
                 'sets the rock TEMPERATURE and not $R_{g,\\mathrm{eff}}$, '
                 'so the field solve is shared down each column.',
                 y=1.02, fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/map_gradient.pdf'); plt.close(fig)


def fig_confound(G):
    """N_lay with the span free, and with it pinned."""
    free = sorted(ok(G['Nlay_free']), key=lambda r: r['N_lay'])
    pin = sorted(ok(G['Nlay_pinned']), key=lambda r: r['N_lay'])
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.0))
    for ax, z, lab in ((axes[0], 'e_prime', r"$e'$ / kWh m$^{-1}$"),
                       (axes[1], 'RTE', r'$\eta_\mathrm{RTE}$'),
                       (axes[2], 'span_K', 'cascade span / K')):
        ax.plot([r['N_lay'] for r in free], [r[z] for r in free], 'o-',
                color=RED, lw=1.8, ms=4, label='span free')
        ax.plot([r['N_lay'] for r in pin], [r[z] for r in pin], 's--',
                color=NAVY, lw=1.8, ms=4, label='span pinned at 48.9 K')
        ax.set_xlabel(r'$N_\mathrm{lay}$'); ax.set_ylabel(lab)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(r'Breaking the $N_\mathrm{lay}$ / span confound, at R4. '
                 'Sweeping the layer count at fixed grading moves the span '
                 'with it; pinning the span isolates the count.',
                 y=1.03, fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/confound_Nlay.pdf'); plt.close(fig)
    return free, pin


def write_csv(G):
    rows = []
    for name, rs in G.items():
        for r in rs:
            rows.append(dict(grid=name, **r))
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(f'{OUT}/par_results.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    return len(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    BP.use_paper_style()
    G = load()
    n = write_csv(G)
    fig_cascade(G); fig_hm(G); fig_geometry(G); fig_grad(G)
    free, pin = fig_confound(G)
    print(f'par pack written to {OUT}/  ({n} rows, 5 figures)')
    ef = np.array([r['e_prime'] for r in free])
    ep = np.array([r['e_prime'] for r in pin])
    print(f"  N_lay, span free   : e' spread "
          f"{100*(ef.max()-ef.min())/ef.mean():5.2f} %")
    print(f"  N_lay, span pinned : e' spread "
          f"{100*(ep.max()-ep.min())/ep.mean():5.2f} %")
    bad = sum(1 for rs in G.values() for r in rs if 'error' in r)
    print(f'  {bad} infeasible points, hatched on the maps')


if __name__ == '__main__':
    main()
