"""Morris elementary-effects screen for THUMS.

WHY MORRIS AND NOT A TWO-LEVEL FACTORIAL. A 2^(k-p) design computes the
SIGNED mean effect of each factor. `psi` has an interior optimum in
DT_3A_4A near 20 K, so the effects are positive below it and negative above
it and the signed mean nearly cancels: a factorial would report the one
factor we already know has an optimum as unimportant. Morris reports mu*
(mean ABSOLUTE effect) and sigma (its spread) separately, so a factor that is
important and non-monotone lands at low mu, high mu*, high sigma -- visibly
different from a factor that is merely unimportant.

WHAT IT DOES NOT DO. Morris says a factor interacts or curves; it never says
with whom. There are no interaction partners here and no variance shares.
That is stage two (LHS + surrogate + Sobol). This is a screen.

RANGES ARE A PROPOSAL. Edit FACTORS below. Every bound is argued in the
comment beside it, and several are coupled -- see `cascade_margin`, which is
sampled as a margin rather than a temperature precisely so that the
span constraint cannot be violated by construction.
"""
import itertools
import json
import math
import pathlib
import numpy as np

import thums as T

# ----------------------------------------------------------------------
# Pipe schedule. r_i and r_e are NOT independent -- they are a schedule, and
# sampling them separately would specify pipes that do not exist. The factor
# is nominal size; both radii follow. Schedule 80 throughout.
#   NPS : (OD mm, wall mm)
# 1 in is NOT in the list. At this duty it needs 2.4 kg/s per tube and the
# pumping reaches 1057 kW against a 1000 kW output -- f_pump = 2.87 and
# eta_RTE goes NEGATIVE. That is not a bad design, it is not a plant, and
# including it would let a single broken corner dominate mu* for every
# efficiency KPI. Smallest size kept is the current one.
SCH80 = {1.25: (42.16, 4.85), 1.50: (48.26, 5.08), 2.00: (60.33, 5.54)}
NPS_LIST = sorted(SCH80)


def pipe_sch80(nps):
    """(r_i, r_e) in metres for a nominal size, schedule 80."""
    od, wall = SCH80[min(SCH80, key=lambda k: abs(k - nps))]
    r_e = od / 2000.0
    return r_e - wall / 1000.0, r_e


# ----------------------------------------------------------------------
# THE FACTORS.  (low, high, kind)
#   kind 'f' continuous, 'i' integer, 'nps' pipe index into NPS_LIST
FACTORS = {
    # --- store hardware ------------------------------------------------
    # Nominal size, Sch 80. Three-way: bigger pipe means LESS PCM (area goes
    # as r^2), a SHORTER conduction path to the merge radius, a WORSE film
    # coefficient (h ~ D^-1.8 at fixed flow) and MUCH cheaper pumping
    # (dp ~ D^-4.75). The sign of the net effect is not obvious, which is
    # the point of including it.
    'nps':            (0, len(NPS_LIST) - 1, 'nps'),
    # The report's sweep crosses from capacity-bound to rate-bound between
    # 12 and 16 fins, so the range must straddle that.
    'num_fins':       (8, 32, 'i'),

    # --- temperature level ---------------------------------------------
    # The top of the cascade. Capped at the current value because
    # T_m_C + DT_4C_M must stay below water saturation at P = 1 MPa
    # (179.9 C); at DT_4C_M = 20 the charge inlet is already 170 C. Going
    # higher needs P raised, which is a separate scenario, not a factor.
    'T_m_C':          (125.0, 150.0, 'f'),

    # --- glides, now genuinely independent ------------------------------
    # DT_3C_2C sets the CHARGE flow only, provided DT_3D_2D is set
    # explicitly -- which it is below. Leaving DT_3D_2D at None would make
    # it fall back to DT_3C_2C and silently re-couple the two flows, which
    # is the defect DN-18 exists to prevent.
    'DT_3C_2C':       (40.0, 70.0, 'f'),
    'DT_3D_2D':       (40.0, 70.0, 'f'),

    # --- the cascade ----------------------------------------------------
    # NOT sampled as a temperature. The span is bounded BELOW by
    # glide - DT_4C_M, and that bound moves with two other factors, so a box
    # in DT_cascade would put a large share of the design infeasible. We
    # sample the fractional margin ABOVE the bound instead, and every draw
    # is feasible by construction. See `build_case`.
    'cascade_margin': (0.05, 0.40, 'f'),
    'N_lay':          (5, 15, 'i'),

    # --- approaches -----------------------------------------------------
    'DT_4C_M':        (5.0, 20.0, 'f'),
    'DT_M_1D':        (5.0, 20.0, 'f'),
    # The binding ORC constraint. Below ~2 K the evaporator explodes; above
    # ~12 K the boiling temperature and eta_ORC collapse.
    'DT_pinch_ORC':   (2.0, 12.0, 'f'),
    # The ONE loose approach in the plant at the design point (v0.14 pinch
    # table: the other three sit hard against their specified value). Worth
    # testing whether tightening it pays.
    'DT_2H_3C':       (4.0, 20.0, 'f'),

    # --- the source side ------------------------------------------------
    # Promoted at v0.13. Sets the evaporating temperature AND the geothermal
    # flow, and has an interior optimum in psi near 20 K. Upper bound keeps
    # reinjection above T_reinject_min_C = 25 C.
    'DT_3A_4A':       (5.0, 30.0, 'f'),

    # --- operation ------------------------------------------------------
    # Sets how much RATE is demanded, which is what makes fins worth
    # anything. Charge and discharge windows moved together.
    't_cycle_h':      (4.0, 16.0, 'f'),
}

# KPIs pulled from performance_indices. The four objectives first.
KPIS = ['eta_RTE', 'psi', 'wells_per_MWe', 'UA_per_MWe',
        'rho_P_kW_m3', 'kW_per_well', 'eps_RTE', 'dE_elec_kWh',
        'rho_E_kWh_m3', 'I_per_MWh', 'f_pump', 'resource_utilisation',
        'N_geo', 'deviation', 'merge_proximity_max', 'regime_pinch_bound',
        'eps_cycled', 'n_warnings']


def build_case(x):
    """Map a unit-cube point to a Case. Feasibility is built in, not screened.

    `cascade_margin` is the only non-obvious one. The span must satisfy
        span = DT_cascade (N-1)/N  >  glide - DT_4C_M
    or the driving difference inverts at the far end of the well. Sampling
    DT_cascade in a box would violate that over much of the range, and the
    violated corner is not a bad design -- it is a run that cannot be
    interpreted. So the sampled quantity is the fractional margin above the
    bound, and DT_cascade is reconstructed from it.
    """
    v = {}
    for (name, (lo, hi, kind)), u in zip(FACTORS.items(), x):
        if kind == 'i':
            v[name] = int(round(lo + u * (hi - lo)))
        elif kind == 'nps':
            v[name] = NPS_LIST[int(round(lo + u * (hi - lo)))]
        else:
            v[name] = lo + u * (hi - lo)

    r_i, r_e = pipe_sch80(v['nps'])
    N_lay = max(2, v['N_lay'])
    glide = v['DT_3C_2C']
    bound = max(1e-6, glide - v['DT_4C_M'])             # span must exceed this
    span_target = bound * (1.0 + v['cascade_margin'])
    DT_cascade = span_target * N_lay / (N_lay - 1)

    # n_segments a MULTIPLE of N_lay, so layer_map never puts a boundary
    # inside a segment. Otherwise N_lay carries a numerical artefact that
    # the screen would read as a physical effect.
    n_seg = int(N_lay * max(4, round(100 / N_lay)))

    return T.Case(
        r_i=r_i, r_e=r_e, num_fins=int(v['num_fins']),
        T_m_C=v['T_m_C'],
        DT_3C_2C=glide, DT_3D_2D=v['DT_3D_2D'], DT_cascade=DT_cascade,
        N_lay=N_lay, n_segments=n_seg,
        DT_4C_M=v['DT_4C_M'], DT_M_1D=v['DT_M_1D'],
        DT_pinch_ORC=v['DT_pinch_ORC'], DT_2H_3C=v['DT_2H_3C'],
        DT_3A_4A=v['DT_3A_4A'],
        t_ch=v['t_cycle_h'], t_dc=v['t_cycle_h'])


def evaluate(x):
    """One design point -> dict of KPIs, or None if it cannot be run."""
    try:
        case = build_case(x)
        if any(k == 'error' for k, _ in T.validate_case(case, verbose=False)):
            return None
        r = T.simulate_css_corrected(case)
        k = T.performance_indices(case, r)
        # A design whose parasitics take half the gross output is not a
        # worse plant, it is not a plant. Left in, the pumping corner --
        # which spans a factor of 1400 in f_pump across this cube, because
        # dp ~ mdot^2.75 / D^4.75 -- would dominate mu* for every efficiency
        # KPI and the screen would report nothing but "pipe size and cycle
        # length are everything".
        if not (k['f_pump'] < 0.5):
            return None
        return {n: float(k[n]) for n in KPIS if n in k}
    except Exception:
        return None


def trajectories(r, p=6, seed=0):
    """Morris trajectories in the unit cube: r paths of k+1 points each.

    Each step moves ONE factor by +/- delta on a p-level grid, so the
    difference between consecutive points is an elementary effect of that
    factor and nothing else. delta = p/(2(p-1)) makes the design balanced.
    """
    rng = np.random.default_rng(seed)
    k = len(FACTORS)
    delta = p / (2.0 * (p - 1.0))
    B = np.tril(np.ones((k + 1, k)), -1)
    out = []
    for _ in range(r):
        xstar = rng.integers(0, p // 2, size=k) / (p - 1.0)
        D = np.diag(rng.choice([-1.0, 1.0], size=k))
        P = np.eye(k)[:, rng.permutation(k)]
        Bs = (np.ones((k + 1, 1)) @ xstar.reshape(1, -1)
              + (delta / 2.0) * ((2.0 * B - np.ones((k + 1, k))) @ D
                                 + np.ones((k + 1, k)))) @ P
        out.append(np.clip(Bs, 0.0, 1.0))
    return out, delta


def run(r=10, p=6, seed=0, out='morris_raw.json', verbose=True):
    """Evaluate every point of every trajectory. ~17 s per point."""
    trajs, delta = trajectories(r, p, seed)
    k = len(FACTORS)
    rows, n_bad = [], 0
    for ti, Bs in enumerate(trajs):
        for si in range(k + 1):
            res = evaluate(Bs[si])
            if res is None:
                n_bad += 1
            rows.append({'traj': ti, 'step': si,
                         'x': Bs[si].tolist(), 'kpi': res})
        if verbose:
            print(f'  trajectory {ti + 1}/{r} done, {n_bad} infeasible so far',
                  flush=True)
    pathlib.Path(out).write_text(json.dumps(
        {'factors': list(FACTORS), 'kpis': KPIS, 'delta': delta,
         'r': r, 'p': p, 'seed': seed, 'rows': rows}))
    if verbose:
        print(f'wrote {out}: {len(rows)} points, {n_bad} infeasible '
              f'({100 * n_bad / len(rows):.1f} %)')
    return rows


def analyse(path='morris_raw.json'):
    """mu, mu*, sigma per factor per KPI, with bootstrap CI on mu*.

    A step whose either end is infeasible contributes NO elementary effect.
    Dropping them is not neutral -- it biases mu* toward the interior -- so
    the count of usable effects is reported per factor and should be read
    beside the ranking.
    """
    d = json.loads(pathlib.Path(path).read_text())
    names, delta, k = d['factors'], d['delta'], len(d['factors'])
    rows = d['rows']
    ee = {n: {f: [] for f in names} for n in d['kpis']}
    for ti in range(d['r']):
        path_rows = [x for x in rows if x['traj'] == ti]
        path_rows.sort(key=lambda z: z['step'])
        for a, b in zip(path_rows[:-1], path_rows[1:]):
            if a['kpi'] is None or b['kpi'] is None:
                continue
            dx = np.array(b['x']) - np.array(a['x'])
            j = int(np.argmax(np.abs(dx)))
            if abs(dx[j]) < 1e-12:
                continue
            for n in d['kpis']:
                if n in a['kpi'] and n in b['kpi']:
                    ee[n][names[j]].append((b['kpi'][n] - a['kpi'][n]) / dx[j])
    rng = np.random.default_rng(1)
    res = {}
    for n in d['kpis']:
        res[n] = {}
        scale = max(1e-30, np.mean([abs(v) for f in names for v in ee[n][f]])
                    or 1e-30)
        for f in names:
            v = np.array(ee[n][f], float)
            if v.size == 0:
                res[n][f] = dict(mu=0., mu_star=0., sigma=0., n=0, ci=0.)
                continue
            bs = [np.mean(np.abs(rng.choice(v, v.size))) for _ in range(400)]
            res[n][f] = dict(mu=float(v.mean() / scale),
                             mu_star=float(np.abs(v).mean() / scale),
                             sigma=float(v.std(ddof=1) / scale) if v.size > 1
                             else 0.0,
                             n=int(v.size),
                             ci=float(np.percentile(bs, 95) / scale
                                      - np.percentile(bs, 5) / scale))
    return res, d


if __name__ == '__main__':
    import sys
    r = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    print(f'Morris screen: {len(FACTORS)} factors, r = {r} trajectories, '
          f'{r * (len(FACTORS) + 1)} runs')
    run(r=r)


# ======================================================================
#  PLOTS.  Four standard Morris figures plus one diagnostic specific to
#  this model.  All take the dict returned by `analyse`.
# ======================================================================
def _mpl():
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'figure.dpi': 110})
    return plt


GRAY, BLUE, ORNG, REDC = '#898781', '#2a78d6', '#eb6834', '#d03b3b'


def _classify(m):
    """Colour by quadrant: negligible / linear / non-linear / non-monotone."""
    if m['mu_star'] < 0.20:
        return GRAY, 'o'
    if abs(m['mu']) < 0.45 * m['mu_star']:      # signed mean cancels
        return REDC, '*'
    if m['sigma'] > m['mu_star']:
        return ORNG, '^'
    return BLUE, 's'


def plot_plane(res, kpi, ax=None, label=True):
    """THE Morris plot: mu* against sigma, one KPI."""
    plt = _mpl()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    hi = max(max(m['mu_star'] for m in res[kpi].values()), 1e-9)
    ax.plot([0, hi], [0, hi], ls='--', lw=0.9, color='#c3c2b7', zorder=1)
    for f, m in res[kpi].items():
        c, mk = _classify(m)
        ax.scatter([m['mu_star']], [m['sigma']], c=c, marker=mk,
                   s=150 if mk == '*' else 45, zorder=4,
                   edgecolors='white', linewidths=0.6)
        if label and m['mu_star'] > 0.12 * hi:
            ax.annotate(f, (m['mu_star'], m['sigma']), fontsize=7.5,
                        xytext=(4, 0), textcoords='offset points',
                        va='center', color='#3a3936')
    ax.set_xlabel(r'$\mu^*$'); ax.set_ylabel(r'$\sigma$')
    ax.set_title(kpi, loc='left', fontsize=10)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    return ax


def plot_planes(res, kpis=('eta_RTE', 'psi', 'wells_per_MWe', 'UA_per_MWe',
                           'rho_P_kW_m3', 'kW_per_well')):
    """The grid that makes the multi-objective argument visually.

    Same axes, one panel per KPI. If the rankings agree there is a single
    objective and this whole discussion was unnecessary; if they disagree,
    ranking designs on one of them recommends the corner that is free in
    that column and expensive in the others.
    """
    plt = _mpl()
    n = len(kpis)
    fig, axes = plt.subplots((n + 2) // 3, 3, figsize=(13, 3.6 * ((n + 2) // 3)))
    for ax, k in zip(np.ravel(axes), kpis):
        plot_plane(res, k, ax)
    for ax in np.ravel(axes)[n:]:
        ax.axis('off')
    fig.tight_layout()
    return fig


def plot_mu_mustar(res, kpi):
    """Signed mean against absolute mean -- the NON-MONOTONE detector.

    On the diagonal a factor pushes one way everywhere. Far below it, the
    effect changes sign across the space, which is the fingerprint of an
    interior optimum. `plot_plane` cannot distinguish that from noise.
    """
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(7, 5))
    hi = max(max(m['mu_star'] for m in res[kpi].values()), 1e-9)
    ax.plot([0, hi], [0, hi], ls='--', lw=0.9, color='#c3c2b7')
    ax.plot([0, hi], [0, -hi], ls='--', lw=0.9, color='#c3c2b7')
    for f, m in res[kpi].items():
        c, mk = _classify(m)
        ax.scatter([m['mu_star']], [m['mu']], c=c, marker=mk,
                   s=150 if mk == '*' else 45, edgecolors='white', linewidths=0.6)
        if m['mu_star'] > 0.12 * hi:
            ax.annotate(f, (m['mu_star'], m['mu']), fontsize=7.5,
                        xytext=(4, 0), textcoords='offset points', va='center')
    ax.axhline(0, lw=0.8, color='#c3c2b7')
    ax.set_xlabel(r'$\mu^*$   (absolute)'); ax.set_ylabel(r'$\mu$   (signed)')
    ax.set_title(f'{kpi} — points far off the diagonals change sign', loc='left',
                 fontsize=10)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return fig


def plot_ranked(res, kpi):
    """Ranked mu* with bootstrap CI, and the count of usable effects.

    The count matters: a step with an infeasible endpoint contributes no
    elementary effect, so a factor whose range strays outside the feasible
    region is estimated from fewer samples AND biased toward the interior.
    """
    plt = _mpl()
    it = sorted(res[kpi].items(), key=lambda kv: kv[1]['mu_star'])
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(it) + 1.4))
    y = np.arange(len(it))
    ax.barh(y, [m['mu_star'] for _, m in it],
            xerr=[m['ci'] / 2 for _, m in it],
            color=[_classify(m)[0] for _, m in it], height=0.62,
            error_kw=dict(lw=0.9, ecolor='#5f5e5a'))
    ax.set_yticks(y)
    ax.set_yticklabels([f'{f}  (n={m["n"]})' for f, m in it], fontsize=8.5)
    ax.set_xlabel(r'$\mu^*$   (normalised)')
    ax.set_title(f'{kpi} — ranked, with 90 % bootstrap interval', loc='left',
                 fontsize=10)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return fig


def plot_rankshift(res, kpis=('eta_RTE', 'psi', 'wells_per_MWe', 'UA_per_MWe')):
    """Where a factor ranks under each KPI, joined by lines.

    Flat lines mean the KPIs agree. Crossings are the whole argument for
    not collapsing the objectives into one number.
    """
    plt = _mpl()
    facs = list(res[kpis[0]])
    order = {k: {f: i for i, f in enumerate(
        sorted(facs, key=lambda g: -res[k][g]['mu_star']))} for k in kpis}
    fig, ax = plt.subplots(figsize=(1.9 * len(kpis) + 3.4, 0.36 * len(facs) + 1.6))
    cmap = plt.get_cmap('tab20')
    for i, f in enumerate(facs):
        ys = [order[k][f] for k in kpis]
        ax.plot(range(len(kpis)), ys, '-o', ms=4, lw=1.4, color=cmap(i % 20))
        ax.annotate(f, (0, ys[0]), xytext=(-8, 0), textcoords='offset points',
                    ha='right', va='center', fontsize=8)
        ax.annotate(f, (len(kpis) - 1, ys[-1]), xytext=(8, 0),
                    textcoords='offset points', ha='left', va='center', fontsize=8)
    ax.set_xticks(range(len(kpis))); ax.set_xticklabels(kpis, fontsize=9)
    ax.set_xlim(-1.4, len(kpis) - 1 + 1.4)
    ax.invert_yaxis(); ax.set_yticks([]); ax.set_ylabel('rank by $\\mu^*$')
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return fig


def plot_ee_strip(path, kpi, max_factors=13):
    """Raw elementary effects per factor -- NOT standard, and worth having.

    Morris compresses non-linearity and regime switching into one sigma and
    cannot tell them apart. A factor whose effects are merely SCATTERED is
    non-linear; one whose effects are BIMODAL has a regime change in its
    range -- the ORC flipping between pinch-bound and hot-end-bound, or the
    fin count crossing from capacity- to rate-bound. We know both switches
    exist, so look at the distribution rather than its summary.
    """
    plt = _mpl()
    d = json.loads(pathlib.Path(path).read_text())
    names, rows = d['factors'], d['rows']
    ee = {f: [] for f in names}
    for ti in range(d['r']):
        pr = sorted([x for x in rows if x['traj'] == ti], key=lambda z: z['step'])
        for a, b in zip(pr[:-1], pr[1:]):
            if a['kpi'] is None or b['kpi'] is None or kpi not in a['kpi']:
                continue
            dx = np.array(b['x']) - np.array(a['x'])
            j = int(np.argmax(np.abs(dx)))
            if abs(dx[j]) > 1e-12:
                ee[names[j]].append((b['kpi'][kpi] - a['kpi'][kpi]) / dx[j])
    fig, ax = plt.subplots(figsize=(8, 0.42 * len(names) + 1.4))
    rng = np.random.default_rng(0)
    for i, f in enumerate(names[:max_factors]):
        v = np.array(ee[f], float)
        if v.size:
            ax.scatter(v, np.full(v.size, i) + rng.normal(0, .07, v.size),
                       s=22, alpha=.75, color=BLUE, edgecolors='none')
    ax.axvline(0, lw=0.8, color='#c3c2b7')
    ax.set_yticks(range(len(names[:max_factors])))
    ax.set_yticklabels(names[:max_factors], fontsize=8.5)
    ax.set_xlabel(f'elementary effect on {kpi}')
    ax.set_title('scattered = non-linear;  BIMODAL = a regime change in range',
                 loc='left', fontsize=10)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return fig


def summary_table(res, kpi):
    """mu, mu*, sigma, n as a DataFrame, ranked."""
    import pandas as pd
    df = pd.DataFrame(res[kpi]).T[['mu', 'mu_star', 'sigma', 'ci', 'n']]
    return df.sort_values('mu_star', ascending=False).round(4)
