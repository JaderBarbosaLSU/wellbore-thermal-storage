"""Study 2: the borehole battery alone, by fractional factorial.

THE SPLIT. The plant is a chain and only one thing crosses between the
machines at the surface and the rock below: hot water in a pipe. Four
temperatures describe it completely --

    how hot the water goes DOWN when charging      T_4c
    how much cooler it comes BACK                  the charge glide
    how cold it goes DOWN when discharging         T_3d
    how much warmer it comes BACK                  the discharge glide

plus how much heat has to move per cycle. Fix those five and the two halves
stop talking to each other: the heat pump and the ORC do not care what the
rock is made of, and the rock does not care what the machines are.

This module holds all five FIXED at the design values and varies only what is
below ground: the PCM itself, the tubes, the fins, the cascade. The surface
study is the mirror image and lives elsewhere.

WHY A FACTORIAL AND NOT MORRIS. The Morris screen of v0.16 was chosen because
the model was expected to be non-linear. It largely is not: for psi and
eta_RTE, eleven of thirteen factors had |mu|/mu* = 1.00, meaning every
elementary effect had the same sign. Morris's advantage did not materialise,
and its cost did -- it cannot give two-factor interactions. For the borehole
those interactions are the interesting part: does the value of fins depend on
the PCM conductivity? does pipe size interact with the cycle length? A
Resolution V factorial answers both in the same 128 runs.

THE DESIGN IS SEARCHED AND VERIFIED, NOT QUOTED. The generators below were
found by exhaustive search over the candidate words of length >= 4, taking the
minimum-aberration set among those with no short defining word. `verify_design`
then forms all 9 main-effect and 36 two-factor-interaction columns and checks
X'X = 128 I exactly. Run it before spending an hour.

At nine factors the search returns something better than it did at eleven.
The defining relation is I = ABCDE = ABCFG = DEFG -- three words, all of
length 5 or 6, none of length 4 -- which is RESOLUTION VI: two-factor
interactions are now aliased only with FOUR-factor interactions, a rung
cleaner than the Resolution V design the eleven-factor version had to accept.
64 runs was checked and cannot do it: no Resolution V design exists for nine
factors in 64 runs, so 128 stands.
"""
import itertools
import json
import pathlib
import numpy as np

import thums as T

# ----------------------------------------------------------------------
# The design. 7 basic factors, 4 generated, 128 runs, Resolution V.
# Found by search; see verify_design().
BASIC = 7
GENERATORS = [(0, 1, 2, 3, 4),     # H = ABCDE
              (0, 1, 2, 5, 6)]     # J = ABCFG
LETTERS = 'ABCDEFGHJ'
NF = BASIC + len(GENERATORS)       # 9 factors


def design_matrix():
    """The 128 x 9 matrix in -1/+1."""
    n = 2 ** BASIC
    B = np.array([[1 if (i >> j) & 1 else -1 for j in range(BASIC)]
                  for i in range(n)])
    G = np.column_stack([np.prod(B[:, list(g)], axis=1) for g in GENERATORS])
    return np.column_stack([B, G])


def verify_design(verbose=True):
    """Check that every main effect and every 2fi is estimable and unconfounded.

    This is the whole justification for the design. If it fails, the
    interaction table the study produces is fiction, so it is checked in code
    rather than asserted in a comment.
    """
    X = design_matrix()
    cols, labs = [], []
    for i in range(NF):
        cols.append(X[:, i]); labs.append(LETTERS[i])
    for i, j in itertools.combinations(range(NF), 2):
        cols.append(X[:, i] * X[:, j]); labs.append(LETTERS[i] + LETTERS[j])
    M = np.column_stack(cols)
    XtX = M.T @ M
    off = np.abs(XtX - np.diag(np.diag(XtX))).max()
    ok = (off == 0) and bool(np.all(np.diag(XtX) == len(X)))
    if verbose:
        print(f"design: {X.shape[0]} runs x {X.shape[1]} factors")
        print(f"balanced: {bool(np.all(X.sum(0) == 0))}")
        print(f"model matrix: {M.shape[1]} columns "
              f"({NF} main + {M.shape[1] - NF} two-factor interactions)")
        print(f"largest off-diagonal of X'X: {off}")
        print("RESOLUTION V OR BETTER CONFIRMED" if ok
              else "*** DESIGN IS CONFOUNDED ***")
        print("\nshortest defining word is length 5, so 2fi alias with 4fi;")
        print("main effects and 2fi are clean.")
    return ok


# ----------------------------------------------------------------------
# The factors. Two levels each: (low, high, kind).
# Everything here is BELOW GROUND. Nothing on this list changes the four
# water temperatures or the duty -- those are the interface and are frozen.
FACTORS = {
    # --- the PCM itself. Absent from the v0.16 screen entirely, which for a
    #     study of a phase-change battery was the wrong omission: latent heat
    #     and conductivity are the two things material selection actually buys.
    'h_m_kJkg':   (180.0, 380.0, 'f'),   # latent heat of fusion
    'k_scale':    (0.5, 2.0, 'f'),       # multiplies BOTH k_s and k_l,
                                         # keeping the solid/liquid ratio
    'rho_scale':  (0.85, 1.15, 'f'),     # multiplies both densities
    'cp_scale':   (0.8, 1.2, 'f'),       # multiplies both specific heats

    # --- the borehole and what is in it
    'L_ft':       (3000.0, 8000.0, 'f'),  # well depth
    'num_tubes':  (1, 3, 'i'),            # HAIRPINS per hole (2 legs each)
    'nps':        (0, 1, 'nps'),          # 1.25 in or 2 in, Sch 80
    'num_fins':   (8, 32, 'i'),        # fin LENGTH is frozen -- see FROZEN
    'N_lay':      (5, 15, 'i'),        # cascade layers; the SPAN is held
                                       # fixed, so this varies granularity
                                       # only -- see build_case
}
NPS_LIST = [1.25, 2.0]
SCH80 = {1.25: (42.16, 4.85), 2.00: (60.33, 5.54)}

# The interface, frozen. These are the five numbers the surface study owns.
BASE = T.Case()
FROZEN = dict(T_m_C=BASE.T_m_C, DT_4C_M=BASE.DT_4C_M, DT_M_1D=BASE.DT_M_1D,
              DT_3C_2C=BASE.DT_3C_2C, DT_3D_2D=55.0,
              t_ch=BASE.t_ch, t_dc=BASE.t_dc,
              W_dot_el_out=BASE.W_dot_el_out,
              # Fin LENGTH fixed at 3.5 mm, down from the 7.5 mm default. Only
              # the fin COUNT varies now (8 -> 32), so the design asks "how
              # many fins" and not "how much fin", which is the cleaner
              # question and the one a manufacturer actually answers. 3.5 mm
              # also fits inside the tightest cell in the design -- three
              # hairpins of 2 in pipe -- so no corner is lost to geometry.
              fin_L=0.0035,
              # Fouling fixed at the design value. In the v0.16 screen it
              # ranked last of thirteen; sampling it here would spend a
              # column of a 9-factor design on a known non-effect.
              Rf_i=BASE.Rf_i)
SPAN_FIXED = 55.0 * (9 - 1) / 9          # the v0.16 design-point span, 48.89 K

# eps_RTE is on this list even though it is exactly eta_RTE * eps_cycled and
# so was always recoverable from the two of them. It was left off the first
# run and had to be reconstructed by hand to answer a question about it --
# which worked, but only because the identity is exact. A KPI that has to be
# rebuilt outside the file is one that does not appear in the plots, the
# summary table or the curvature test, and this one turned out to carry the
# study's most interesting disagreement with eta_RTE (conductivity and fins
# rank 9th and 8th on efficiency, 2nd and 5th on utilisation). Cheap to store,
# expensive to notice late.
KPIS = ['eta_RTE', 'eps_RTE', 'psi', 'wells_per_MWe', 'kW_per_well',
        'rho_E_kWh_m3', 'rho_P_kW_m3', 'dE_elec_kWh', 'eps_cycled', 'f_pump',
        'UA_per_MWe', 'deviation', 'merge_proximity_max',
        'regime_pinch_bound', 'n_warnings', 'N_geo', 'I_per_MWh']


def build_case(row):
    """Map one -1/+1 design row to a Case. The interface stays frozen."""
    v = {}
    for (name, (lo, hi, kind)), s in zip(FACTORS.items(), row):
        u = 0.5 * (s + 1.0)                       # -1/+1 -> 0/1
        if kind == 'i':
            v[name] = int(round(lo + u * (hi - lo)))
        elif kind == 'nps':
            v[name] = NPS_LIST[int(round(u))]
        else:
            v[name] = lo + u * (hi - lo)

    od, wall = SCH80[v['nps']]
    r_e = od / 2000.0
    r_i = r_e - wall / 1000.0
    N = max(2, v['N_lay'])

    # The span is pinned so that N_lay cannot move the bottom of the cascade.
    # T_m_bottom is an INTERFACE temperature -- the ORC sees it -- so letting
    # N_lay change it would leak a surface effect into a borehole study.
    DT_cascade = SPAN_FIXED * N / (N - 1)

    return T.Case(
        h_m=v['h_m_kJkg'] * 1000.0,
        k_s=0.60 * v['k_scale'], k_l=0.45 * v['k_scale'],
        rho_s=1550.0 * v['rho_scale'], rho_l=1450.0 * v['rho_scale'],
        cp_s=1280.0 * v['cp_scale'], cp_l=1800.0 * v['cp_scale'],
        L_ft=v['L_ft'], num_tubes=int(v['num_tubes']),
        r_i=r_i, r_e=r_e,
        num_fins=int(v['num_fins']),
        N_lay=N, n_segments=int(N * max(4, round(100 / N))),
        DT_cascade=DT_cascade,
        **FROZEN)


def geometry_ok(case):
    """The fins must fit inside the cell the leg owns.

    r_cell = R / sqrt(2 n_t), so packing three hairpins into a 7 in hole with
    2 in pipe and 15 mm fins puts the fin tips outside the cell. That is not a
    worse design, it is not a geometry.
    """
    return case.r_e + case.fin_L < case.r_cell


def evaluate(row):
    try:
        case = build_case(row)
        if not geometry_ok(case):
            return None
        if any(k == 'error' for k, _ in T.validate_case(case, verbose=False)):
            return None
        r = T.simulate_css_corrected(case)
        k = T.performance_indices(case, r)
        if not (k['f_pump'] < 0.5):
            return None
        return {n: float(k[n]) for n in KPIS if n in k}
    except Exception:
        return None


def run(centre_points=4, out='borehole_raw.json', verbose=True):
    """128 design points plus a few centre points to test for curvature.

    The centre points are the insurance. A two-level design ASSUMES the
    response is planar between the corners; running the midpoint a few times
    and comparing it against the plane the corners predict turns that
    assumption into a measurement. Four extra runs, under two minutes.
    """
    X = design_matrix()
    rows = []
    bar = None
    try:
        bar = Progress(len(X) + centre_points, 'borehole factorial')
    except Exception:
        pass
    for i, row in enumerate(X):
        res = evaluate(row)
        rows.append({'run': i, 'x': row.tolist(), 'centre': False, 'kpi': res})
        if bar:
            bar.tick(res is not None, f'corner {i + 1}/{len(X)}')
        if (i + 1) % 16 == 0:
            pathlib.Path(out).write_text(json.dumps(
                {'factors': list(FACTORS), 'kpis': KPIS, 'rows': rows,
                 'frozen': {k: v for k, v in FROZEN.items()},
                 'span_fixed': SPAN_FIXED}))
    for j in range(centre_points):
        res = evaluate(np.zeros(NF))
        rows.append({'run': len(X) + j, 'x': [0.0] * NF,
                     'centre': True, 'kpi': res})
        if bar:
            bar.tick(res is not None, f'centre {j + 1}/{centre_points}')
    if bar:
        bar.close()
    pathlib.Path(out).write_text(json.dumps(
        {'factors': list(FACTORS), 'kpis': KPIS, 'rows': rows,
         'frozen': {k: v for k, v in FROZEN.items()},
         'span_fixed': SPAN_FIXED}))
    nb = sum(1 for r in rows if r['kpi'] is None)
    if verbose:
        print(f'\n{len(rows)} runs, {nb} infeasible ({100 * nb / len(rows):.1f} %)')
    return rows


def analyse(path='borehole_raw.json'):
    """Main effects, two-factor interactions, and a curvature test.

    ESTIMATED BY LEAST SQUARES, NOT BY A DIFFERENCE OF AVERAGES. The first
    version of this function used the textbook shortcut

        effect_i = 2 * mean(x_i * y)

    which is exact -- and much easier to explain -- but ONLY when the design
    is balanced. The first real run lost 10 of 128 corners, and not at random:
    every one of them was a deep well with a single small pipe, where the
    pumping screen f_pump < 0.5 bites. A structured loss is precisely the case
    the shortcut cannot survive, because the surviving runs no longer carry
    each factor at its high level exactly half the time, so every effect picks
    up a share of every other effect.

    It was not a small error. On psi and eta_RTE the shortcut disagreed with
    least squares by more than 100 % OF THE LARGEST EFFECT IN THE TABLE -- it
    was not perturbing the ranking, it was inventing it.

    Least squares has no such requirement: it conditions on the runs that
    exist. The surviving design is still well behaved (largest correlation
    between any two of the 45 effect columns is 0.09, worst variance inflation
    1.18x), so the estimates stay sharp. The shortcut is kept alongside and
    reported as `naive_gap` so the discrepancy is visible rather than silent.
    """
    d = json.loads(pathlib.Path(path).read_text())
    names = d['factors']
    nf = len(names)
    corner = [r for r in d['rows'] if not r['centre'] and r['kpi']]
    centre = [r for r in d['rows'] if r['centre'] and r['kpi']]
    X = np.array([r['x'] for r in corner])
    pairs = list(itertools.combinations(range(nf), 2))
    labs = list(names) + [f'{names[i]} x {names[j]}' for i, j in pairs]
    M = np.column_stack([np.ones(len(X))]
                        + [X[:, i] for i in range(nf)]
                        + [X[:, i] * X[:, j] for i, j in pairs])
    inv = np.linalg.inv(M.T @ M)
    out = {}
    for kpi in d['kpis']:
        y = np.array([r['kpi'][kpi] for r in corner if kpi in r['kpi']])
        if len(y) != len(X):
            continue
        beta, *_ = np.linalg.lstsq(M, y, rcond=None)
        eff = 2.0 * beta[1:]
        mu = float(beta[0])
        res = y - M @ beta
        dof = max(1, len(y) - M.shape[1])
        sig = float(np.sqrt(res @ res / dof))
        se = 2.0 * sig * np.sqrt(np.diag(inv)[1:])
        naive = np.array([2.0 * np.mean(X[:, i] * y) for i in range(nf)])
        big = float(np.abs(eff).max()) or 1.0

        # CURVATURE. Compared against the fitted INTERCEPT, which is the
        # model's own prediction at the centre, not against the mean of the
        # corners -- with runs missing those are not the same number.
        #
        # Note what the centre points can and cannot do here. The model is
        # deterministic, so repeating the centre gives the identical answer
        # every time: the replicates measure no noise and four of them carry
        # exactly as much information as one. All they buy is the curvature
        # offset. Budget one next time.
        cy = [r['kpi'][kpi] for r in centre if kpi in r['kpi']]
        curv = (float(np.mean(cy)) - mu) if cy else float('nan')
        out[kpi] = {
            'effects': dict(zip(labs, map(float, eff))),
            'se': dict(zip(labs, map(float, se))),
            'mean': mu, 'curvature': curv,
            'curv_rel': curv / (float(np.std(y)) or 1.0),
            'curv_pct': 100.0 * curv / (abs(mu) or 1.0),
            'naive_gap': float(100.0 * np.abs(naive - eff[:nf]).max() / big),
            'n_corner': len(y), 'n_centre': len(cy)}
    return out, d


def top(res, kpi, n=12, kind='both'):
    """Ranked effects for one KPI, with the t-ratio that says which are real.

    `t` is the effect divided by its own standard error. Above about 2 the
    effect is larger than the model's own residual scatter can account for;
    below it, the row is telling you the factor did nothing.
    """
    import pandas as pd
    e, s_e = res[kpi]['effects'], res[kpi]['se']
    if kind == 'main':
        e = {k: v for k, v in e.items() if ' x ' not in k}
    elif kind == 'inter':
        e = {k: v for k, v in e.items() if ' x ' in k}
    s = pd.Series(e).sort_values(key=np.abs, ascending=False).head(n)
    return pd.DataFrame({
        'effect': s.round(5),
        '% of mean': (100 * s / res[kpi]['mean']).round(1),
        't': (s / pd.Series({k: s_e[k] for k in s.index})).abs().round(1)})


try:
    from doe_morris import Progress            # reuse the live progress bar
except Exception:                              # pragma: no cover
    Progress = None


if __name__ == '__main__':
    if verify_design():
        run()
