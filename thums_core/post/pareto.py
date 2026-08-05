"""Pareto fronts -- the analysis the legacy `Pareto_for_THUMS_DoE` notebook is named
after but does not contain.

That notebook has no dominance test anywhere in its 181 lines; "Pareto" there refers
to a Pareto *bar chart* of effect magnitudes. The scalarisation it uses instead --
`weighted_RTE = RTE_real * m_m_eff` -- collapses a genuine two-objective trade-off
into one number with an unstated weighting. These functions do it properly, using
the `direction` field of the KPI registry.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import kpi as kpi_mod


def _signed(df: pd.DataFrame, kpis) -> np.ndarray:
    """Return a matrix where larger is always better."""
    cols = []
    for n in kpis:
        s = df[n].to_numpy(dtype=float)
        cols.append(s if kpi_mod.REGISTRY[n].direction == kpi_mod.MAX else -s)
    return np.column_stack(cols)


def is_dominated(df: pd.DataFrame, kpis) -> np.ndarray:
    """Boolean mask: True where a row is dominated by at least one other row.

    Rows containing NaN in any objective are marked dominated (they cannot be
    shown to be non-dominated, and silently keeping them would put failed runs
    on the front).
    """
    Y = _signed(df, kpis)
    n = len(Y)
    bad = ~np.isfinite(Y).all(axis=1)
    dom = bad.copy()
    for i in range(n):
        if bad[i]:
            continue
        better_or_equal = np.all(Y >= Y[i], axis=1)
        strictly_better = np.any(Y > Y[i], axis=1)
        cand = better_or_equal & strictly_better & ~bad
        dom[i] = bool(cand.any())
    return dom


def front(df: pd.DataFrame, kpis) -> pd.DataFrame:
    """Non-dominated subset, sorted by the first objective."""
    kpis = list(kpis)
    missing = [k for k in kpis if k not in kpi_mod.REGISTRY]
    if missing:
        raise KeyError(f"not in the KPI registry: {missing}")
    kpi_mod.check_independent(kpis)      # refuse fronts over dependent objectives
    out = df[~is_dominated(df, kpis)].copy()
    return out.sort_values(kpis[0], ascending=kpi_mod.REGISTRY[kpis[0]].direction == kpi_mod.MIN)


def plot_tradeoff(df, kx, ky, ax=None, colour_by=None, annotate=None):
    """Scatter of all runs with the non-dominated front drawn through them."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(5.5, 4.2))
    f = front(df, [kx, ky])
    sc = ax.scatter(df[kx], df[ky], s=14, alpha=0.35,
                    c=None if colour_by is None else df[colour_by],
                    label="dominated", zorder=1)
    ax.plot(f[kx], f[ky], "-o", ms=5, lw=1.4, color="crimson",
            label="Pareto front", zorder=3)
    if annotate:
        for _, r in f.iterrows():
            ax.annotate(f"{r[annotate]:g}", (r[kx], r[ky]),
                        fontsize=6, xytext=(3, 3), textcoords="offset points")
    R = kpi_mod.REGISTRY
    ax.set_xlabel(f"{R[kx].label}  [{R[kx].unit}]")
    ax.set_ylabel(f"{R[ky].label}  [{R[ky].unit}]")
    ax.legend(frameon=False, fontsize=8)
    if colour_by is not None:
        plt.colorbar(sc, ax=ax, label=colour_by)
    return ax
