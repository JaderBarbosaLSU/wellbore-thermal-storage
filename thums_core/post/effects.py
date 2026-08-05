"""Factor effects computed from the design, not from group means.

The legacy code estimates effects as `df[df[f]==high][y].mean() - df[df[f]==low][y].mean()`.
That is correct only if the design is perfectly balanced and every run converged.
Because failures are written as NaN and `pandas.mean()` skips NaN silently, a single
failed run makes the contrast non-orthogonal without any message.

Here the contrast is computed from the coded design matrix, and an incomplete
design raises. Significance uses Lenth's pseudo standard error, which is the
standard tool for an unreplicated factorial -- the legacy code ranks effects by
`abs()` with no threshold at all.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import kpi as kpi_mod


def main_effects(df: pd.DataFrame, factors, response: str,
                 require_complete: bool = True) -> pd.DataFrame:
    """Orthogonal main-effect estimates from coded columns.

    Effect = 2 * beta, i.e. the change in response from low to high level.
    """
    coded = ["coded_" + f for f in factors]
    missing = [c for c in coded if c not in df.columns]
    if missing:
        raise KeyError(f"coded columns absent: {missing}. Effects must be computed "
                       "from the design, not inferred from physical values.")
    d = df[df.get("is_centre", False) == False] if "is_centre" in df else df  # noqa: E712
    if "converged" in d.columns:
        n_bad = int((~d["converged"]).sum())
        if n_bad:
            if require_complete:
                raise ValueError(
                    f"{n_bad} non-converged run(s) in the factorial portion. The "
                    "contrasts are no longer orthogonal. Re-run the failures, or "
                    "call with require_complete=False and report the bias."
                )
            print(f"[thums_core] WARNING: {n_bad} non-converged run(s) included; "
                  "contrasts are not orthogonal")
    X = d[coded].to_numpy(float)
    y = d[response].to_numpy(float)
    if not np.isfinite(y).all():
        raise ValueError(f"{response} contains non-finite values in the factorial runs")
    n = len(y)
    G = X.T @ X
    if not np.allclose(G, n * np.eye(len(factors))):
        raise ValueError("coded design is not orthogonal; cannot use simple contrasts")
    beta = X.T @ y / n
    return (pd.DataFrame({"factor": list(factors), "effect": 2.0 * beta})
            .assign(abs_effect=lambda t: t.effect.abs())
            .sort_values("abs_effect", ascending=False)
            .reset_index(drop=True))


def lenth(effects: np.ndarray, alpha: float = 0.05):
    """Lenth's pseudo standard error and the margin of error.

    Returns (PSE, ME, SME). |effect| > ME is significant at alpha individually;
    > SME is significant simultaneously across the family.
    """
    e = np.asarray(effects, float)
    m = len(e)
    s0 = 1.5 * np.median(np.abs(e))
    keep = np.abs(e) < 2.5 * s0
    pse = 1.5 * np.median(np.abs(e[keep])) if keep.any() else s0
    d = m / 3.0
    from scipy import stats
    me = stats.t.ppf(1 - alpha / 2, d) * pse
    gamma = (1 + 0.95 ** (1 / m)) / 2
    sme = stats.t.ppf(gamma, d) * pse
    return pse, me, sme


def screen(df: pd.DataFrame, factors, response: str, alpha: float = 0.05,
           require_complete: bool = True) -> pd.DataFrame:
    eff = main_effects(df, factors, response, require_complete)
    pse, me, sme = lenth(eff.effect.to_numpy())
    eff["pse"] = pse
    eff["ME"] = me
    eff["SME"] = sme
    eff["significant"] = eff.abs_effect > me
    eff["significant_simultaneous"] = eff.abs_effect > sme
    return eff


def curvature(df: pd.DataFrame, response: str) -> dict:
    """Factorial-vs-centre contrast.

    Note: for a deterministic simulator, replicated centre points give exactly zero
    pure error, so no valid F test can be built from them. This returns the contrast
    and says so, rather than manufacturing a p-value.
    """
    if "is_centre" not in df:
        raise KeyError("no is_centre column; the design carried no centre points")
    f = df[~df.is_centre][response].astype(float)
    c = df[df.is_centre][response].astype(float)
    nf, nc = len(f), len(c)
    diff = f.mean() - c.mean()
    ss = (nf * nc / (nf + nc)) * diff ** 2 if nf and nc else np.nan
    return {"corner_mean": f.mean(), "centre_mean": c.mean(),
            "difference": diff, "SS_curvature": ss, "n_corner": nf, "n_centre": nc,
            "note": "Deterministic model: centre replicates give zero pure error. "
                    "Assess against a lack-of-fit criterion, not an F test."}


def half_normal(eff: pd.DataFrame, ax=None):
    """Daniel half-normal plot -- the visual companion to Lenth."""
    import matplotlib.pyplot as plt
    from scipy import stats
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    a = np.sort(eff.abs_effect.to_numpy())
    m = len(a)
    q = stats.norm.ppf(0.5 + 0.5 * (np.arange(1, m + 1) - 0.5) / m)
    ax.plot(q, a, "o", ms=4)
    order = eff.sort_values("abs_effect").factor.tolist()
    for x, y, name in list(zip(q, a, order))[-5:]:
        ax.annotate(name, (x, y), fontsize=7, xytext=(4, -2), textcoords="offset points")
    if "ME" in eff:
        ax.axhline(eff.ME.iloc[0], ls="--", lw=1, color="crimson", label="Lenth ME")
        ax.axhline(eff.SME.iloc[0], ls=":", lw=1, color="crimson", label="Lenth SME")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("half-normal quantile")
    ax.set_ylabel("|effect|")
    return ax
