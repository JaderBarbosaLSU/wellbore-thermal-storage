"""The results schema — one tidy row per run, identical for every study type.

Three column blocks:

    inputs   physical factor settings, plus coded_<factor> when part of a design
    kpis     exactly the names in thums_core.kpi.REGISTRY
    status   converged, n_iter, warnings, version, git_commit, study_id, run_id

Two rules make the schema worth having:

1. `converged` is mandatory. The legacy bisections return a bracket bound, a
   collapsed midpoint, or the last iterate after 100 iterations -- all
   indistinguishable from success in the output table.
2. Failed runs are WRITTEN, not dropped. A failure is a row with
   `converged=False` and NaN KPIs. Post-processing must then either exclude them
   explicitly and report the count, or refuse to compute a balanced contrast.
   `pandas.mean()` skipping NaN is how an unbalanced, non-orthogonal effect
   estimate reaches a figure without a warning.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np
import pandas as pd

from . import kpi as kpi_mod

STATUS_COLUMNS = ("converged", "n_iter", "residual", "warnings",
                  "version", "git_commit", "study_id", "run_id", "wall_time_s")


def git_commit(short: bool = True) -> str:
    """The commit this package was built from.

    Order matters, and the old implementation had it wrong. It ran `git` in the
    *current working directory*, which in Colab is `/content` -- not a repository,
    so the stamp read "unknown"; and in any other checkout it would have
    cheerfully reported that unrelated repository's HEAD. The version stamp is
    the thing that tells you whether Colab fetched the new notebook
    (Colab_GitHub_Token_Workflow section 8), so it has to describe the installed
    package, not the shell's location.
    """
    # 1. Installed from a git URL. pip records the exact resolved commit, which
    #    is what `pip install git+https://...` in a Colab cell produces.
    try:
        import json
        from importlib.metadata import distribution
        raw = distribution("thums_core").read_text("direct_url.json")
        if raw:
            cid = json.loads(raw).get("vcs_info", {}).get("commit_id")
            if cid:
                return cid[:7] if short else cid
    except Exception:
        pass
    # 2. Running from a source checkout: ask git about THIS file's repository,
    #    not the caller's working directory.
    try:
        import pathlib
        out = subprocess.run(
            ["git", "-C", str(pathlib.Path(__file__).resolve().parent),
             "rev-parse", "--short" if short else "HEAD", "HEAD"],
            capture_output=True, text=True, timeout=5)
        return out.stdout.strip().splitlines()[0] if out.stdout.strip() else "unknown"
    except Exception:
        pass
    return "unknown"


def git_ref() -> str:
    """The branch or tag that was installed, when pip recorded one.

    `pip install git+https://...@v0.2` stores "v0.2" here. Printing it next to
    the commit is what makes a stale Colab tab obvious: the tab shows the ref you
    asked for AND the commit it resolved to, so "v0.2" against an unexpected
    commit is visible immediately.
    """
    try:
        import json
        from importlib.metadata import distribution
        raw = distribution("thums_core").read_text("direct_url.json")
        if raw:
            return json.loads(raw).get("vcs_info", {}).get(
                "requested_revision") or ""
    except Exception:
        pass
    return ""


@dataclass
class RunRecord:
    inputs: dict
    kpis: dict = field(default_factory=dict)
    converged: bool = False
    n_iter: int = 0
    residual: float = float("nan")
    warnings: tuple = ()
    run_id: int = -1
    wall_time_s: float = float("nan")

    def to_row(self, version: str, commit: str, study_id: str) -> dict:
        row = dict(self.inputs)
        for name in kpi_mod.REGISTRY:
            row[name] = self.kpis.get(name, np.nan)
        unknown = set(self.kpis) - set(kpi_mod.REGISTRY)
        if unknown:
            raise KeyError(f"KPIs not in the registry: {sorted(unknown)}. "
                           "Declare them in thums_core/kpi.py first.")
        row.update(converged=self.converged, n_iter=self.n_iter,
                   residual=self.residual, warnings="; ".join(self.warnings),
                   version=version, git_commit=commit, study_id=study_id,
                   run_id=self.run_id, wall_time_s=self.wall_time_s)
        return row


def write(df: pd.DataFrame, manifest: dict, stem: str | Path) -> tuple[Path, Path]:
    """Write results and their manifest together. Never one without the other."""
    stem = Path(stem)
    csv, man = stem.with_suffix(".csv"), stem.with_suffix(".manifest.json")
    df.to_csv(csv, index=False)
    manifest = dict(manifest)
    manifest["n_rows"] = int(len(df))
    manifest["n_converged"] = int(df["converged"].sum())
    manifest["n_failed"] = int((~df["converged"]).sum())
    man.write_text(json.dumps(manifest, indent=2, default=str))
    return csv, man


def read(stem: str | Path) -> tuple[pd.DataFrame, dict]:
    stem = Path(stem)
    df = pd.read_csv(stem.with_suffix(".csv"))
    manifest = json.loads(stem.with_suffix(".manifest.json").read_text())
    return df, manifest


def converged_only(df: pd.DataFrame, *, require_all: bool = False) -> pd.DataFrame:
    """Filter to converged runs, loudly.

    require_all=True raises if anything failed -- use before computing a balanced
    contrast, where a missing cell destroys orthogonality.
    """
    n_bad = int((~df["converged"]).sum())
    if n_bad and require_all:
        raise ValueError(
            f"{n_bad} of {len(df)} runs did not converge. A balanced contrast "
            "cannot be computed from an incomplete design; re-run the failures "
            "or drop the affected factor."
        )
    if n_bad:
        print(f"[thums_core] excluding {n_bad} non-converged run(s) of {len(df)}")
    return df[df["converged"]].copy()
