"""Freeze the v0.1 reference results BEFORE any physics change.

Step 0 of the migration. Every later port step re-runs this and compares; a
number that moves is attributable to the one module that was swapped.

The fixture is deliberately generated with `charge_uses_wall_conductivity=False`
and `strict=False`, i.e. with the legacy defects intact. It is a record of what
the code did, not of what is right.

Usage:
    python tools/freeze_baseline.py [-o verification/reference/thums_v0_baseline.json]
"""

import argparse
import json
import platform
import sys

import numpy as np

import thums_core
from thums_core.config import BASELINE
from thums_core.errors import ThumsError
from thums_core.system import run_cycle

# (DT_3C_2C, N_lay) points that converge under the legacy physics.
POINTS = [(55.0, 9), (40.0, 12), (30.0, 20)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out",
                    default="verification/reference/thums_v0_baseline.json")
    a = ap.parse_args()

    fixture = {
        "schema": 1,
        "generated_by": "tools/freeze_baseline.py",
        "thums_version": thums_core.__version__,
        "git_commit": thums_core.results.git_commit(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "note": ("Legacy physics, defects intact: charging uses PCM liquid "
                 "conductivity as the wall/fin conductivity, delta_max exceeds "
                 "the borehole radius, melt fronts are unbounded. This fixture "
                 "records behaviour, not correctness."),
        "points": [],
    }

    for DT, N_lay in POINTS:
        case = BASELINE.with_(DT_3C_2C=DT, N_lay=N_lay,
                              N_wells_bracket=(1.0, 400.0))
        entry = {"inputs": {"DT_3C_2C": DT, "N_lay": N_lay},
                 "case_key": case.key()}
        try:
            r = run_cycle(case)
            entry.update(converged=True, kpis=r.kpis,
                         detail={k: v for k, v in r.detail.items()
                                 if isinstance(v, (int, float))},
                         warnings=r.warnings)
            print(f"DT={DT:5.1f} N_lay={N_lay:3d}  N_wells={r.kpis['N_wells']:8.4f}  "
                  f"eta_rte={r.kpis['eta_rte']:.6f}  eps_pcm={r.kpis['eps_pcm']:.6f}")
        except ThumsError as e:
            entry.update(converged=False, error=f"{type(e).__name__}: {e}")
            print(f"DT={DT:5.1f} N_lay={N_lay:3d}  FAILED: {e}")
        fixture["points"].append(entry)

    with open(a.out, "w") as f:
        json.dump(fixture, f, indent=2, default=float)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
