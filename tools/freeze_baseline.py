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
from thums_core.config import BASELINE, LEGACY
from thums_core.errors import ThumsError
from thums_core.system import run_cycle

NOTE_LEGACY = (
    "Legacy physics, defects intact: charging uses PCM liquid conductivity as "
    "the wall/fin conductivity, delta_max exceeds the borehole radius, melt "
    "fronts are unbounded. This fixture records behaviour, not correctness.")

NOTE_MARCHED = (
    "v0.2 marched front. The melt front advances from the heat the resistance "
    "network delivers (rho*h_m*dA/dt = q'), so energy closure is structural; "
    "the melt state carries from charging into discharging; and the well count "
    "is max(heat transfer, inventory). Still present and NOT fixed here: the "
    "charging path uses PCM liquid conductivity as k_wall "
    "(charge_uses_wall_conductivity=False), delta_max exceeds the borehole "
    "radius, and there is no formation heat loss, no natural convection in the "
    "melt, and no sensible heat in the PCM. Compare against tag v0.2, which "
    "holds the legacy values this file replaced.")

# (DT_3C_2C, N_lay) points that converge under the legacy physics.
POINTS = [(55.0, 9), (40.0, 12), (30.0, 20)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out",
                    default="verification/reference/thums_v0_baseline.json")
    ap.add_argument("--front", choices=("marched", "legacy"), default="marched",
                    help="melt-front formulation; 'legacy' reproduces v0.1")
    a = ap.parse_args()
    base = BASELINE if a.front == "marched" else LEGACY

    fixture = {
        "schema": 1,
        "generated_by": "tools/freeze_baseline.py",
        "thums_version": thums_core.__version__,
        "git_commit": thums_core.results.git_commit(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "front": a.front,
        "note": (NOTE_MARCHED if a.front == "marched" else NOTE_LEGACY),
        "points": [],
    }

    for DT, N_lay in POINTS:
        case = base.with_(DT_3C_2C=DT, N_lay=N_lay,
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
