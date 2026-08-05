"""Regenerate thums_core/_legacy_physics.py from the v0.1 notebook.

The legacy physics is the frozen reference the unified package is validated
against. It is extracted mechanically rather than copied by hand so that
"the reference" cannot quietly drift from what the conference paper ran.

Usage:
    python tools/extract_legacy.py [notebook.ipynb] [--last-cell 19]

Only definition cells are taken (up to --last-cell); driver cells are excluded,
because the drivers are what the package replaces.
"""

import argparse
import json
from pathlib import Path

DEFAULT_NB = "thums_core/THUMS_Multilayer_BB_param_Jan_12b_converg_improved.ipynb"
OUT = Path("thums_core/_legacy_physics.py")

HEADER = '''"""Auto-generated from {nb} (v0.1).

DO NOT EDIT. Regenerate with tools/extract_legacy.py.
This is the frozen reference physics: the unified package is validated against it.

Known defects retained deliberately (see docs/THUMS_code_audit_2026-08-05.md):
  * the charging call chain has no steel conductivity (k_m_l passed where k_w belongs)
  * h_e = 1e9 on the delta -> 0 branch, commented as "large resistance" (it is the opposite)
  * melt-front and bisection failures return 0.0 / a bracket bound silently
They are reproduced here so that step 1 of the migration is a pure refactor.
"""
# flake8: noqa
'''


def main():
    p = argparse.ArgumentParser()
    p.add_argument("notebook", nargs="?", default=DEFAULT_NB)
    p.add_argument("--last-cell", type=int, default=19)
    p.add_argument("--out", default=str(OUT))
    a = p.parse_args()

    nb = json.loads(Path(a.notebook).read_text())
    parts = [HEADER.format(nb=Path(a.notebook).name)]
    n = 0
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] != "code":
            continue
        src = "".join(c["source"])
        if src.lstrip().startswith("!"):
            continue
        if i > a.last_cell:
            break
        parts.append(f"# ---------- notebook cell {i} ----------")
        parts.append(src)
        parts.append("")
        n += 1
    Path(a.out).write_text("\n".join(parts))
    print(f"wrote {a.out}: {n} cells, {len(''.join(parts).splitlines())} lines")


if __name__ == "__main__":
    main()
