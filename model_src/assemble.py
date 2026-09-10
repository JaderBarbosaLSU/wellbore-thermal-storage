"""Extract the model functions from the repo into one flat, self-contained module.

Extraction rather than retyping, so the numerics cannot shift. Cleaning is limited
to removing stray `print(` warnings and the duplicated mid-file imports; no
expression is rewritten.

Writes /tmp/build/thums_model_core.py
"""
import pathlib
import re

SRC = pathlib.Path('/tmp/src/thums_core/_legacy_physics.py').read_text()
OUT = pathlib.Path('/tmp/build')
OUT.mkdir(exist_ok=True)

WANT = [
    # (function name, section label)
    ('double_stage_rankine', 'cycles'),
    ('two_stage_htheatpump_2regs', 'cycles'),
    ('calculate_pressure_drop', 'hydraulics'),
    ('evaluate_Q_ratio_ch', 'solvers'),
    ('evaluate_Q_ratio_dc', 'solvers'),
    ('find_N_wells_for_Q_ratio_ch_fast', 'solvers'),
    ('find_m_dot_d_well1_for_Q_ratio_dc_fast', 'solvers'),
    ('compute_delta2_fast', 'front_closed_form'),
    ('compute_T_re', 'network'),
    ('compute_U_i', 'network'),
    ('make_cp_state', 'properties'),
    ('get_props_state', 'properties'),
    ('compute_h_i_from_state', 'properties'),
    ('temperature_profile_melt', 'march_legacy'),
    ('time_profiles_melt', 'march_legacy'),
]

# index every top-level def
starts = [(m.start(), m.group(1)) for m in re.finditer(r'^def (\w+)\(', SRC, re.M)]
bounds = {}
for i, (pos, name) in enumerate(starts):
    end = starts[i + 1][0] if i + 1 < len(starts) else len(SRC)
    bounds[name] = (pos, end)


def grab(name):
    a, b = bounds[name]
    body = SRC[a:b]
    lines = []
    for ln in body.splitlines():
        s = ln.strip()
        # drop the notebook cell markers and the stray re-imports
        if s.startswith('# ---------- notebook cell'):
            continue
        if re.match(r'^(import|from)\s', s) and 'CoolProp' not in s or \
           s in ('import numpy as np', 'import pandas as pd',
                 'import matplotlib.pyplot as plt'):
            continue
        # drop bare warning prints, keep everything computational
        if re.match(r'^print\(["\']Warning', s) or re.match(r'^#\s*print\(', s):
            continue
        lines.append(ln)
    return '\n'.join(lines).rstrip() + '\n'


chunks = {}
for name, sec in WANT:
    chunks.setdefault(sec, []).append((name, grab(name)))

HEADER = '''"""THUMS model -- flat, self-contained.

Assembled from the repository source by docs/assemble.py. The functions in the
`cycles`, `network`, `front_closed_form`, `march_legacy` and `hydraulics`
sections are the ORIGINAL v0.1 notebook functions, extracted verbatim; only
warning prints and duplicated imports were removed. Nothing was rewritten, so
this module reproduces the conference-paper numbers exactly.

The `front_energy_balance`, `march`, `sizing` and `run` sections are the v0.2
corrections and are new code.
"""
import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from math import log, pi
import pandas as pd
'''

parts = [HEADER]
order = ['properties', 'cycles', 'network', 'front_closed_form',
         'march_legacy', 'hydraulics', 'solvers']
for sec in order:
    parts.append(f'\n\n# {"=" * 74}\n# {sec}\n# {"=" * 74}\n')
    for name, body in chunks[sec]:
        parts.append('\n' + body)

path = OUT / 'thums_model_core.py'
path.write_text('\n'.join(parts))
n = len(path.read_text().splitlines())
print(f'wrote {path}  ({n} lines)')
for sec in order:
    print(f'  {sec:20s} {[n for n, _ in chunks[sec]]}')
