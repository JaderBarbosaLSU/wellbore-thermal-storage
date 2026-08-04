"""Freeze the notebook's configuration constants, and fail if they drift.

Part I of the report and its tables are generated: every number is read out of
the notebook or the fixtures at build time, so they cannot go stale.

Part II is different. The model sections quote configuration constants inline,
inside prose and derivations -- "the damper floor is 0.05", "the coil is sized
at 3000 W/K per kg/s" -- and prose cannot be templated without becoming
unreadable. Those numbers are therefore typed, which means they can rot.

This script is the guard. It snapshots every scalar assignment in the
configuration cell into config_fixture.json and, on every later build, compares
the live notebook against that snapshot. A changed constant does not silently
contradict Part II: the build stops and names the symbol, so the prose can be
re-checked and the fixture re-frozen deliberately.

    python3 check_config.py           # compare (exit 1 on drift)
    python3 check_config.py --freeze  # re-record after a deliberate change
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NB = os.path.join(ROOT, 'geothermal_cooling', 'dbhe_master.ipynb')
FIX = os.path.join(HERE, 'config_fixture.json')

# scalar assignment at column 0: NAME = <number>, optional trailing comment
ASSIGN = re.compile(r'^([A-Za-z_]\w*)\s*=\s*'
                    r'(-?(?:\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?))'
                    r'\s*(?:#.*)?$', re.M)


def config_scalars():
    """Every top-level numeric constant in the configuration cell."""
    nb = json.load(open(NB))
    ver = 'unknown'
    for c in nb['cells']:
        m = re.search(r'NB_VERSION\s*=\s*"([^"]+)"', ''.join(c['source']))
        if m:
            ver = m.group(1)
            break
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        src = ''.join(c['source'])
        # the configuration cell is the one carrying the lettered sections
        if 'SECTION A' not in src.upper() and 'SECT' not in src[:400].upper():
            continue
        vals = {m.group(1): float(m.group(2)) for m in ASSIGN.finditer(src)}
        if len(vals) > 20:
            return ver, vals
    raise SystemExit('check_config: configuration cell not found')


def freeze():
    ver, vals = config_scalars()
    json.dump({'notebook_version': ver, 'constants': vals},
              open(FIX, 'w'), indent=1, sort_keys=True)
    print(f'check_config: froze {len(vals)} constants at v{ver} -> '
          f'{os.path.basename(FIX)}')


def compare():
    if not os.path.exists(FIX):
        print('check_config: no fixture yet; run with --freeze', file=sys.stderr)
        return 1
    ref = json.load(open(FIX))
    ver, vals = config_scalars()
    old = ref['constants']

    changed = {k: (old[k], vals[k]) for k in old
               if k in vals and old[k] != vals[k]}
    removed = sorted(set(old) - set(vals))
    added = sorted(set(vals) - set(old))

    if not (changed or removed):
        note = f'  (+{len(added)} new: {", ".join(added)})' if added else ''
        print(f'check_config: {len(old)} constants unchanged since '
              f'v{ref["notebook_version"]}{note}')
        return 0

    print(f'check_config: CONFIGURATION DRIFT  (fixture v'
          f'{ref["notebook_version"]} -> notebook v{ver})', file=sys.stderr)
    for k, (a, b) in sorted(changed.items()):
        print(f'    {k:24s} {a!r:>16}  ->  {b!r}', file=sys.stderr)
    for k in removed:
        print(f'    {k:24s} {old[k]!r:>16}  ->  REMOVED', file=sys.stderr)
    print('\n  Part II quotes these constants in prose. Re-read the model '
          'sections that\n  mention them, then re-freeze with '
          '`python3 check_config.py --freeze`.', file=sys.stderr)
    return 1


if __name__ == '__main__':
    if '--freeze' in sys.argv:
        freeze()
    else:
        sys.exit(compare())
