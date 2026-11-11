"""Guard for the student snapshot: every code cell must parse, and every
GITREF assignment in a notebook must name the SAME ref.

Born from two failures in the v1.0 snapshot (2026-10-31):

  * BOREAS_model.ipynb carried a dangling string in its pack-fetch cell and
    was a SyntaxError.  V18 hashes the EMBEDDED MODULES and said nothing
    about the notebook's own cells, so the freeze captured the break.
  * the model notebook assigns GITREF in TWO cells; the pin to the tag
    reached only one of them, so a snapshot advertised as frozen fetched
    half its inputs from the live branch.

Run it on the student branch before any push.  Usage:

    python3 student_check.py notebooks/*.ipynb [--expect REF]
"""
import ast
import json
import re
import sys

ASSIGN = re.compile(r"^\s*(?:[A-Za-z_][\w.]*\.)?GITREF\s*=\s*['\"]([^'\"]+)['\"]")
MAGIC = re.compile(r'^\s*[!%]')


def check(path, expect=None):
    nb = json.load(open(path))
    bad, refs = [], {}
    for i, c in enumerate(nb.get('cells', [])):
        if c.get('cell_type') != 'code':
            continue
        lines = ''.join(c['source']).split('\n')
        for ln in lines:
            m = ASSIGN.match(ln)
            if m:
                refs.setdefault(m.group(1), []).append(i)
        if any(MAGIC.match(l) for l in lines):
            continue            # IPython magic: not Python, not ours to parse
        try:
            ast.parse('\n'.join(lines))
        except SyntaxError as e:
            bad.append((i, e.lineno, (e.text or '').strip()))

    ok = True
    name = path.split('/')[-1]
    for i, ln, txt in bad:
        print(f'  {name}: cell {i} line {ln} SyntaxError: {txt}')
        ok = False
    if len(refs) > 1:
        print(f'  {name}: GITREF disagrees across cells -> '
              + ', '.join(f'{r!r} in cells {c}' for r, c in refs.items()))
        ok = False
    elif not refs:
        print(f'  {name}: no GITREF assignment found')
        ok = False
    elif expect and list(refs) != [expect]:
        print(f'  {name}: GITREF is {list(refs)[0]!r}, expected {expect!r}')
        ok = False
    if ok:
        print(f'  {name}: {len(nb["cells"])} cells parse, '
              f'GITREF = {list(refs)[0]!r} in cells {refs[list(refs)[0]]}')
    return ok


def main(argv):
    expect = None
    if '--expect' in argv:
        k = argv.index('--expect')
        expect = argv[k + 1]
        argv = argv[:k] + argv[k + 2:]
    print('student snapshot check')
    good = all([check(p, expect) for p in argv])
    print('PASS' if good else 'FAIL')
    return 0 if good else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
