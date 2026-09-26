"""Build everything, in dependency order, from one set of sources.

    model_part2.py + model_part3.py + thums_model_core.py
        |
        +--> thums_model.py        the engineering build (everything)
                |
                +--> thums.py      the live model only        (make_thums.py)
                |       |
                |       +--> P2H2P_model.ipynb                (build_student.py)
                |
                +--> P2H2P_verification.ipynb                 (build_notebook.py)

Both notebooks are generated from the same model sources, so a change cannot
land in one and miss the other. Run this, not the individual scripts.
"""
import pathlib
import re
import subprocess
import sys

# The model version, in ONE place. It has now gone stale twice in three
# commits -- notebook headers and module docstrings carrying 0.6 and 0.7
# after the physics had moved on. A version string that lags the code is
# the first thing anyone checks when two runs disagree, so the build
# refuses to proceed if the generators disagree with this.
VERSION = "0.10"

HERE = pathlib.Path(__file__).parent
BUILD = pathlib.Path('/tmp/build')


def assemble():
    core = (BUILD / 'thums_model_core.py').read_text()
    p2 = (HERE / 'model_part2.py').read_text()
    p3 = (HERE / 'model_part3.py').read_text()
    out = core + p2 + p3
    (BUILD / 'thums_model.py').write_text(out)
    print(f'thums_model.py  {len(out.splitlines())} lines')


def run(script):
    print(f'\n--- {script}')
    r = subprocess.run([sys.executable, str(HERE / script)],
                       capture_output=True, text=True)
    print(r.stdout.rstrip())
    if r.returncode:
        print(r.stderr[-3000:])
        raise SystemExit(f'{script} failed')


def check_version():
    """Every generator must agree with VERSION, or the build stops."""
    want = {
        'build_student.py': f'*Model version {VERSION} \u00b7 notebook built',
        'make_thums.py': f'The live model, v{VERSION}.',
    }
    bad = []
    for name, needle in want.items():
        txt = (HERE / name).read_text()
        if needle not in txt:
            found = re.search(r'v?(?:ersion )?0\.\d+[a-z]?', txt)
            bad.append(f'  {name}: expected {needle!r}, found '
                       f'{found.group(0) if found else "nothing"}')
    if bad:
        raise SystemExit('version strings are stale:\n' + '\n'.join(bad))
    print(f'version {VERSION}: generators agree')


def check_names(path):
    """Every name a notebook USES must be defined by an EARLIER cell.

    The verification notebook shipped v0.9 with a NameError in it:
    `cycle_state_points` gained a call to `feasible_rankine`, and the cell
    that emits it was never told to emit the new function too. Nothing in
    the build executed the notebook, so nothing noticed.

    Executing both notebooks in the build is too slow to be a gate. This is
    the cheap 95 %: walk the code cells in order, collect what each defines,
    and flag any global load that no earlier cell defined. It catches the
    whole class -- a `src(...)` list that forgot a dependency -- in about a
    tenth of a second.
    """
    import builtins
    import json
    import symtable

    def scan(tab, defined, needed):
        """Walk a symbol table and its nested scopes.

        `symtable` does the scope resolution, which is the whole point: a
        loop variable inside a function is local and must NOT be reported,
        while a function calling a name it never binds IS a global reference
        and must be. Hand-rolling this with ast.walk gets it wrong in both
        directions.
        """
        top = tab.get_type() == 'module'
        for s in tab.get_symbols():
            name = s.get_name()
            if top:
                if s.is_assigned() or s.is_imported():
                    defined.add(name)
                elif s.is_referenced():
                    needed.add(name)
            elif s.is_global() and s.is_referenced():
                needed.add(name)
        for child in tab.get_children():
            if child.get_type() == 'function':
                defined.add(child.get_name())
            scan(child, defined, needed)

    nb = json.loads(pathlib.Path(path).read_text())
    known = set(dir(builtins)) | {'__name__', '__file__', 'get_ipython',
                                  'display', 'In', 'Out'}
    missing = []

    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] != 'code':
            continue
        text = ''.join(cell['source'])
        if text.lstrip().startswith('%%'):
            continue                      # %%writefile and friends
        text = '\n'.join('' if ln.lstrip().startswith(('%', '!')) else ln
                         for ln in text.splitlines())
        try:
            tab = symtable.symtable(text, f'cell{i}', 'exec')
        except SyntaxError:
            continue
        defined, needed = set(), set()
        scan(tab, defined, needed)
        missing += [(i, n) for n in sorted(needed - known - defined)]
        known |= defined

    if missing:
        lines = [f'  cell {i}: {n}' for i, n in missing]
        raise SystemExit(f'{pathlib.Path(path).name}: {len(missing)} name(s) '
                         f'used before being defined:\n' + '\n'.join(lines))
    print(f'  {pathlib.Path(path).name}: names resolve')


if __name__ == '__main__':
    check_version()
    assemble()
    run('make_thums.py')
    run('build_student.py')
    run('build_notebook.py')
    print('\nchecking that every name resolves')
    check_names(BUILD / 'P2H2P_model.ipynb')
    check_names(BUILD / 'P2H2P_verification.ipynb')
    print('\nall built.')
