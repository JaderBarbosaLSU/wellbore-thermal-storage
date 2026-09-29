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
import hashlib
import os
import pathlib
import re
import subprocess
import sys

# The model version, in ONE place. It has now gone stale twice in three
# commits -- notebook headers and module docstrings carrying 0.6 and 0.7
# after the physics had moved on. A version string that lags the code is
# the first thing anyone checks when two runs disagree, so the build
# refuses to proceed if the generators disagree with this.
VERSION = "0.14"

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
    """Every generator must agree with VERSION -- AND VERSION must be current.

    The first half of this check has existed since v0.7 and it has now failed
    to catch the same thing THREE times: v0.9 physics shipped stamped 0.8,
    and v0.14 and v0.14a both shipped stamped 0.13. Each time the check
    passed, because all the generators agreed -- on a stale value.

    That is a guard reading only its own inputs, which is the failure mode
    DN-27 describes for the Lorenz reference. Agreement between three strings
    I edit by hand can only catch disagreement; it can never catch the case
    where I changed the physics and forgot to touch any of them.

    So the second half compares a HASH of the physics sources against the
    hash recorded when VERSION was last set. If the sources moved and VERSION
    did not, the build stops. That is an independent input, and it is what
    makes this a gate rather than a consistency check.

    To cut a version: bump VERSION, run with THUMS_STAMP_VERSION=1 to record
    the new hash, and commit `.version_hash` alongside.
    """
    src = b''.join(sorted(
        (HERE / f).read_bytes() for f in ('model_part2.py', 'model_part3.py')))
    digest = hashlib.sha256(src).hexdigest()[:16]
    stamp = HERE / '.version_hash'
    if os.environ.get('THUMS_STAMP_VERSION'):
        stamp.write_text(f'{VERSION} {digest}\n')
        print(f'recorded {VERSION} at physics hash {digest}')
    elif stamp.exists():
        was_v, was_d = stamp.read_text().split()
        if was_d != digest and was_v == VERSION:
            raise SystemExit(
                f'the physics sources changed but VERSION is still '
                f'{VERSION!r}.\n'
                f'  recorded hash {was_d}, current {digest}\n'
                f'  Bump VERSION (and the two generator strings), then '
                f're-run with THUMS_STAMP_VERSION=1.\n'
                f'  This is the check that v0.9, v0.14 and v0.14a all got '
                f'past by agreeing on a stale value.')
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
    import ast
    import builtins
    import json
    import symtable

    def loop_only_names(text):
        """Module-level for/comprehension targets never bound by `=`.

        These are SCRATCH, not API. A plotting cell writing

            for c_, k in zip(cmap, ks):

        leaks `k` into the notebook namespace as an int, and a later cell
        writing k['UA_total_kW_K'] then RESOLVES statically and dies at
        run time with "'int' object is not subscriptable". That is exactly
        what shipped in v0.11. Static name resolution cannot see the type,
        but it can decline to treat a loop variable as an export.
        """
        try:
            tree = ast.parse(text)
        except SyntaxError:
            return set()
        targets, assigned = set(), set()
        SCOPED = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef,
                  ast.Lambda)

        def names(node):
            return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}

        def visit(node):
            """Recurse through MODULE scope, not descending into functions.

            The loop that leaked `k` was nested inside another loop, so
            looking only at tree.body missed it -- which is why the first
            version of this check still passed the bug.
            """
            nonlocal targets, assigned
            for child in ast.iter_child_nodes(node):
                if isinstance(child, SCOPED):
                    continue                        # its own scope
                if isinstance(child, (ast.For, ast.AsyncFor)):
                    targets |= names(child.target)
                elif isinstance(child, ast.comprehension):
                    targets |= names(child.target)
                elif isinstance(child, (ast.Assign, ast.AnnAssign,
                                        ast.AugAssign)):
                    tgt = child.targets if isinstance(child, ast.Assign) \
                        else [child.target]
                    for t in tgt:
                        assigned |= names(t)
                visit(child)

        visit(tree)
        return targets - assigned

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
        known |= defined - loop_only_names(text)

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
