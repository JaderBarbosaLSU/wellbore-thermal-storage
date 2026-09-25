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
VERSION = "0.9"

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


if __name__ == '__main__':
    check_version()
    assemble()
    run('make_thums.py')
    run('build_student.py')
    run('build_notebook.py')
    print('\nall built.')
