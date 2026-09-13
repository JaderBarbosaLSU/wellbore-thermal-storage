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
import subprocess
import sys

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


if __name__ == '__main__':
    assemble()
    run('make_thums.py')
    run('build_student.py')
    run('build_notebook.py')
    print('\nall built.')
