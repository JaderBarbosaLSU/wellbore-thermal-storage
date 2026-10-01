"""Build the BOREAS series. Separate from build_all.py ON PURPOSE.

P2H2P/THUMS is frozen at v0.16 and must stay runnable forever -- that was the
condition on starting a new series rather than cutting v0.17. Two build
scripts, two version stamps, two hash gates, no shared mutable state. If this
script can break `thums.py`, the separation has failed.

    boreas_part2.py + boreas_part3.py + thums_model_core.py
        |
        +--> boreas_model.py      the engineering build
                |
                +--> boreas.py    the live model only   (make_boreas.py)
"""
import hashlib
import os
import pathlib
import subprocess
import sys

VERSION = "1.0"
SERIES = "BOREAS"

HERE = pathlib.Path(__file__).parent
BUILD = pathlib.Path('/tmp/build')


def assemble():
    core = (BUILD / 'thums_model_core.py').read_text()
    p2 = (HERE / 'boreas_part2.py').read_text()
    p3 = (HERE / 'boreas_part3.py').read_text()
    out = core + p2 + p3
    (BUILD / 'boreas_model.py').write_text(out)
    print(f'boreas_model.py  {len(out.splitlines())} lines')


def check_version():
    """Same physics-hash gate as THUMS, on its own sources and its own stamp.

    The gate exists because a version string that lags the code got past a
    pure consistency check three times in the THUMS line: all the generators
    agreed, on a stale value. Hashing the sources gives it an independent
    input. Starting a new series does not earn an exemption from the lesson.
    """
    src = b''.join(sorted((HERE / f).read_bytes()
                          for f in ('boreas_part2.py', 'boreas_part3.py')))
    digest = hashlib.sha256(src).hexdigest()[:16]
    stamp = HERE / '.boreas_version_hash'
    if os.environ.get('BOREAS_STAMP_VERSION'):
        stamp.write_text(f'{VERSION} {digest}\n')
        print(f'recorded {SERIES} {VERSION} at physics hash {digest}')
        return
    if stamp.exists():
        was_v, was_d = stamp.read_text().split()
        if was_d != digest and was_v == VERSION:
            raise SystemExit(
                f'the BOREAS sources changed but VERSION is still {VERSION!r}.'
                f'\n  recorded {was_d}, current {digest}'
                f'\n  Bump VERSION, then re-run with BOREAS_STAMP_VERSION=1.')
    print(f'{SERIES} v{VERSION}: physics hash {digest}')


def check_thums_untouched():
    """THE POINT OF THE WHOLE EXERCISE, so it is checked rather than assumed.

    The new series exists so that P2H2P_model keeps running. That promise is
    worth exactly as much as the test behind it, and `cp` plus a careless sed
    is how it would quietly be broken.
    """
    stamp = HERE / '.version_hash'
    if not stamp.exists():
        print('  (no THUMS stamp to compare against)')
        return
    want_v, want_d = stamp.read_text().split()
    src = b''.join(sorted((HERE / f).read_bytes()
                          for f in ('model_part2.py', 'model_part3.py')))
    got = hashlib.sha256(src).hexdigest()[:16]
    if got != want_d:
        raise SystemExit(
            f'THUMS sources have MOVED: stamped {want_d} at v{want_v}, '
            f'now {got}.\n  BOREAS must not edit the THUMS line. '
            f'Revert model_part2.py / model_part3.py.')
    print(f'  THUMS v{want_v} untouched ({got})')


def run(script):
    print(f'\n--- {script}')
    r = subprocess.run([sys.executable, str(HERE / script)],
                       capture_output=True, text=True)
    print(r.stdout.rstrip())
    if r.returncode:
        print(r.stderr[-3000:])
        raise SystemExit(f'{script} failed')


if __name__ == '__main__':
    check_version()
    print('\nchecking the THUMS line is untouched')
    check_thums_untouched()
    assemble()
    run('make_boreas.py')
    print('\nbuilt.')
