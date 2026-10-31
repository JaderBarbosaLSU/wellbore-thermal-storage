"""Provenance stamps for generated files (handoff 6.1).

v2.2 shipped a tab_ladder.tex that was byte-for-byte the v2.1 file, while
the return note said every derived table had been recomputed and the
README said nothing in the pack is typed by hand. Both statements were
false and nothing caught it: a stale LaTeX table is still valid LaTeX.

Every generator now writes the SHA-256 of the pack_runs.json it was built
from into the file, and V16 checks every stamp against the current cache.
A file that was not regenerated now fails a numbered check instead of
quietly disagreeing with the CSV in the third digit.
"""
import hashlib
import os
import re

CACHE = 'pack_runs.json'
MARK = 'boreas-pack-sha256'


def pack_hash(path=CACHE):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]


def stamp(comment='%', path=CACHE):
    """One comment line naming the cache this file was built from."""
    return f'{comment} {MARK}: {pack_hash(path)}\n'


def read_stamp(path):
    try:
        with open(path, 'r', errors='ignore') as f:
            head = f.read(4096)
    except OSError:
        return None
    m = re.search(MARK + r':\s*([0-9a-f]+)', head)
    return m.group(1) if m else None


def check(dirname, patterns=('.tex',), path=CACHE):
    """(stale, unstamped, ok) lists for every generated file in `dirname`."""
    want = pack_hash(path)
    stale, unstamped, ok = [], [], []
    for fn in sorted(os.listdir(dirname)):
        if not fn.endswith(patterns):
            continue
        got = read_stamp(os.path.join(dirname, fn))
        if got is None:
            unstamped.append(fn)
        elif got != want:
            stale.append((fn, got))
        else:
            ok.append(fn)
    return stale, unstamped, ok
