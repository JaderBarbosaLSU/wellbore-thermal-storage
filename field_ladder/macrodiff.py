"""Every macro whose value moved, against the previous commit.

Each return note has claimed a macro diff and v2.4/v2.5 got it wrong --
six values had moved and the note said none had. Claiming it by hand is
the problem; this computes it.

    python3 macrodiff.py                 # against HEAD
    python3 macrodiff.py <git-ref>       # against anything else
"""
import re
import subprocess
import sys
import pathlib

FILES = ['results_pack_v2/boreas_numbers.tex',
         'results_pack_v2/boreas_profile_numbers.tex',
         'results_pack_v2/boreas_orange_numbers.tex',
         'results_pack_v2/boreas_water_numbers.tex',
         'results_pack_v2/boreas_ts_numbers.tex',
         'results_pack_v2/boreas_video_numbers.tex',
         'par_pack/par_numbers.tex']
REF = sys.argv[1] if len(sys.argv) > 1 else 'HEAD'
REPO = sys.argv[2] if len(sys.argv) > 2 else '/tmp/w/repo'
PAT = re.compile(r'\\newcommand\{\\(\w+)\}\{(.*?)\}')


def parse(text):
    return dict(PAT.findall(text or ''))


def old(path):
    r = subprocess.run(['git', '-C', REPO, 'show', f'{REF}:{path}'],
                       capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else ''


def live(path):
    # the working copy FIRST, by full relative path -- matching on the
    # bare filename found the repository's copy and reported no change
    for root in ('.', '/tmp/par', REPO):
        p = pathlib.Path(root, path)
        if p.exists():
            return p.read_text()
        p = pathlib.Path(root, pathlib.Path(path).name)
        if p.exists():
            return p.read_text()
    return ''


tot_c = tot_a = tot_r = 0
for f in FILES:
    o, n = parse(old(f)), parse(live(f))
    if not o and not n:
        continue
    chg = [(k, o[k], n[k]) for k in sorted(set(o) & set(n)) if o[k] != n[k]]
    add = sorted(set(n) - set(o))
    rem = sorted(set(o) - set(n))
    tot_c += len(chg); tot_a += len(add); tot_r += len(rem)
    print(f'\n{f}  {len(o)} -> {len(n)}: '
          f'{len(chg)} changed, {len(add)} added, {len(rem)} removed')
    for k, a, b in chg:
        print(f'    \\{k:30s}{a:>14s} -> {b}')
    if rem:
        print('    REMOVED: ' + ', '.join(rem))
    if add:
        print('    added: ' + ', '.join(add))
print(f'\nTOTAL against {REF}: {tot_c} changed, {tot_a} added, '
      f'{tot_r} removed')
if tot_r:
    print('A REMOVED macro will break the paper build -- check each one.')
