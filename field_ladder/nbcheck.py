"""Is the notebook's embedded code the code the pack was built from?

BOREAS_model.ipynb carries boreas.py, boreas_plots.py and
boreas_ground.py inside %%writefile cells so Colab needs nothing from
the repository. The cost of that is a second copy of every module, and a
copy goes stale the moment a module is edited without rebuilding the
notebook.

This has now happened four times in this project. Each time the
notebook still ran, still produced plausible figures, and still
disagreed with the pack -- a stale copy is valid Python, exactly as a
stale table is valid LaTeX. Only a comparison catches it, so V18 makes
one.
"""
import hashlib
import json
import pathlib

NB = 'BOREAS_model.ipynb'
MODULES = ('boreas.py', 'boreas_plots.py', 'boreas_ground.py')


def embedded(nb_path=NB):
    nb = json.load(open(nb_path))
    out = {}
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        s = ''.join(c['source'])
        if s.startswith('%%writefile'):
            out[s.splitlines()[0].split()[1]] = '\n'.join(s.splitlines()[1:])
    return out


def digest(text):
    return hashlib.sha256(text.rstrip().encode()).hexdigest()[:12]


def check(nb_path=NB, root='.'):
    """(stale, missing, ok) -- stale entries are (name, nb_hash, live_hash)."""
    emb = embedded(nb_path)
    stale, missing, ok = [], [], []
    for fn in MODULES:
        live = pathlib.Path(root, fn)
        if not live.exists():
            continue
        if fn not in emb:
            missing.append(fn)
            continue
        a, b = digest(emb[fn]), digest(live.read_text())
        (ok if a == b else stale).append(fn if a == b else (fn, a, b))
    return stale, missing, ok


if __name__ == '__main__':
    st, ms, okk = check()
    for fn, a, b in st:
        print(f'STALE    {fn}: notebook {a}, live {b}')
    for fn in ms:
        print(f'MISSING  {fn} is not embedded')
    for fn in okk:
        print(f'current  {fn}')
    raise SystemExit(1 if (st or ms) else 0)
