"""Build notebooks/P2H2P_doe.ipynb -- the Morris screen, for Colab.

A THIRD notebook, and the reason it is separate from the other two: the run
takes about an hour. The student notebook is for changing a parameter and
seeing what moves in twenty seconds; the verification notebook is for reading
every line of the model. Neither should have an hour-long cell in it.

The model itself is FETCHED from the repository rather than embedded, because
a screen has to say which version produced it and a 3000-line `%%writefile`
cell would be copied and then drift. The fetched version is printed and
recorded in the output file.
"""
import datetime as _dt
import json
import pathlib

STAMP = (_dt.datetime.utcnow() - _dt.timedelta(hours=3)).strftime(
    '%Y-%m-%d %H:%M BRT')
RUNNER = pathlib.Path('/tmp/build/doe_morris.py')
OUT = pathlib.Path('/tmp/build/P2H2P_doe.ipynb')
RAW = ('https://raw.githubusercontent.com/JaderBarbosaLSU/'
       'wellbore-thermal-storage/notebook-model')


def md(t):
    return {'cell_type': 'markdown', 'metadata': {},
            'source': t.strip('\n').splitlines(True)}


def code(t):
    return {'cell_type': 'code', 'metadata': {}, 'execution_count': None,
            'outputs': [], 'source': t.strip('\n').splitlines(True)}


cells = []

cells.append(md(rf"""
# P2H2P — Morris elementary-effects screen

Which of thirteen design factors actually move the plant, and which of them do
so **non-linearly**. This is a *screen*: it decides what the next study should
spend its runs on. It does not find an optimum and it does not give you
interaction partners.

*Runner v1 · built {STAMP} · model version is fetched and printed below*

---

### Why Morris and not a two-level factorial

A $2^{{k-p}}$ design estimates the **signed** mean effect of each factor. We
already know that $\psi$ has an interior optimum in `DT_3A_4A` near 20 K — the
elementary effects are positive below it and negative above it, so the signed
mean nearly cancels. **A factorial would report the one factor we know has an
optimum as unimportant.**

Morris separates the two things that a signed mean confuses:

| | |
|---|---|
| $\mu$ | signed mean effect — the direction |
| $\mu^*$ | mean **absolute** effect — the importance |
| $\sigma$ | spread of the effects — non-linearity **or** interaction |

A factor with small $\mu$, large $\mu^*$ and large $\sigma$ is important and
non-monotone. A factorial cannot see the difference between that and noise.

### What it costs

13 factors, $r$ trajectories, $r(k+1)$ runs at about 26 s each.

| $r$ | runs | wall time |
|---|---|---|
| 3 | 42 | ~18 min |
| 6 | 84 | ~36 min |
| **10** | **140** | **~62 min** |

Start with $r=3$ to see the shape, then run $r=10$ for anything you will quote.
**The run cell checkpoints after every trajectory**, so a Colab disconnect
costs you one trajectory, not the whole screen — re-run the cell and it
resumes.
"""))

# ---------------------------------------------------------------- setup
cells.append(md("""
---
## 1. Setup

`CoolProp` is not in a stock Colab image. The model and the runner are pulled
from the repository, so this notebook always screens a *named* version rather
than a copy that has quietly drifted.
"""))

cells.append(code(r"""
!pip install -q CoolProp
"""))

cells.append(code(rf"""
import urllib.request, pathlib, sys, re

BRANCH = 'notebook-model'          # pin a tag here for a quotable run
RAW = ('https://raw.githubusercontent.com/JaderBarbosaLSU/'
       'wellbore-thermal-storage/' + BRANCH)

for name, url in [('thums.py',      RAW + '/model_src/thums.py'),
                  ('doe_morris.py', RAW + '/model_src/doe_morris.py')]:
    urllib.request.urlretrieve(url, name)
    print(f'{{name:>14}}  {{len(open(name).read().splitlines()):>5}} lines')

import thums as T, doe_morris as dm
MODEL_VERSION = re.search(r'live model, v([\d.]+?)\.\s',
                          open('thums.py').read()).group(1)
print(f'\nmodel v{{MODEL_VERSION}} · {{len(dm.FACTORS)}} factors · '
      f'{{len(dm.KPIS)}} KPIs')
"""))

# ---------------------------------------------------------------- factors
cells.append(md(r"""
---
## 2. The factors and their ranges

**Edit `dm.FACTORS` here, not in the file.** Every bound is argued in a comment
beside it in `doe_morris.py`; three of the choices are more than range-setting
and are worth knowing about before you change them.

**`cascade_margin` replaces `DT_cascade` as the sampled quantity.** The cascade
span must exceed `DT_3C_2C − DT_4C_M`, and that bound moves with two other
factors. A box in `DT_cascade` would put much of the design infeasible — and
those corners are not bad designs, they are runs you cannot interpret. Sampling
the fractional margin above the bound makes every draw feasible by
construction.

**`n_segments` is derived** as a multiple of `N_lay`. Otherwise varying `N_lay`
drags the layer-boundary-inside-a-segment artefact along with it, and the
screen reads a numerical effect as a physical one.

**1″ pipe is excluded, and `f_pump < 0.5` is a hard screen.** At 1″ this duty
needs 2.4 kg/s per tube and the pumping reaches 1057 kW against a 1000 kW
output — $\eta_{\rm RTE}$ goes *negative*. That is not a worse design, it is
not a plant. Left in, one broken corner would dominate $\mu^*$ for every
efficiency KPI and the screen would report nothing but "pipe size and cycle
length are everything".
"""))

cells.append(code(r"""
import pandas as pd
rows = []
for f, (lo, hi, kind) in dm.FACTORS.items():
    if kind == 'nps':
        rows.append((f, f'{dm.NPS_LIST[0]}"', f'{dm.NPS_LIST[-1]}"', 'Sch 80'))
    else:
        rows.append((f, lo, hi, 'integer' if kind == 'i' else ''))
pd.DataFrame(rows, columns=['factor', 'low', 'high', 'note'])
"""))

# ---------------------------------------------------------------- sanity
cells.append(md("""
---
## 3. Sanity check before spending an hour

Two things that would waste the run if wrong: the trajectories must move
**one factor at a time**, and a mid-cube point must build a `Case` that
`validate_case` accepts.
"""))

cells.append(code(r"""
import numpy as np
k = len(dm.FACTORS)
tr, delta = dm.trajectories(2, seed=0)
bad = sum(int(np.sum(np.abs(b - a) > 1e-12) != 1)
          for Bs in tr for a, b in zip(Bs[:-1], Bs[1:]))
moved = [len({int(np.argmax(np.abs(b - a))) for a, b in zip(Bs[:-1], Bs[1:])})
         for Bs in tr]
print(f'one factor per step : {"OK" if bad == 0 else f"FAIL ({bad} steps)"}')
print(f'factors per path    : {moved}  (want {k} each)')

c = dm.build_case(np.full(k, 0.5))
print(f'\nmid-cube case:  OD {c.r_e*2000:.1f} mm · {c.num_fins} fins · '
      f'T_m {c.T_m_C:.0f} C · N_lay {c.N_lay} · n_seg {c.n_segments} '
      f'(multiple: {c.n_segments % c.N_lay == 0})')
print(f'  span {c.cascade_span:.1f} K must exceed '
      f'glide - DT_4C_M = {c.DT_3C_2C - c.DT_4C_M:.1f} K')
errs = [m for lvl, m in T.validate_case(c, verbose=False) if lvl == 'error']
print(f'  validate_case: {"clean" if not errs else errs[0][:80]}')
"""))

# ---------------------------------------------------------------- run
cells.append(md("""
---
## 4. Run the screen

**Checkpointed.** Results are written after every trajectory. If Colab drops
the session, re-run this cell — it reads what is already on disk and continues
from there. To start over, delete `morris_raw.json`.

Set `R = 3` for a first look (~18 min); `R = 10` for anything you will quote.
"""))

cells.append(code(r"""
import json, pathlib, time

R, SEED, CKPT = 10, 0, 'morris_raw.json'

trajs, delta = dm.trajectories(R, seed=SEED)
k = len(dm.FACTORS)
done = {}
if pathlib.Path(CKPT).exists():
    old = json.loads(pathlib.Path(CKPT).read_text())
    if old.get('r') == R and old.get('seed') == SEED:
        done = {(x['traj'], x['step']): x for x in old['rows']}
        print(f'resuming: {len(done)} of {R*(k+1)} points already done')

rows, t0 = [], time.time()
for ti, Bs in enumerate(trajs):
    for si in range(k + 1):
        if (ti, si) in done:
            rows.append(done[(ti, si)]); continue
        rows.append({'traj': ti, 'step': si, 'x': Bs[si].tolist(),
                     'kpi': dm.evaluate(Bs[si])})
    nb = sum(1 for x in rows if x['kpi'] is None)
    pathlib.Path(CKPT).write_text(json.dumps(
        {'factors': list(dm.FACTORS), 'kpis': dm.KPIS, 'delta': delta,
         'r': R, 'p': 6, 'seed': SEED, 'model': MODEL_VERSION, 'rows': rows}))
    el = time.time() - t0
    print(f'  trajectory {ti+1}/{R}  ·  {len(rows)} points  ·  '
          f'{nb} infeasible  ·  {el/60:.1f} min elapsed', flush=True)

nb = sum(1 for x in rows if x['kpi'] is None)
print(f'\ndone: {len(rows)} points, {nb} infeasible ({100*nb/len(rows):.1f} %)')
"""))

cells.append(md("""
**Read the infeasible fraction before anything else.** A step whose either end
is infeasible contributes no elementary effect, and dropping those steps is not
neutral — it biases $\\mu^*$ toward the interior of the cube. Above roughly
15 %, narrow the ranges of whichever factors are pushing runs out of bounds
rather than interpreting the ranking.
"""))

# ---------------------------------------------------------------- analyse
cells.append(md("""
---
## 5. Results

`analyse` returns $\\mu$, $\\mu^*$, $\\sigma$, a 90 % bootstrap interval on
$\\mu^*$, and **n** — the count of usable elementary effects. Read n beside the
ranking.
"""))

cells.append(code(r"""
res, raw = dm.analyse('morris_raw.json')
dm.summary_table(res, 'psi')
"""))

cells.append(code(r"""
dm.summary_table(res, 'eta_RTE')
"""))

cells.append(md("""
### 5.1 The Morris plane, one panel per objective

Same axes throughout. **If the rankings agree there is a single objective and
this whole discussion was unnecessary; if they disagree, ranking designs on one
of them recommends the corner that is free in that column and expensive in the
others.**

Colour is by quadrant: grey negligible, blue important and near-linear, orange
non-linear or interacting, red star = the signed mean cancels.
"""))

cells.append(code(r"""
dm.plot_planes(res);
"""))

cells.append(md(r"""
### 5.2 Signed against absolute — the non-monotone detector

On the diagonals a factor pushes one way everywhere. **Far off them, the effect
changes sign across the space**, which is the fingerprint of an interior
optimum. §5.1 cannot tell that apart from noise; this can.

`DT_3A_4A` is the prediction to check: it should sit well below the diagonal
for $\psi$ and on it for $\eta_{\rm RTE}$.
"""))

cells.append(code(r"""
dm.plot_mu_mustar(res, 'psi');
"""))

cells.append(code(r"""
dm.plot_mu_mustar(res, 'eta_RTE');
"""))

cells.append(md("""
### 5.3 Ranked, with uncertainty

With ten trajectories the estimates are noisy. An unbracketed ranking invites
over-reading, so the bars carry a 90 % bootstrap interval and the label carries
n.
"""))

cells.append(code(r"""
dm.plot_ranked(res, 'psi');
"""))

cells.append(md("""
### 5.4 Does the ranking survive changing the objective?

Flat lines mean the KPIs agree. **Crossings are the argument against
collapsing the objectives into one number** without explicit cost weights.
"""))

cells.append(code(r"""
dm.plot_rankshift(res);
"""))

cells.append(md("""
### 5.5 Raw elementary effects — a diagnostic, not a standard Morris plot

Morris compresses non-linearity and *regime switching* into one $\\sigma$ and
cannot tell them apart. A factor whose effects are merely **scattered** is
non-linear; one whose effects are **bimodal** has a regime change inside its
range — the ORC flipping between pinch-bound and hot-end-bound, or the fin
count crossing from capacity- to rate-bound.

We know both switches exist in this model, so look at the distribution rather
than its summary.
"""))

cells.append(code(r"""
dm.plot_ee_strip('morris_raw.json', 'psi');
"""))

cells.append(md("""
### 5.6 How often did the ORC change regime?

`regime_pinch_bound` is 1 when the evaporator pinch set the boiling
temperature and 0 when the hot-end rule did. If this is neither 0 nor 1 across
the screen, some of the spread in §5.1 is a regime change rather than physics.
"""))

cells.append(code(r"""
import numpy as np
v = [x['kpi']['regime_pinch_bound'] for x in raw['rows']
     if x['kpi'] and 'regime_pinch_bound' in x['kpi']]
w = [x['kpi']['n_warnings'] for x in raw['rows']
     if x['kpi'] and 'n_warnings' in x['kpi']]
print(f'pinch-bound in {100*np.mean(v):.0f} % of feasible runs '
      f'({int(sum(v))} of {len(v)})')
print(f'validate_case warnings per run: mean {np.mean(w):.2f}, max {max(w):.0f}')
"""))

# ---------------------------------------------------------------- export
cells.append(md("""
---
## 6. Export

`morris_raw.json` holds every design point and every KPI, with the model
version recorded. Keep it: the screen is only reproducible if the version that
produced it is known.
"""))

cells.append(code(r"""
import pandas as pd
tab = pd.concat({k: dm.summary_table(res, k) for k in
                 ('eta_RTE', 'psi', 'wells_per_MWe', 'UA_per_MWe',
                  'rho_P_kW_m3', 'kW_per_well')}, axis=0)
tab.to_csv('morris_summary.csv')
print(f"model {raw.get('model','?')} · r={raw['r']} · "
      f"{len(raw['rows'])} points -> morris_summary.csv")
try:
    from google.colab import files
    files.download('morris_summary.csv'); files.download('morris_raw.json')
except Exception:
    pass
"""))

cells.append(md("""
---
## 7. What this does not tell you

- **No interaction partners.** $\\sigma$ says a factor interacts or curves; it
  never says with whom. That needs Sobol $S_{ij}$.
- **No variance shares.** $\\mu^*$ is a ranking, not a decomposition — you
  cannot say "this explains 40 % of the spread".
- **No optimum, no response surface, no Pareto front.**

The screen decides which factors the next 300–500 LHS points should vary. That
study gives the trade-offs; this one gives the shortlist.
"""))

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3",
                                  "language": "python", "name": "python3"},
                   "language_info": {"name": "python"},
                   "colab": {"provenance": []}},
      "nbformat": 4, "nbformat_minor": 0}
OUT.write_text(json.dumps(nb, indent=1))
n_code = sum(1 for c in cells if c['cell_type'] == 'code')
print(f'wrote {OUT}')
print(f'  {len(cells)} cells: {n_code} code / {len(cells)-n_code} markdown')
