"""Build notebooks/BOREAS_parametric.ipynb.

Separate from BOREAS_model.ipynb on purpose. That notebook is about ONE
case and runs in a couple of minutes. This one is about 291 cyclic
steady-state evaluations, which is roughly 2.4 hours of compute, so it is
built to REPRODUCE FROM THE COMMITTED CACHE and only recompute when asked.
"""
import datetime as _dt, json, pathlib

STAMP = (_dt.datetime.utcnow() - _dt.timedelta(hours=3)).strftime('%Y-%m-%d %H:%M BRT')
HERE = pathlib.Path('.')
OUT = pathlib.Path('BOREAS_parametric.ipynb')
RAW = ('https://raw.githubusercontent.com/JaderBarbosaLSU/'
       'wellbore-thermal-storage/notebook-model/')

md = lambda t: {'cell_type': 'markdown', 'metadata': {},
                'source': t.strip('\n').splitlines(True)}
code = lambda t: {'cell_type': 'code', 'metadata': {}, 'execution_count': None,
                  'outputs': [], 'source': t.strip('\n').splitlines(True)}
cells = []

cells.append(md(rf"""
# BOREAS — parametric studies
*generated {STAMP}, boreas.py v1.3*

`BOREAS_model.ipynb` evaluates **one** case. This notebook covers the
**sweeps**: the cascade design map, the geometry maps, the field maps, the
storage-temperature study and the layer-count confound — every parametric
figure in the paper.

### Read this before running anything

The maps behind the paper are **291 cyclic steady-state evaluations**. One
costs about 30 s, so recomputing them from scratch is **roughly 2.4 hours**.

This notebook therefore **loads the committed cache** and redraws every
figure in seconds. Recomputation is opt-in, resumable, and never silently
triggered. The cache is keyed by the case parameters, so adding a point to
a grid computes only that point.
"""))

cells.append(md("""
---
## 0. Setup
"""))
cells.append(code("!pip install -q CoolProp"))
for name in ('boreas.py', 'boreas_plots.py', 'boreas_sweep.py', 'arr.py'):
    src = (HERE / name).read_text()
    cells.append(code(f'%%writefile {name}\n' + src))

cells.append(code(rf"""
import urllib.request, os, json
CACHE = 'boreas_sweep_cache.json'
if not os.path.exists(CACHE):
    url = '{RAW}data/' + CACHE
    try:
        urllib.request.urlretrieve(url, CACHE)
        print('cache downloaded:', len(json.load(open(CACHE))), 'entries')
    except Exception as e:
        print('COULD NOT FETCH THE CACHE:', e)
        print('Every figure below will be empty until the grids are '
              'recomputed (section 6, ~2.4 h). The cache is data, not code: '
              'it is the only thing here that is expensive to reproduce.')
else:
    print('cache present:', len(json.load(open(CACHE))), 'entries')
"""))

cells.append(code(r"""
import matplotlib.pyplot as plt, numpy as np
import boreas as B, boreas_plots as BP, boreas_sweep as S, arr
BP.use_paper_style()
print('serif:', BP.active_serif()[0])

# Coverage: what the cache can redraw without computing anything
GRIDS = [('cascade  T_m x glide', arr.CASCADE,  arr.ARRAY),
         ('geometry L x N_lay',   arr.G_LN,     arr.ARRAY_NOSEG),
         ('geometry n_ft x L',    arr.G_FT,     arr.ARRAY),
         ('field    N x theta',   arr.G_WF,     arr.ARRAY_NOFAN),
         ('field    G x L',       arr.G_GL,     arr.ARRAY),
         ('N_lay, span pinned',   arr.G_NLPIN,  arr.ARRAY_NOSEG)]
for name, g, fx in GRIDS:
    rows = S.collect('css', g, fixed=fx)
    print(f'  {name:22s} {len(rows):3d}/{len(g):3d} cached, '
          f'{sum(1 for r in rows if r.get("feasible")):3d} feasible')
"""))

cells.append(md(r"""
---
## 1. How the sweep driver works

Three things in `boreas_sweep.py` are there because of something that went
wrong, and they are worth knowing before you extend a grid.

**The cache is keyed by the sorted parameter dict, with a version tag.**
Adding a swept variable invalidates only the entries that actually change.
Bumping `VERSION` invalidates everything — which is what you want when the
metrics or the feasibility gate change, and what stopped a set of poisoned
results from being served silently.

**Every point is saved the moment it lands, not at the end of the run.** A
long sweep can be killed by a tool or session timeout with no chance to
clean up; a driver that persists only on completion loses the whole call's
work. That happened once, for two minutes of runs.

**Every point passes `validate_case` before it is evaluated.** Section 2
shows what that gate catches.

`run_grid_parallel` adds processes — three workers on four cores roughly
doubles throughput — and a wall-clock budget, so a sweep can be advanced in
bounded chunks and resumed.
"""))

cells.append(md(r"""
---
## 2. Three traps in parametric work with this model

All three produced wrong numbers at some point. None announced itself.
""" "\n\n"
r"""### 2.1 Past 180 °C the water flashes to steam

The secondary loop is water at `case.P` = 10 bar, which saturates at
179.9 °C. The hottest water in the plant is the charging inlet
$T_\mathrm{4c} = T_m + \Delta T_\mathrm{4C,M}$. Past that, CoolProp does
not fail — it returns a **vapour** enthalpy about ten times the liquid one,
and every pinch computed from it is meaningless. A $T_m$ sweep run without
the gate produced a false $\eta_\mathrm{ORC}$ optimum at 185 °C.
"""))

cells.append(code(r"""
import CoolProp.CoolProp as CP
c = B.Case()
T_sat = CP.PropsSI('T', 'P', c.P, 'Q', 0, 'Water') - 273.15
print(f'water at {c.P/1e5:.0f} bar saturates at {T_sat:.2f} C\n')
print(f"{'T_4c [C]':>10}{'h [kJ/kg]':>12}   phase")
for Tm in (150., 165., 175., 185., 200.):
    T4c = Tm + c.DT_4C_M
    h = CP.PropsSI('H', 'T', T4c + 273.15, 'P', c.P, 'Water') / 1e3
    print(f'{T4c:>10.1f}{h:>12.1f}   '
          + ('liquid' if T4c < T_sat else '*** SUPERHEATED STEAM ***'))
print()
for Tm in (150., 185.):
    ok, why = S.feasibility(B, B.Case(T_m_C=Tm))
    print(f'T_m = {Tm:5.1f} C -> feasible={ok}  {why[:78]}')
"""))

cells.append(md(r"""
### 2.2 The model defaults describe an *isolated* well

`Case()` ships with `n_wells_field=None`, `packing='disc'` and
`fan_angle_deg=45`. That is the conservative bound, not the reference case.
A sweep that forgets to override all three evaluates isolated wells and
returns $\eta_\mathrm{RTE} \approx 0.06$ where the array gives 0.294 — a
plausible-looking number that is simply a different machine.
"""))

cells.append(code(r"""
d = B.Case()
print('defaults :', {k: getattr(d, k) for k in
                     ('n_wells_field', 'packing', 'n_islands', 'fan_angle_deg')})
print('reference:', {k: v for k, v in arr.ARRAY.items() if k != 'n_segments'})
print('\nalways pass the field explicitly, e.g. fixed=arr.ARRAY')
"""))

cells.append(md(r"""
### 2.3 $N_\mathrm{lay}$ is confounded with the cascade span

The span of the melting temperatures is one layer width short of the
grading parameter, $\mathrm{span} = \Delta T_\mathrm{cascade}(N-1)/N$.
Sweeping $N_\mathrm{lay}$ at fixed $\Delta T_\mathrm{3C,2C}$ therefore
widens the span as well, and the apparent "fewer layers is better" is the
span, not the count. Section 5 separates them.
"""))

cells.append(code(r"""
print(f"{'N_lay':>6}{'span, free [K]':>16}{'DT_cascade to pin it':>24}")
for n in arr.NL:
    print(f'{n:>6}{55.0*(n-1)/n:>16.2f}{arr.SPAN0*n/(n-1):>24.2f}')
"""))

cells.append(md(r"""
---
## 3. The configuration table

The same fourteen field scenarios that `BOREAS_model.ipynb` computes in
section 5, read here from cache in a second rather than recomputed in a
quarter of an hour. Use whichever suits: that notebook is self-contained
and slow, this one is fast and depends on the cache.

**array** and **fan** describe the trajectories --- vertical throughout,
against deviating at half-angle $\theta$ below the kickoff. **disc**,
**rows** and **ring** describe the plan view of the wellheads: a filled
circle, two straight lines side by side, and the same double line closed
into a loop. Every field case is one of each.
"""))

cells.append(code(r"""
import pandas as pd
cfg_grid = [{k: v for k, v in d.items() if k != 'label'} for d in arr.CONFIGS]
cfg = S.collect('config', cfg_grid, fixed=arr.CONFIG_FIXED)
tab = {}
for d, r in zip(arr.CONFIGS, cfg):
    if not r.get('feasible'):
        continue
    tab[d['label']] = {'loss %': r['loss_pct'], 'eta_storage': r['eta_storage'],
                       'eta_RTE': r['eta_RTE'], 'eps_cycled': r['eps_cycled'],
                       'T_2d [C]': r['T_2d_realised_C'],
                       'wells/MWe': r['wells_per_MWe']}
print(f'{len(tab)} of {len(arr.CONFIGS)} configurations cached')
display(pd.DataFrame(tab).T.round(4))
"""))

cells.append(md(r"""
**Closing the rows into a ring buys almost nothing.** The ring is exposed
on both faces, so its perimeter is the rim length twice over; against the
double row it saves only the two end-caps, $2w_\mathrm{row}$ per island.
That is 1.6 % of the perimeter at 1000 wells and 14 % at 100, where the
caps are a larger share of a shorter row. The arrangement that matters is
whether the wells are strung out at all, not how the line is closed.
"""))

cells.append(code(r"""
import numpy as np
m = {d['label']: r for d, r in zip(arr.CONFIGS, cfg) if r.get('feasible')}
print(f"{'case':24s}{'rows':>9}{'ring':>9}{'disc':>9}   (loss %)")
for stem in ('fan, 100 wells', 'fan, 1000 wells', 'fan, 1000, 70 deg',
             'vertical array'):
    vals = []
    for pk in ('rows', 'ring', 'disc'):
        r = m.get(f'{stem}, {pk}')
        vals.append(f"{r['loss_pct']:9.1f}" if r else f"{'--':>9}")
    print(f'{stem:24s}' + ''.join(vals))
print()
for N in (100., 1000.):
    for pk, ni in (('disc', 1), ('rows', 4), ('ring', 4)):
        cc = B.Case(n_wells_field=N, fan_angle_deg=0., packing=pk, n_islands=ni)
        P0 = float(np.atleast_1d(B.field_perimeter(cc, 0.))[0])
        print(f'  N={int(N):5d}  {pk:5s}  P0 = {P0:8.1f} m  '
              f'-> {P0/N:.3f} m per well')
"""))

cells.append(md(r"""
---
## 3. The maps, redrawn from the cache

Nothing below computes. Each call is a dictionary lookup plus a plot.

Highlighting follows the paper: the **reference** cell has its value in
bold with a red outline, the **best** cell has its value in bold at a
larger size. Values are printed on the cells because several of these maps
are deliberately flat, and flatness cannot be judged from colour.
"""))

cells.append(code(r"""
cas = S.collect('css', arr.CASCADE, fixed=arr.ARRAY)
fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.3))
for ax, (z, cl, cm, bst) in zip(axes, [
        ('eta_RTE', r'$\eta_\mathrm{RTE}$', 'viridis', 'max'),
        ('kWh_per_m_well', r"$e'$ / kWh m$^{-1}$", 'magma', 'max'),
        ('loss_frac', 'loss fraction', 'cividis', 'min')]):
    _, Z = BP.heatmap(cas, 'DT_3C_2C', 'T_m_C', z, arr.DG, arr.TM, ax=ax,
                      cmap=cm, cbar_label=cl, mark=(55.0, 150.0), best=bst)
    ax.set_xlabel(r'$\Delta T_\mathrm{3C,2C}$ / K')
    print(f'{z:16s} spread {BP.spread(Z):6.1f} %')
axes[0].set_ylabel(r'$T_\mathrm{m}$ / $^\circ$C')
fig.tight_layout(); fig.show()
"""))

cells.append(md(r"""
**$\eta_\mathrm{RTE}$ and $e'$ optimise in opposite corners.** Efficiency
peaks at a cool store with a narrow glide, capacity at a hot one. The
reference case sits at neither. A narrower glide helps *both* at every
storage temperature — there is no cell where 55 K beats 25 K — so if the
55 K is set by the PCM selection rather than by choice, that belongs in the
text.
"""))

cells.append(code(r"""
ln = S.collect('css', arr.G_LN, fixed=arr.ARRAY_NOSEG)
ft = S.collect('css', arr.G_FT, fixed=arr.ARRAY)
fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.3))
_, Z1 = BP.heatmap(ln, 'L_ft', 'N_lay', 'kWh_per_m_well', arr.LW, arr.NL,
                   ax=axes[0], cmap='magma', cbar_label=r"$e'$ / kWh m$^{-1}$",
                   mark=(5000.0, 9))
axes[0].set_ylabel(r'$N_\mathrm{lay}$')
_, Z2 = BP.heatmap(ln, 'L_ft', 'N_lay', 'eta_RTE', arr.LW, arr.NL, ax=axes[1],
                   cmap='viridis', cbar_label=r'$\eta_\mathrm{RTE}$',
                   mark=(5000.0, 9))
_, Z3 = BP.heatmap(ft, 'L_ft', 'n_ft', 'eta_RTE', arr.LW, arr.NFT, ax=axes[2],
                   cmap='viridis', cbar_label=r'$\eta_\mathrm{RTE}$',
                   mark=(5000.0, 3))
axes[2].set_ylabel(r'$n_\mathrm{ft}$')
for a in axes: a.set_xlabel(r'$L$ / ft')
fig.tight_layout(); fig.show()

# the separability result, stated as numbers rather than eyeballed
for i, n in enumerate(arr.NL):
    r_ = Z1[i][np.isfinite(Z1[i])]
    print(f"  N_lay={n:>2}: e' varies {100*(r_.max()-r_.min())/r_.mean():5.2f} % across L")
print(f"  across N_lay the same map spans {BP.spread(Z1):.1f} %")
"""))

cells.append(code(r"""
wf = S.collect('css', arr.G_WF, fixed=arr.ARRAY_NOFAN)
gl = S.collect('css', arr.G_GL, fixed=arr.ARRAY)
fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.3))
BP.heatmap(wf, 'n_wells_field', 'fan_angle_deg', 'loss_frac', arr.NW, arr.FA,
           ax=axes[0], cmap='cividis', cbar_label='loss fraction',
           mark=(1000.0, 0.0), best='min')
axes[0].set_ylabel(r'$\theta$ / deg')
BP.heatmap(wf, 'n_wells_field', 'fan_angle_deg', 'eta_RTE', arr.NW, arr.FA,
           ax=axes[1], cmap='viridis', cbar_label=r'$\eta_\mathrm{RTE}$',
           mark=(1000.0, 0.0))
for a in axes[:2]: a.set_xlabel(r'$N$ wells')
BP.heatmap(gl, 'L_ft', 'grad_K_per_km', 'loss_frac', arr.LW, arr.GR,
           ax=axes[2], cmap='cividis', cbar_label='loss fraction',
           mark=(5000.0, 56.5), best='min')
axes[2].set_ylabel(r'$G$ / K km$^{-1}$'); axes[2].set_xlabel(r'$L$ / ft')
fig.tight_layout(); fig.show()
"""))

cells.append(md(r"""
**These panels use the rows packing throughout**, so the 45° column is *not*
comparable with the disc-packing fan rows of the paper's formation table.
Mixing the two packings is how an earlier draft concluded that the array and
the fan were "the same machine".
"""))

cells.append(md(r"""
---
## 4. Storage temperature

Cheap: these are cycle-only evaluations, no borehole march. `auto_P` raises
the loop pressure to keep the water liquid; the fixed-pressure branch goes
infeasible past 165 °C and that is the point of plotting both.
"""))

cells.append(code(r"""
Ts = [float(t) for t in range(105, 216, 5)]
a_P = S.collect('cycle', [{'T_m_C': t, 'auto_P': True} for t in Ts])
fixP = S.collect('cycle', [{'T_m_C': t} for t in Ts])
fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.0))
_, info = BP.cycle_vs_Tm(a_P, fixP, B.Case(), axes=axes)
fig.tight_layout(); fig.show()
for k, v in info.items(): print(f'  {k:22s} {v}')
"""))

cells.append(md(r"""
---
## 5. Breaking the $N_\mathrm{lay}$ / span confound

Pinning $\Delta T_\mathrm{cascade} = \mathrm{span}_0 N/(N-1)$ holds the span
at the reference 48.9 K while the layer count varies. Setting
`DT_cascade = DT_3C_2C` is a **no-op** — the model falls back to
`DT_3C_2C` when `DT_cascade` is None, so passing the same value changes
nothing. That mistake is why an earlier check reported no effect.
"""))

cells.append(code(r"""
pin = sorted(S.collect('css', arr.G_NLPIN, fixed=arr.ARRAY_NOSEG),
             key=lambda r: r['N_lay'])
fm = {(r['L_ft'], r['N_lay']): r for r in ln}
free = [fm[(5000.0, n)] for n in arr.NL]
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.8))
N = np.array(arr.NL)
for a, key, lab in ((ax[0], 'eta_RTE', r'$\eta_\mathrm{RTE}$'),
                    (ax[1], 'kWh_per_m_well', r"$e'$ / kWh m$^{-1}$")):
    a.plot(N, [r[key] for r in free], 'o-', color=BP.RED, lw=1.8, ms=5,
           label=r'$\Delta T_\mathrm{3C,2C}$ fixed (span varies)')
    a.plot(N, [r[key] for r in pin], 's--', color=BP.NAVY, lw=1.8, ms=5,
           label='span pinned at 48.9 K')
    a.set_xlabel(r'$N_\mathrm{lay}$'); a.set_ylabel(lab); a.legend()
fig.tight_layout(); fig.show()

e_free = np.array([r['eta_RTE'] for r in free])
e_pin = np.array([r['eta_RTE'] for r in pin])
hi = [i for i, n in enumerate(arr.NL) if n >= 7]
print(f'  free   spread over all N_lay : {100*(e_free.max()-e_free.min())/e_free.mean():5.2f} %')
print(f'  pinned spread for N_lay >= 7 : '
      f'{100*(e_pin[hi].max()-e_pin[hi].min())/e_pin[hi].mean():5.2f} %')
"""))

cells.append(md(r"""
---
## 6. Extending or recomputing a grid

**Cost first.** A cyclic steady-state point is ~30 s. The 180 points behind
the maps are ~1.5 h serially, ~45 min on three workers. Cycle-only points
are ~6 s.

`run_grid` and `run_grid_parallel` are resumable and budgeted: they skip
anything already cached, stop when `budget_s` is spent, and return how many
remain. Call again to continue. Nothing is recomputed twice.

The cell below is a **small** extension — four new glide values at the
reference temperature, about two minutes — to show the mechanism without
committing you to hours. Raise `budget_s` and widen the grid when you mean
it.
"""))

cells.append(code(r"""
EXTRA = [{'T_m_C': 150.0, 'DT_3C_2C': float(d), 'auto_P': True}
         for d in (30, 35, 45, 50)]
cache, done, left = S.run_grid_parallel('boreas', 'css', EXTRA,
                                        fixed=arr.ARRAY, budget_s=150,
                                        workers=3, verbose=True)
rows = S.collect('css', EXTRA, fixed=arr.ARRAY)
for r in sorted(rows, key=lambda r: r['DT_3C_2C']):
    print(f"  glide {r['DT_3C_2C']:5.1f} K -> eta_RTE {r['eta_RTE']:.4f}, "
          f"e' {r['kWh_per_m_well']:.4f}")
print(f'\n{done} computed this call, {left} still to do '
      f'(re-run this cell to continue)')
"""))

cells.append(md(r"""
---
### What this notebook will not tell you

- **Whether the optima are real optima.** The grids are 6×6. A maximum on a
  coarse grid is a cell, not a stationary point.
- **Anything about the arrangement of a diverging field.** $P(z)$ scales a
  row perimeter with the fan radius, and wells spread over 6500 acres are
  not in rows at depth. The rows and disc results bound the arrangement;
  neither is the answer.
- **Multi-year behaviour.** Every point is evaluated at one field age,
  `t_operation_yr`. The loss falls as the rock warms, and that sweep has
  not been run.
"""))

nb = {'cells': cells,
      'metadata': {'colab': {'provenance': [], 'toc_visible': True},
                   'kernelspec': {'name': 'python3', 'display_name': 'Python 3'},
                   'language_info': {'name': 'python'}},
      'nbformat': 4, 'nbformat_minor': 0}
OUT.write_text(json.dumps(nb, indent=1))
n_code = sum(1 for c in cells if c['cell_type'] == 'code')
print(f'{OUT.name}: {len(cells)} cells, {n_code} code / {len(cells)-n_code} markdown')
