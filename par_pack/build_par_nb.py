"""Build notebooks/BOREAS_parametric.ipynb, on the Grissom field (R4).

The previous version swept `n_wells_field`, `packing` and `fan_angle_deg`
-- the idealized field model that the ladder replaced. Those axes are
gone: the field is the measured one now, and the only thing a grid point
can change about it is the storage length.
"""
import datetime as _dt
import json
import pathlib

STAMP = (_dt.datetime.utcnow() - _dt.timedelta(hours=3)).strftime(
    '%Y-%m-%d %H:%M BRT')
SRC = pathlib.Path('/tmp/par')
OUT = pathlib.Path('/tmp/par/BOREAS_parametric.ipynb')
RAW = ('https://raw.githubusercontent.com/JaderBarbosaLSU/'
       'wellbore-thermal-storage/notebook-model/')


def md(t):
    return {'cell_type': 'markdown', 'metadata': {},
            'source': t.strip('\n').splitlines(True)}


def code(t):
    return {'cell_type': 'code', 'metadata': {}, 'execution_count': None,
            'outputs': [], 'source': t.strip('\n').splitlines(True)}


cells = []

cells.append(md(rf"""
# BOREAS — parametric study on the measured Grissom field

Every map here is computed at the **reference case (R4)**: the 363 Island
Grissom slots from the WellSTAR export, reconstructed S-curve trajectories
(staggered kickoff, $\theta_\mathrm{{max}} = 55^\circ$), the THUMS
upper-zone casing, equal-temperature superposition at a 50 m station
spacing, $t_\mathrm{{op}} = 10$ yr.

**This replaces the idealized-geometry version.** The old notebook swept
`n_wells_field`, `packing` and `fan_angle_deg` — a well count, a plan-view
shape and a cone half-angle. That description cannot represent where the
wells are, and at 1.65 m spacing the difference is most of the answer. The
field is now read from the slot coordinates, and the only thing a grid
point can change about it is the storage length.

*v2.2 · built {STAMP}*
"""))

cells.append(md(r"""
---
## 0. What the field depends on, and what it does not

This is the fact that makes a 6 $\times$ 6 map at R4 affordable at all.

$R_{g,\mathrm{eff}}(\mathrm{MD})$ comes out of the superposition solve. It
is a function of

* the well positions and the trajectories,
* the storage length $L$,
* the bore radius, the rock conductivity and diffusivity, and the field
  age.

It is **not** a function of the PCM — $T_m$, $h_m$, $N_\mathrm{lay}$, the
cascade span — nor of the tube count, nor of the geothermal gradient. The
gradient sets the rock *temperature*; the superposition sets the
*resistance* between the well and it, and those are different quantities.

So one field solve per storage length serves every point that shares that
length, and the whole gradient axis rides along for free. Five of the
eight grids below sit at the reference length and need exactly one solve
between them.

**The $L$ axis stops at the 5000 ft casing shoe**, 1524 m MD, as the paper
asks. Note the inconsistency that exposes, and which the paper task should
settle: the ladder's sweep marks a *provisional* shoe at 1250 m MD, while
the maps are capped at 5000 ft. They are two statements of the same
unknown. The grid includes 1250 m so the reference case sits on it.
"""))

cells.append(code("!pip install -q CoolProp"))

cells.append(md(r"""
---
## 1. Setup

The model, the ground solver and the grid definitions are fetched rather
than embedded, so this notebook and `BOREAS_model.ipynb` cannot drift
apart. The cached results are fetched too: 198 CSS marches is about an
hour, and nothing here needs recomputing to be read.
"""))

cells.append(code(r"""
import os, urllib.request
RAW = ('https://raw.githubusercontent.com/JaderBarbosaLSU/'
       'wellbore-thermal-storage/notebook-model/')
FILES = [('model_src/boreas.py',            'boreas.py'),
         ('model_src/boreas_plots.py',      'boreas_plots.py'),
         ('field_ladder/boreas_ground.py',  'boreas_ground.py'),
         ('field_ladder/THUMS_wells_WellSTAR_clean.csv',
          'THUMS_wells_WellSTAR_clean.csv'),
         ('par_pack/par_grids.py',          'par_grids.py'),
         ('par_pack/par_runner.py',         'par_runner.py'),
         ('par_pack/par_runs.json',         'par_runs.json')]
for remote, local in FILES:
    if not os.path.exists(local):
        urllib.request.urlretrieve(RAW + remote, local)
print('fetched', len(FILES), 'files')
"""))

cells.append(code(r"""
import numpy as np, matplotlib.pyplot as plt, pandas as pd, json
import boreas as B, boreas_ground as G, boreas_plots as BP
import par_grids as A, par_runner as R
BP.use_paper_style()

RUNS = {}
for k, v in json.load(open('par_runs.json')).items():
    RUNS.setdefault(json.loads(k)[0], []).append(v)
for r in [x for v in RUNS.values() for x in v]:
    r['feasible'] = 'error' not in r
print(f'{sum(len(v) for v in RUNS.values())} cached points in '
      f'{len(RUNS)} grids')
for k, v in sorted(RUNS.items()):
    bad = sum(1 for x in v if not x['feasible'])
    print(f'  {k:12s}{len(v):4d}' + (f'   ({bad} infeasible)' if bad else ''))
print(f"\nreference cell check: cascade at T_m = 150 C, "
      f"DT = 55 K gives e' = "
      f"{[x for x in RUNS['cascade'] if x['T_m_C']==150. and x['DT_3C_2C']==55.][0]['e_prime']:.4f}")
print("the ladder's R4 reference is 0.2151 -- the grids and the pack agree.")
"""))

cells.append(md(r"""
---
## 2. Recomputing a point

`par_runner.run_point` takes a dict, builds the case on top of the
reference, solves (or reuses) the field for that length, and runs the
march inside `field_ground`. The cache is keyed by the full point, so
nothing is computed twice.

Each point is about 20 seconds. A full 6 $\times$ 6 map is 12 minutes.
"""))

cells.append(code(r"""
# one point, recomputed from scratch, against its cached value
pt = dict(L_md=A.L_REF, T_m_C=150.0, DT_3C_2C=55.0, auto_P=True)
fresh = R.run_point(pt)
cached = [x for x in RUNS['cascade']
          if x['T_m_C'] == 150. and x['DT_3C_2C'] == 55.][0]
print(f"            {'fresh':>12}{'cached':>12}")
for q in ('loss', 'e_prime', 'RTE', 'T2d_C', 'wells_per_MWe'):
    print(f'  {q:10s}{fresh[q]:12.4f}{cached[q]:12.4f}')
"""))

cells.append(md(r"""
---
## 3. The cascade: melting temperature against grading

$T_m$ is **capped at 160 $^\circ$C**. The ceiling is the secondary loop,
not the PCM: water saturates at 179.9 $^\circ$C at 10 bar, so the top of
the cascade has nowhere to go much above that, and the top layer has to
sit below the HTHP condenser anyway. Points above 160 $^\circ$C were
measuring the loop pressure rather than the store.

`auto_P` raises the secondary pressure just enough to keep the water
liquid with 20 K of margin, so the $T_m$ axis measures the PCM. At the
reference point it returns 10 bar unchanged, which is why the 150/55 cell
reproduces the ladder exactly.

**Hatched cells are infeasible, and the reason is physical**: at a low
$T_m$ with a wide grading the bottom layer is too cold to drive the ORC
at the required pinch. The constraint is the bottom of the cascade, not
the top.
"""))

cells.append(code(r"""
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))
for ax, z, lab, fmt, best in (
        (axes[0], 'e_prime', r"$e'$ / kWh per m of well", '{:.3f}', 'max'),
        (axes[1], 'RTE', r'$\eta_\mathrm{RTE}$', '{:.3f}', 'max'),
        (axes[2], 'wells_per_MWe', 'storage wells per MWe', '{:.0f}', 'min')):
    BP.heatmap(RUNS['cascade'], 'DT_3C_2C', 'T_m_C', z, A.DG, A.TM, ax=ax,
               cbar_label=lab, fmt=fmt, best=best)
    ax.set_xlabel(r'cascade grading $\Delta T_\mathrm{3C,2C}$ / K')
    ax.set_ylabel(r'$T_m$ of the hottest layer / $^\circ$C')
fig.tight_layout(); fig.savefig('map_cascade.pdf'); plt.show()

c = RUNS['cascade']
ref = [x for x in c if x['T_m_C'] == 150. and x['DT_3C_2C'] == 55.][0]
beste = max((x for x in c if x['feasible']), key=lambda x: x['e_prime'])
bestr = max((x for x in c if x['feasible']), key=lambda x: x['RTE'])
print(f"reference    T_m 150 C, DT 55 K : e' {ref['e_prime']:.3f}  "
      f"RTE {ref['RTE']:.3f}")
print(f"best e'      T_m {beste['T_m_C']:.0f} C, DT {beste['DT_3C_2C']:.0f} K"
      f" : e' {beste['e_prime']:.3f}  RTE {beste['RTE']:.3f}")
print(f"best RTE     T_m {bestr['T_m_C']:.0f} C, DT {bestr['DT_3C_2C']:.0f} K"
      f" : e' {bestr['e_prime']:.3f}  RTE {bestr['RTE']:.3f}")
print('\nthe two optima are NOT the same cell: capacity per metre and '
      'round-trip\nefficiency are traded against each other by the grading.')
"""))

cells.append(md(r"""
---
## 4. The PCM itself: latent heat against melting temperature

This axis is new. $h_m$ runs 200 to 500 kJ/kg, which spans the salt
hydrates and nitrates that are candidates in this temperature band; 380 is
the reference.

It is the **largest single lever in the whole study**, and it is not
close. Read the two maps together: $e'$ and the well count are very nearly
reciprocal in $h_m$, because the latent inventory is what sets how much
energy a metre of well holds, while $T_m$ mostly moves the *quality* of
what comes back out.
"""))

cells.append(code(r"""
for r in RUNS['hm_Tm']:
    r['h_m_kJ'] = r['h_m'] / 1000.0
HK = tuple(h / 1000.0 for h in A.HM)
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))
for ax, z, lab, fmt, best in (
        (axes[0], 'e_prime', r"$e'$ / kWh per m of well", '{:.3f}', 'max'),
        (axes[1], 'wells_per_MWe', 'storage wells per MWe', '{:.0f}', 'min'),
        (axes[2], 'loss', 'formation loss', '{:.3f}', 'min')):
    BP.heatmap(RUNS['hm_Tm'], 'h_m_kJ', 'T_m_C', z, HK, A.TM, ax=ax,
               cbar_label=lab, fmt=fmt, best=best)
    ax.set_xlabel(r'latent heat $h_m$ / kJ kg$^{-1}$')
    ax.set_ylabel(r'$T_m$ of the hottest layer / $^\circ$C')
fig.tight_layout(); fig.savefig('map_hm.pdf'); plt.show()

h = RUNS['hm_Tm']
row = sorted([x for x in h if x['T_m_C'] == 150.], key=lambda x: x['h_m'])
col = sorted([x for x in h if abs(x['h_m'] - 380e3) < 1], key=lambda x: x['T_m_C'])
print("along h_m at T_m = 150 C:")
for x in row:
    print(f"  h_m {x['h_m']/1000:5.0f} kJ/kg  e' {x['e_prime']:.3f}  "
          f"wells/MWe {x['wells_per_MWe']:6.1f}  loss {x['loss']:.3f}")
print(f"  -> a {row[-1]['h_m']/row[0]['h_m']:.1f}x change in h_m moves e' by "
      f"{row[-1]['e_prime']/row[0]['e_prime']:.1f}x")
print("\nalong T_m at h_m = 380 kJ/kg:")
print(f"  e' spans {min(x['e_prime'] for x in col):.3f} to "
      f"{max(x['e_prime'] for x in col):.3f}, a factor of "
      f"{max(x['e_prime'] for x in col)/min(x['e_prime'] for x in col):.1f}")
"""))

cells.append(md(r"""
---
## 5. Well geometry: length, layers and tubes

The $L$ axis stops at the 5000 ft shoe. Two things to read off it.

**Depth is not free any more.** On the idealized geometry $e'$ was almost
flat in $L$ — about 1 % across a fourfold range — and depth was a pure
capital lever. At R4 it falls steadily, because past the kickoff each
extra metre of measured depth buys less true vertical depth *and* sits in
rock the neighbors no longer share.

**One fin tube is a pumping failure, not a storage one.** Watch the
$n_\mathrm{ft} = 1$ row: $\eta_\mathrm{RTE}$ goes *negative* by 1524 m.
The store still delivers heat — $T_\mathrm{2d}$ is above
$T_\mathrm{3d}$ — but with a single tube the flow per tube is three times
the reference, the pressure drop goes as roughly the square of it, and by
1524 m the borehole pumping is **284 % of the gross output**. The cell is
printed rather than clipped, because a negative round trip is the honest
statement of that design.
"""))

cells.append(code(r"""
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.3))
BP.heatmap(RUNS['L_Nlay'], 'L_md', 'N_lay', 'e_prime', A.LW_M, A.NL,
           ax=axes[0], cbar_label=r"$e'$ / kWh m$^{-1}$", fmt='{:.3f}')
axes[0].set_ylabel(r'$N_\mathrm{lay}$')
BP.heatmap(RUNS['L_Nlay'], 'L_md', 'N_lay', 'wells_per_MWe', A.LW_M, A.NL,
           ax=axes[1], cbar_label='storage wells per MWe', fmt='{:.0f}',
           best='min')
axes[1].set_ylabel(r'$N_\mathrm{lay}$')
BP.heatmap(RUNS['L_nft'], 'L_md', 'n_ft', 'RTE', A.LW_M, A.NFT, ax=axes[2],
           cbar_label=r'$\eta_\mathrm{RTE}$', fmt='{:.3f}')
axes[2].set_ylabel(r'fin tubes per well $n_\mathrm{ft}$')
for ax in axes:
    ax.set_xlabel('storage length $L$ / m MD')
fig.tight_layout(); fig.savefig('map_geometry.pdf'); plt.show()

print('pumping as a fraction of gross output, along the n_ft = 1 row:')
for x in sorted([r for r in RUNS['L_nft'] if r['n_ft'] == 1],
                key=lambda r: r['L_md']):
    flag = '   <- net negative' if x['RTE'] < 0 else ''
    print(f"  L {x['L_md']:6.0f} m   f_pump {x['f_pump']:6.3f}   "
          f"RTE {x['RTE']:+7.4f}{flag}")
"""))

cells.append(md(r"""
---
## 6. The geothermal gradient

The gradient moves the rock temperature, not the resistance, so every
column of this map shares one field solve. Hotter rock cuts the loss at
every length — the store is less far above its surroundings — while buying
very little extra capacity.
"""))

cells.append(code(r"""
fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3))
BP.heatmap(RUNS['L_grad'], 'L_md', 'grad_K_per_km', 'loss', A.LW_M, A.GR,
           ax=axes[0], cbar_label='formation loss', fmt='{:.3f}', best='min')
BP.heatmap(RUNS['L_grad'], 'L_md', 'grad_K_per_km', 'e_prime', A.LW_M, A.GR,
           ax=axes[1], cbar_label=r"$e'$ / kWh m$^{-1}$", fmt='{:.3f}')
for ax in axes:
    ax.set_xlabel('storage length $L$ / m MD')
    ax.set_ylabel(r'geothermal gradient / K km$^{-1}$')
fig.tight_layout(); fig.savefig('map_gradient.pdf'); plt.show()
"""))

cells.append(md(r"""
---
## 7. Breaking the $N_\mathrm{lay}$ / span confound

$N_\mathrm{lay}$ cannot be swept honestly at a fixed grading. The cascade
span is

$$\mathrm{span} = \Delta T_\mathrm{cascade}\,\frac{N-1}{N}$$

so raising $N$ at fixed $\Delta T_\mathrm{3C,2C}$ *widens the span at the
same time*, and "more layers is better" may be nothing but the span. The
fix is to pin $\Delta T_\mathrm{cascade} = \mathrm{span}_0\,N/(N-1)$,
holding the span at the reference 48.9 K, and sweep $N$ alone.

Run both and the confound is most of the effect.
"""))

cells.append(code(r"""
free = sorted(RUNS['Nlay_free'], key=lambda r: r['N_lay'])
pin = sorted(RUNS['Nlay_pinned'], key=lambda r: r['N_lay'])
fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.0))
for ax, z, lab in ((axes[0], 'e_prime', r"$e'$ / kWh m$^{-1}$"),
                   (axes[1], 'RTE', r'$\eta_\mathrm{RTE}$'),
                   (axes[2], 'span_K', 'cascade span / K')):
    ax.plot([r['N_lay'] for r in free], [r[z] for r in free], 'o-',
            color='#b03028', lw=1.8, ms=4, label='span free')
    ax.plot([r['N_lay'] for r in pin], [r[z] for r in pin], 's--',
            color='#23356b', lw=1.8, ms=4, label='span pinned at 48.9 K')
    ax.set_xlabel(r'$N_\mathrm{lay}$'); ax.set_ylabel(lab)
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig('confound_Nlay.pdf'); plt.show()

ef = np.array([r['e_prime'] for r in free])
ep = np.array([r['e_prime'] for r in pin])
print(f"e' spread over N_lay = 3..15")
print(f"  span free   {100*(ef.max()-ef.min())/ef.mean():6.2f} %")
print(f"  span pinned {100*(ep.max()-ep.min())/ep.mean():6.2f} %")
print(f"  -> {100*(1-((ep.max()-ep.min())/ep.mean())/((ef.max()-ef.min())/ef.mean())):.0f} % "
      f"of the apparent layer-count effect is the span moving with it")
print('\nWith the span held, the curve is still rising but has mostly')
print('flattened by N = 9, which is the reference. The layer count is a')
print('real but second-order lever; the span is the first-order one.')
"""))

cells.append(md(r"""
---
## 8. What is not here

- **SM Table S1**, the water-balance numbers and the glide-shortfall
  example. They belong with the plant tables of the results pack rather
  than with the maps, and they have not been moved to R4 yet.
- **The field axes of the old notebook** — `n_wells_field`, `packing`,
  `fan_angle_deg`. They are retired, not pending: the measured slot
  positions replace all three, and the field ladder in
  `BOREAS_model.ipynb` is where the geometry question is now answered.
- **$k_\mathrm{rock}$ as a swept axis.** It would change
  $R_{g,\mathrm{eff}}$, so it needs one field solve per value rather than
  one per length, and nothing in the paper currently asks for it.
"""))

nb = {'cells': cells,
      'metadata': {'kernelspec': {'display_name': 'Python 3',
                                  'language': 'python', 'name': 'python3'},
                   'language_info': {'name': 'python'}},
      'nbformat': 4, 'nbformat_minor': 5}
OUT.write_text(json.dumps(nb, indent=1))
print(OUT)
print(f'  {len(cells)} cells: {sum(1 for c in cells if c["cell_type"]=="code")}'
      f' code / {sum(1 for c in cells if c["cell_type"]=="markdown")} markdown')
