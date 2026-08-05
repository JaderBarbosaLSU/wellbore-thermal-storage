# THUMS — unifying the four notebooks around one KPI schema

Companion to `THUMS_unified_design_spec.md`. That document covers the solver; this one
covers everything downstream of it — how the four notebooks collapse into one code, and
how post-processing is built so that adding a KPI or a study type does not mean writing
another notebook.

---

## 1. The observation the design rests on

The four notebooks are not four models. They are **one model driven by four different
experiment types**, each of which grew its own private copy of the physics and its own
private definition of the metrics.

| Notebook | What it actually is | Post-processing it produces |
|---|---|---|
| `THUMS_..._Feb_22_IHTC` | a single design point, in full detail | state-point tables, T–s diagrams, borehole profiles |
| `THUMS_..._BB_param_Jan_12b` | a 2-factor grid sweep (`DT_3C_2C` × `N_lay`) | contour maps |
| `THUMS_..._DoE_Dec_5` | a 16-factor fractional factorial | main effects |
| `Pareto_for_THUMS_DoE` | consumer of the DoE CSV | effect bar charts, multi-KPI table |

So the unification has three parts, and only the first is about the solver:

1. **One solver** — covered in the design spec.
2. **One results schema**, produced identically by every study type.
3. **One post-processing layer**, which consumes only that schema and therefore works for
   any study, present or future.

The current arrangement fails at (2) and (3), and the consequences are not cosmetic.
`RTE` is computed in three places and means different things in each; the CSV that the
Pareto notebook reads cannot be traced to the code that produced it; and effects are
averaged over runs that silently failed.

---

## 2. The KPI registry — the centre of the design

Every metric is declared once, in `thums_core/kpi.py`, and computed in exactly one place.

```python
@dataclass(frozen=True)
class KPI:
    name:       str                      # 'eta_rte'
    label:      str                      # r'$\eta_{\mathrm{RTE}}$'  — used by every figure
    unit:       str                      # '-' , 'MWh', 'wells/MW'
    direction:  str                      # 'max' | 'min'
    group:      str                      # thermodynamic | storage | sizing | economic
    depends_on: tuple[str, ...]          # other KPIs this is algebraically built from
    fn:         Callable[[CycleResult], float]

REGISTRY: dict[str, KPI] = {...}
```

Three fields carry real weight:

**`direction`** — a Pareto front cannot be computed without knowing which way is better.
The notebook named `Pareto_for_THUMS_DoE` contains no dominance test anywhere; "Pareto"
there refers only to a Pareto *bar chart*. With `direction` declared, a genuine
non-dominated sort is four lines.

**`depends_on`** — the audit found the multi-KPI "consensus" averaging seven metrics of
which at least four are algebraic products of the others, so the same physical mechanism is
counted two or three times. With dependencies declared, the post-processor can select a
maximal independent subset automatically, or refuse and say why. A statistical error becomes
structurally impossible rather than a thing to remember.

**`label`** — one LaTeX string per KPI, used by every axis, legend and table in the paper.
Renaming a quantity updates every figure.

### Proposed canonical KPI set

Deliberately different from the current one in one respect: **no KPI is baked into
another.**

**Thermodynamic**

| name | definition | dir |
|---|---|---|
| `eta_rte` | `(W_out − W_pump,dc)·t_dc / ((W_in + W_pump,ch)·t_ch)` — the actual round trip | max |
| `cop_hp` | HTHP COP, isentropic efficiencies applied *inside* the cycle | max |
| `eta_orc` | ORC thermal efficiency, turbine and generator losses included | max |
| `eta_ex` | exergy efficiency of the whole PtHtP chain | max |
| `X_dest_hthp`, `X_dest_well`, `X_dest_orc` | exergy destruction by component | min |

**Storage**

| name | definition | dir |
|---|---|---|
| `eps_pcm` | melted PCM mass / total PCM mass — utilization, **reported on its own** | max |
| `E_well` | energy stored per well [MWh] | max |
| `rho_E` | energy density [kWh per m³ of well] | max |
| `self_discharge` | loss rate [% per day] — meaningful only once formation losses exist | min |

**Sizing and parasitics**

| name | definition | dir |
|---|---|---|
| `N_wells` | wells for the stated duty | min |
| `wells_per_MW` | specific well count | min |
| `f_pump` | pumping work / gross work output | min |

**Economic** (you already hard-code `PCM_cost_per_kWh = 20`)

| name | definition | dir |
|---|---|---|
| `c_pcm` | PCM cost per kWh stored | min |
| `LCOS` | if you choose to go there | min |

### On `weighted_RTE`

The reason that quantity was invented is clear: you wanted to trade efficiency against PCM
utilization in a single number. But `RTE × ε_PCM` is not an efficiency of any control
volume — it falls when utilization falls even though electricity in and out are unchanged,
and a reader will read 0.30 as 30 % round-trip.

**A Pareto front is exactly the tool for that trade-off.** Keep `eta_rte` and `eps_pcm`
orthogonal, plot the front, and the trade-off is shown rather than asserted with an
arbitrary weighting. You named a notebook after the right idea and then implemented a
scalarisation instead.

---

## 3. The results schema

One tidy row per run, three blocks of columns, identical for every study type.

```
inputs      T_m_C, rho_m_s, ..., N_lay, T_sink_C, T_source_C
            + coded_<factor> when the run belongs to a design
kpis        eta_rte, eps_pcm, N_wells, ... (exactly the registry names)
status      converged (bool), n_iter, residual, warnings (list),
            version, git_commit, study_id, run_id, wall_time_s
```

Two rules make this worth having:

- **`status.converged` is mandatory.** Nothing writes a row without it. The audit found
  three bisections that return a bracket bound, a collapsed midpoint, or the last iterate
  after 100 iterations, all indistinguishable from success.
- **Failed runs are written, not dropped.** A failed run is a row with `converged=False`
  and NaN KPIs. Post-processing must then either exclude them explicitly and report the
  count, or refuse to compute a balanced contrast. `pandas.mean()` silently skipping NaN
  is how an unbalanced, non-orthogonal effect estimate gets published without a warning.

Alongside every results file, a `run_manifest.json`:

```json
{"study_id": "...", "study_type": "FractionalFactorial",
 "version": "1.0.0", "git_commit": "abc1234", "timestamp_brt": "...",
 "design": {"n_factors": 16, "n_runs": 128, "resolution": "IV",
            "generators": ["H=ABC", "..."], "defining_relation": "...",
            "alias_table": {...}},
 "factor_levels": {...},
 "n_converged": 128, "n_failed": 0}
```

That manifest is the answer to "which model version produced this CSV" — currently
unanswerable, with two notebooks defining `run_parametric_simulation` with identical
signatures and different physics.

---

## 4. Studies as configuration, not code

```python
class Study(Protocol):
    def design(self) -> pd.DataFrame: ...     # rows of factor settings
    def metadata(self) -> dict: ...           # goes verbatim into the manifest

CasePoint(config)                              # 1 run   → Feb_22's role
FactorialGrid({'DT_3C_2C': [...], 'N_lay': [...]})   # grid → Jan_12b's role
FractionalFactorial(factors, levels, resolution=4)   # 128  → Dec_5's role
LatinHypercube(factors, n)                     # optional, for surrogates
```

One runner for all of them:

```python
df, manifest = run_study(study, base_config, n_jobs=-1)
```

`run_study` owns parallelism, the failure counter, progress reporting, and the cache. It is
**the only place** allowed to catch an exception from the solver — it records
`converged=False` and continues. Everywhere else, errors propagate.

**Caching.** Key on `hash(config) + version`. A 128-run design over an expensive solve
should not be recomputed because a plotting function changed. `Jan_12b` already has an
ad-hoc `_mass_eff_cache`; formalise it, make it content-addressed, and invalidate it on
version change so a physics fix can never be served from stale cache.

---

## 5. Post-processing

`thums_core/post/` consumes the schema and nothing else, so every function works for every study.

```
post/
  case.py       report_case(result)        state-point tables, T–s diagrams, profiles
  fields.py     contour(df, x, y, kpi)     the Jan_12b maps, any KPI, any factor pair
  effects.py    main_effects(df, design)   orthogonal contrasts from the design metadata,
                                           Lenth PSE, half-normal plot, curvature contrast
  pareto.py     pareto_front(df, kpis)     real non-dominated sorting via KPI.direction
                tradeoff(df, kx, ky)       2-D fronts with the dominated cloud behind
  exergy.py     waterfall(result)          destruction by component
  style.py      one rcParams, one palette, stamp_figure() on every figure
  tables.py     LaTeX export straight into the manuscript
```

Two things follow from building it this way. `effects.py` computes contrasts from the
design metadata rather than from `df[df[f]==low].mean()`, so an unbalanced subset raises
instead of quietly producing a biased number. And every figure in the paper is generated by
one script from one results file, so "which run made Figure 7" always has an answer.

---

## 6. What each notebook becomes

| Now | Becomes |
|---|---|
| `Feb_22` cells 3–19 (physics) | `thums_core/` package — deleted from the notebook |
| `Feb_22` cell 21 (baseline) | `notebooks/01_case_study.ipynb` — ~20 lines |
| `Feb_22` cell 22 (profile plots) | `post.case.report_case()` |
| `Feb_22` cell 24 `run_parametric_simulation` | `thums_core.system.run_cycle()` |
| `Jan_12b` cells 21–23 (grid + contours) | `notebooks/02_parameter_maps.ipynb` + `post.fields` |
| `Dec_5` cell 22 (design generation) | `thums_core.doe.FractionalFactorial` (Resolution IV) |
| `Dec_5` cells 24–33 (effects, RSM) | `post.effects` — the rank-deficient RSM cells dropped |
| `Pareto_for_THUMS_DoE` | `post.pareto` + `post.effects`, with an actual front |
| all four | archived under `legacy/`, tagged `v0.1`, never edited again |

Keeping the originals under `legacy/` at a tag matters: they are the record of what the
conference paper ran, and step 1 of the migration is proving the new code reproduces them.

---

## 7. Build order

Follows the design spec's migration table; this is the downstream half.

1. `kpi.py` registry with the current (defective) definitions, so the port is a pure
   refactor and the frozen fixture still matches.
2. `results.py` schema + `run_manifest`. Immediately re-emit the existing DoE CSV through
   it, so provenance exists before anything else changes.
3. `study.py` + `run_study` with cache and failure counting. Reproduce `Jan_12b`'s contour
   maps through the new path and diff the figures.
4. Fix the KPI definitions — split `weighted_RTE` into `eta_rte` and `eps_pcm`, apply
   component efficiencies consistently. Record the delta in `docs/changes_from_IHTC.md`.
5. `post.effects` with correct contrasts, Lenth PSE, half-normal plot; regenerate the DoE
   on the Resolution IV design.
6. `post.pareto` — the front you have been describing but not computing.
7. `post.exergy` once the exergy KPIs exist.

Steps 1–3 change no physics and no numbers. That is the point: they are what make step 4
onwards attributable.
