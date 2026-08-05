# THUMS unified model — design specification

**Repository:** `JaderBarbosaLSU/wellbore-thermal-storage`
**Scope decided:** THUMS only (PCM-in-repurposed-well storage). The `cold_tes`, `hot_tes`
and `geothermal_cooling` notebooks stay as they are.
**Results policy decided:** full recompute on the unified core, with every change from the
IHTC conference paper documented.
**Companion documents:** `THUMS_code_audit_2026-08-05.md` (findings),
`THUMS_handover_brief_v2.md` (repo state and the `Jan_12b` vs `Feb_22` question).

---

## 1. The design principle

The current code has four notebooks each containing a copy of the physics, and — more
importantly — charge and discharge implemented as **two independent boundary-value
problems** joined by a bisection on cumulative energy. Discharge restarts its front from
zero in PCM assumed uniformly solid; the liquid annulus left by charging is discarded.

Everything else follows from fixing that. The design rule is:

> **One state object is marched through the cycle. Charging produces it, discharging
> consumes it. Nothing is reset between them.**

Consequences that come free, rather than by discipline:

- Round-trip efficiency is bounded by the storage inventory, because discharge can only
  draw on the state charge left behind.
- `m_m_eff ≤ 1` by construction.
- An asymmetry like the `k_w` defect becomes unrepresentable — there is one
  `Geometry.k_wall`, read by one resistance function, used in both directions.
- The energy balance is a property of the marching, not a post-hoc ratio check.

Everything in §2–§5 is machinery in service of that sentence.

---

## 2. Package layout

```
thums_core/
  __init__.py        __version__, VERSION_STAMP (version + git short hash + UTC time)
  errors.py          ThumsPropertyError, ThumsConvergenceError, ThumsGeometryError,
                     ThumsValidityError
  config.py          frozen dataclasses (§3)
  props.py           the ONLY place CoolProp is called; cached state objects; raises
  network.py         thermal resistance network: h_e, fin efficiency, U_i, T_re
  stefan.py          moving-boundary advance, delta(t) or d(delta)/dt
  well.py            segment marching along the tube -> WellState
  cycle.py           HTHP and ORC thermodynamics
  system.py          charge(state, ...) -> state; discharge(state, ...) -> state;
                     run_cycle(...) -> CycleResult
  sizing.py          well-count and flow solvers, returning ConvergenceResult
  doe.py             Resolution IV design generation, coding, effects, Lenth PSE
  post.py            figures, tables, genuine Pareto front

notebooks/           one per case study: setup + call + plot, 10-30 lines
tests/               regression and verification suite (§6)
verification/        frozen reference fixtures (extend the existing folder)
```

Notebooks remain runnable in Colab. First cell:

```python
!pip install -q git+https://github.com/JaderBarbosaLSU/wellbore-thermal-storage@main#subdirectory=.
import thums_core; print(thums_core.VERSION_STAMP)
```

Printing the stamp on every run satisfies the "printed when the code runs" row of the
version-stamping table in the workflow note, without anyone having to remember.

---

## 3. Configuration objects

Frozen dataclasses. No function in the package takes more than a config object, a state,
and a time. The sixteen-positional-argument call signature is what allowed `k_m_l` to be
passed where `k_w` belonged; removing it removes the whole class of defect.

```python
@dataclass(frozen=True)
class Geometry:
    L_well: float; L_tube: float; D_well: float
    r_i: float; r_e: float
    fin_t: float; fin_L: float; num_fins: int; num_tubes: int
    k_wall: float          # steel — used in BOTH directions, no exceptions
    Rf_i: float            # fouling resistance
    geo_grad: float; T_geo_surface: float

@dataclass(frozen=True)
class PCM:
    T_m: float; h_m: float
    rho_s: float; rho_l: float
    cp_s: float; cp_l: float
    k_s: float;  k_l: float
    def k(self, phase): ...      # explicit phase selection, never positional

@dataclass(frozen=True)
class Numerics:
    n_segments: int; delta_tol: float; delta_max: float
    max_iter: int; time_grid: TimeGrid
```

`delta_max` must be validated at construction against `D_well/2`. The current value of
0.5 m against a 0.089 m borehole radius is not a convergence parameter, it is a licence to
melt PCM that does not exist.

---

## 4. State and results

```python
@dataclass
class StorageState:
    delta: np.ndarray          # melt-front position per segment, per layer  [m]
    liquid_fraction: float
    energy_stored: float       # [J] — tracked, not inferred
    t: float

@dataclass
class ConvergenceResult:
    converged: bool
    iterations: int
    residual: float
    message: str = ""
```

**Rule:** every solver returns `(value, ConvergenceResult)`, and nothing may write a row to
a results table or draw a figure without recording `converged`. The audit found three
places where a bisection returns a bracket bound, a collapsed midpoint, or the last iterate
after 100 iterations, all indistinguishable from success. A required field makes that
impossible.

---

## 5. Fail-loud policy

Only `props.py` touches CoolProp, and it raises `ThumsPropertyError`. Nowhere else in the
package is there a `try/except` that substitutes a value.

Specifically, these current behaviours are removed:

| Current | Replace with |
|---|---|
| `h_e = 1e9` when δ ≈ 0 (comment claims "large resistance" — it is the opposite) | analytic δ→0 limit, or raise |
| `except Exception: delta_new = 0.0` | `ThumsConvergenceError` |
| `compute_delta2` returning `0.0` when the root exceeds `delta_max` | `ThumsGeometryError` — the front left the domain |
| `hp_cop = float('inf')`, `rank_eff = float('inf')` | `ThumsConvergenceError` |
| `except ValueError: return 1e-9, 0, 0, 0` (segment silently adiabatic) | `ThumsPropertyError` |
| bare `except Exception: pass` around cycle calls | let it propagate; the DoE driver catches once, records `converged=False`, and counts |

The DoE driver is the single place allowed to catch. It records the failure, increments a
counter, and **reports the converged-run count**. A run that fails must be visible in the
output, not absorbed by `pandas.mean()` skipping NaN.

Add a `Ste = cp·ΔT/h_m` check in `stefan.py` that raises `ThumsValidityError` above a
stated threshold. The quasi-steady assumption is currently never tested, and the DoE sweeps
`DT_3C_2C` across a range where it could stop holding.

---

## 6. Verification without external data

No external validation set has been identified yet. That is a real gap for a journal
submission, and the honest response is a verification ladder that a referee will accept in
its place — plus one item that would turn the gap into a contribution.

**Tier A — must exist before any result is published**

1. **Energy closure.** `∫Q dt` against `ρ·h_m·V_melt + sensible terms`, per segment and
   globally, to a stated tolerance. This single test would have caught the geometric
   inconsistency between the finned fluid-side perimeter and the bare-cylinder Stefan front.
2. **Inventory bounds.** `V_melt ≤ V_well`, `m_m_eff ≤ 1`, front merging between adjacent
   tubes detected and raised.
3. **Zero-loss round trip.** With losses and parasitics disabled, `RTE = 1` to tolerance.
   Impossible to pass while charge and discharge are disconnected — which is the point.
4. **Grid and tolerance independence.** `n_segments`, the time grid, `delta_tol`,
   `max_iter`. None has ever been studied; a referee will ask for exactly this table.

**Tier B — the one that is worth doing properly**

5. **An independent reference solver.** A small 1-D enthalpy-method finite-volume
   solidification model in `tests/`, written deliberately in a different formulation from
   the production code. Compare δ(t) and Q(t) against the quasi-steady model across the
   parameter ranges the DoE actually sweeps.

   This is not housekeeping. It answers the question the audit could not — *where does the
   quasi-steady assumption stop being valid?* — and the output is a validity map in Stefan
   number and dimensionless time. That is a publishable figure, and it converts "we did not
   validate" into "we established the validity domain of the reduced model." For a journal
   paper it is worth more than agreement with someone else's code.

**Tier C — cheap, do it anyway**

6. `Huang_2020_validation_no_axial_conduction.ipynb` already exists in your Drive (7 June
   2026). Whatever state it is in, push it to the repo. Comparison against a published
   model is the first thing a referee asks for, and having *something* beats having nothing.

---

## 7. Migration order

Every step ends with a recorded numerical comparison. Refactoring and bug-fixing in the
same step destroys the ability to attribute a change to a cause.

| # | Step | Exit condition |
|---|---|---|
| 0 | Freeze `Jan_12b` baseline outputs to `verification/reference/thums_v0_baseline.json` | fixture committed, before any edit |
| 1 | Port to `thums_core/` package, physics **unchanged**, including current defects | reproduces the frozen fixture to machine precision |
| 2 | Fix `k_w` → `Geometry.k_wall` = 45 W/m·K in both directions | Δ recorded; both variants run once and compared |
| 3 | Unify charge/discharge on `StorageState` | zero-loss `RTE = 1` passes |
| 4 | Convert fallbacks to exceptions | list of newly-failing cases produced and triaged |
| 5 | Tier A tests green | energy closure within tolerance; grid-independence table |
| 6 | Regenerate DoE at Resolution IV, port onto unified core | defining relation and alias table published |
| 7 | Tier B reference solver, validity map | validity domain figure |
| 8 | Regenerate every paper figure from one script | figures carry version stamps |

Step 0 is not optional. It is the same discipline the `verification/` folder already
applies to the DBHE work, and it is what lets step 2 answer "how much did `k_w` matter"
with a number instead of an opinion.

---

## 8. The conference-to-journal change log

Decided policy: full recompute, every change documented. Maintain
`docs/changes_from_IHTC.md` from step 1, one row per change:

| Change | Audit ref | Effect on `N_wells` | on `RTE` | on `m_m_eff` |
|---|---|---|---|---|

Two expected entries are already known, and both push the same way — they suppressed
charging heat transfer, so correcting them should *improve* the reported performance:

- `k_w` = 0.45 → 45 W/m·K on the charge side (fin efficiency ≈0.33 → ≈0.97)
- melt-layer conductance, plane-slab → cylindrical in `compute_U_i` (~34 % low at
  δ ≈ 15.6 mm; `Feb_22` regression, already correct in `Jan_12b`)

Being able to say "the conference results were conservative, here is exactly why and by how
much" is a stronger position than quietly publishing different numbers. It is also the
answer to the referee who has the IHTC paper open next to your manuscript.

---

## 9. One note on the scope decision

Scope is THUMS only, which is right for reaching a submission. One cheap hedge: keep
`network.py` and `well.py` free of any PCM-specific assumption — they should describe a
finned coaxial well exchanging heat with *some* surrounding medium. The PCM enters through
the `stefan.py` boundary condition. That costs nothing now and is what would let the
`coaxial_ctes` and DBHE work reuse the core later without a rewrite.
