# Verification fixtures

Frozen records of what each notebook version *produces*, so that any later change
can be shown to move only what it was meant to move.

## Why this exists

`dbhe_master.ipynb` stores its outputs, but those are overwritten on every run.
Once the code changes there is nothing left to compare against, and "the results
look reasonable" quietly replaces "the results are unchanged". These files are the
missing half.

## Naming rule — fixtures are never overwritten

Every fixture carries the version that produced it:

```
reference/v10.2_baseline.json
reference/v10.3_baseline.json      <- a NEW file; v10.2 is left untouched
```

**A fixture is immutable once committed.** Overwriting it in place would make the
reference a moving target, which defeats the entire purpose — the comparison would
always pass because the thing being compared against was updated at the same time.
Keeping every version also gives a history of what moved and when.

(The notebook itself *is* overwritten in place, since it is one file. Past states
are recovered from git tags: `github.com/JaderBarbosaLSU/wellbore-thermal-storage/tree/v10.2`.)

## Usage

```bash
cd verification

# record a fixture for the current notebook
python harness.py baseline reference/v10.3_baseline.json

# check the current notebook against a frozen fixture
python harness.py check reference/v10.2_baseline.json

# depth x loop-flow sweep (2 workers; ~25 s per 20-year case)
python harness.py sweep 4500,4677 3,4,5,6 out.json
```

`check` exits 0 and prints PASS only if every recorded value is reproduced to
1e-9. Otherwise it lists each value that moved, with the year and the magnitude,
and exits 1.

## How to read a failure

A failure is not automatically a bug — it depends on what the change was meant to do.

| Change | Expected result |
|---|---|
| Renaming variables, restructuring code, adding diagnostics | **PASS.** Any moved value means the refactor altered the physics. |
| Adding a component, fixing a correlation, changing a parameter | **FAIL, and that is correct.** Read the diff, confirm every moved value is explainable, then record a new fixture under the new version. |

The second row is the one that requires discipline: record the new fixture only
*after* the diff has been explained, never before.

## `harness.py`

Loads every code cell of the notebook except the demo/dashboard/chart cells, then
calls the notebook's own `run_facility`. It is not a re-implementation — if the
harness and the notebook ever disagree, the harness is wrong by definition.

The notebook path is resolved at call time rather than bound as a default
argument. That was a real bug during development: binding it at import time made
overriding the path silently ineffective, so a deliberately perturbed notebook
still reported PASS. A regression check that cannot fail is worse than no check,
so the harness is validated both ways — an untouched notebook must PASS, and a
notebook perturbed by +0.45 % in formation conductivity must FAIL (it moves 242
recorded values, the largest being T_ret at year 1 by 0.109 K).

## Files

| File | Contents |
|---|---|
| `harness.py` | Loader, runner, sweep driver, comparator |
| `reference/v10.2_baseline.json` | Default case (4500 m, 10 m³/h, insulated tubing): 16 series x 20 annual samples, plus wheel-regime counters |
| `reference/v10.2_flow_depth_sweep.json` | 21 coupled 20-year runs: depth 4000/4500/4677 m x flow 2–10 m³/h |
| `reference/v10.2_flow_depth_sweep.png` | The figure drawn from that sweep |

### Sweep fixture schema

The sweep file predates `harness.py sweep` and uses a flatter, per-case schema
(`T_ret_20`, `T_wa_20`, `T_ro_20`, `Q_well_20`, `T_ret_yr`, `T_wa_yr`,
`wheel_reversed`, …). It is kept as the record of the chiller-feasibility study;
regenerate future sweeps with `harness.py sweep`.

### What that sweep established

- The cement reliability limit (`T_cement_limit = 165 °C`) caps depth at **4677 m**;
  at 4500 m the design is already at 96 % of it, so deepening is not a usable lever.
- Reducing loop flow is. At 4500 m, 10 -> 4 m³/h lifts year-20 return from
  **72.8 to 86.4 °C** and the accumulator from **70.4 to 78.9 °C**, clearing the
  78 °C LiBr generator minimum for all 20 years. Cost: 41 % of the well duty
  (336 -> 198 kW).
- The response is **non-monotonic**, peaking near 4–5 m³/h: below that, delivered
  power collapses and the accumulator can no longer hold temperature.
- Margin at 4500 m / 4 m³/h is only **0.9 K** at year 20; 4677 m / 4 m³/h gives 3.2 K.
- **The knob conflicts.** Lower flow reduces coil duty, dropping regeneration air
  from 50.8 to 44.3 °C — below the ~48 °C point where the desiccant wheel reverses
  and begins humidifying the process air. The v10.2 wheel-regime counter recorded
  reversal in 0.1 % of evaluations at 4 m³/h. Resolving this needs independent
  allocation of hot water between generator and regeneration coil, i.e. the
  expanded valve network of the upgraded plant.
