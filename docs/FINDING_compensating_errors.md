# Finding: the `k_w` defect and the unbounded melt front were compensating

**Date:** 2026-08-05 · **Code:** `thums_core` v1.0.0.dev0 driving `_legacy_physics` (v0.1)
**Status:** measured, not inferred — every number below came from executing the model.

---

## What was run

`run_cycle` (the strangler boundary) was run at three matched operating points, twice
each: once reproducing the legacy charging conductivity, once with the physical value.
Nothing else differs between the two columns.

| DT_3C_2C | N_lay | | `k_wall = 0.45` (legacy) | `k_wall = 45` (physical) |
|---|---|---|---|---|
| 55 | 9 | `N_wells` | **24.12** | **1.76** |
| | | `eta_rte` | 0.4303 | **−0.162** |
| | | `eps_pcm` | 0.520 | **1.470** |
| 40 | 12 | `N_wells` | 19.98 | 1.39 |
| | | `eta_rte` | 0.4198 | **−1.005** |
| | | `eps_pcm` | 0.598 | **1.421** |
| 30 | 20 | `N_wells` | 17.97 | 1.24 |
| | | `eta_rte` | 0.3945 | **−2.636** |
| | | `eps_pcm` | 0.643 | **1.346** |

## What it means

Fixing `k_w` alone does not improve the model. It destroys it.

With the physical steel conductivity, charging heat transfer is roughly an order of
magnitude faster, so the sizing bisection drives the well count from ~20 down to ~1.4.
Two things then break simultaneously:

**`eps_pcm` exceeds 1.** The model melts more PCM than exists. With one or two wells the
melt front is far outside the borehole — `delta_max = 0.5 m` against a borehole radius of
0.089 m (audit 2.5). Nothing bounds it.

**`eta_rte` goes negative.** All the flow now passes through one or two wells, so the
pressure drop and pumping power explode until parasitic power exceeds gross output.

**The legacy `k_w` error was holding the model together.** It suppressed charging, which
kept the well count high, which kept the flow per well low, which kept pumping power small
and kept the melt front inside the borehole. Two errors in opposite directions, and the
answer looked reasonable.

## The structural omission underneath

Sizing solves `Q_ratio_ch = 1` — *can the wells absorb the energy in the charging window?*
It is purely a heat-transfer criterion. There is no constraint requiring the wells to
**contain** enough PCM to hold that energy.

The inventory number is computed. `N_wells_ideal = 12.2` is printed by the baseline cell as
"ideal number of wells" and then never used again. It should be a lower bound:

```
N_wells = max(N_heat_transfer, N_inventory)
```

With the crippled conductivity the heat-transfer requirement happened to be the larger of
the two, so the missing constraint never bound and the omission stayed invisible for the
whole conference-paper campaign.

## Consequence for the migration order

The plan said: fix `k_w`, re-run, record the delta. That is now known to be wrong. The
three must move together:

1. Bound the melt front — `delta_max ≤ D_well/2`, detect front merging between the four
   legs in the borehole, raise on `V_melt > V_well`.
2. Add the inventory constraint to the sizing.
3. Then fix `k_w`.

Doing 3 before 1 and 2 produces `eps_pcm = 1.47` and negative round-trip efficiency, which
is what the table above shows.

I would also expect the pumping model to need the buoyancy term before the corrected case
gives a trustworthy `eta_rte` at low well counts: at 1,524 m with 60 K between the legs the
static head is ±7.5 bar, comparable to or larger than friction, and it is currently absent
entirely.

## What is now frozen

`verification/reference/thums_v0_baseline.json` records the three legacy points with
defects intact — behaviour, not correctness. The cycle-level numbers reproduce the
notebook's own printed output exactly:

| | package | notebook cell 21 |
|---|---|---|
| `cop_hp` | 2.956330409531823 | 2.956330409531823 |
| `eta_orc` | 0.23685132027873143 | 0.23685132027873143 |
| `eta_rte_nopump` | 0.4348350503020928 | 0.4348350503020928 |

So the strangler boundary is faithful, and any number that moves from here is attributable
to the change that moved it.

## One incident worth recording

The first attempt at this experiment reported "did not converge, last iterate = 120.0" at
every point, in zero seconds. The real cause was a bug of mine: I passed a 13-element
geometry vector where v0.1 expects 11 (the `Feb_22` branch appended `geo_grad` and
`T_geo_baseline_K`, which it computes and never uses).

The legacy solver's `except Exception: val = 1e9` swallowed that `ValueError`, both bracket
endpoints became `1e9`, the root appeared unbracketed, and the function returned its upper
bound as an answer. Had I not probed `Q_ratio_ch` directly, the honest-looking conclusion
would have been "the model cannot size this case."

A hard type error became a plausible engineering result. That is the exact failure mode
`docs/__Colab_GitHub_Token_Workflow.docx` warns about, encountered here in the space of one
afternoon.
