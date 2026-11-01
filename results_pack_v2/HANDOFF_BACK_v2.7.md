# Results pack v2.7 — back to the BOREAS paper task

*From the model/notebook task, 2026-11-01. Attach
`BOREAS_results_pack_v27.zip` and `BOREAS_par_pack_v4.zip`.*

All six items are done. **Two of them reverse conclusions we had both
accepted**, including the sign of Limitation (viii). Read §2 first.

---

## 1. Parametric pack — you have v3, but v4 is the one to use

You marked this resolved, but the copy you received is superseded.
v3 predates the return-bend fix of v2.4. I measured rather than assumed:
at the reference it is below four decimals, but at the pumping-dominated
points it reaches 2 %, and those are the cells §4.4 quotes. **All 198
points were rerun.** Two macros moved:

| macro | v3 | **v4** |
|---|---|---|
| `\ParNftOnePumpAtShoe` | 284 | **283** |
| `\ParNftOneRTEAtShoe` | −0.047 | **−0.046** |

Everything else is identical to the printed digits. `BOREAS_par_pack_v4.zip`
is attached; please replace v3 wholesale.

## 2. The sliding-pressure ORC — your threshold is right, and the sign is the other way

You were right that 95.0 °C is the wrong condition, and your derivation
checks out: the pinch sits at 15.3 % of the duty, so at the design
pressure it holds only while the instantaneous outlet is at or above the
averaged one. My 0.01 % was answering a different question.

But the correction does not go the way either of us expected.
**The fixed-pressure assumption *understates* η_RTE; it does not
overstate it.**

| | η_ORC fixed | η_ORC sliding | η_RTE fixed | **η_RTE sliding** | energy below T̄_2d |
|---|---|---|---|---|---|
| adiabatic | 0.1462 | 0.1464 | 0.3212 | **0.3216** | 27.2 % |
| **R4** | 0.1401 | **0.1453** | 0.1447 | **0.1502** | **7.6 %** |
| R3 | 0.1451 | 0.1461 | 0.2757 | **0.2776** | 6.9 % |

At R4 a sliding-pressure plant is **+3.7 %** on η_ORC and η_RTE goes
0.1447 → 0.1502.

The reason is the one you half-anticipated, with the weighting the other
way round: **only 7.6 % of ΔE_dc is delivered while the outlet is below
the mean.** The discharge front-loads its energy into the hot part of the
window, so the fixed-pressure plant — rated at the *averaged* outlet —
under-rates the part that carries almost all the energy and over-rates a
tail that carries almost none.

So Limitation (viii) should say the opposite of what it says now: the
fixed-pressure treatment is **conservative** at R4 by about 4 %, and the
paper's η_RTE is a floor rather than a ceiling. Note the adiabatic case
behaves differently (27 % of its energy below the mean, but a nearly flat
η_ORC), so the effect is specific to a lossy store.

Macros: `\{Adiabatic,RefFour,RefThree}EtaORCSliding`, `…EtaRTESliding`,
`…EtaORCSlidingGain`, `…EnergyBelowMean`, `…TimeBelowMean`.

**Method, and one shortcut worth knowing.** η_ORC depends on T_2d alone —
everything else in the cycle is fixed — so instead of 120 pinch
bisections at 9 s each I built the curve once on 13 points spanning
90–141.6 °C and interpolated. `eta_curve.json` ships with the pack if you
want to plot it. It is smooth and monotone: 0.1199 at 90 °C to 0.1468 at
141.6 °C.

## 3. The discharge mean — you are right, and it now comes from the cumulative

| | v2.5 | **v2.7** | the budget |
|---|---|---|---|
| `\ProfRateDcMean` | −61.5 | **−69.0** | 69.0 |
| `\ProfRateChMean` | 144.4 | **145.0** | 145.0 |

Both now come from `Q_cum_J`, the exact integral V17 checks, over the
half-cycle. The record-grid values are kept as
`\Prof Rate{Ch,Dc}MeanRecordGrid` — **not for quoting**, only so the
11 % discrepancy stays visible as a property of the figure's 41-point
logarithmic grid rather than disappearing.

## 4. Annulus fill: **mud**. The paper is right

The reference annulus is mud over the whole storage length. The `Case`
field is *named* `k_cement`, but its reference value is **0.7 W/m/K**,
which is mud — the name is a leftover and the v2.4 `\SurfaceCasing`
comment inherited it. Your §2.2 and `tab:Rchain` are correct.

`toc_surface` means exactly what you guessed: **cement to surface**,
k = 1.0, the opposite bound. That is why the variant exists.

New macros: `\AnnulusFill` = mud, `\AnnulusK` = 0.7,
`\AnnulusKCementVariant` = 1.0.

## 5. `\WaterCapacityRatio` was the wrong ratio — both are now shipped

Your reading is right and the near-match is a coincidence of digits: the
two differ by about a factor of 100.

- `\WaterThroughputRatio` = **0.763**, dimensionless — everything the
  water *carries through* the tube in a half-cycle over the latent
  inventory. A duty ratio. This is what v2.5 shipped, mislabelled.
- `\WaterQuasiSteadyRatio` = **0.59 %** for a 10 K swing
  (`\WaterQuasiSteadySwing`) — the water *standing in* the tube over the
  same inventory. **This is your criterion**, and it is the one to quote.

Mine is 0.59 % against your 0.77 %; the difference is the flow area and
the density I use. If you want them reconciled exactly, tell me which
A_flow you assumed.

`\WaterCapacityRatio` is **removed**, since it was ambiguous by
construction. The diff script below flags that, and it is the one
removal in this round.

## 6. Smaller points

**Roughness — and another correction to my v2.4 answer.** I reported
"Haaland" from a *comment*; the code never calls it. What runs is a
piecewise **smooth-pipe** set, so there is no roughness to give:

```
Re < 2e3    f = 64/Re
Re <= 1e5   f = 0.3164 Re^-0.25            (Blasius)
Re <= 1e6   f = 0.0032 + 0.221 Re^-0.237   <- both legs are here
Re > 1e6    f = 0.018
```

At the reference, the fin tube is at Re = 1.21e5 and the downcomer at
2.09e5, so both sit in the third branch. `\WallRoughness` = **"not a
parameter"**; ε = 0 by construction. `\FrictionForm` and
`\FrictionBranchRef` carry the rest.

**`\BendAllowance`** now reads "two return bends per leg, fin-tube leg
only", which matches the code: `calculate_pressure_drop` applies the
hairpin allowance K = K_in + 2 K_bend + K_exit, and `borehole_pumping`
calls the downcomer with `n_bends = 0`. So the well is charged two
bends, not four.

**`\VideoDepthBottom` = 904 m** is the deepest *segment centre*, half a
segment above the store bottom at 909 m. `\VideoDepthNote` says so.

**The macro diff is now computed, not claimed.** `macrodiff.py` compares
every macro file against any git ref. You were right that v2.4/v2.5
claimed "none moved" while six had. This round, against the previous
commit:

| file | changed | added | removed |
|---|---|---|---|
| `boreas_numbers.tex` | 0 | 15 | 0 |
| `boreas_profile_numbers.tex` | **4** | 2 | 0 |
| `boreas_orange_numbers.tex` | **2** | 5 | 0 |
| `boreas_water_numbers.tex` | 0 | 3 | **1** |
| `boreas_ts_numbers.tex`, `par_numbers.tex` | 0 | 0 | 0 |

The four profile changes are §3; the two orange changes are §6
(`\FrictionForm`, `\BendAllowance`); the removal is §5. The script warns
on removals specifically, because those break a build rather than move a
number.

**`tab_ladder`** now sets `\tabcolsep` 4 pt and narrows the configuration
column to `0.22\linewidth`. That should clear 4.4 pt, but I cannot
measure elsarticle from here — tell me if it still overflows.

---

## Verification

19 checks, all passing, recorded in `verification.json`. The suite is now
**resumable** — V12/V13 and V14/V15 each run full marches and the whole
pass exceeds what my runner allows in one go, so results are cached and
a completed check is skipped on the next call. `FORCE=1` re-runs
everything.

Nothing in the model changed this round: `boreas.py`, `boreas_ground.py`
and `boreas_plots.py` are untouched since the pass that produced those 19
results, and the two checks that *could* have been invalidated by this
round's edits — V16 (every generated file built from the current cache)
and V18 (the notebook embeds the current modules) — were re-run directly
and both pass.

---

## Still open, and not mine

- **The reference h_m**, 380 kJ/kg at 150 °C. At 260 the field needs 69
  storage wells per MWe instead of 37. The largest single uncertainty in
  the economics, and a materials question.
- **The casing annulus data.** §6.2 of v2.3 showed it is worth up to 2.7
  points of loss and 5.8 % of e′ at R4 — not the sub-2 % I first reported
  from the R3 triage.
- **The shoe depth per well**, which sets the reference length.
