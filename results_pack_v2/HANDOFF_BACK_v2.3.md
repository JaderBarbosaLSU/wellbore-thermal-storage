# Results pack v2.3 — back to the BOREAS paper task

*From the model/notebook task, 2026-10-30. Attach `BOREAS_results_pack_v23.zip`
and `BOREAS_par_pack_v2.zip`.*

All of §6 is done. **No existing macro changed value: 161 → 181, twenty
added, none changed, none removed.** v2.3 fixes a stale file, adds what
was missing and changes presentation; it moves no number the paper is
already printing.

---

## 6.1 The stale table — my error, and the guard

You are right on both counts: `tab_ladder.tex` was the v2.1 file, and the
return note's claim that every derived table had been recomputed was
false. Regenerated; it now reads **97.6** MWh, matching the CSV and
`\RefFieldMWh`.

The guard is **V16**. Every generator (`mktab_ladder`, `mktables`,
`mknumbers`, `mkprofiles`) stamps the SHA-256 of the `pack_runs.json` it
read into the file it writes, and `verify.py` compares every stamp with
the current cache. Stale LaTeX is still valid LaTeX, which is exactly why
nothing caught this; a hash does catch it. All **17 checks pass**, with
V16 covering 8 generated files.

## 6.2 The three missing items

**`tab_cycle_states.tex`**, label `tab:cycle-states`, at R4: COP 2.473,
y = 0.8200, x_8h = 0.4365, both IHX effectivenesses, and all six
pressures. Macros `\RefCOP`, `\RefBleedY`, `\RefFlashX`, `\RefEpsRone`,
`\RefEpsRtwo`. They are not per MWe and never were — a COP and a pressure
level describe the machine, not its size.

**The deviation 𝒟** is a row of `tab_ref_results` and four macros:
adiabatic −12.07 %, R4 **−53.89 %**, R3 −19.89 %, R0v −63.11 %.

**The casing sensitivities — and the placeholder is wrong.** The paper's
"less than 0.2 points or 2.2 %" does not hold at the reference case:

| unknown | largest \|Δloss\| | largest \|Δe′\| |
|---|---|---|
| top of cement (surface ↔ shoe) | 1.49 pt | 3.20 % |
| k_mud, 0.6–1.5 W/(m K) | **2.71 pt** | **5.81 %** |
| casing weight, 32 vs 36 lb/ft | 1.06 pt | 4.34 % |
| **all three** | **2.71 pt** | **5.81 %** |

Macros as requested, plus `\SensCasingMaxLoss` and `\SensCasingMaxEp`.

**This is a correction to my own v2.0 advice**, where I told you the
casing unknowns were "all sub-2 %" and not worth waiting for. That triage
was done at **R3**, where the field is well shielded and the near-well
annulus is a small share of the total resistance. At **R4** the shielding
is three and a half times weaker — R_g,eff 0.915 against 3.230 m K/W — so
the same annulus change matters proportionally more. At R3 the numbers
are still small (0.09–0.27 pt); at the reference case they are not. The
WellSTAR annulus data is worth having after all.

## 6.3 Widths

`tab_T3e` now uses the short labels (adiabatic / R4 / R3 / R0v) with the
definitions in the caption and `\footnotesize`; `tab_UA` carries siunitx
units in the headers; `tab_ladder`'s configuration column is shortened
("isolated, vertical", "isolated, S-curve TVD"). The adjustbox wrappers
should not be needed.

## 6.4 `off_station_weight` — withdrawn as evidence, defined as a quantity

You are right that it is not monotone and not a property of the field. The
README now defines it exactly — *of the total kernel weight assembled at a
given spacing, the fraction whose target and source sit at different
station indices* — prints your non-monotonicity table, and states that the
v2.2 reading of the 5.2 % is **withdrawn**. Δs = 50 m is justified by
convergence alone.

The quantity is kept for the one case where it is exact and structural:
for parallel vertical wells it is **0.00 %** at every spacing, which is a
proof, not a measurement.

## 6.5 The cutoff is robust

Reran the R4 sweep at Δs = 50, 25 and 12.5 m. `cutoff_robustness.csv`:

| Δs | first 100 m | threshold | 300–400 | 400–500 | cutoff |
|---|---|---|---|---|---|
| 50 m | 0.2617 | 0.1309 | 0.155 | 0.124 | **400–500 m** |
| 25 m | 0.2619 | 0.1310 | 0.155 | 0.122 | **400–500 m** |
| 12.5 m | 0.2617 | 0.1309 | 0.154 | 0.122 | **400–500 m** |

`\CutoffRobust` = **yes**. The marginal at 400–500 m moves by 0.002 across
a fourfold refinement, so the 5 % margin is real and not a discretization
artifact.

`cutoffs.py` now reports the **first crossing only** and flags any return
above the threshold in a separate `recrosses_at` column, with a
`beyond_shoe` boolean. R4 re-crosses at 1250–1600 m (`\CutoffRecrossLo`,
`\CutoffRecrossHi`), which is beyond the shoe. R4_th45 re-crosses at
1000–1250 and R4_top1 at 800–1000; those are inside the shoe and the paper
may want to say so.

## 6.6 Section 4.3 at R4

The four profile figures are regenerated at R4, inside `field_ground`, one
PDF per panel, with the integrals axis **per fin tube** (q_prime is per
unit length of the whole well, so the per-tube rate is that over n_ft).

`boreas_profile_numbers.tex`, 19 macros, all prefixed `\Prof`:

| | |
|---|---|
| T_2c over the charge | 94.1 → 111.4 °C |
| T_2d over the discharge | 141.6 → 92.0 °C |
| subcooling at the end of discharge | 59.2 K max, 32.2 K mean |
| ε̄_ch / ε̄_dc | 0.5885 / **0.0072** |
| U_i range | 64.3 → 798.8 W m⁻² K⁻¹ |
| mean / surface-end U_i | 138.4 / 77.4, ratio 1.79 |
| charge rate per tube, start / end / mean | 60.3 / 42.9 / 57.3 kW |
| discharge rate per tube, start / end / mean | −46.2 / 0.0 / −41.8 kW |

**ε̄_dc = 0.0072 is worth a look before the paper leans on it.** The
discharge cycles essentially none of the PCM, and the discharge rate
reaches exactly zero by the end of the window. That is consistent with the
52 % loss at R4 — there is little left to give back — but it is an extreme
number and I would rather flag it than let it pass as routine.

## 6.7 One file per panel

Every multi-panel figure is now written one PDF per panel, no panel letter
or suptitle in the image, same size for panels in a row, legend in the
first only:

`comp_{a,b,c}{1,2}`, `BOREAS_grissom_sweep_{a,b,c}`,
`BOREAS_grissom_Rg_{a,b}`, `BOREAS_grissom_geometry_{a,b,c,d}`,
`BOREAS_grissom_ladder_{a,b}`, `BOREAS_prof_{Tfluid,melt,Ui}_{a,b}`,
`BOREAS_prof_integrals_{a,b}`, and in the parametric pack
`par_{cascade,hm,geom,grad,nlay}_{a,b,c}`.

A note on how: the two axes are created on **different figures** and both
handed to the plotting helper, so each panel is drawn natively. My first
attempt copied artists out of a two-panel figure and silently lost the
pinch annotation and the guide lines — most of what panel 2 is for.

`BOREAS_grissom_geometry` gained a fourth panel, the radial distribution
of the slots, which is the evidence for "ring, not disc" and was only ever
in prose.

## 6.8 The parametric pack

**Macro names fixed.** `\ParHm200Eprime` and its five siblings are now
`\ParHmTwoHundredEprime` … `\ParHmFiveHundredEprime`, and
`\ParParNPoints` is `\ParNPoints`. Both `par_numbers.py` and
`mknumbers.py` now **raise** on any name that is not `^[A-Za-z]+$` — the
failure belongs here, not in your preamble.

**Your mechanism correction is right and I have removed mine.** My own
CSV confirms it: Q_loss per well *rises* 2.04 → 2.34 MWh, +14 %, while
ΔE_ch roughly doubles. New macros `\ParHmQlossLow`, `\ParHmQlossHigh`,
`\ParHmQlossRatio` (1.14), `\ParHmTtwodLow`/`High` (105.8 → 118.2 °C) so
the text can state it from generated numbers.

**Also done:** `feasible` column in `par_results.csv`; the reference cell
outlined in red in every map; one PDF per panel.

**New macros** as requested: the Q_loss pair and ratio; the gradient-map
corners (`\ParGrad{Low,High}{Short,Shoe}{Eprime,Loss}`); e′ at both ends
of the L axis at the reference PCM (`\ParRefEprimeAtShort`,
`\ParRefEprimeAtShoe`); the T_m that maximizes e′ at each end of the h_m
axis (`\ParHmBestTmTwoHundred` = 130 °C, `\ParHmBestTmFiveHundred` =
160 °C — your qualitative observation, now a number).

And the one you asked for last: **`\ParNlaySpreadPinnedSeven` = 0.96 %**.
With the span pinned and N ≥ 7 the layer count does essentially nothing.
The 16.4 % pinned spread is almost entirely the N = 3 point.

---

## Every macro whose value changed

**None.** 161 → 181 in `boreas_numbers.tex`: twenty added, zero changed,
zero removed. The twenty are the four `*Deviation`, the five cycle
indicators, the eight casing sensitivities, `\CutoffRobust` and the two
`\CutoffRecross*`.

Separately, `boreas_profile_numbers.tex` (19, new file) and
`par_numbers.tex` (63, six renamed and sixteen added).

---

## Still open

- **SM Table S1, the water-balance numbers and the glide-shortfall
  example.** Still not moved to R4. They are the last §5.6 item.
- **ε̄_dc = 0.0072**, flagged above.
- **The casing unknowns are worth more than I said**, §6.2. If the
  WellSTAR annulus data is reachable, it is now worth the wait.
- **The reference h_m.** Your \chk stands: 380 kJ/kg at 150 °C is on the
  high side, and at 260 the field needs 69 wells per MWe instead of 37.
  That is the largest single uncertainty in the economics and it is a
  materials question, not a modelling one.
