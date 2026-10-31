# Results pack v2.4 — back to the BOREAS paper task

*From the model/notebook task, 2026-10-31. Attach
`BOREAS_results_pack_v24.zip` and `BOREAS_par_pack_v3.zip`.*

All of §1–§7 and §10 is done, and §8 is partly done. **You were right on
item 1 and I had it backwards in two ways at once.** 18 checks pass.

---

## 1. The rates — confirmed, and worse than you said

You are right that `q'_j` is already per fin tube. I verified it against
the budget before changing anything: **n_ft × Q_cum_J reproduces ΔE_ch and
ΔE_dc per well exactly** (2×10⁻¹⁶). So v2.3 divided by n_ft on top of a
quantity that was already per tube, and the plotted rates were a factor of
three low. Corrected:

| | start | end | mean |
|---|---|---|---|
| charge, per tube | 180.9 | 128.8 | **144.4 kW** |
| discharge, per tube | −138.7 | −0.1 | **−61.5 kW** |
| charge, per well | | | 433 kW |

Your cross-check wanted 145 kW per tube on charge; it is 144.4.

**1b, the means, also confirmed.** They were arithmetic means over a
logarithmic grid whose steps run 0.3 s to 8491 s, so the first minute
weighed as much as the last two hours. Now step-weighted (trapezoidal).
For the discharge that is −125 → −61.5 kW.

**The V-check you asked for is V17**, and it is exact rather than
approximate, because I used the solver's own cumulative rather than the
record grid: `n_ft × Q_cum_J` against ΔE per well, 2×10⁻¹⁶ on both
half-cycles. One honest note: the *plotted* curve integrates to within
0.4 % (charge) and 11 % (discharge) of that, because the 41-point
logarithmic record grid under-resolves a rate that collapses. That is the
figure's grid, not the physics, and V17 deliberately does not depend on
it.

## 2. ε̄_dc — your reading is right, mine was wrong

Relabelled: `\ProfEpsCh` is "melted fraction at the END OF CHARGE",
`\ProfEpsDc` is "residual melted fraction at the END OF DISCHARGE", and I
added `\ProfEpsCycled` = **0.5811** so the quantity the text wants exists
rather than being a subtraction the reader must do.

**The below-pinch share, for all three cases:**

| | T_3e + ΔT_pinch | outlet at the end | **energy below** | time below |
|---|---|---|---|---|
| adiabatic | 98.9 °C | 131.6 °C | 0.00 % | 0.0 % |
| **R4** | 95.0 °C | 92.0 °C | **0.01 %** | **11.8 %** |
| R3 | 98.2 °C | 118.9 °C | 0.00 % | 0.0 % |

So the effect is real but immaterial: at R4 the outlet spends 11.8 % of
the window below what the ORC can use, and delivers **0.01 %** of ΔE_dc
there — because by the time the outlet falls, the rate has already
collapsed. The overstatement of η_RTE is bounded at 0.01 %.

**I did not do (b), and I would argue against it.** A time-resolved ORC
with a sliding evaporating pressure is a real piece of work, and (a)
already bounds what it would refine at one part in ten thousand.
Limitation (viii) can state the bound instead of the caveat.

## 3. Profile figures

Fixed: temperatures in °C with "temperature / °C"; the depth axis
labelled **"TVD / m"** on *both* panels (in per-panel mode there is no
shared left axis to carry it); a colorbar, "time from the start of the
half-cycle / h", on panel (b); the integrals in the corrected units with
the dashed curve labelled "melted fraction of the PCM".

`BOREAS_grissom_Rg_b` now reads "reference (R4)", "KOP 300 ft (all
wells)", "hold 45°", "hold 70°", with the R4 curve added.

## 4. The ring, measured

| | |
|---|---|
| `\GrissomExtentLong` | **202 m** |
| `\GrissomExtentShort` | **139 m** |
| `\GrissomRminOverRmax` | **0.51** |

Along the principal axes of the 363 slots. Your revised "about 200 by
140" and "half" are both right, and your diagnosis of the 290 m is
consistent with what I see — the 363-slot set has no such outlier.

## 5. Widths

`tab_T3e` now sets `\footnotesize` **and** `\tabcolsep` 4 pt;
`tab_ladder`'s configuration column is `p{0.24\linewidth}`. Please check
against your text width — I cannot measure elsarticle from here.

## 6. The T–s figure at R4

Delivered. You were right that it was simply absent — the notebook drew
it, the pack never carried it. `BOREAS_Ts_combined_water.pdf`, at the
**realised** outlet (116.54 °C) rather than the provisional 146.11 °C,
plus `BOREAS_Ts_cycles.pdf` and `boreas_ts_numbers.tex`. README §7 is
rewritten and no longer claims the parametric notebook is idealized.

## 7. R3 casing maxima

`\SensCasingMaxLossRthree` = **0.27**, `\SensCasingMaxEpRthree` =
**2.20**. Against 2.71 and 5.81 at R4, so "an order of magnitude less" is
right for the loss and a factor of 2.6 for e′.

## 10. The section 10 checks — `boreas_orange_numbers.tex`, 27 macros

| | |
|---|---|
| `\RefTfTop` / `\RefTfBottom` | **18.3** / **58.7 °C** (TVD 891 m) — your text is right |
| `\RefFourCorrPassDT` | **0.0** — exactly zero in double precision, at all three cases, better than the claimed 1e-4 |
| `\ConvNsBase/Fine`, `\ConvNtFine` | loss 52.415 → 52.583 → 52.596 %; e′ 0.2151 → 0.2137 → 0.2130 |
| `\IntegratorCheck` | **5×10⁻¹¹** |
| `\PumpEta` | **none** — hydraulic power only, no pump or motor efficiency |
| `\FrictionForm` | **Haaland** turbulent, laminar below transition |
| `\SurfaceCasing` | **not represented** — one string over the whole storage length with annulus fill `k_cement`; no second shallower string |

**Three things in §10 need more than a number.**

**(a) §3.1 is out of date about the scheme.** BOREAS does not step
explicitly. `advance_segment` solves dE/dt = K(T₀ − T_pcm(E)) **in closed
form inside each enthalpy branch**, so there is no stability limit and
Δt/τ is not a stability number. I report it — `\TauRatioLast` = **61.6**
(Δt 8491 s, τ 138 s) — but the paper should not present it as a CFL-like
constraint. The 7.7 from v0.12 was from a scheme that no longer exists.
`\IntegratorScheme` carries the one-line description. The meaningful
check is the one below.

**(b) The integrator check is a different test from v0.12's.** Error
lives only at branch crossings, so I compared the closed-form step with
2×10⁵ substeps *of itself*, over four branch states and two step lengths:
worst relative enthalpy error **5×10⁻¹¹**. The 3×10⁻¹³ in the text is not
this number and should not be quoted as if it were.

**(c) §2.5 — the return-bend double count is real, and immaterial.**
`calculate_pressure_drop` applies two 180° bends, correct for a hairpin,
and `borehole_pumping` called it twice per well — four bends where the
geometry has one shared bottom turn. Fixed: the downcomer leg is now
called with `n_bends = 0`. Measured effect at R4: **f_pump 4.86 % →
4.84 %**, with η_RTE and e′ unchanged to four decimals. Removing the
minor losses *entirely* only reaches 4.80 %, because friction over
L/D ≈ 3×10⁴ is three orders larger. Worth fixing, not worth a sentence.

**(d) §3.2 confirmed.** The report no longer says "RTE is independent of
the deviation" or "η_storage must be 1"; both were replaced in v2.2 and
I re-checked the current source.

## 8. Partly done

**Done:** the water-balance numbers and the glide-shortfall example, at
R4, in `boreas_water_numbers.tex` — velocity 0.811 m/s per fin tube,
one-way transit 25.7 min (4.3 % of the half-cycle), water-to-latent
capacity ratio 0.76, and the shortfall as **ΔE_dc/ΔE_ch = 0.4758** against
**(T_2d − T_3d)/span = 0.5202**.

**Not done:** SM Table S1, and Video S1 at R4.

---

## Parametric pack, and the last §10 row

`\ParRefAtCutoff*`, at the reference PCM:

| | L = 500 m (the cutoff) | L = 1250 m (reference) |
|---|---|---|
| e′ | **0.2817** | 0.2151 |
| η_RTE | **0.1958** | 0.1447 |
| storage wells per MWe | 71.0 | **37.2** |
| island output | **5.1 MWe** | **9.8 MWe** |

**This is a trade, not a recommendation, and the direction may surprise
the text.** At the cutoff the wells are 31 % better per metre and the
round trip is 35 % better — but you need nearly twice as many of them,
and Grissom has a fixed 363 slots, so going to the cutoff **halves the
island's output**. The cutoff marks where extra depth stops paying *per
metre drilled*; it does not mark where the field should stop.

---

## Macros changed

`boreas_numbers.tex`: **181, unchanged** — no value moved. The twelve
`\Prof*` rate and eps macros changed as described in §1 and §2; they were
never quoted, by your own note. New files: `boreas_orange_numbers.tex`
(27), `boreas_water_numbers.tex` (10), `boreas_ts_numbers.tex` (5).
`par_numbers.tex` 63 → 75.

## Still open

- SM Table S1 and Video S1.
- The grid is good to a few tenths of a point of loss and about 1 % of
  e′, not better (§10 convergence). The README now says so.
- The reference h_m at 380 kJ/kg remains the largest single uncertainty
  in the economics, and it is a materials question.
