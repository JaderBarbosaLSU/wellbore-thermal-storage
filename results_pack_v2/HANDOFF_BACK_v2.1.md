# Results pack v2.1 — back to the BOREAS paper task

*From the model/notebook task, 2026-10-30. Paste this whole file as the
prompt and attach `BOREAS_results_pack_v2.zip`.*

Everything in §1, §2, §2b and §2c of your review is implemented, except
where noted in §4 below. **Every number has been recomputed from scratch at
one discretization.** The ladder, the sweep, the sensitivities and the
plant tables now all come from the same 111 runs.

---

## 0. The rule, unchanged

**Do not retype a single number from this file into the paper.**
`tab_*.tex` are machine-generated; `boreas_numbers.tex` now carries **151
`\newcommand` macros** covering every number in §1–§3 here, including all
the §2b intensive quantities. `\input` both.

---

## 1. Your blocking items

### B3 / 2c.3 — one discretization. I did not do what you asked, and here is why

You asked for a fixed station **count**, default 50. A count cannot be the
rule. 50 stations is 25 m apart in a 1250 m hole and 2 m apart in a 100 m
one, against a well spacing of 1.65 m — so the short runs would fail V11
*by construction*, which is the B1 defect itself.

The rule is now a fixed station **spacing**, `ds = 50 m`; the count follows
as `max(3, round(L/ds))`. It is an input in cell 2.1 of the notebook and it
is part of the cache key. V11 is an assertion, with no fallback search: a
negative strength fails the run and says so.

This buys two things a fixed count could not:

- the validity condition is identical at every length;
- the stations sit at the **same measured depths** whatever *L* is, so
  R_g,eff at a given MD is comparable across the sweep. **V10 now agrees to
  0.00 %** between L = 200 m and L = 2000 m, rather than approximately.

**The grid study says something you should see** (`station_convergence.csv`).
R4 at L = 1250 m:

| stations | Δs | negatives | off-station weight | loss | e′ |
|---|---|---|---|---|---|
| 25 | 50.0 m | **0** | 5.2 % | 52.39 % | 0.2153 |
| 50 | 25.0 m | 130 | 25.8 % | 52.43 % | 0.2151 |
| 100 | 12.5 m | 67 | 27.6 % | 52.51 % | 0.2147 |
| 200 | 6.25 m | 490 | 32.7 % | 52.68 % | 0.2142 |

Refining does **not** converge to a better answer. The loss moves 0.29
points and e′ 0.51 % across an eightfold refinement, while the negative
strengths — the signature of the line-source assumption breaking — appear
from 25 m spacing on. The coarse grid is the well-posed one. 50 m gives
zero negatives at every length from 100 to 2000 m, for R1, R3 and R4.

**This removes the artifact I flagged last round.** The cutoff is now
400–500 m for R4, R4_th45 **and** R4_th70, where v2.0 put the two
inclination variants one interval later in an order that was not monotonic
in angle. Below the kickoff the three inclinations are identical to four
decimals, as they must be. A deep kickoff removes the cutoff altogether:
`R4_kop300` has none within the sweep.

### B1 — fixed, and now agreeing with your prototype

| L (m MD) | v2.0 | **v2.1** | your prototype + THUMS |
|---|---|---|---|
| 100 | 83.2 % | **18.6 %** | 19.3 % |
| 200 | 38.7 % | **21.2 %** | 21.9 % |
| 300 | 27.0 % | **26.5 %** | 26.8 % |
| 1250 | 52.4 % | **52.4 %** | — |

Within 0.7 points everywhere. `q_spread` is now positive by construction
and is defined in the README: the coefficient of variation of the solved
strengths at a station, averaged over stations. The negative values in v2.0
were a symptom of negative strengths, not a quantity.

### B2 — the cutoff, applied as defined

Marginal RTE below half its value over the **first 100 m**, which is
RTE(100) because E(0) = 0. **The 400–500 m cutoff stands**, from a
first-100 m value of 0.262 and a threshold of 0.131. The 300–400 m claim in
the v2.0 note is withdrawn.

### B3 — one configuration, one row

The 1250 m rows that looked duplicated are `superposition` and `flux` — two
ground models at the same geometry, one row each, and the CSV has a
`ground` column saying so. With the discretization fixed they are cleanly
separated: 3.05 / 10.18 / 13.37 % against 2.62 / 10.02 / 12.97 %.

---

## 2. Section 2c.2 — one normalization basis

Every per-MWe quantity is now divided by **W_net = W\*_el,out − W_f,dc**,
the numerator of η_RTE, and by nothing else. Duties come from the realized
energies (ΔE_dc/t_dc, ΔE_ch/t_ch); UA is integrated at those duties on the
realized composite curves; one helper `per_MWe()` is the only divisor in
the code. Target-based values survive only as `*_per_MWe_target`.

- **V12** — storage wells per MWe = t_dc/(e′L) — holds to 1.1e-16.
- **V13** — per-MWe outputs identical at 1 MW and 10 MW — holds to 1.2e-11.

**Your own sanity check was right.** R4 gives **37.2 storage wells per
MWe**, so the 363 Grissom slots give **9.8 MWe**. You estimated 37 and
about 10.

Storage wells, geothermal wells and the total are now three separate
quantities. At R4: 37.2 + 19.0 = 56.2 per MWe.

---

## 3. Section 2b — the intensive tables

All at L = 1250 m MD, t_op = 10 yr, THUMS casing.

### The headline table (`tab_T3e.tex`)

| case | T̄_2d | T_3e | η_ORC | wells/MWe |
|---|---|---|---|---|
| adiabatic (k_r → 0) | 139.5 °C | 93.9 °C | 0.1462 | 18.5 |
| **reference case (R4)** | **116.6 °C** | **90.0 °C** | **0.1401** | **37.2** |
| vertical wells, same slots (R3) | 135.3 °C | 93.2 °C | 0.1451 | 20.5 |
| isolated vertical well (R0v) | 111.4 °C | 89.2 °C | 0.1388 | 47.2 |

Read the last column rather than the η_ORC one. **The formation loss barely
touches the cycle efficiency** — 0.1462 to 0.1401, a 4 % relative spread —
because T_3e is set by the evaporator pinch and the pinch moves little.
What it does instead is **double the well count**: R4 needs 37.2 wells per
MWe against 18.5 adiabatic. The loss is paid in drilling, not in
thermodynamics, and that is the sentence the paper should carry.

### Verified numerically, as you asked

`chk_scale.py` runs R4 at 1 MW and 10 MW. Per-tube flow, e′, RTE, loss,
η_storage and T_2d agree to 2e-16; every per-MWe quantity to 1.2e-11;
N_wells scales 10.000000×. **The 1 MW runs in this pack are valid at any
plant size** — the reporting basis was the only thing that needed changing.

### What is in the pack

| file | what |
|---|---|
| `tab_cycles.tex` | duties as Q/W_net, HTHP input 6.885 MW per MWe, water flows per MWe, per-tube flows |
| `tab_UA.tex` | Q/W_net, UA in kW/K per MWe, ΔT_eff and pinch unchanged |
| `tab_T3e.tex` | the four ladder steps above |
| `tab_ref_results.tex` | 24 rows per storage well at R4, with R3, R0v and the adiabatic limit as columns. No N_w, no field totals |
| `comp_a/b/c.pdf` | composite curves, adiabatic / R4 / R3, on one shared temperature range |
| `intensive.json`, `composites.json` | the same numbers as data |
| `boreas_numbers.tex` | 151 macros |

The pinch is interior in all three: 5.01 K at 16.0 % of duty (adiabatic),
5.00 K at 15.3 % (R4), 5.00 K at 15.9 % (R3).

---

## 4. Still open, and one thing I did not do

**Not done: the parametric maps, SM Table S1, the water-balance numbers,
the glide-shortfall example.** Those live in `BOREAS_parametric.ipynb`,
which is still on the idealized geometry. That is the next piece of work
here and it is the largest remaining item. The T–s figure, the profiles and
the animation **are** already at R4 — they were moved last round.

**Open, and worth a line in the limitations:**

- Shoe depth per well, which sets the reference length. Everything else
  about the casing is settled: top of cement moves the loss by 0.09 points,
  k_mud by 0.17, casing grade by 2.2 % of e′.
- Surveyed trajectories. The S-curve is a build-and-hold reconstruction
  from SPE-2562, SPE-11360 and SPE-15101, not survey data.
- Ground-surface and well-end effects, unrepresented as in R_g.
- The cutoff is resolved to one interval of the L grid and no finer: the
  marginal curve is 0.123 against a threshold of 0.131 over 400–500 m.

**V1 is no longer the v1.3 "7.93 % / 0.2983" reference.** That was the
shielded case under the packing/perimeter field model, and that model no
longer exists — it could not represent where the wells are, which is the
whole reason for the ladder. V7 replaces it: with the casing set back to
the legacy string, R3 and R4 reproduce the prototype at 17.39 % and
70.59 % against 17.4 and 70.7.

**"No non-delivering runs" still holds.** Short isolated wells were never
affected by B1. V6 confirms an isolated 300 m well on the legacy casing
still reports `delivers = False`, e′ = 0, T_2d = 85.5 °C, and none of the
111 runs with the THUMS casing fails to deliver.
