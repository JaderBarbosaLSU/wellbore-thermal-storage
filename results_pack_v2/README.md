# BOREAS results pack v2.4

Generated from `pack_runs.json` by `pack_writer.py`, `cutoffs.py`,
`mktab_ladder.py`, `mknumbers.py` and `verify.py`. Nothing in this pack is
typed by hand.

**Do not retype a number from here into the paper.** `tab_ladder.tex` is
machine-generated; `boreas_numbers.tex` carries a `\newcommand` for every
number the prose quotes. `\input` both.

---

## 1. The one thing that changed everything: a fixed station spacing

v2.0 chose the superposition station count **per run**, as the first count
that solved with no negative strengths. Two runs of the same sweep were
therefore discretized differently — at L = 500 m the reference case landed
on 25 stations and the 45° variant on 6 — and that moved e′ by 2.8 %
against the 0.4 % the inclination was worth. The cutoff came out
non-monotonic in angle as a result.

v2.1 fixes **one station spacing, `ds` = 50 m**, for every run of the
ladder, the L sweep and every sensitivity. The count follows as
`max(3, round(L/ds))`. It is an input, in cell 2.1 of the notebook, and it
is part of the cache key: a run at a different discretization is a
different run, not a replacement for one already there.

### Why a spacing and not a count

Fix the **count** and the station spacing changes with L: 50 stations is
25 m apart in a 1250 m hole and 2 m apart in a 100 m one. Fix the
**spacing** and the stations sit at the same measured depths whatever L is,
which is what makes V10 a real check rather than an approximate one:
R_g,eff at MD = 50 m agrees to **0.00 %** between L = 200 m and L = 2000 m.

### The matrix is not symmetric, and v2.1 symmetrized it

Each row says *the rise at my location, on well i at station s, due to
well k*, and the nearest point on k seen from there is at some other
station s′. The reverse entry is about a different target point at a
different distance, so **G ≠ Gᵀ**: measured, ‖G−Gᵀ‖/‖G‖ is 0.25 at 50 m
spacing and 0.90 at 25 m.

v2.1 solved (G + Gᵀ)/2, which averages rows that are not about the same
place. It produced 130 negative well strengths at 25 m spacing and 569 at
8 m, and v2.1 read those as the line-source assumption breaking under
refinement. **They were an artifact of the symmetrization.** Solved as
assembled there are no negative strengths at any spacing tested, and the
v2.1 claim that "the coarse grid is the well-posed one" is **withdrawn**.

Two solver consequences, both in `Ground.profiles`:

- the system is solved as assembled, `A = G + (R_near + G_ii) I`;
- above `_DIRECT_MAX` it uses **BiCGSTAB** with a Jacobi preconditioner,
  not CG. CG needs a symmetric positive-definite matrix, and this one is
  neither symmetric nor diagonally dominant — the diagonal is 0.567 m K/W
  against a largest off-diagonal row sum of 3.73, because near the surface
  all 363 wells are inside the 80 m cutoff. The v2.1 code comment claiming
  diagonal dominance was wrong and is corrected.

### The convergence test, redone (`station_convergence.csv`)

R4, L = 1250 m, no symmetrization:

| Δs | stations | negatives | off-station weight | loss | e′ |
|---|---|---|---|---|---|
| 104 m | 12 | 0 | 25.9 % | 51.75 % | 0.2177 |
| **50 m** | **25** | **0** | **5.2 %** | **52.42 %** | **0.2151** |
| 25 m | 50 | 0 | 25.8 % | 52.42 % | 0.2150 |
| 12.5 m | 100 | 0 | 27.6 % | 52.50 % | 0.2147 |
| 6.25 m | 200 | 0 | 32.7 % | 52.53 % | 0.2146 |

Δs = 50 m is adopted on **convergence**: the loss moves 0.001 points and e′
0.05 % to the next finer grid, and 0.12 points across the whole eightfold
refinement. 104 m is too coarse (0.67 points). At L = 400 m the same
refinement leaves about 0.4 points of residual drift, so the honest
statement is that the discretization is worth up to a few tenths of a point
of loss — not that it is exact.

Two structural results, both measured:

- **Vertical fields are exactly block diagonal in station.** The
  off-station kernel weight is 0.00 % and every block is the same matrix,
  so R_g,eff is independent of both the spacing and L: R1 gives 15.783 m
  K/W and R3 gives 3.230 everywhere. G is also symmetric there, which is
  why the symmetrization went unnoticed. The solver exploits this, which is
  what makes the disc rung computable at all.
- **The S-curve is the only case where the grid does anything**, because
  there a neighbor's nearest point can land in an adjacent station.

### `off_station_weight`, defined — and not evidence of anything

The quantity is: of the total kernel weight assembled at a given station
spacing, Σ G_ik over all interaction entries, the fraction whose target
and source sit at **different** station indices.

**It is not monotone in the spacing and must not be read as a property of
the field.** Measured, R4:

| | 104 m | 50 m | 25 m | 12.5 m | 6.25 m |
|---|---|---|---|---|---|
| L = 1250 m | 25.9 % | 5.2 % | 25.8 % | 27.6 % | 32.7 % |
| L = 400 m | 0.1 % | 0.6 % | 22.0 % | 22.7 % | 25.8 % |

The low values at 50 m come from where the station grid happens to fall
relative to the shortest-distance points, not from the field being nearly
block diagonal there. The v2.2 README used the 5.2 % in that sense and
that reading is **withdrawn**; Δs = 50 m is justified by convergence
alone. The quantity is kept only for the one case where it is exact and
structural: for parallel vertical wells it is **0.00 %** at every spacing,
which is a proof and not a measurement.

`q_spread` is the coefficient of variation of the solved well strengths at
a station, `std(q)/mean(q)` over the field, averaged over stations. It
measures how unequally the field sheds heat — edge wells carry more — and
is positive by construction once V11 holds. A negative value in v2.0 was a
symptom of negative strengths, not a quantity.

---

## 2. What moved from v2.0, and why

| | v2.0 | v2.1 | why |
|---|---|---|---|
| R4 loss, L = 100 m | 83.2 % | **18.6 %** | B1. The prototype gives 19.3 % |
| R4 loss, L = 200 m | 38.7 % | **21.2 %** | B1. Prototype 21.9 % |
| R4 loss, L = 1250 m | 52.4 % | **52.4 %** | unchanged |
| R1 / R2 / R3, L = 1250 m | 3.1 / 10.2 / 13.4 % | **3.0 / 10.2 / 13.4 %** | unchanged |
| cutoff, R4 | 400–500 m | **400–500 m** | unchanged |
| cutoff, θ = 45° / 70° | 500–600 m | **400–500 m** | the station count, not physics |
| wells per MWe | storage + geothermal, gross basis | **storage only, net basis** | §2c.2 |

**The cutoff is 400–500 m**, and it no longer depends on the inclination:
R4, R4_th45 and R4_th70 all give 400–500 now, where v2.0 had the two
variants one interval later. That flip was the per-run station count and
nothing else. The criterion is applied as defined — marginal RTE below half
its value over the **first 100 m**, which is RTE(100) because E(0) = 0.

**"No non-delivering runs" still holds.** Short isolated wells were never
affected by B1: V6 confirms an isolated 300 m well on the legacy casing
still reports `delivers = False` with e′ = 0 and T_2d = 85.5 °C, and no run
in the 111 with the THUMS casing fails to deliver.

---

## 3. One normalization basis

Every per-MWe quantity is divided by **W_net**, the realized net output of
the field, and by nothing else:

    W_net = W*_el,out − W_f,dc

which is the numerator of η_RTE. v1.4 used three bases at once — UA from
the *target* budget duties over the realized *gross* output, wells over the
realized gross, N_geo over the *target* — so no two per-MWe numbers in the
same table were comparable. At the reference case this is not a detail: R4
loses 52 % of the charge, so anything per *target* MWe understates the real
plant by about a factor of two.

Implemented as: duties from the realized energies (ΔE_dc/t_dc, ΔE_ch/t_ch);
UA integrated at those duties on the realized composite curves; one helper
`per_MWe(x) = x / (W_net/1000)` through which every per-MWe output passes.
Target-based values survive only as `*_per_MWe_target`, for traceability.

---

## 4. Verification, V1–V13

All pass. `verify.py` regenerates `verification.json`.

| # | check | result |
|---|---|---|
| V1 | charge closes, Clausius gate positive | closure 3.4e-16, +1.55e+11 |
| V2 | G_ii degenerates to the shipped isolated R_g | 5.7e-7 m K/W apart |
| V3 | one well in the field solver = isolated well | 1.7e-6 relative |
| V4 | equal flux is the optimistic model, R_flux > R_temp | 3.343 > 3.230 |
| V5 | k_rock → 0: zero loss, η_storage = 1 | exact |
| V6 | short isolated legacy well reports non-delivery | delivers False, e′ = 0 |
| V7 | legacy-casing rungs regress to the prototype | 17.39 % (17.4), 70.59 % (70.7) |
| V8 | R4 short lengths against the prototype wrapper | 18.6 % (19.3), 26.5 % (26.8) |
| V9 | disc < double row < real slots, all vertical | 3.05 < 10.18 < 13.37 % |
| V10 | R_g,eff at MD = 50 m independent of L | 0.00 % between L = 200 and 2000 m |
| V10b | the same at L = 100 m, where the 3-station floor bites | 1.28 % |
| V11 | no negative strengths in any run | 111 runs, 0 |
| V12 | wells per MWe = t_dc/(e′L) | 1.1e-16 relative |
| V13 | per-MWe outputs identical at 1 MW and 10 MW | 1.3e-10 worst |

**V1 is not the v1.3 "7.93 % / 0.2983" reference.** That number was the
shielded case under the packing/perimeter field model, and that model no
longer exists — it cannot represent where the wells are, which is the whole
reason for the ladder. V7 is its replacement as a regression: with the
casing set back to the legacy string, R3 and R4 reproduce the prototype.

**V10b is reported rather than hidden.** With ds = 50 m the stations sit at
25, 75, 125 … m for every length, so L = 200 and L = 2000 sample the same
depths exactly. L = 100 m is the one length where the floor of three
stations bites (ds = 33.3 m there), so its stations are offset from every
other run's, and the 1.28 % is that offset.

---

## 5. Contents

| file | what |
|---|---|
| `ladder_results.csv` | 111 runs, one row per (configuration, ground model), now carrying `n_stations`, `ds_station`, `off_station_weight` and `n_negative` |
| `tab_ladder.tex` | `\input` this |
| `boreas_numbers.tex` | 71 `\newcommand` macros for the prose |
| `cutoffs.csv` | the criterion, its threshold and the bracketing interval, per step |
| `station_convergence.csv` | the grid study of §1 |
| `verification.json` | V1–V13 with margins |
| `profiles.npz` | R_g,eff(MD), the flux variant, TVD and q_spread for **all eleven** steps and sensitivities, plus the spacing used |
| `tab_cycle_states.tex` | COP, bleed fraction, flash quality, IHX effectivenesses, six pressures |
| `boreas_profile_numbers.tex` | 19 macros for section 4.3, at R4 |
| `cutoff_robustness.csv` | the cutoff at three station spacings |
| `casing_sensitivity.csv` | the three annulus unknowns at R3 and R4 |

**Figures are one PDF per panel** (no panel letter or suptitle in the
image, same size across a row, legend in the first only):
`BOREAS_grissom_geometry_{a,b,c,d}`, `BOREAS_grissom_Rg_{a,b}`,
`BOREAS_grissom_ladder_{a,b}`, `BOREAS_grissom_sweep_{a,b,c}`,
`comp_{a,b,c}{1,2}`, `BOREAS_prof_{Tfluid,melt,Ui}_{a,b}` and
`BOREAS_prof_integrals_{a,b}`. The integrals axis is **per fin tube**.

Every generated `.tex` carries the SHA-256 of the `pack_runs.json` it was
built from, and **V16** checks every stamp. v2.2 shipped a `tab_ladder.tex`
that was the v2.1 file while the note claimed otherwise; stale LaTeX is
still valid LaTeX, so only a hash catches it.

The `field` column for R0v and R0i is **363 isolated wells** — the same
count, with the shielding switched off. It is a bound, not an arrangement.

Rows with `beyond_shoe = True` (L = 1600 and 2000 m) are kept in the CSV
and shaded in `grissom_sweep.pdf`: they extrapolate the geometry past the
provisional casing shoe, which is the one genuinely blocking unknown.

---

## 6. The intensive tables (section 2b)

All at L = 1250 m MD, t_op = 10 yr, THUMS casing. Four cases: the adiabatic
limit (k_rock → 0), the reference case R4, vertical wells on the same slots
R3, and one isolated vertical well R0v.

| file | what |
|---|---|
| `tab_cycles.tex` | duties as Q/W_net, HTHP input per MWe, water flows per MWe, per-tube flows |
| `tab_UA.tex` | Q/W_net, UA in kW/K per MWe, ΔT_eff and pinch |
| `tab_T3e.tex` | mean T_2d, T_3e and η_ORC at the four steps |
| `tab_ref_results.tex` | 24 rows per storage well, four columns, nothing extensive |
| `comp_a/b/c.pdf` | composite curves, adiabatic / R4 / R3, shared temperature range |
| `intensive.json`, `composites.json` | the same numbers as data |

The result worth reading off them: **only the ORC conversion is
insensitive to the loss.** η_ORC goes 0.1462 (adiabatic) → 0.1401 (R4), a
4 % relative spread, because T_3e is set by the evaporator pinch and the
pinch hardly moves. Everything else does move, and the loss is paid
**twice**:

- in wells — 18.5 → 37.2 storage wells per MWe;
- in charged electricity — the HTHP draws 6.9 MW per MWe against about
  3.1 adiabatic, so 2.2× the electricity per MWh delivered, and
  η_RTE falls 0.321 → 0.145.

An earlier version of this note said the loss is "paid in drilling, not in
thermodynamics". That is half right and the half it drops is the larger
one: the round trip more than halves.

`chk_scale.py` verifies the model is exactly linear in plant size: R4 at
1 MW and 10 MW gives identical per-tube flows, e′, RTE, loss and T_2d to
2e-16, identical per-MWe quantities to 1.2e-11, and N_wells scaling
10.000000×. The 1 MW runs in this pack are valid at any plant size.

## 7. What is here, and what is still owed

**Delivered in v2.4:** the T–s figure at R4
(`BOREAS_Ts_combined_water.pdf`, drawn at the *realised* outlet, 116.54 °C
against the provisional 146.11 °C); the four section 4.3 profile figures
at R4, one PDF per panel, with the rates **per fin tube** and the depth
axis labelled TVD; the section 3.4 water balance and glide shortfall
(`boreas_water_numbers.tex`); and the section 10 checks
(`boreas_orange_numbers.tex`).

`BOREAS_parametric.ipynb` is at R4 — that moved in the parametric pack
and this README said otherwise until now.

**Still owed:** SM Table S1, and Video S1 at R4. The animation helper
exists in `boreas_plots` but no animation is shipped in this pack.

### A correction worth reading before the rates are used

v2.3 divided the per-tube heat rate by `n_ft` on top of `q_prime`, which
is **already per unit length of fin tube**. The plotted rates were a
factor of three low and the half-cycle means were arithmetic means over a
logarithmic grid, which weights the first minute as heavily as the last
two hours. Both are fixed, and **V17** gives the rates the independent
reference they lacked: `n_ft` times the exact per-tube integral must equal
ΔE per well, and it does to 2×10⁻¹⁶.

### On the integrator

`advance_segment` solves dE/dt = K(T₀ − T_pcm(E)) **in closed form inside
each enthalpy branch**, so there is no explicit-stability limit and
Δt/τ = 61.6 at the last step is not a stability number. The error lives
at branch crossings only; measured against 2×10⁵ substeps of itself, the
worst relative enthalpy error is 5×10⁻¹¹.

### On the grid

Doubling the segment count moves the loss by 0.17 points and e′ by
−0.7 %; doubling the time levels moves them by 0.18 points and −1.0 %. The
shipped grid is therefore good to a few tenths of a point of loss and
about 1 % of e′ — not better, and the pack says so rather than implying
convergence it does not have.
