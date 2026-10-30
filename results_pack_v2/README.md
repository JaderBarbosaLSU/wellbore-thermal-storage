# BOREAS results pack v2.1

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

The method treats each station as an infinite line source in the plane
normal to the trajectory, so a station must represent a segment long
compared with the distance to the neighboring wells. Fix the **count** and
that condition changes with L: 50 stations is 25 m apart in a 1250 m hole
and 2 m apart in a 100 m one, against a well spacing of 1.65 m. The short
case then returns negative strengths and R_g,eff goes negative with them —
which is exactly the B1 defect.

Fix the **spacing** and the condition is the same at every length, the
stations sit at the same measured depths whatever L is, and V10 becomes
meaningful rather than approximate: R_g,eff at MD = 50 m now agrees to
**0.00 %** between L = 200 m and L = 2000 m.

### The convergence test (`station_convergence.csv`)

R4, L = 1250 m:

| stations | ds (m) | ds/B₀ | negative strengths | off-station weight | loss | e′ |
|---|---|---|---|---|---|---|
| 25 | 50.0 | 30.3 | **0** | 5.2 % | 52.39 % | 0.2153 |
| 50 | 25.0 | 15.2 | 130 | 25.8 % | 52.43 % | 0.2151 |
| 100 | 12.5 | 7.6 | 67 | 27.6 % | 52.51 % | 0.2147 |
| 200 | 6.25 | 3.8 | 490 | 32.7 % | 52.68 % | 0.2142 |

Read it the other way round from a normal grid study. **Refining does not
converge to a better answer here**: the loss moves 0.29 points and e′ 0.51 %
across an eightfold refinement, which is small, while the negative
strengths — the signature of the line-source assumption breaking — appear
from 25 m spacing on. The coarse grid is the well-posed one. 50 m gives
**zero negatives at every length from 100 to 2000 m**, for R1, R3 and R4.

Two structural results make most of the grid moot, and both are measured,
not assumed:

- **Vertical fields are exactly block diagonal in station.** For parallel
  vertical wells the nearest point on a neighbor is at the same depth, so
  the off-station kernel weight is 0.00 % and every station block is the
  same matrix. R_g,eff is then independent of both the count and L: R1
  gives 15.783 m K/W and R3 gives 3.230 at every length and every count
  tested. The solver now exploits this, which is what made the disc rung
  tractable at all.
- **The S-curve is block tridiagonal.** Neighbors' nearest points land at
  most one station away, carrying 5.2 % of the kernel weight at the
  reference. This is the only case where the count does anything.

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
| `grissom_geometry.pdf` | slots by KOP group, trajectories, spacing vs MD |
| `grissom_Rg.pdf` | the shielding profiles |
| `grissom_ladder.pdf` | the bars |
| `grissom_sweep.pdf` | e′, RTE and marginal RTE against L, with the beyond-shoe band shaded |

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

The result worth reading off them: **the formation loss barely touches the
cycle efficiency** (η_ORC 0.1462 adiabatic against 0.1401 at R4, a 4 %
relative spread, because T_3e is set by the evaporator pinch) and instead
**doubles the well count**, 18.5 → 37.2 storage wells per MWe. The loss is
paid in drilling, not in thermodynamics.

`chk_scale.py` verifies the model is exactly linear in plant size: R4 at
1 MW and 10 MW gives identical per-tube flows, e′, RTE, loss and T_2d to
2e-16, identical per-MWe quantities to 1.2e-11, and N_wells scaling
10.000000×. The 1 MW runs in this pack are valid at any plant size.

## 7. Still owed

`BOREAS_parametric.ipynb` — the maps, the T_m study (to be capped at
160 °C), the h_m study and the N_lay confound — is still on the idealized
geometry, as are SM Table S1, the water-balance numbers and the
glide-shortfall example. The T–s figure, the profiles and the animation are
already at R4.
