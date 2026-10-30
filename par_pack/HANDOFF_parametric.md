# Parametric study on the Grissom field (R4) — to the BOREAS paper task

*From the model/notebook task, 2026-10-30. Attach `BOREAS_par_pack.zip`.*

This is the §5.6 item: `BOREAS_parametric.ipynb` moved to the reference
case. **198 points, all at R4** — measured slot positions, S-curve
trajectories, THUMS casing, 50 m station spacing, t_op = 10 yr. The
idealized `n_wells_field` / `packing` / `fan_angle_deg` axes are retired,
not pending.

Same rule as the results pack: numbers are generated, never retyped.
`par_results.csv` is the source and `par_numbers.tex` carries the macros.

---

## 0. The one thing that made this affordable

R_g,eff(MD) depends on the well positions, the trajectories, L, the bore
radius, the rock properties and the field age. It does **not** depend on
the PCM (T_m, h_m, N_lay, span), on the tube count, or on the geothermal
gradient — the gradient sets the rock *temperature*, not the resistance
between the well and it.

So one field solve per storage length serves every point that shares it,
and the whole gradient axis rides for free. Five of the eight grids sit at
the reference length and share a single solve.

**Consistency check, in the notebook:** the cascade cell at T_m = 150 °C,
ΔT = 55 K reproduces the ladder's R4 reference, e′ = 0.2151, exactly.

---

## 1. The result that dominates everything else: h_m

This axis is new and it is not close.

| h_m (kJ/kg) | e′ | storage wells per MWe | loss |
|---|---|---|---|
| 200 | 0.065 | 123.6 | 0.765 |
| 260 | 0.116 | 69.0 | 0.658 |
| 320 | 0.167 | 47.8 | 0.581 |
| **380** | **0.215** | **37.2** | **0.524** |
| 440 | 0.259 | 30.9 | 0.481 |
| 500 | 0.298 | 26.9 | 0.447 |

At T_m = 150 °C, a **2.5× change in h_m moves e′ by 4.6×** and the well
count by the same factor in reverse. Along T_m at fixed h_m, e′ spans only
1.6×.

The amplification is worth a sentence in the paper: e′ is *more* than
proportional to h_m, because a bigger latent inventory does not just store
more — it holds the store nearer its melting plateau for longer, so the
driving difference to the rock is smaller and the loss falls with it,
0.765 → 0.447. **Latent heat buys capacity and cuts losses at the same
time**, which is not true of any other parameter here.

---

## 2. The cascade, with T_m capped at 160 °C

Capped as asked. The ceiling is the secondary loop, not the PCM: water
saturates at 179.9 °C at 10 bar. `auto_P` raises the pressure just enough
to keep the water liquid with 20 K of margin, so the T_m axis measures the
store; at the reference point it returns 10 bar unchanged.

**The two optima are different cells**, and the paper should not quote one
as "the optimum":

| | T_m | ΔT_3C,2C | e′ | η_RTE |
|---|---|---|---|---|
| reference | 150 °C | 55 K | 0.215 | 0.145 |
| best e′ | 140 °C | 40 K | **0.223** | 0.162 |
| best η_RTE | 120 °C | 40 K | 0.205 | **0.182** |

Capacity per metre and round-trip efficiency are traded against each other
by the grading. Both optima sit at ΔT = 40 K rather than the reference's
55 K, so the grading is the one cascade parameter clearly off its best.

**Six cells are infeasible and hatched.** The reason is physical and is
the bottom of the cascade, not the top: at a low T_m with a wide grading
the bottom layer is too cold to drive the ORC at the required pinch.

---

## 3. Well geometry, L capped at the 5000 ft shoe

**A unit inconsistency the paper task should settle.** The ladder's sweep
marks a *provisional* shoe at 1250 m MD; this review asked for the maps to
stop at the 5000 ft shoe, 1524 m MD. They are two statements of the same
unknown. The map grid includes 1250 m so the reference sits on it and ends
at 1524 m. Pick one for the paper.

**Depth is no longer free.** On the idealized geometry e′ was flat in L to
about 1 % and depth was a pure capital lever. At R4 it falls steadily,
because past the kickoff each extra metre of measured depth buys less true
vertical depth and sits in rock the neighbors no longer share.

**One fin tube is a pumping failure, not a storage one** — and this is the
sharpest new finding after h_m:

| L (m MD) | pumping / gross | η_RTE |
|---|---|---|
| 250 | 0.014 | +0.194 |
| 750 | 0.390 | +0.107 |
| 1250 | 1.671 | +0.019 |
| 1524 | **2.842** | **−0.047** |

The store still delivers heat — T_2d stays above T_3d — but with a single
tube the flow per tube is three times the reference, the pressure drop
goes as roughly its square, and by 1524 m the borehole pumping is 284 % of
the gross output. The cell is printed, not clipped: a negative round trip
is the honest statement of that design, and it is a different failure mode
from the non-delivery flag.

---

## 4. The N_lay / span confound, at R4

N_lay cannot be swept at fixed grading: span = ΔT_cascade (N−1)/N, so
raising N widens the span at the same time.

| | e′ spread over N_lay = 3…15 |
|---|---|
| span free | 46.19 % |
| span pinned at 48.9 K | 16.36 % |

**65 % of the apparent layer-count effect is the span moving with it.**
With the span held the curve still rises but has mostly flattened by
N = 9, the reference. The layer count is a real but second-order lever;
the span is the first-order one.

---

## 5. What is in the pack

| file | what |
|---|---|
| `par_results.csv` | all 198 points, one row each, with `feasible` |
| `par_numbers.tex` | macros for every number above |
| `map_cascade.pdf` | T_m × grading: e′, η_RTE, wells per MWe |
| `map_hm.pdf` | h_m × T_m: e′, wells per MWe, loss |
| `map_geometry.pdf` | L × N_lay and L × n_ft |
| `map_gradient.pdf` | L × geothermal gradient |
| `confound_Nlay.pdf` | free against pinned span |
| `par_grids.py`, `par_runner.py`, `par_runs.json` | the definitions and the cache, fetched by the notebook |

---

## 6. Still owed

SM Table S1, the water-balance numbers and the glide-shortfall example.
They belong with the plant tables of the results pack rather than with the
maps, and they have not been moved to R4.

k_rock is not swept. It would change R_g,eff, so it needs one field solve
per value rather than one per length; nothing in the paper currently asks
for it.

---

## 7. One request

These maps have not been checked against anything. Every round so far has
caught a real method bug, and in both cases it was because this task
recomputed something independently rather than reading my code — B1 from
your prototype wrapper, the symmetrization from your `symtest.py`. Nothing
has played that role for the parametric results.

If you can spare it, the cheapest useful check is the **h_m row at
T_m = 150 °C**: six points, no field solve needed if you hold R_g,eff
fixed at 0.915 m K/W, and it is the result the paper will lean on hardest.
