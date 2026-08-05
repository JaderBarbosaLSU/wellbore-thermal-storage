# Fix: energy closure in the wellbore model

**Date:** 2026-08-05 · measured, not argued — every number below came from running both formulations.

---

## The inconsistency, located exactly

v0.1 computes the heat leaving the fluid and the position of the melt front from two
**independent** models:

| | surface it uses | perimeter at the design geometry |
|---|---|---|
| fluid side (`compute_U_i`, `compute_T_re`) | finned: `2πr_e + 2·n_fin·fin_L`, with fin efficiency | **0.4924 m** |
| melt front (`compute_delta2_fast`) | bare concentric cylinder of radius `r_e` | **0.1324 m** |

A factor of **3.72**. The fluid gives up heat through a finned surface; the front absorbs
it through a bare tube. They are not the same energy, and nothing in the code ever compares
them.

The Stefan solver is not *wrong* — `term(x) = ½x²ln x − ¼x² + ¼` with `x = 1+δ/r_e` is the
correct quasi-steady solution **for a bare cylinder**. The fins simply do not appear in it.

## Measured closure error

`Q_in` = ∫ of the heat given up by the fluid. `Q_latent` = `ρ·h_m·V_melt` from the front
position. These must agree to within the neglected sensible heat of the melt layer, which
is 4.7 % here.

| N_wells | k_wall (charge) | `Q_in / Q_latent` |
|---|---|---|
| 24 | 0.45 (legacy defect) | **1.03** |
| 24 | 45 (correct steel) | **2.02** |
| 12 | 0.45 | 1.18 |
| 12 | 45 | 2.16 |

**This is why fixing `k_w` alone broke the model.** With `k_w = 0.45` the fin efficiency
collapses to ≈0.33, so the fluid-side conductance is throttled to roughly what a bare
cylinder can absorb — accidentally matching the front. Two errors of opposite sign, and the
energy balance closes to 3 %.

Restore the steel conductivity and the fins start working on the fluid side while the front
still cannot see them. The gap opens to a factor of two.

## The correction

Only one of those two is a free model. The other follows from conservation:

```
rho_m * h_m * dA_melt/dt = q'(t)          [W per m of tube]
delta = sqrt(r_e**2 + A_melt/pi) - r_e
```

where `q'` is whatever the resistance network actually delivers. Integrating this makes
closure exact by construction instead of something to verify afterwards, and the fins enter
the front advance automatically because `q'` comes from the finned network.

`thums_core/well_marched.py` implements it. **The only change is where δ comes from** — internal
convection, wall and fouling resistance, fin efficiency, the melt-layer conduction
resistance and the NTU segment march are the legacy functions, called unchanged. That
isolates the effect.

## Result

| N | k_wall | legacy ε_PCM | legacy closure error | marched ε_PCM | marched closure error | change in Q_in |
|---|---|---|---|---|---|---|
| 12 | 0.45 | 0.785 | 15.5 % | 0.936 | **0** | +0.7 % |
| 12 | 45 | 0.527 | 53.7 % | 1.124 | **0** | −1.2 % |
| 24 | 0.45 | 0.522 | 3.2 % | 0.542 | **0** | +0.6 % |
| 24 | 45 | 0.286 | 50.5 % | 0.578 | **0** | 0.0 % |

Three things to read from that table.

**Closure is exact** — 10⁻¹⁶, i.e. machine precision, because it is now structural.

**The heat delivered barely moves** (≤1.2 %). The fluid side was never the problem. What was
wrong was the bookkeeping of where that heat went.

**ε_PCM becomes almost insensitive to the `k_w` defect.** Legacy: 0.522 → 0.286, a 45 %
swing. Marched: 0.542 → 0.578, a 7 % swing. Once the melted volume is set by the delivered
energy rather than by an independent front solve, the conductivity error largely stops
mattering for the storage KPIs. The `k_w` argument should still be corrected — it is
unambiguously wrong — but it is no longer capable of wrecking the answer.

## Three further things this formulation gives for free

1. **`T_re` no longer has to be constant since t = 0.** The closed-form Stefan solution
   assumes it; the marched state does not. A varying inlet temperature is handled correctly.
2. **State carries across the cycle.** `WellState.A_melt` is what charging produces and
   discharging consumes. The legacy formulation cannot do this at all — discharge restarts
   from a bare tube in fully solid PCM, which is why round-trip efficiency is not bounded by
   the storage inventory.
3. **Bounds are checkable.** `stefan.check_bounds` raises when the front leaves the borehole
   or the melted volume exceeds the PCM inventory. At N = 12 with correct conductivity it
   fires immediately: ε = 1.12.

## Correction to my earlier audit

I previously listed `h_e = 1e9` at δ → 0 as a defect. **The value is right** — a
zero-thickness melt layer has zero conduction resistance, so a very large `h_e` is the
correct limit. The comment ("Assign a small value, leading to large resistance") is
inverted, and I read the comment rather than the physics.

The real defect is narrower and still stands: *error paths route into that branch*, so a
failed solve is reported as maximum heat transfer rather than as a failure.

## What this changes about sizing

At N = 12 with correct conductivity, ε = 1.12 — the store is over-melted, so **more** wells
are needed, not fewer. Combined with the missing inventory constraint, the corrected sizing
rule is

```
N_wells = max(N_heat_transfer, N_inventory)
```

with `N_inventory` solving `ε(N) = 1`. On the numbers above that lands between 12 and 24
wells for the baseline case — the same order as the published figure, not a different
regime.

## Next

1. Wire `well_marched.march` into `system.run_cycle` behind the existing interface and
   re-run the frozen fixture; the KPIs will move, and the change log records why.
2. Add the inventory constraint to sizing and re-solve.
3. Correct the `k_w` argument. Now safe to do, and now a small effect.
4. Then the closure test becomes a permanent regression test rather than a diagnostic.
