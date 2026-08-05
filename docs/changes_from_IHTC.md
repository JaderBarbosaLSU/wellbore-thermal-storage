# Changes from the IHTC model

One row per change, with the mechanism and the measured effect. Every number
below came from running both formulations on this machine; none is quoted from a
previous session.

**Reference point.** `DT_3C_2C = 55 K`, `N_lay = 9`, the design geometry, 10 h
charge / 10 h discharge, 1 MW ORC. Two other points (`40/12`, `30/20`) are
carried through the tables to show the trends are not particular to one case.

**Where the "before" lives.** `verification/reference/thums_v0_baseline.json` has
been regenerated with the marched front, as instructed. The legacy values it
replaced are preserved in git at **tag `v0.2`**, and transcribed in full in
§6 below. `git show v0.2:verification/reference/thums_v0_baseline.json` recovers
the file byte for byte.

**Reproducing any row.** `Case(front="legacy")` runs the v0.1 physics unchanged;
`Case(front="marched")` is the v0.2 default. Both go through the same
`run_cycle`, so a row is a two-line script.

---

## 1. Headline numbers

| | v0.1 (IHTC) | v0.2 marched | v0.2 marched + correct `k_w` |
|---|---|---|---|
| `N_wells` | 24.123 | 21.651 | **12.787** |
| binding constraint | *(none — only one existed)* | heat transfer | **inventory** |
| `eps_pcm` | 0.5200 | 0.5597 | **0.9992** |
| `eta_rte` | 0.43027 | 0.42924 | 0.42047 |
| `eta_rte_nopump` | 0.43484 | 0.43484 | 0.43484 |
| `E_well` [MWh] | 2.1674 | 2.4149 | 4.0890 |
| `rho_E` [kWh/m³] | 78.29 | 87.23 | 147.7 |
| `f_pump` | 0.01508 | 0.01843 | 0.04733 |
| `cop_hp`, `eta_orc` | 2.95633, 0.236851 | unchanged | unchanged |
| energy closure error | 3–54 % | 2e-16 | 2e-16 |

**The one-sentence version.** The energy-closure fix moves the *sizing* a long
way and the *efficiency* hardly at all: `N_wells` falls by 47 % once the
conductivity defect is also corrected, while `eta_rte` moves by 2.3 %.

That asymmetry is the useful result. Round-trip efficiency is set by the HTHP and
ORC cycles, which none of this touches — `eta_rte_nopump` is identical to twelve
significant figures across all three columns. What the storage model determines
is how much hardware the cycle needs, and that is where v0.1 was wrong.

---

## 2. Change log

### C1 — melt front advanced from the delivered heat

`thums_core/stefan.py`, `thums_core/well_marched.py`.

v0.1 took the heat leaving the fluid from a **finned** outer surface (perimeter
0.4924 m, with fin efficiency) and the melt-front position from an **independent**
bare-cylinder Stefan solution (perimeter 0.1324 m). Nothing compared the two.
Replaced by `rho*h_m*dA/dt = q'`, integrated from whatever the resistance network
delivers.

| effect | before | after |
|---|---|---|
| `Q_in / Q_latent` (charge, correct `k_w`) | 2.021 | 1.000 |
| `N_wells` | 24.123 | 21.651 |
| `eps_pcm` | 0.5200 | 0.5597 |
| `eta_rte` | 0.43027 | 0.42924 |

Note what did *not* change: `Q_in` itself moves by ≤1.2 %. The fluid side was
never the problem. The bookkeeping of where that heat went was.

### C1a — the closed-form solve has a *second*, independent defect

**This was not in `FIX_energy_closure.md`, and it changes the argument.**

That document attributes the closure failure to the fin/bare surface mismatch,
and reports the marched closure of 1e-16 as the evidence that the fix works. The
1e-16 is not evidence of anything: `advance` integrates the same `q'` that
`closure_error` then differences, so closure is exact by construction. It is an
arithmetic check, not a physics check.

The real test is to remove the fins. With `num_fins = 0` the fluid-side and
front-side surfaces become identical, so on the stated diagnosis the discrepancy
should collapse to the ~5 % neglected sensible heat. It does not:

| `num_fins` | `k_w` | legacy `Q_in / Q_latent` |
|---|---|---|
| 24 | 0.45 | 1.032 |
| 24 | 45 | 2.021 |
| **0** | **0.45** | **0.642** |
| **0** | **45** | **0.642** |

With the fins removed the legacy front *over*-melts by 56 %, in the opposite
direction, and it does so essentially independently of `k_w`.

The cause is the other assumption in the closed form: it applies the *current*
`T_re` as though it had held since `t = 0`. Measured at mid-well over the
charging window, `T_re − T_m` runs **0 → 8.62 K**. It is not constant, it starts
at zero, and the closed form is being handed its largest value and told to apply
it retroactively to the whole history.

The algebra of the solve is fine — given a genuinely constant `T_re` it returns
the exact root, residual 1e-15 against
`½x²ln x − ¼x² + ¼ = k_m ΔT t / (ρ h_m r_e²)`. It is the *premise* that fails.

So v0.1 contains **three** errors in this chain, not two:

1. finned fluid surface vs bare front surface — front under-melts;
2. frozen `T_re` applied retroactively — front over-melts, ~1.56× ;
3. `k_w = 0.45` collapsing fin efficiency to ≈0.33 — throttles the fluid side.

The reported 3 % closure at the legacy design point is a **triple**
cancellation, not a double one. This makes the case for the marched formulation
stronger rather than weaker: there is no version of the closed-form solve that is
correct here, because the condition it is derived under never holds.

*Evidence: `verification/verify_fins_off.py`, `verification/verify_frozen_Tre.py`.*

### C2 — the well count now has two constraints

`thums_core/sizing.py`. v0.1 sized only on *can the wells absorb the energy in
time?* The inventory number existed in the notebook as a printed
`N_wells_ideal = 12.2` and was never used.

`N_wells = max(N_heat_transfer, N_inventory)`:

| point | `N_heat` | `N_inventory` | binding | `N_wells` |
|---|---|---|---|---|
| 55/9, legacy `k_w` | 21.651 | 9.986 | heat transfer | 21.651 |
| 40/12, legacy `k_w` | 18.242 | 9.597 | heat transfer | 18.242 |
| 30/20, legacy `k_w` | 16.683 | 9.134 | heat transfer | 16.683 |
| **55/9, correct `k_w`** | 8.014 | **12.787** | **inventory** | **12.787** |
| **40/12, correct `k_w`** | 6.991 | **13.664** | **inventory** | **13.664** |
| **30/20, correct `k_w`** | 6.504 | **14.759** | **inventory** | **14.759** |

The constraint that binds *flips* when the conductivity is corrected. This is why
the missing constraint stayed invisible for so long: with `k_w = 0.45` it never
bound. It is also why correcting `k_w` alone made the model worse — the heat
transfer criterion collapses to ~8 wells, at which point the store melts more PCM
than it contains.

### C3 — melt state carries across the cycle

`thums_core/system.py`. v0.1 discharges from a bare tube in fully solid PCM
regardless of what charging produced, so the recovered energy is not bounded by
the stored energy. The marched discharge starts from `r_ch.state`.

New measurable: `eta_storage = E_discharged / E_stored`.

| point | `eta_storage`, legacy `k_w` | `eta_storage`, correct `k_w` | `eps_pcm` at end of discharge |
|---|---|---|---|
| 55/9 | 0.9531 | 0.9033 | 0.026 → 0.097 |
| 40/12 | 0.9515 | 0.8454 | 0.032 → 0.155 |
| 30/20 | 0.9516 | 0.7815 | 0.035 → 0.219 |

With the legacy conductivity `eta_storage` lands near 0.953, which is close to
the `loss_surplus = 0.05` that v0.1 *assumed* — the assumption was roughly right
for the wrong reasons. With the correct conductivity it falls to 0.78–0.90 and
the assumed flat 5 % is no longer defensible. `loss_surplus` should become an
output of the model, not an input. **Not done in v0.2** — it inverts the energy
bookkeeping chain, and is the natural next change.

### C4 — latent inventory uses one density in both directions

`thums_core/config.py`, `PCM.rho_latent`.

Exposed by C3, and unreachable before it. `march` selected `rho_l` when melting
and `rho_s` when freezing, mirroring the (correct) directional choice for `k_m`.
With state carried across the cycle that returns `rho_s/rho_l` = **1.069** times
more latent energy than went in — a 6.9 % round-trip gain from bookkeeping alone.

The front is now carried as a *solid-equivalent* area, so `V_melt / V_well` is
the melted **mass** fraction and stored energy equals available energy. Confirmed:
a discharge run to saturation returns `E_dc / E_stored = 1.0000`.

Not modelled, and now explicit: the 6.9 % volume change on melting is real. Its
consequence is that the liquid layer is ~2.3 % thicker in radius than
`delta_from_area` reports, so the melt-layer conduction resistance is mildly
under-predicted.

### C5 — segments that exhaust locally no longer warm the fluid for free

`thums_core/well_marched.py`.

Found while wiring C3. When a segment's melt is exhausted mid-discharge, the
front was clamped but the fluid still left at the temperature the *unclamped*
heat flow implied. Downstream segments then saw an inlet temperature no energy
balance supported. The limit is now applied inside the segment loop and the fluid
temperature is taken from the limited heat flow.

The effect is not small:

| | before C5 | after C5 |
|---|---|---|
| discharge flow ratio | 2.398 | 0.965 |
| `pumping_dc` [kW] | 110.7 | 8.70 |
| `f_pump` | 0.1205 | 0.01843 |
| `eta_rte` | 0.38505 | 0.42924 |

Spuriously warmed fluid starves the downstream segments, so the solver compensates
with flow, and pumping power goes as roughly the cube of it. An intermediate
version of this work concluded that `ratio_bracket = (0.2, 2.0)` was too narrow
and widened it; with C5 in place the root sits at 0.965 and the original bracket
is correct. It has been restored.

Still not modelled: sensible heat in the PCM. A segment whose latent store is
exhausted stops contributing rather than cooling further, so `Q_rejected` is a
lower bound on what a sensible-heat model would recover.

### C6 — failures are failures

Carried over from v0.1 and now exercised. The legacy bisections return the
bracket bound and report it as an answer; the driver cell `continue`s. Both
sizing loops now raise `ThumsConvergenceError` with the last iterate. The
discharge solve additionally reports the stored energy alongside the shortfall,
so an unreachable target reads as a statement about the store rather than about
the solver.

---

## 3. What has *not* changed

Worth stating plainly, because it bounds what the fix can be blamed for.

- **The cycles.** `cop_hp = 2.9563304095`, `eta_orc = 0.2368513203`,
  `eta_rte_nopump = 0.4348350503` are identical across every row. The HTHP and
  ORC state points are untouched.
- **`Case(front="legacy")` reproduces the frozen v0.1 fixture** to 4e-11 on every
  KPI. The residual is CoolProp/NumPy version drift (this machine runs Python
  3.10 / NumPy 2.2.6; the fixture was generated on 3.11 / 2.4.4), not the
  refactor. This is a permanent regression test.
- **The `k_w` defect is still the default.** `charge_uses_wall_conductivity=False`.
  The "correct `k_w`" column is a switch, not the shipped configuration.
- **`delta_max = 0.5 m` against an 0.0889 m borehole radius** still warns rather
  than raises, because `Case.strict = False`.

## 4. Still absent from the model

Unchanged from the audit, repeated because the numbers above should not be read
as final: no formation heat loss; no natural convection in the melt; no buoyancy
head in the 1,524 m loop (±7.5 bar, and asymmetric between charge and discharge);
no sensible heat in the PCM; no volume change on melting.

## 5. Open, in order

1. Make `loss_surplus` an output (C3). It is now measurable and the assumed value
   is wrong by 5–22 % once `k_w` is corrected.
2. Correct the `k_w` default. Safe now — the two-constraint sizing absorbs it —
   and it is the single largest mover in the table.
3. Set `Case.strict = True` and fix `delta_max`.
4. Add PCM sensible heat, which bounds C5's residual.
5. Port the DoE onto the unified core using the Resolution IV design.

## 6. The v0.1 values this file replaced

Recoverable with `git show v0.2:verification/reference/thums_v0_baseline.json`.

| | 55/9 | 40/12 | 30/20 |
|---|---|---|---|
| `eta_rte` | 0.43027101049765343 | 0.419756645691732 | 0.3945257609555505 |
| `eta_rte_nopump` | 0.4348350503020928 | 0.4348350503020928 | 0.4348350503020928 |
| `cop_hp` | 2.956330409531823 | 2.956330409531823 | 2.956330409531823 |
| `eta_orc` | 0.23685132027873143 | 0.23685132027873143 | 0.23685132027873143 |
| `eps_pcm` | 0.5200473147142692 | 0.5983257937144157 | 0.6429616691151382 |
| `E_well` | 2.167417633092469 | 2.6166121373585245 | 2.9098255357678453 |
| `rho_E` | 78.29462518690391 | 94.52108510425177 | 105.11296771041152 |
| `N_wells` | 24.123429921001843 | 19.982154265412923 | 17.96861933430409 |
| `f_pump` | 0.015081686346236986 | 0.05011372132110686 | 0.1361728651415322 |
