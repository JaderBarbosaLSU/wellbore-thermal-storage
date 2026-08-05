# THUMS code audit — repurposed-well PCM storage, power-to-heat-to-power

**Repository:** `JaderBarbosaLSU/wellbore-thermal-storage`
**Files audited:** `THUMS_Multilayer_DoE_Feb_22_IHTC.ipynb`, `THUMS_Multilayer_DoE_Dec_5.ipynb`, `Pareto_for_THUMS_DoE.ipynb`
**Date:** 2026-08-05 · commit `7cb1357`
**Purpose:** assess readiness for upgrading the IHTC conference paper to a journal submission.

---

## How to read this

Line numbers refer to a flattened export of each notebook (code cells concatenated in order).
Cell numbers are given alongside so you can find each item in Colab.

Findings are graded by what they do to the *paper*, not by how hard they are to fix:

- **Tier 1 — invalidating.** A referee who finds this asks for the results to be recomputed.
- **Tier 2 — must fix.** Wrong or unjustifiable, but the conclusions might survive correction.
- **Tier 3 — will be asked about.** Defensible with an added sentence, indefensible in silence.

Every Tier 1 item below was verified directly against the source, independently of the
automated pass that first surfaced it. Tier 2 and 3 items come from a systematic read;
they are specific and line-referenced, but treat them as claims to confirm, not verdicts.

**The overall picture.** The code runs and produces plausible numbers, and the physics
scaffolding is sound in structure. But "it runs" is doing a lot of work here. The model
contains several paths that respond to failure by *silently returning the most optimistic
possible answer* — which is precisely the failure mode your own workflow document warns
about, appearing in your own code. Combined with the two Tier 1 items, my assessment is
that the numbers currently in the conference paper cannot be carried into a journal
submission unchanged.

---

## Tier 1 — invalidating

### 1.1 Charging and discharging use wall conductivities that differ by a factor of 100

**Where:** Cell 17, `evaluate_Q_ratio_ch`, line 1669 (charge) vs `evaluate_Q_ratio_dc`, line 1707 (discharge). Same defect at lines 3527 and 4206.

`time_profiles_melt` takes the tube-wall and fin conductivity as its sixth positional
argument, `k_w`. The discharge path passes `k_w` correctly. The charge path passes
`k_m_l` — the *liquid PCM* conductivity — into that slot:

```python
# charging (line 1669)
time_profiles_melt(times_ch, geom_par_vector,
    T_4c_kelvin, N_lay, T_m_lay, k_m_l, Rf_i_prime, ...)
#                                 ^^^^^ this is k_w's position

# discharging (line 1707)
time_profiles_melt(times_dc, geom_par_vector,
    T_3d_kelvin, N_lay, T_m_lay_dc, k_w, Rf_i_prime, ...)
```

With `k_w = 45 W/m·K` (steel) and `k_m_l = 0.45 W/m·K`, the charging model runs with a
steel wall and steel fins that conduct like liquid PCM. Two consequences compound:

- wall resistance `(r_i/k_w)·ln(r_e/r_i)` inflates by 100×;
- fin efficiency, through `fin_m = sqrt(2·h_e/(k_w·fin_t))`, collapses from ≈0.97 to ≈0.33,
  dragging overall surface effectiveness `eta_o` down with it.

Net effect: charging heat transfer is throttled by close to an order of magnitude relative
to discharging, in a model whose entire purpose is to compare the two.

**Why it invalidates rather than merely biases.** `find_N_wells_for_Q_ratio_ch` sizes the
well field against this crippled charging conductance, and the discharge mass flow is then
derived from that well count. So the error does not cancel — it propagates into
`optimal_N_wells_ch`, into pumping power on both sides, into `m_m_eff`, and into every
round-trip efficiency in the DoE table. Meanwhile the borehole profile *figures* in Cell 22
(line 3607 onward) recompute with `k_w = 45`, correctly. **The published figures and the
published table therefore come from two different models.**

**Fix:** pass `k_w` in both. Then re-run everything. There is no way to patch the results.

---

### 1.2 The 128-run design is Resolution III — every main effect is aliased with two-factor interactions

**Where:** Dec_5, Cell 22, `build_128run_fractional_16factor`, lines 3096–3143.

Nine of the sixteen factors are assigned to *two-factor* interaction columns:

```python
H = A*B;  I = A*C;  J = A*D;  K = A*E;  L = A*F;  M = A*G;  N = B*C;  O = B*D;  P = B*E
```

I enumerated the full defining relation (all 2⁹−1 = 511 words). The word-length
distribution is `{3: 12, 4: 27, 5: 30, 6: 48, 7: 88, 8: 95}`. **Twelve words of length
three**, so the design is 2⁽¹⁶⁻⁹⁾ of **Resolution III**, not IV. The twelve are:

```
ABH  ACI  ADJ  AEK  AFL  AGM  BCN  BDO  BEP  HIN  HJO  HKP
```

What this means for the factor ranking, in your own variable names:

| Reported "main effect" | Is inseparable from |
|---|---|
| `T_m_C` | `rho_m_s·h_m`, `rho_m_l·r_e`, `cp_m_s·fin_t`, `cp_m_l·fin_L`, `k_m_s·t_ch`, `k_m_l·DT_3C_2C` |
| `h_m` | `T_m_C·rho_m_s`, `r_e·N_lay`, `fin_t·T_sink_C`, `fin_L·T_source_C` |
| `r_e` | `T_m_C·rho_m_l`, `h_m·N_lay` |
| `N_lay` | `rho_m_s·rho_m_l`, `h_m·r_e` |
| `t_ch` | `T_m_C·k_m_s` |
| `DT_3C_2C` | `T_m_C·k_m_l` |
| `T_sink_C` | `rho_m_s·cp_m_s`, `h_m·fin_t` |
| `T_source_C` | `rho_m_s·cp_m_l`, `h_m·fin_L` |

Note the trap in the top rows: `T_m_C` is aliased with `rho_m_s·h_m` while `h_m` is
aliased with `T_m_C·rho_m_s` — the same interaction. If those are your two largest
effects, they cannot be separated from each other by any argument about effect sparsity.

The only acknowledgement in the code is a comment at line 3141: *"This is just one
possible alias structure."* The generators are never written down as a defining relation,
no alias table is printed, and the words *resolution*, *generator* and *alias* appear
nowhere in either notebook. A reviewer who knows DoE will ask for the defining relation in
the first round, and the answer is currently unfavourable.

**The good news: this is free to fix.** I verified that assigning the same nine factors to
*odd-order* columns — `ABC, ABD, ABE, ABF, ABG, ACD, ACE, ACF, ACG` — gives a Resolution IV
design (shortest word length 4, 43 words of length 4) in exactly the same 128 runs. Same
cost, main effects clear of all two-factor interactions. Any standard tool
(`pyDOE.fracfact`, JMP, Minitab) would have produced Resolution IV by default.

**Fix:** regenerate the design with odd-order generators and re-run the 128 cases, or fold
over the existing design to reach Resolution IV at 256 runs. Publish the defining relation
and the alias table either way.

---

### 1.3 The column labelled `RTE` is not a round-trip efficiency

**Where:** Cell 21, line 3580 — and again at line 4224 as `weighted_RTE`.

```python
'RTE': RTE_real * m_m_eff   # Store this as a <weighted RTE>
```

The actual round-trip efficiency is computed at lines 3564–3570 and then **discarded**.
What is stored, tabulated, exported to CSV and plotted as `ε_RTE` in the Pareto notebook
is an efficiency multiplied by a PCM *utilization fraction*.

That product is not an efficiency of any control volume. It falls when utilization falls
even though electricity in and electricity out are unchanged. Any reader seeing 0.30 in
that column will read it as 30 % round-trip, which it is not.

**Fix:** store both, and label them separately. This one is cheap — but it means every
reported efficiency number in the conference paper needs to be checked against which
quantity it actually was.

---

## Tier 2 — must fix

### 2.1 Two functions read their energy target out of the global namespace

**Where:** lines 1687 and 1717.

```python
D_E_out_HP = find_N_wells_for_Q_ratio_ch.__globals__['D_E_out_HP']  # Temporary access, should be passed as argument
D_E_in_ORC = find_m_dot_d_well1_for_Q_ratio_dc.__globals__['D_E_in_ORC']
```

The values that *are* passed as arguments (call sites 3453, 3477, 4157, 4177) are never
used. Inside `run_parametric_simulation`, `D_E_out_HP` and `D_E_in_ORC` are **locals** —
so each DoE point computes its own storage duty and then sizes the well field against the
*baseline* cell's duty instead. Every KPI in the DoE CSV is converged against the wrong
target. If Cell 21 was never run, this raises `KeyError`, which is swallowed by the
`except` at 1789 → `return 1e9` → the bracket test fails → the solver returns the upper
bound, 100 wells, with only a printed warning.

This is also what makes the notebooks impossible to refactor safely, and it is the single
strongest argument for extracting the solver into an importable module.

### 2.2 Failure paths return the *most optimistic* answer available

Three places, all silent:

- **line 2145 / 2236:** when `delta` is negligible, `h_e = 1e9`. The comment reads
  *"Assign a small value, leading to large resistance"* — the comment is inverted.
  1e9 W/m²·K is essentially **zero** resistance, i.e. maximum heat transfer.
- **lines 2450–2457:** `except Exception: delta_new = 0.0`, with the diagnostic `print`
  commented out. δ = 0 routes straight into the `h_e = 1e9` branch above. So *any*
  exception in the melt-front solve makes that segment report the best possible
  performance, undetectably.
- **`compute_delta2`, lines 2031–2049 and 2095–2099:** when the root lies beyond
  `delta_max`, or the iteration fails to converge, the function returns `0.0` — again with
  the warning commented out. These are the segments where the PCM is most nearly
  exhausted, and they report as if the front had not moved at all.

Your own workflow document argues that a guard which fails loudly is worth far more than
one which carries on with an approximation. That principle applies here more sharply than
anywhere else in the codebase.

### 2.3 `compute_T_re` and `compute_U_i` use different melt-layer conductances

**Where:** line 2141 vs line 2233.

```python
h_e = k_m / r_e / log_term     # compute_T_re  — cylindrical shell
h_e = k_m / delta              # compute_U_i   — plane slab
```

So the interface temperature that drives the Stefan front is computed from one thermal
network, and the heat drawn from the fluid from a different one. They diverge as melting
proceeds — roughly 22 % at δ = 10 mm, and worse beyond. A converged `Q_ratio = 1.00`
therefore does not certify an energy balance; it certifies agreement between two
inconsistent models.

Related: at line 2230–2236 the `if log_term > 1e-9` branch assigns the plane form in
*both* arms, and if the condition is false `h_e` is never bound at all — a latent
`UnboundLocalError` which finding 2.2 would then swallow.

### 2.4 The geothermal gradient is computed, stored, and never used

**Where:** `T_geo` computed at line 2414, written to the results DataFrame at line 2531,
and referenced nowhere else in the file.

The formation is effectively adiabatic and the PCM is implicitly initialised at its own
melting temperature everywhere. With `T_geo_baseline = 18 °C` and `geo_grad = 0.05 K/m`
over a 1524 m well, the PCM actually starts between 18 and 56 °C, not at `T_m`. The
sensible preheat is on the order of 30 % of the latent capacity and is neither charged for
nor credited.

The only storage loss anywhere in the model is the flat 5 % at line 3179
(`D_E_out_HP = 1.05 * D_E_in_ORC`) — independent of duration, geometry, and ΔT. Note also
the commented-out alternative on the next line quoting a 13.6 % Alamooti correction: two
mutually inconsistent loss figures, one of them arbitrary.

Since `geo_grad` is a declared factor in the DoE, this means one of your design factors
**cannot influence any KPI** — its estimated effect is structurally zero.

### 2.5 The melt front is geometrically unbounded

**Where:** `delta_max = 0.5` m at line 2656, against a borehole radius of 0.089 m and
`r_e = 0.021` m. `A_melt = π((r_e+δ)² − r_e²)` at line 2514 has no cap, and the per-well
volume multiplies by `num_tubes` (line 3532) without accounting for overlap between the
hairpins inside a single borehole.

So the model can melt PCM that does not exist, `m_m_eff` can exceed 1, and the "weighted
RTE" of 1.3 can exceed the unweighted one. Nothing checks `V_melt ≤ V_well`.

### 2.6 `float('inf')` fallbacks pass every downstream guard

`hp_cop = float('inf')` (lines 93, 378) and `rank_eff = float('inf')` (491, 592) are
returned when a denominator is small **or when any enthalpy is NaN**. Every downstream
check tests only `np.isnan` (3240, 3408, 3508, 4106–4116, 4194), and `np.isnan(inf)` is
`False`.

Trace it through: one CoolProp failure → `hp_cop = inf` → `W_dot_C_HP = Q/inf = 0` →
`W_dot_el_in = 0` → `RTE_real` becomes a large finite number, written to the results table
as if converged.

### 2.7 Bisection failures are indistinguishable from convergence

`find_N_wells_for_Q_ratio_ch` returns a bracket bound when the root is not bracketed
(1804–1807), the collapsed midpoint (1828–1829), or the last midpoint after 100 iterations
(1848–1850). Same in the discharge solver (1921–1929, 1960–1962). No convergence flag
reaches the results dict. A design needing 250 wells silently returns 100 — the upper
bracket — and its RTE is tabulated as valid.

### 2.8 Runs that fail are dropped from the effect estimates without a count

`run_parametric_simulation` wraps every thermodynamic call in bare `except Exception: pass`
(eight occurrences in Dec_5 between lines 2643 and 2688) and returns all-NaN KPIs. The DoE
loop writes that row to CSV indistinguishably from a converged run. Both effect routines
then use pandas `.mean()`, which skips NaN silently — so the low and high groups stop
having 64 runs each, the contrast stops being orthogonal, and nothing is printed.

Nowhere does either notebook report `df.isna().sum()` or a converged-run count.

---

## Tier 3 — will be asked about

- **The response-surface cell is rank-deficient.** Dec_5 line 3454: with all corner runs at
  coded ±1, `X1**2` and `X2**2` are identical columns. `lstsq` returns the minimum-norm
  solution, splitting one estimable curvature contrast arbitrarily between the two
  quadratic coefficients. Worse, `top2 = ["T_m_C", "r_e"]` and `r_e` is column A·C, so the
  "interaction" column `X1*X2` **is exactly the `rho_m_l` main-effect column** — I verified
  this. The contour and 3-D surface figures should be withdrawn, not repaired.
- **`top2` is a hard-coded placeholder.** Line 3429 carries your own comment,
  *"# example – replace with what your table shows"*.
- **Curvature is never tested.** The markdown promises the test; the code prints two means
  and leaves "far off" to the eye. And since the model is deterministic (no `np.random`
  anywhere), the 8 centre points are bit-identical replicates giving exactly zero pure
  error — no valid F test can be built from them. One centre run carries the same
  information; the other seven are wasted.
- **Effect magnitudes are not comparable across factors.** Coded half-ranges span 5.4×,
  from ±5.26 % (`t_ch`) to ±28.57 % (`fin_L`). A factor given a wider window automatically
  ranks higher. The three temperature factors are coded in **Celsius**, whose zero is
  arbitrary, so their percentage ranges are physically meaningless.
- **No significance screening.** No standard errors, no half-normal plot, no Lenth PSE.
  Factors are ranked by `abs()` and separated by an unstated threshold.
- **The Pareto notebook contains no Pareto front.** No dominance test, no non-dominated
  sorting anywhere in its 181 lines — "Pareto" refers only to a Pareto bar chart of effect
  magnitudes. If the manuscript describes a Pareto front of designs, that analysis does not
  exist in this code.
- **The multi-KPI consensus double-counts.** `Pareto_for_THUMS_DoE.py` lines 53–61 average
  seven KPIs of which at least four are algebraic products of the others.
- **Provenance is not recorded.** The producer writes to `/content/...` (Colab scratch,
  wiped on disconnect); the consumer reads a bare relative filename. No version string, no
  timestamp, no commit hash, no seed, no model-identifier column. Two notebooks define a
  `run_parametric_simulation` with identical 16-argument signatures and different physics,
  so nothing in the CSV records which one produced it. **This is the finding I would fix
  first**, because it costs almost nothing and because every other fix on this list makes
  it more urgent.
- **Pumping power carries no pump or motor efficiency** (line 1640 is hydraulic power),
  while the compressor and turbine both do. Understates parasitic load by roughly 45 % on
  both sides of the RTE ratio.
- **Charging pumping power uses the wrong temperature pair** — lines 3503–3504 evaluate
  properties between `T_4c` (borehole inlet) and `T_4a` (source-side loop), not `T_2c`
  (borehole return). Discharge uses the matched pair correctly.
- **IHX effectivenesses are outputs, never checked for admissibility.** ε > 1 or ε < 0
  produce a finite, high COP that feeds straight into the RTE with no flag.
- **Cycle efficiencies are ideal.** No isentropic efficiency step in either the HTHP or the
  ORC; `Turb_eff`/`Comp_eff` are applied outside, but `rank_eff` is exported raw as "ORC
  efficiency" and used in `DE_per_well_eff_kWh` with no turbine or generator loss at all.
- **Cumulative energy is integrated on a log grid starting at t = 1 s**, over an integrand
  that peaks at t → 0 (where δ → 0 and `h_e = 1e9`). The interval [0, 1 s] is omitted and
  trapezoid over a convex decaying integrand overestimates. No grid-independence study
  exists for `n_segments = 100`, the 40-point time grid, `delta_maxiter = 20`, or
  `delta_tol = 1e-5`.
- **`delta_tol` is a length passed as a dimensionless residual tolerance** (line 2453 vs
  the comparison inside `compute_delta2`), so its effective strictness varies by orders of
  magnitude across the time sweep.
- **PCM layers are stratified along the tube coordinate, not depth.** `z_lay` spans
  `0…L_tube = 2·L_well`, but the hairpin descends and returns — and the code knows this,
  computing `depth = z if z <= L/2 else L − z` ten lines later for `T_geo`. The return leg
  is assigned a fictitious continuation of the layer sequence instead of re-traversing the
  same physical layers in reverse. For the multilayer-PCM concept that is the paper's
  central variable, this is worth checking carefully.
- **`from re import T`** at line 2591 binds the global name `T` to a regex flag.
- **`np.trapezoid`** (line 2581) requires NumPy ≥ 2.0 and fails on 1.x — exactly the kind
  of Colab-vs-elsewhere difference your workflow note anticipates.

---

## What I would do, in order

**Before any new physics.**

1. **Fix 1.1** (`k_w`) and re-run the baseline. Everything downstream depends on it.
2. **Fix 1.3** (`RTE` labelling) and check which quantity each number in the conference
   paper actually was.
3. **Regenerate the design at Resolution IV** (1.2) — same 128 runs, and publish the
   defining relation.
4. **Turn the silent fallbacks (2.2) into exceptions.** Then re-run and see what breaks.
   Expect breakage; that is the point. Every case that now fails loudly was previously
   reporting an optimistic number.
5. **Add provenance**: a version string and commit hash stamped into every CSV and every
   figure caption. Cheap, and it is what lets us trust any of the re-runs.

**Then the refactor.** Extract the solver into `thums_core/` — cycles, borehole, DoE, post —
with the notebooks reduced to case setup plus plotting. Two reasons beyond tidiness: the
`__globals__` back-doors (2.1) cannot be removed while everything lives in one namespace,
and finding 1.1 exists *only* because the same call is written out twice with different
arguments. A single `charge()`/`discharge()` pair sharing one parameter object makes that
class of bug impossible rather than merely absent.

**Then the tests.** Three that would have caught most of the above:

- energy closure: `∫Q dt` against `ρ·h_m·V_melt`, per segment and globally;
- `V_melt ≤ V_well` and `m_m_eff ≤ 1`;
- charge/discharge symmetry: a zero-loss case must return `RTE = 1` to within tolerance.

You already have `verification/harness.py` and frozen reference fixtures in this repo from
the DBHE work. The same pattern applied to THUMS is the single highest-value structural
addition, and it is the thing a journal referee is most likely to reward.

---

## A note on what this audit is not

I have not run the code. Everything above is from reading it, with the design-resolution
and aliasing claims verified by re-executing the design generator in isolation. Once the
repository is set up so I can execute and commit, the natural next step is to reproduce
the baseline, then work down the list — fixing, re-running, and recording what each fix
does to the headline numbers. That record is itself a useful thing to have when the
referees ask what changed between the conference paper and the journal version.
