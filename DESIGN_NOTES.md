# Design notes

Why the model is the way it is, and not the way it looks like it should be.

Each note records one decision: what was done first, why it was wrong, what
replaced it, and **what the change was worth** — because a correction whose
effect is never measured is indistinguishable from a preference.

These are referenced from the code as `see DN-n`, and in full in the project
report. The superseded code is in version control; it is not in `thums.py`,
and it is not in the student notebook.

| | note | worth |
|---|---|---|
| DN-1 | The state variable must be enthalpy, not melted area | 8.2 % of cycle energy |
| DN-2 | The closing conductance is not $U_iP_i$ | 4–32 % on the segment conductance |
| DN-3 | Explicit stepping is inadmissible off the plateau | divergence past 8×10¹⁰ K |
| DN-4 | The history is stamped where the state is exact | 2.36 h of label error |
| DN-5 | The discharge march is mirrored; results are not | 0.677 in ε |
| DN-6 | The state mirrors with the cascade | 48.9 K of enthalpy datum |
| DN-7 | Sizing at CSS is a simulation, not a solve | the problem has no root |
| DN-8 | Only inlets may be prescribed | 6.9 percentage points of deviation |
| DN-9 | The capacity criterion must resolve the cascade | N 11.573 → 9.574 |
| DN-10 | The melt front is bounded by the cell, not the borehole | factor 2.9 in δ |
| DN-11 | The deviation cancels out of round-trip efficiency | 1.2 % of apparent η |
| DN-12 | `layer_map` uses a node convention *(open)* | ~0.1 pp on the deviation |
| DN-13 | The conduction path was on the wrong side of the front while freezing | a spurious factor of 7 in $U_i$ |
| DN-14 | A single-phase cell has no front, so it cannot have a side | $U_i(t\to0)$ 823 → 260 |

---

## DN-1 — The state variable must be enthalpy, not melted area

**First approach.** Carry the melt front as a layer thickness δ, or equivalently
a melted cross-sectional area, and advance it with the heat the resistance
network delivers. Two variants were implemented: a closed-form Stefan solve
(*Formulation A*, through v0.1 and the conference paper) and an energy balance
on the melted area (*Formulation B*, v0.2–v0.3).

**Why it fails.** Ask such a model for the PCM temperature in a control volume
and the only available answer is $T_m$ — that is all the state encodes. The
moment a control volume completes its phase change, the network still reports a
positive heat flow, and the model must either grow the melted area past the PCM
the volume contains (melting material that is not there, ε > 1) or clip and
discard the surplus (breaking its own energy balance). No refinement of the
marching scheme repairs either: the deficiency is in the state.

**Replacement.** One scalar per segment: the enthalpy per unit tube length,
measured from a datum of fully solid material at the local $T_m$. Subcooled
solid, two-phase and superheated liquid are three branches of one single-valued
curve, and ε ∈ [0,1] becomes structural rather than enforced.

**Worth.** At CSS, 82 of the 100 segments finish the charge fully molten —
mean ε = 0.9706, peak superheat 15.6 K — and 8.2 % of the energy passing through
the store each cycle traverses the sensible branches (4.3 % superheat, 3.9 %
subcooling). That is the energy an area-based state must fabricate or throw
away, and it is being clipped in most of the exchanger at the moment the store
is fullest. A secondary and larger effect is on the *distribution*: sensible
capacity is self-levelling, and the spread in local melt fraction at equal
exchanger size falls from 0.407 to 0.000 on a charge from cold.

**Note on a retired number.** Earlier drafts quoted "the discarded heat reached
148 % of the heat delivered during discharge". That figure was produced by
Formulation B under the flow-reversal defect of DN-6, cannot be reproduced now
that B is deleted, and has been withdrawn.

---

## DN-2 — The closing conductance is not $U_iP_i$

**First approach.** Close the segment equation with the overall coefficient
times the perimeter, $K = U_i P_i$ — the obvious choice, and the one the
fluid-side and PCM-side balances appear to license.

**Why it fails.** Deriving the segment equation from both control volumes
instead gives the effectiveness form. The two agree only as NTU → 0. More
importantly, $U_iP_i$ carries no capacity-rate limit: as the wall conductance
improves it permits a segment to extract more heat than the stream can supply,
and in the limit to drive heat from cold to hot.

**Replacement.**

    K = ṁ c_p (1 − e^−NTU) / Δz        NTU = 2π r_i U_i Δz / (ṁ c_p)

**Worth.** Over the charge, NTU runs 0.081 → 0.588, so the ratio
(1 − e^−NTU)/NTU runs 0.96 → 0.76: using $U_iP_i$ overstates the segment
conductance by 4 % to 32 %, with the largest error early, when the melt layer is
thin and the exchanger is at its best.

---

## DN-3 — Explicit stepping is inadmissible once sensible capacity exists

**Why it was not caught earlier.** On the latent plateau the equation is
$dE'/dt = K(T_f - T_m)$ with $T_m$ fixed — the capacity is *infinite* and
explicit Euler is exact for any step. Area-based formulations live entirely on
the plateau, so they never met the constraint. It appears only when the PCM is
allowed to leave $T_m$, i.e. only after DN-1.

**The constraint.** With sensible capacity the branch becomes a relaxation with
time constant $C/K$, and explicit Euler needs Δt < C/K. The time grid here is
logarithmic — the natural choice for a diffusion-limited process — and its final
step exceeds the limit by a factor of 7.7.

**Replacement.** Integrate each branch in closed form and split the step at
branch crossings. Unconditionally stable, monotone, and energy closure becomes
structural rather than checked.

**Worth.** At Δt/τ = 1.9 explicit Euler oscillates; at 7.7 — the ratio the grid
actually reaches — it diverges, alternating in sign and passing 8×10¹⁰ K within
twelve steps. The branchwise-exact scheme agrees with a 200 000-substep explicit
reference to 3×10⁻¹³.

**Related trap.** The branch test must not let a state come to rest exactly on a
boundary. If the test uses closed intervals, a state landing on $E'_{\rm lat}$ is
still classified two-phase, the crossing time evaluates to zero, and the segment
is pinned there forever. The symptom is not a crash but an apparent success:
every segment reporting ε = 1 and zero superheat. The state is nudged past the
boundary by 10⁻⁹ max(E'_lat, 1) for this reason.

---

## DN-4 — Frames are stamped where the state is exact

**What was wrong.** `march_h` computed the state at the *start* of a step, used
it to advance the segment, then appended that pre-step state to the history —
stamped with the time at the *end* of the step. It also appended `E` after the
update while `A_melt`, `T_pcm` and `delta` in the same frame came from before
it, so one frame did not describe one instant.

**Why it mattered here and not on a uniform grid.** On a logarithmic grid the
last step spans 8 491 s. The final frame of the record was therefore the store
as it stood 2.36 h — almost a quarter of the charge — before its label, and the
true end state was never recorded at all.

**How it surfaced.** A figure and the text beneath it disagreed: mean ε at the
end of discharge read 0.356 from the map and 0.1545 from the march.

**Replacement.** `t_hist` stamps the instants at which the *state* arrays are
exact; the profile arrays in the same frame are those evaluated *from* that
state. There are `len(t)+1` frames — one at t = 0 and one at the end of every
accepted step — the last closed by one extra read-only sweep
(`segment_profile`).

**Worth.** No reported quantity, because nothing in the model reads the history.
Every plotted profile in the results sections, because everything does.

---

## DN-5 — The discharge march is mirrored; results returned to the caller are not

**What was wrong.** The discharge marches from the opposite end of the tube, so
its segment index runs backwards relative to the charge. The sizing and CSS
loops un-mirror when they close the cycle — which is why the well count and the
deviation were right — but nothing un-mirrored before plotting, so every
discharge figure was reflected in depth.

**Replacement.** `unmirror_march` returns every discharge result in depth
indexing. The consequence worth stating is conceptual: **there is one cascade,
not two.** $T_m$ belongs to the PCM at a physical position and does not care
which way the fluid is flowing; the mirrored copy is an implementation detail of
the march and no longer leaves the routine that needs it.

**Worth.** Compared without un-mirroring, the two edges of the CSS melt map
differ by 0.677 in ε; compared correctly, by nothing. The map closing on itself
is now printed as a check.

---

## DN-6 — The state mirrors with the cascade

**What was wrong.** The discharge cascade was correctly mirrored, but the melt
state handed across the half-cycle boundary was not. Since $E'$ is measured from
a datum of *solid at the local $T_m$*, and the local $T_m$ had just been
reversed, the datum moved by up to 48.9 K.

**How it surfaced.** Only once $T_{\rm pcm}$ was plotted rather than ε: the model
reported PCM at 176.6 °C against a charging inlet of 160.0 °C — hotter than any
fluid that had ever touched it.

**Replacement.** Mirror both. Building the discharge cascade as the exact
reverse of the charge cascade, rather than by a separate `layer_map` call, also
removes a one-layer (6.111 K) mismatch that arises whenever `n_segments` is not
a multiple of `N_lay`.

**Worth.** $N_{\rm wells}$ unchanged — the field is sized on the charge, which was
never affected. Flow ratio +2.06 %, η_RTE −0.20 %, and the CSS deviation
−3.34 % → −5.74 %: much larger on the cycle than on a single pass, as one would
expect of an error that compounds across repeated half-cycle boundaries.

---

## DN-7 — At cyclic steady state, sizing is a simulation and not a solve

**First approach.** Size the field by requiring it to absorb the specified
energy within the charging window, and solve for $N_{\rm wells}$.

**Why it fails at CSS.** The state returns to itself and the store is adiabatic,
so η_storage ≡ 1 *identically*, for every well count and every flow — computed
values agree to eight decimal places. The charge and discharge requirements
therefore collapse into one equation, which fixes the ratio of the two flow
rates and says nothing whatever about the exchanger size. There is no root, and
the degeneracy is structural rather than numerical.

Two corollaries. The assumed loss surplus λ is unattainable: the CSS charge pins
at 1/(1+λ) of the requirement for every well count, because the model has no
loss path to absorb it — on the first cycle the surplus is simply *parked in the
store* as residual melt, which works only because that cycle starts empty.
And setting λ = 0 removes the contradiction without restoring the ability to
size.

**Replacement.** Specify the hardware — $N$ from latent heat alone, both flows
pinned by their glides — march to CSS, and report the deviation of delivered
energy from target as a result rather than driving it to zero.

**Worth.** It converts an unsatisfiable constraint into a reported output, and
it dissolves an older objection: that the discharge flow could not be pinned to
its glide because delivered energy then saturates below the requirement for any
well count. Under this framing that saturation is not a failure — it is the
answer.

---

## DN-8 — Only the inlets may be prescribed

**What was wrong.** The discharge side was closed on the **outlet**: $T_{2d}$
pinned at the top layer's melting temperature, and the inlet obtained by
subtracting the glide.

**Why it fails.** Only the two exchanger inlets are boundary conditions on the
march; both outlets are results of the heat transfer. An area-based model can
blur this, because a saturated control volume holds its surface at $T_m$ and the
outlet is bounded by the cascade whatever happens. After DN-1 it cannot: the
store finishes each charge superheated, so at the start of the following
discharge the water leaves at 157.93 °C — nearly 8 K *above* the temperature it
had been pinned to. The charge side was no better: $T_{2c}$ prescribed at
105.00 °C against a realised flow-mixed mean of 104.38 °C, ranging over 20 K
within the window.

The implied inlet approach was also set by nobody:
ΔT = ΔT_3C,2C / N_lay = 6.111 K, exactly one layer width, moving silently with
the number of layers. Against the 10 K the charge receives, the discharge was
given 61 % of the driving difference at every layer — and in the linear-glide
idealisation the water arrives at each layer's exit face *at* that layer's
melting temperature, pinching to zero nine times over, while the charge never
comes closer than 3.889 K.

**Replacement.** Prescribe the two inlets symmetrically:

    charge:    T_4c = T_m,top    + ΔT_4C,M     (superheat, drives melting)
    discharge: T_3d = T_m,bottom − ΔT_M,1D     (subcooling, drives freezing)

$T_{2d}$ is demoted to a provisional estimate. Since η_ORC follows it, one
correction pass re-evaluates the plant at the outlet the store actually
produced.

**Worth.** At the symmetric ΔT_M,1D = 10 K the field moves from 5.69 % below the
1 MWe target to 1.23 % above it — **most of what had been reported as a
shortfall of the store was a property of the closure.** Verified to be a pure
reparameterisation: at ΔT_M,1D = 6.111 K the v0.5a sources taken from version
control and run beside the current build agree to the last bit.

**Not adopted: 9.30 K**, where the deviation crosses zero. Choosing an approach
because it zeroes the deviation converts the one honestly reported output back
into a satisfied constraint. 10 K is adopted because it is symmetric with the
charge.

**One pass suffices, and that is checkable.** Across a sweep in which the inlet
moves 9 K and the deviation moves twelve percentage points, the realised outlet
moves less than 0.3 K and the second pass shifts it by at most 0.03 K. The shift
is returned as `dT_2d_pass` so it is reported rather than assumed.

---

## DN-9 — The capacity criterion must resolve the cascade

**What was wrong.** The energy one well can absorb was computed with the *top*
layer's melting temperature for the whole store, and with an arbitrary
discharge-side reference. That is the no-cascade limit, $N_{\rm lay}=1$, of the
correct expression — and it was not a bound: the model exceeded it by 1.6 %,
4.822 against 4.744 MWh, with a mean superheat of 13.81 K rather than the 10 K
assumed.

**Replacement.** Use the volume-weighted mean melting temperature over the
cascade, ⟨T_m⟩ = T_m,top − (N_lay−1)/(2 N_lay) ΔT_glide, and include the solid
subcooling the discharge actually reaches.

**Worth.** N_capacity 11.573 → 9.574, so the rate criterion became the binding
one and N_wells moved 11.573 → 11.255. Two conclusions reversed with it: the fin
optimum moved from 16 to about 20 fins and became shallow, and the melt-fraction
spread previously quoted as 1.321 → 0.078 turned out not to be a like-for-like
comparison.

**Reporting rule that came out of it.** `N_capacity` is a *lower bound*, not a
result, and must never be read across a parameter sweep on its own, because
`N_wells` is a maximum and which criterion binds can change mid-sweep.

*(Superseded in practice by DN-7: at CSS neither criterion is used.)*

---

## DN-10 — The melt front is bounded by the cell, not the borehole

**What was wrong.** The geometric check on the melt-layer thickness compared it
against the borehole wall. But each tube owns only its share of the annulus, so
the fronts of neighbouring tubes meet at the **cell** radius,
$r_{\rm cell} = (D_{\rm well}/2)/\sqrt{2 n_t}$, long before either reaches the
wall.

**Worth.** The old check was a factor 2.9 too permissive in δ and could
therefore never fire. Under the enthalpy formulation the corrected limit cannot
be *violated* — the enthalpy cap holds A_melt ≤ A_avail, and the two conditions
turn out to be the same line — but it is being **saturated** over 99 % of the
well. That is not a violation; it is the solution sitting exactly at the limit
of validity of the annular conduction resistance, and it is now the leading open
item.

---

## DN-11 — The deviation cancels out of the round-trip efficiency

**The trap.** At CSS the field absorbs 1.0123 × the target charge and delivers
1.0123 × the target discharge — identically the same ratio, forced by
η_storage ≡ 1. Reporting the delivered surplus against the *target* electrical
input gives η_RTE = 0.4543 and the impression that the design has gained 1.2 %
in round-trip efficiency.

**It has not.** Since both the output and the required compressor input scale
with (1 + d),

    η_RTE  ∝  (1+d) E_out,target / ((1+d) E_in,target)

is independent of d exactly. Computed: 0.448825 with the realised input, which
is *identically* the no-deviation value.

**Consequence for reading any sweep.** The deviation is a statement about
**plant size**, not performance: this field is 1.23 % larger than a 1 MWe plant,
in both directions at once. A sweep that moves the deviation and leaves the
residual melt fraction alone is changing the size; one that moves the residual
melt fraction is changing the physics.

---

## DN-12 — `layer_map` uses a node convention *(open)*

**What is wrong.** `layer_map` assigns melting temperatures on the v0.1 **node**
grid — $n_s+1$ points — then drops the first and returns the rest as if they
were segment values. The effect is systematic: the first layer is always one
segment short and the last one segment long, whatever the resolution. At
$n_s=99$, $N_{\rm lay}=9$ the layers get 10, 11 × 7, 12 segments instead of 11
each.

**The fix** is to assign by segment midpoint, which gives exactly equal layers
whenever $N_{\rm lay}$ divides $n_s$:

```python
mid = (np.arange(n_segments) + 0.5) / n_segments
j   = np.minimum((mid * N_lay).astype(int), N_lay - 1)
return np.asarray(T_m_lay)[j]
```

**Worth, measured.** At the design point the deviation moves 1.228 % → 1.304 %,
about 0.08 percentage points; $N_{\rm wells}$ moves −0.003. Across
$n_s \in \{99, 100, 108, 117, 180\}$ the deviation under the current convention
scatters between +0.89 % and +1.23 % non-monotonically, which is the signature
of boundary alignment rather than of resolution.

**Not yet applied**, because it requires re-freezing the fixture and updating
both documents for a 0.08 pp change. `validate_case` warns when `n_segments` is
not a multiple of `N_lay`, which is the case where it is largest.

---

## DN-13 — The conduction path was on the wrong side of the front while freezing

**What was wrong.** The PCM-side resistance in `compute_U_i` is an annulus
growing **outward from the tube wall**,

    R' = ln(1 + delta/r_e) / (2 pi k),

so the model must say which material that annulus is made of. It used the
**melted** thickness, with `k_l` whenever any melt existed, in *both*
half-cycles:

```python
delta = delta_from_area(A_melt, ...)
k_m   = case.k_l if E[i] > 0.0 else case.k_s
```

That is right while melting. Melting begins at the tube wall and the front moves
outward, so the liquid shell genuinely lies between the tube and the remaining
solid, and heat genuinely has to cross it.

It is inverted while freezing. Solidification also begins at the tube wall — the
tube is the driven boundary in both directions — so the shell against the tube
is the **frozen** material, of thickness corresponding to `A_avail - A_melt` and
conductivity `k_s`. The liquid is displaced outward, beyond the front, and is
not in the conduction path at all.

**How bad.** The two are complementary, so the error is a time reversal rather
than a scale factor. Over a discharge:

| ε | δ used | δ wanted | R used | R wanted | ratio |
|---|---|---|---|---|---|
| 1.00 | 23.37 mm | 0.00 mm | 0.2639 | 0.0000 | ∞ |
| 0.50 | 14.32 | 14.32 | 0.1833 | 0.1375 | 1.33 |
| 0.08 | 3.48 | 22.01 | 0.0541 | 0.1897 | 0.28 |
| 0.00 | 0.00 | 23.37 | 0.0000 | 0.1979 | 0 |

(m·K/W per unit tube length.) The discharge runs top to bottom. **The modelled
resistance fell by a factor of five as the store froze; it should rise from
nothing to its maximum.** A segment that had frozen solid was given *zero*
PCM-side resistance, exactly where the physical resistance is largest.

**How it survived.** Every result that mattered for a long time was a *melting*
result — the first-cycle sizing, the self-levelling comparison, the stability
and order-of-accuracy studies — and the melting direction was correct
throughout. The discharge produced a plausible-looking `U_i` that rose through
the half-cycle, and that rise was rationalised in §14 of the verification
notebook as "the annulus refrozen, and solid PCM conducts better than liquid".
The conductivity ratio is only k_s/k_l = 1.33; the observed factor was **7.1**,
and nobody asked where the other 5.3 came from.

Asking the code one question would have found it: *what resistance do you give a
segment at ε = 0 during a discharge?*

**The fix.** `conduction_shell` selects the shell by the sign of the driving
temperature difference, which is known before the conductance and has the same
sign as q'. All four limits then come out right with no special cases: fully
solid and melting → δ = 0; fully melted and freezing → δ = 0; and the two
saturated cases give the full annulus with the appropriate conductivity.
`Case.front_geometry = 'melt_side'` reproduces the retired behaviour exactly,
which is how the change was separated from everything else.

**Worth, measured at the design point.**

| | melt_side | directional |
|---|---|---|
| U_i mid-well, t→0, charge / discharge | 822 / 117 | 823 / 829 |
| U_i mid-well, end, charge / discharge | 117 / 821 | 117 / 148 |
| deviation | +1.23 % | **+4.73 %** |
| residual ε at end of discharge | 0.0823 | **0.0094** |
| ε at end of charge | 1.0000 | 0.9706 |
| N_wells | 11.743 | 11.626 |
| η_RTE | 0.4320 | 0.4370 |
| cycles to CSS | 32 | 73 |

The first two rows are the real result. The two half-cycles are now near mirror
images: at t → 0 they agree to 0.5 % because the shell has zero thickness at the
start of *both*, and at the end the ratio is 1.265, close to k_s/k_l = 1.333 as
it should be. The spurious factor of seven is gone.

**A prediction I got wrong, recorded because the reasoning is instructive.** I
expected the fix to make the store perform *worse*, on the grounds that the
model was 3.5× too optimistic at the end of the discharge, which is where the
rate limit bites. The deviation went the other way, +1.23 % → +4.73 %. The error
was to weigh only the end of the half-cycle: at the *start* the retired
treatment was infinitely too pessimistic, and that is when the driving
difference is largest and most of the energy moves. The early over-penalty
dominates the integral.

**Consequence for the physical reading.** The store now very nearly refreezes
(residual ε = 0.0094 against 0.0823), so the discharge is much closer to
inventory-limited at the design point than it appeared. The rate-limit
*diagnostic* survives — residual ε still rises with well count, 0.0000 at N = 9
to 0.0826 at N = 16 — but the design point now sits near the bottom of that
curve rather than in the middle of it.

**Not fixed: the three-region geometry.** One thickness describes one front. At
cyclic steady state the charge begins with residual liquid in the *outer* part
of the cell, so melting produces liquid at the tube, solid in the middle and
liquid outside. Neither treatment can represent that. `conduction_shell` is the
better approximation, not a correct one, and that is what H5 actually assumes.

**Found by** a reader asking which material the model places next to the tube
during discharge — the right question, asked of the right function.

---

## DN-14 — A single-phase cell has no front, so it cannot have a side

**What was wrong.** DN-13 made the conduction shell follow the direction of the
phase change. Asked what to do for a cell that is *entirely* one phase — all
subcooled solid, or all superheated liquid — the same rule returned a resistance
that depended on the direction of the heat flow:

| state | direction | δ | R′ [m·K/W] |
|---|---|---|---|
| subcooled solid | heated | 0 | **0** |
| subcooled solid | cooled | δ_merge | **0.1979** |
| superheated liquid | cooled | 0 | **0** |
| superheated liquid | heated | δ_merge | **0.2639** |

Same body, same material, two answers. There is no front in a single-phase
cell, so "which side of it" is not a question that has an answer.

**What actually changes at a branch boundary is the meaning of the state
variable.** In two-phase, T_pcm = T_m is the temperature *of the front*, and the
resistance that belongs with it is tube-to-front. In single phase there is no
front; T_pcm is the *volume mean* of the cell, and the resistance that belongs
with a mean is mean-to-surface.

**The derivation.** Take the cell as an annulus r_e ≤ r ≤ r_o, adiabatic at r_o
— its neighbour is identical, so that boundary is a symmetry plane — storing
sensible heat at a uniform volumetric rate s = ρc dT/dt. Quasi-steady, the heat
crossing radius r is what the material beyond it stores:

    −k 2πr dT/dr = s π(r_o² − r²)

    T(r_e) − T(r) = (s/2k)[ r_o² ln(r/r_e) − (r² − r_e²)/2 ]

Averaging over the volume and dividing by Q′ = s π(r_o² − r_e²) gives
R′_bulk = S/(2πk) with, for β = r_o/r_e,

    S = [ β⁴ lnβ − β⁴/2 + β²/2 − (β²−1)²/4 ] / (β²−1)²

**S depends on geometry alone** — k cancels — so one shape factor serves both
phases and only the conductivity changes. For the default cell β = 2.1086 and
S = 0.34672, against lnβ = 0.74604 for the surface-to-surface annulus: smaller
by a factor 2.152, the cylindrical analogue of the familiar 1/3 for a slab with
uniform generation. Checked against numerical quadrature of the same integral to
six decimals.

Fed back as an equivalent thickness, ln(1 + δ_eq/r_e) = S gives
δ_eq = r_e(e^S − 1) = **8.74 mm**, so the existing fin-efficiency and
finned-area machinery handles the single-phase branches unchanged instead of
needing a second path. Like S, it is a property of the geometry alone.

**The two retired values bracketed the answer** from opposite sides: zero, and
2.152× too large — and which one you got depended on the direction of the heat
flow.

**Worth.**

| | v0.7 | v0.8 |
|---|---|---|
| deviation | +4.73 % | **+4.53 %** |
| residual ε at end of discharge | 0.0094 | 0.0116 |
| ε at end of charge | 0.9706 | 0.9715 |
| N_wells | 11.6255 | 11.6322 |
| η_RTE | 0.4370 | 0.4367 |
| U_i mid-well at t→0, charge | 823 | **260** |
| U_i mid-well at end, charge | 117 | 214 |

The integrated effect is small — 0.20 percentage points on the deviation —
because the two errors partly cancelled over a half-cycle: too little resistance
at one end, too much at the other. **The qualitative change is larger.** The
spurious infinite conductance at ε = 0 and ε = 1 is gone, and with it the claim
that the exchanger is at its best at t → 0. It is not:

| t [h] | ε | U_i | which resistance |
|---|---|---|---|
| 0.00 | 0.000 | 260.2 | bulk, subcooled solid |
| 1.16 | 0.039 | **479.4** | shell — thin melt, the best moment of the run |
| 3.41 | 0.245 | 214.8 | shell — thicker, degrading |
| 10.00 | 1.000 | 213.7 | bulk, superheated liquid |

The exchanger is best *shortly after melting starts*, when there is a front at
the tube and almost nothing between them. Before that, heat is warming a solid
body and must reach its mean; after it, the melt layer grows and gets in the
way. The two half-cycles are now near mirror images, 260 → 214 on charge and
214 → 260 on discharge, differing only through k_s/k_l.

**The jump at a branch boundary is real, and it is not a bug.** As ε → 1⁻ while
melting, R′ = lnβ/2πk_l; the instant the cell becomes superheated it drops to
S/2πk_l. Nothing physical moved — T_pcm stopped being the front temperature and
became the volume mean, and the resistance follows the reference. The same
happens at ε → 0, where the mechanism switches from *warming a body* to *moving
a front*. The driving difference is bounded by the capacity rate through
K = ṁc_p(1−e^−NTU)/Δz, so the jump never produces an unbounded flux.

**Still not resolved.** All of this is a correction bolted onto a body that
formally has no internal profile (H1, H5). R′_bulk is the best available within
one scalar per segment; the honest treatment of the sensible branches is radial
discretisation — an apparent-heat-capacity scheme on a 1-D radial grid — which
would also dissolve the three-region problem left open in DN-13.

**Raised by** a reader asking whether assuming zero resistance in the
single-phase limits was better than computing one. It was not.

**One thing this exposed.** The equivalence test in §16.5 classified
`eta_storage` as a closed-form quantity. It is identically 1 at CSS as a matter
of algebra, but the *computed* value depends on which cycle the drift test
stopped at, so it inherits the 1e-9 tolerance like any marched quantity, and it
reported 1.0000000012. An identity in the mathematics is still only converged-to
in the arithmetic. Moved to the iteration-terminated group.

---

## DN-15 — The ORC evaporator ran heat uphill

**What was wrong.** The ORC evaporating temperature was set from the hot end of
the water alone:

    T_3e = T_2d - DT_2D_3E

and nothing ever looked at the rest of the exchanger. That is safe only if the
two streams have similar heat-capacity curves. They do not. The water glides
55 K, 146.1 → 91.1 °C, while the cyclopentane takes 76 % of its heat boiling at
one constant temperature, 134.1 °C. Composite curves of that shape cross, and
these crossed badly: the pinch was **−29.8 K**, worst at 23 % of the duty where
boiling begins, and the cold end was inverted by 12.7 K as well. Only the top
fifth of the exchanger had a positive driving difference at all.

So eta_ORC = 0.2320 was not optimistic, it was **unreachable** — no isentropic
efficiency, no exchanger area and no amount of money buys it, because the heat
was flowing from cold to hot over most of the surface.

**Why the hot-end rule looked reasonable.** DT_2D_3E = 12 K is a perfectly
sensible approach *at the hot end*, and the hot end is where anyone would check.
The failure is a full 30 K away from the place the parameter names.

**The fix.** `orc_pinch` builds both composite curves properly — real enthalpy
on both sides, the latent plateau as a vertical segment rather than a straight
line in T, an even grid in Q that includes every break point — and returns the
minimum gap wherever it falls. `feasible_rankine` keeps the hot-end rule when it
already clears `DT_pinch_ORC`, and otherwise bisects T_3e downward until the
pinch is met. The new field is `Case.DT_pinch_ORC = 5.0`.

    T_3e      134.11 °C  ->  95.43 °C      (first pass; 95.89 at CSS)
    eta_ORC   0.2320     ->  0.1738
    N_wells   11.66      ->  15.66
    eta_RTE   0.4372     ->  0.3140

**What this costs.** Roughly a third of the round-trip efficiency, and a third
more wells for the same duty. That is not a regression; it is the previous
number being a fiction. The honest statement is that a single-pressure ORC is a
poor match for a 55 K sensible glide, and the model now says so instead of
hiding it.

**The design lesson, which is the interesting part.** The glide DT_3C_2C is not
a free parameter to be raised for storage density. A large glide buys PCM
inventory per well and simultaneously destroys the ORC, because a boiling fluid
cannot follow a gliding one. That trade-off was invisible while the pinch went
unchecked, and it is now the most important thing in a parametric study. A
zeotropic mixture, or a second pressure level, is the way to buy it back.

**Raised by** a student reading the cycle code, who noticed that the water
reached the evaporator colder than the fluid boiling inside it.

---

## DN-16 — The heat pump created energy in IHX-2

**What was wrong.** Inside `two_stage_htheatpump_2regs`, the second internal
heat exchanger was closed with

    T_7h = T_12h - epsilon_IHX_2 * x_4h * cpr_2 * (T_12h - T_13h)

and `x_4h` is the **vapour quality at the evaporator inlet**. It is not a flow
split, and a flow split is what belongs there. The stream passing through IHX-2
on the vapour side is the suction stream — the part that left the flash
separator as *liquid*, went through the evaporator, and is now on its way to the
low-stage compressor. That is `1 - x_8h`. (IHX-1, one screen earlier, correctly
uses `x_8h`, the flashed vapour.)

The two are numerically similar enough to look plausible and different enough to
break the first law:

    IHX-2 liquid gave      9.19 kJ/kg
    IHX-2 vapour received 21.76 kJ/kg
    condenser             271.43 kJ/kg
    work + evaporator     258.33 kJ/kg      <- 4.8 % of energy from nowhere

**The fix.** The whole effectiveness-times-cp-ratio loop is gone, replaced by
the enthalpy balance it was approximating:

    h_12h = h_3h  - x_8h       * (h_11h - h_10h)      IHX-1, vapour = x_8h
    h_7h  = h_12h - (1 - x_8h) * (h_1h  - h_13h)      IHX-2, vapour = 1 - x_8h
    x_8h  = (h_7h - h_9h) / (h_10h - h_9h)            separator, h_8h = h_7h

`x_8h` appears on both sides, so this is a fixed point; it is linear and
strongly contracting (gain 0.013), and converges to 1e-12 in a few steps. State
1h is fixed independently by the isentropic compression back from 2h, so
`epsilon_IHX_2` is now an **output** — the effectiveness the regenerator would
need — rather than an input. It is reported, and a value above 1 would mean the
component cannot exist.

    x_8h      0.4160  ->  0.3778
    hp_cop    2.9563  ->  2.8868
    residual  13.11 kJ/kg  ->  1.7e-12 kJ/kg

**Why it survived so long.** Nothing ever summed the cycle. The COP expression
at the end of the function is written in terms of `x_8h` and is correct as
written; it simply consumed an `x_8h` that had been produced by an inconsistent
network. An energy balance over the whole cycle is one line and would have
caught it on day one. It is now an acceptance test.

**Raised by** the same student, from the same reading.

**Postscript: the version guard has a blind spot.** `check_version` was added at
v0.8 precisely to stop a stale version string shipping, and it did not stop this
one. It passed the v0.9 build cleanly, because it checks only that the three
generators *agree with each other* — and they agreed, on `0.8`. Agreement is not
currency. The notebook went out with correct v0.9 physics under a `0.8` stamp,
which is exactly the confusion the guard exists to prevent, arriving by the one
route it does not watch. Bumping `VERSION` is still a manual act, and nothing
ties it to "the physics changed". A guard that compares the build against the
last *released* stamp, and refuses when the model sources have changed but the
version has not, would close it. Not yet written.

---

## DN-17 — The heat-pump evaporator had no approach at all

**What was wrong.** Two independent fields set the two ends of the same
exchanger:

    T_4a  = T_source_C - DT_3A_4A      # source outlet
    T_13h = T_source_C - DT_3A_13H     # evaporating temperature

so the approach between them was their *difference*, `DT_3A_13H - DT_3A_4A`,
and nothing checked its sign. The defaults, `10.0` and `10.0`, made it exactly
**zero**: the source left the evaporator at precisely the evaporating
temperature, which needs infinite area. Swap them — `DT_3A_4A = 12`,
`DT_3A_13H = 10` — and the approach is **−2 K**, a second-law violation.

The model reported the same COP in all three cases.

**Why nothing could have caught it.** `T_4a` was written into the state dict at
one line and **never read again**. The source stream's energy balance was never
closed, so there was no exchanger there at all: the evaporating temperature was
declared, the source outlet was declared, and the two were never required to be
consistent. A crossing could not produce a symptom because nothing downstream
depended on the quantity that would have crossed. The v0.8 guard added for this
was worse than useless — it warned when the two were *equal*, the benign case,
and said nothing when they were reversed.

That is the more general lesson here, and it is worth more than the fix: **a
parameter that influences no output cannot be validated by any amount of
checking the outputs.** `DT_3A_4A` was inert, so every consistency test the
project has — energy closure, the CSS identity, the fixture, the equivalence
tests — passed with it set to a value that violates the second law.

**The fix.** The evaporating temperature is now *derived* from the source
outlet and a stated minimum approach:

    T_13h = T_source_C - DT_3A_4A - DT_pinch_HPE      # DT_pinch_HPE = 5 K

`DT_3A_13H` is retired and `validate_case` raises if a Case still names it,
rather than ignoring it silently. The approach is now positive by construction
and `DT_3A_4A` is no longer inert: it moves the evaporating temperature, and
therefore the COP.

**Why an end approach is sufficient here, when it was not for the ORC.** The
cold stream in this evaporator is isothermal with **no preheat section** —
state 4h enters two-phase and 13h leaves saturated, so the composite is flat
with no kink. The hot stream is monotonic. The gap is therefore monotonic in Q
and its minimum is necessarily at the cold end, which is what `DT_pinch_HPE`
names. No composite-curve search is needed. The ORC evaporator needs one
precisely because its liquid preheat puts a kink before the plateau, which is
what moves the pinch into the interior (DN-15).

**What it costs.**

    T_evap    50.0 -> 45.0 C
    COP     2.8868 -> 2.7594      -4.4 %
    eta_RTE 0.3140 -> 0.3003      -4.4 %

Everything on the discharge side is untouched to the last digit — `eta_T`, net
electric per well, `N_wells`, the CSS deviation, thermal energy per well,
storage density — because the heat pump appears only in the charging branch of
the budget. It buys electricity in, not electricity out.

**Still open.** The source stream is specified but not sized: nothing tracks the
source mass flow, so the cost of cooling it by 10 K rather than 5 K is still
invisible. Closing that would make `DT_3A_4A` a genuine trade rather than a
free choice.

**Raised by** the observation that a zero approach was holding only by accident.
It was not even an accident: it was the absence of a mechanism.

---

## DN-18 — Four temperature differences that were doing more than one job

Not defects, this time: an audit. Each of these was a single field standing in
for two or three independent physical statements. All four changes are exact
no-ops at their defaults, verified to \(10^{-6}\) on every reported quantity.

### 1. `DT_E_sink` referenced the wrong end

`T_1e = T_sink_C + DT_E_sink` measures the ORC condenser approach to the sink
**inlet**. That is the same wrong-end pattern that made `DT_3A_13H` unsafe
(DN-17); it was merely masked, because the sink is modelled as an infinite
reservoir and inlet = outlet. Now `T_1e = T_sink_C + DT_sink_glide +
DT_E_sink`, with `DT_sink_glide = 0`.

**And the ORC condenser is not structurally safe the way the HTHP condenser
is.** State 4e is superheated by 22.5 K, so the hot composite has a
desuperheating kink at the HOT end — 7 % of the duty. Give the sink any glide
and the pinch appears at that corner, **interior, at 93.4 % of the duty**, and
crosses at about 5.4 K:

    sink glide  0 K -> +5.000 K      4 K -> +1.263 K
                5 K -> +0.328 K      6 K -> -0.606 K   crosses

The HTHP condenser has no such kink: state 2h is *saturated vapour*, because
the model forces isentropic compression to land exactly on the dew line, so its
only kink is the 2 K subcooling at the COLD end where the gap is widest. The
rule is not "isothermal hot stream is safe" — it is **where the kink sits**.

### 2. `DT_2D_3E` is inert but not obsolete

`T_3e` is now written as what it always was:

    T_3e = min( T_2d - DT_2D_3E ,  T_3e_pinch )

Neither constraint subsumes the other, and `rank['T_3e_binding']` names the
active one. At the default `DT_2D_3E = 12 K` the pinch binds across the whole
useful range of glides (15–70 K), so the hot-end rule never becomes active —
but it would if raised, and it remains the bisection's upper bracket. The
change is presentational: the same arithmetic, stated as a constraint rather
than as "start from the old rule and correct it".

### 3. The cascade span was tied to the glide, and the bound runs the other way

`DT_cascade` now grades the cascade independently, defaulting to `DT_3C_2C`.

That default is not arbitrary. With span = glide·(N−1)/N, the approach between
the water and the layer it is working on is **exactly `DT_4C_M` at the leading
face of every layer**, decaying to `DT_4C_M - span/N` at each trailing face —
10.000 K and 3.889 K at the design point. That sawtooth *is* the cascade, and
matching the span to the glide is the unique choice that makes it uniform.

The important correction is to the direction of the bound. The span must be
**larger** than glide minus the approach, not smaller:

    charging : span > DT_3C_2C  - DT_4C_M = 45.0 K
    discharge: span > glide_dc  - DT_M_1D = 45.0 K

At the default, span = 48.889 K — only **3.889 K of margin**. Below the bound
the driving difference inverts at the far end of the well, and `validate_case`
now raises. Sweeping the cascade narrower is the dangerous direction, and it is
the one that looks attractive because it raises `T_m,bottom` and hence the ORC
evaporating temperature:

    span 26.7 K -> eta_RTE 0.3508, but residual melt 0.196 and -32.6 % delivery
    span 48.9 K -> eta_RTE 0.3003, residual melt 0.012,        +4.5 % delivery

### 4. `DT_3D_2D`: the two half-cycles need not glide alike

Two of the three couplings were never assumptions. `T_3c = T_4c`, so the HTHP
condenser glide **is** the borehole charging glide — same water, same flow,
closed loop, energy conservation. Likewise the ORC evaporator glide is the
borehole discharging glide. Neither can be relaxed.

What *was* assumed is that the charging and discharging glides are equal: both
mass flows were divided by `DT_3C_2C`. `DT_3D_2D` now separates them, and there
is an interior optimum:

    DT_3D_2D 40 K -> deviation +6.4 %, eta_RTE 0.2764, store fully refreezes
    DT_3D_2D 55 K -> deviation +4.5 %, eta_RTE 0.3003
    DT_3D_2D 65 K -> deviation -8.5 %, eta_RTE 0.3041, residual melt 0.227

A wider discharge glide raises `T_2d`, lets the ORC boil hotter and lifts
`eta_RTE` — until the store cannot keep up and the plant misses its target.

**Why the three were originally tied together.** The cascade is *shared
hardware*: graded once, and both half-cycles must live with it. Matching both
water glides to the span is the only way to hold a uniform approach in both
directions. That is a defensible design choice — but it is a choice, and the
tables above show it is not obviously the best one.

**Still open.** Two interior optima are now visible (cascade span, discharge
glide) and neither has been searched, nor has the glide itself (DN-15). That is
a three-parameter optimisation the model can now express and has never been
asked.
