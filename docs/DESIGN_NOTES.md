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
| DN-1 | The state variable must be enthalpy, not melted area | 9.3 % of cycle energy |
| DN-2 | The closing conductance is not $U_iP_i$ | 4–33 % on the segment conductance |
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

**Worth.** At CSS the store finishes every charge fully molten in all 100
segments with a peak superheat of 18.5 K, and 9.3 % of the energy passing
through the store each cycle traverses the sensible branches (7.1 % superheat,
2.3 % subcooling). That is the energy an area-based state must fabricate or
throw away. A secondary and larger effect is on the *distribution*: sensible
capacity is self-levelling, and the spread in local melt fraction at equal
exchanger size falls from 0.405 to 0.000.

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

**Worth.** Over the charge, NTU runs 0.083 → 0.602, so the ratio
(1 − e^−NTU)/NTU runs 0.96 → 0.75: using $U_iP_i$ overstates the segment
conductance by 4 % to 33 %, with the largest error early, when the melt layer is
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
