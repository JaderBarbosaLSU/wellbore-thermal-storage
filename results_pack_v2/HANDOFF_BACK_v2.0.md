# Results pack v2.0 — back to the BOREAS paper task

*From the model/notebook task, 2026-10-29. Paste this whole file as the
prompt and attach `BOREAS_results_pack_v2.zip`.*

The field ladder is implemented natively and run. **Every number in your
draft that came from the prototype has moved**, because three things
changed and all three were necessary. Nothing here is tuned to match
anything; where this disagrees with the prototype, this is the one to use
and the reason is stated.

---

## 0. The rule that matters most

**Do not retype a single number from this file into the paper.**

`tab_ladder.tex` is machine-generated from `ladder_results.csv`. `\input`
it. Every instance of drift we have found between these two tasks — the
pinned-span spread quoted as 0.17 % when it is 0.19 %, the disc-against-rows
comparison that mixed two packings, the $T$–$s$ figure drawn at a different
case from its own state-point table — was a number restated by hand on one
side after being recomputed on the other. The CSV is the source of truth;
the table and figures are views of it; prose is yours.

If you need a number that is not in the CSV, ask for it rather than reading
it off a figure.

---

## 1. What changed, and why every number moved

### 1.1 The signed-discharge bug (§1.4 of your handoff) — confirmed, fixed

`simulate_css` stored `abs(Q_cum_J)` for the discharge. Confirmed and
reproduced: an isolated 300 m well integrates to $+8.7\times10^{7}$ J/m
with $\overline{T}_{2d} = 85.5$ °C against $T_{3d} = 91.1$ °C — the store
*absorbs*, and `abs()` reported it as delivered. Model is now v1.4: the
value is signed, `delivers` is reported, and $e'$ and $\eta_\mathrm{RTE}$
are zeroed when it is False while the loss is reported as computed.

**Nothing in the paper as it stands was affected.** Every configuration in
your current draft is shielded or long enough to deliver, so the integral
was already negative and `abs()` was inert. The v1.3 reference is
bit-identical: 7.93 % / 0.2983 / 0.3176 / 19.14 / 142.62.

### 1.2 The idealised rungs were at the wrong spacing

The WellSTAR slots give a **median nearest-neighbour of 1.65 m**, 10 %
tighter than the nominal 6 ft that `disc` and `double_row` used. So the
prototype's R2 → R3 step mixed a layout change with a spacing change.

Rebuilt at 1.65 m, the layout step is **+4.2 points, not +3.0** — the real
layout penalty was **understated by 40 %**. At the nominal spacing my port
reproduces the prototype exactly (R1 4.4 / R2 14.4 / R3 17.4), which is
what validates it.

### 1.3 The casing is worth more than the ground model

The THUMS upper-zone string (SPE-2562) needs no new code — it is three
existing `Case` fields. Its 201.2 mm ID carries **46 % more PCM per tube**
than the legacy 177.8 mm, and that is worth **+36 % on $e'$**.

**The pending WellSTAR unknowns are not worth waiting for.** At R3,
L = 1250 m: top of cement anywhere between surface and shoe moves the loss
by 0.09 points; `k_mud` across its whole plausible range 0.6–1.5 moves it
by 0.17; 32# against 36# costs 2.2 % of $e'$. They are all sub-2 %, because
the near-term resistance they move is itself only ~7 % of the total. The
one genuinely blocking unknown is the **shoe depth**, which is a constraint
on $L$ rather than a resistance.

---

## 2. The result

`tab_ladder.tex`, at L = 1250 m MD, N = 363, superposition, THUMS casing.
Prototype values in brackets for orientation only — **do not print them**:

| step | | loss % | $e'$ |
|---|---|---|---|
| R0v | isolated, vertical | 64.1 (85.9) | 0.169 |
| R0i | isolated, S-curve | 68.9 (92.4) | 0.146 |
| R1 | disc, vertical | 3.1 (4.4) | 0.425 |
| R2 | double row, vertical | 10.2 (14.4) | 0.403 |
| R3 | Grissom slots, vertical | 13.4 (17.4) | 0.391 |
| R4 | Grissom slots, S-curve | 52.4 (70.7) | 0.215 |

**No run is non-delivering.** The prototype's 12 were an artefact of the
legacy casing, not a physical regime: with the THUMS string an isolated
100 m well goes from absorbing (81.9 °C, loss 114.7 %) to delivering
(99.2 °C, loss 85.6 %). Both fixes were needed — the signed fix zeroes the
legacy case correctly, the casing removes the regime. **Any sentence about
short isolated wells failing to deliver must come out.**

**The cutoff criterion survives**: marginal RTE is 0.369 over the first
interval and first drops below half of it at **300–400 m MD**, exactly as
the prototype found. Vertical wells still have no cutoff.

**Sensitivities at L = 1250 m:** KOP all at 300 ft moves the loss
52.4 → 35.1 %, the largest single lever, so "the cutoff is set by the
kickoffs" holds. $\theta_\mathrm{max}$ = 45/70° gives 50.0/54.9 %. And
$t_\mathrm{op}$ = 1/30 yr gives **71.7/39.2 %** — the widest spread of
anything tested. Field age matters more than trajectory, and it is new
since the prototype.

---

## 3. A geometry correction you should propagate

I earlier told this project that Grissom was "a blob, not two rows", from a
convex-hull aspect test. **That was wrong.** The hull of a ring contains
its empty middle, so "inside the hull" does not mean "off the rim".

The radial distribution settles it: **no slot lies inside 0.3 R_max and
72 % lie beyond 0.7 R_max**, against 0.09 and 0.51 for a filled disc.
Grissom is a **thick closed loop — a double row bent into a ring**, which
is what Jader described originally. Panel (a) of `grissom_geometry.pdf`
shows it directly.

Consequences: the straight double row (R2) is the right idealised
comparator, and any text calling the layout irregular or cluster-like
should say *ring* instead.

---

## 4. What is in the zip

- `ladder_results.csv` — 111 runs, one row per run, the §6 columns.
- `tab_ladder.tex` — `\input` this.
- `profiles.npz` — $R_{g,\mathrm{eff}}$(MD) and TVD(MD).
- `grissom_geometry.pdf` — slots by KOP group, trajectories, spacing vs MD.
- `grissom_Rg.pdf` — the shielding profiles.
- `grissom_ladder.pdf` — the bars.
- `grissom_sweep.pdf` — $e'$, RTE and marginal RTE against L.

Figures are serif with TrueType embedded (`pdf.fonttype=42`), authored for
a single-column layout.

---

## 5. Status against your verification table

| # | check | result |
|---|---|---|
| V1 | perimeter, v1.3 reference | **pass**, bit-identical |
| V5 | $k_r \to 0$ | **pass**, loss 0, $\eta_\mathrm{storage}$ 1.000000000 |
| V6 | isolated, L = 300 m | **pass**, delivers False, $e'$ = 0 |
| V8 | R3/R4, THUMS casing | **pass**, 13.3 % and $e'$ = 0.391 against your estimate of 13.3 and 0.391 |
| V2–V4 | flux / superposition $R_{g,\mathrm{eff}}$ point checks | not yet run as standalone assertions |
| V7 | ladder against `ladder_results.csv` | superseded: the prototype CSV is at the legacy casing and nominal spacing, so it is no longer the right target. R1–R3 reproduce it exactly when both are set back. |

The V2–V4 and V7 gaps, the README with the full V-table, and extending
`profiles.npz` to every step and sensitivity are the next items here.

---

## 6. Open, and worth a line in the limitations

- Shoe depth per well, which sets the reference length. Everything else
  about the casing is settled.
- Surveyed trajectories; the S-curve is a build-and-hold reconstruction
  from SPE-2562, SPE-11360 and SPE-15101, not survey data.
- Ground-surface and well-end effects, unrepresented as in $R_g$.
- The 80 m superposition cutoff and the 50-station sampling: the profile is
  insensitive to the latter for vertical wells (identical to 4 dp at 15
  stations) but not for the S-curve, where 25 stations move it by up to
  0.35 m K/W locally.
