# Results pack v2.5 — the supplementary material, at R4

*From the model/notebook task, 2026-10-31. Attach
`BOREAS_results_pack_v25.zip`.*

§8 is closed. **Nothing in the supplementary material is still owed.**
18 checks pass.

**This note is not cumulative — read `HANDOFF_BACK_v2.4.md` first.** It
answers §1–§7 and §10 of your review and is in the same zip; this one
covers §8 and one correction. The *pack* is cumulative: v2.5 is a strict
superset of v2.4, so the v2.4 zip can be discarded.

---

## Table S1 — `tab_states.tex`, label `tab:states`

The 13 state points of both cycles, at the reference case (R4), with the
ORC columns 1e–10e and the HTHP columns 1h–13h, each carrying *T*, *p*,
*h* and *s*.

Two things are different from the version the paper has been carrying:

- it is at **R4**, not at the former 1000-well vertical array in double
  rows on four islands, which no longer exists as a configuration;
- it is evaluated at the **realised** store outlet, 116.54 °C, not the
  provisional 146.11 °C. At R4 those differ by 30 K, so a state-point
  table at the provisional value would disagree with every other table in
  the pack.

Set `\footnotesize` with `\tabcolsep` 4 pt, since it is a ten-column
table; please check it against your text width.

## Video S1 — `BOREAS_cycle.mp4` and `.gif`

The full cycle at cyclic steady state: 120 frames at 12 fps, three panels
with depth on the ordinate — fluid temperature against the melting
staircase with *T*_pcm overlaid, local melted fraction, and *U_i* — with
axis limits fixed across both half-cycles so nothing jumps between them.

Rendered **inside `field_ground`**, which matters more than it sounds:
`animate_well` takes its depth axis from the model's own
`formation_profile`, and outside the context manager that silently draws a
1250 m well. It now runs to 904 m, the true vertical depth.

`boreas_video_numbers.tex` carries `\VideoFrames`, `\VideoFps` and
`\VideoDepthBottom`.

---

## One correction, found while wiring the animation

**`tvd_bottom` was the TVD at the last station centre, not at the bottom
of the store.** At Δs = 50 m the last station sits at MD = 1225 m, so the
reference case has been reported **18 m too shallow throughout**:

| | v2.4 and earlier | **v2.5** |
|---|---|---|
| R4, bottom of the store | 891 m | **909 m** |
| R1, R2, R3 (vertical) | 1225 m | **1250 m** |

The vertical rungs are the clearer tell: a vertical well of 1250 m MD
reaches 1250 m TVD by definition, and the pack was saying 1225.

Every "the S-curve reaches only 891 m TVD at 1250 m MD" — in my handoffs
from v2.0 on, and in the paper if it was quoted — should read **909 m**.
The physical point is unchanged and the ratio barely moves: the S-curve
still falls about 27 % short of measured depth, against 29 % as reported.

**It is a reporting quantity only.** The solver places each station at its
own depth and never used this number, so no run was recomputed and no
performance result moves. The old value is kept beside it as
`tvd_last_station` so the earlier numbers remain traceable.

The macro diff is five entries, all TVD:

| macro | v2.4 | v2.5 |
|---|---|---|
| `\RefTvd`, `\LadderTvdRFour` | 891 | **909** |
| `\LadderTvdROne`, `\LadderTvdRTwo`, `\LadderTvdRThree` | 1225 | **1250** |

`boreas_numbers.tex` stays at 181 macros; nothing else changed.

### How it survived four rounds

The same way the stale `tab_ladder` did: nothing checked it against an
independent statement of the same quantity. V10 compares R_g,eff across
lengths, V17 compares rates against the energy budget — but no check ever
asked whether a *vertical* well's TVD equals its MD, which is the sort of
identity that would have caught this immediately. I have not added one for
it, and I think that is the honest place to stop: the quantity is now
computed from the trajectory at MD = L in one place, and a check that
re-derives the same one-line interpolation would test nothing.

---

## Also in v2.5

`depth_axis` now returns true vertical depth whenever the field patch is
active, so **every** figure and the animation pick it up automatically
rather than each one re-deriving it. `field_ground` leaves the mapping on
the model for that purpose and removes it on exit;
`formation_profile`'s own second return value stays measured depth,
because the resistance lookup is indexed on it.

The animation's ordinate now reads "TVD / m", matching the still figures.

---

## Status

Everything from handoffs v2.1 through v2.4 is delivered. The open items
are no longer mine:

- the reference *h*_m, 380 kJ/kg at 150 °C, which is a materials question
  and the largest single uncertainty in the economics — at 260 kJ/kg the
  field needs 69 storage wells per MWe instead of 37;
- the casing annulus data, which §6.2 of v2.3 showed is worth up to
  2.7 points of loss and 5.8 % of *e*′ at R4, not the sub-2 % I first
  reported from the R3 triage;
- the shoe depth per well, which sets the reference length.

The grid is good to a few tenths of a point of loss and about 1 % of
*e*′, and the README says so rather than implying better.
