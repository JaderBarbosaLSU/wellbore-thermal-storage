# BOREAS v1.3 results pack for `P2H2P_enthalpy_paper`

**Revision 6 (2026-10-27).** *Consistency fix.* The T–s diagram was
generated from `Case()` — the **assumed** cascade, T_2d = 146.11 °C —
while the state-point, UA and T_3e tables come from the array reference,
where the **realised** T_2d is 142.62 °C. Figure and table therefore
disagreed (T_3e 95.09 vs 94.46 °C, eta_ORC 0.1481 vs 0.1471, p_evap 3.704
vs 3.647 bar). The figure is now built from the same realised cycle as the
tables. New `tex/tab_reference.tex` states the reference case explicitly.

**Revision 5 (2026-10-27).** Figures re-authored for a **single-column**
layout, so `width=\linewidth` reproduces them near 1:1. The four
three-panel figures (cascade map, geometry maps, field maps, and the T_m
study) are now **independent panels in `subfigure` environments** — they
were unreadable squeezed into one single-column width. Heat-map changes:
contour trendlines removed; the reference and best cells are called out by
**styling the value** (bold, and bold-larger respectively) plus a red cell
outline, instead of a star or ring drawn on top of the number it points at;
both are now explained in the captions. T–s diagram: every annotation
carries a white mask so no label sits on a curve.

**Requires `\usepackage{subcaption}`.** Panel PDFs live in
`figures/panels/`.

**Revision 4 (2026-10-26).** *Coordinate correction.* The profile figures in
revision 3 plotted the solver's `z` and labelled it depth. It is not: `z = 0`
is the **bottom** of the well on charge, so `depth = L_well - z`, and the
figures had the well upside down. They now use `formation_profile`'s own
depth array. The consequence for the text is that the hottest layer
(150 °C) is at the **bottom**, charge enters the fin tubes at the bottom and
rises, and discharge enters at the top and descends — the opposite of what
revision 3 said. PCM temperature is now overlaid on the fluid profile and in
the animation.

**Revision 3 (2026-10-26).** Profile figures redrawn with depth on the ordinate, increasing downward, and an animation of the full cycle added (`figures/BOREAS_cycle.mp4`, `.gif`).

**Revision 2 (2026-10-26).** Review comments applied: the composite-curve figure and the hot-end-rule discussion are removed, the state-point table is now supplementary, and the hot-end column is gone from `tab:T3e`. A new subsection on thermal behaviour along the well has been added.

Everything here was generated from `boreas.py` v1.3 on 2026-10-25 at the
reference case of `tab:conditions`, `n_s = 99`, `t_op = 10 yr`. Paste this
prompt and the attached files into the enthalpy-paper task.

---

## What this closes

| `\todo{}` in `results.tex` | closed by |
|---|---|
| §4.1 *T–s diagram of both cycles … table of state points, pressures, y, x_8h, ε_R1, ε_R2, η_ORC, COP, duties and flow rates* | `tex/res_cycles.tex` → `fig:ts`, `tab:states`, `tab:cycles` |
| §4.1 *ORC evaporator composite curves showing the interior pinch … (v1.0 values; regenerate at v1.3)* | `tex/res_cycles.tex` → `fig:composite`, `fig:composite-stores`, `tab:T3e` |
| §4.1 *UA and pinch of the four plant exchangers* | `tex/res_cycles.tex` → `tab:UA` |
| §4.4 *(partially)* | `tex/res_param.tex` → five new paragraphs and figures |
| §4.3 *Profiles along the well at CSS …* | `tex/res_thermal.tex` → `fig:prof-T`, `fig:prof-melt`, `fig:prof-Ui`, `fig:prof-int` |
| *(new)* one sentence for the model section | `tex/model_hotend_sentence.tex` |
| *(supplementary)* cycle animation | `figures/BOREAS_cycle.mp4`, `figures/BOREAS_cycle.gif` |

**Not** closed, and deliberately left alone: the §4.4 items on storage length
against `z_KO`, the `t_op` sweep, `k_r`, the insulation sweep, and
`N_geo`; the §4.3 profile and convergence plots; the `\chk` on moving the
reference to 4000 ft.

---

## How to install

1. Copy `figures/*.pdf` into the paper's `figures/`.
2. Copy `tex/*.tex` into a `tex/` subdirectory (or flatten and fix the
   `\input` paths — `res_cycles.tex` does `\input{tex/tab_states}` etc.).
3. Replace the whole of `\subsection{Power-to-heat-to-power cycle analysis}`
   with `tex/res_cycles.tex`.
3b. Insert `tex/res_thermal.tex` after *Loss to the formation and the role of
   the well field*. It closes the `\todo{}` on profiles along the well.
3c. Paste the single sentence in `tex/model_hotend_sentence.tex` into the
   model section, immediately after the equation defining `T_3e`. Without it
   the `min()` of two constraints has an unexplained second argument.
3d. Put `tex/tab_states.tex` in the supplementary material (its caption now
   says so).
4. Insert `tex/res_param.tex` into `\subsection{Parametric assessment}`,
   before the existing `\todo{}`s — it adds to that subsection rather than
   replacing it.

### Macros and packages assumed to exist

`\etaORC`, `\COP`, `\Nlay`, `\nft`, `\md`, `\rev{}`, `\chk{}`, `\todo{}`,
plus `siunitx`, `graphicx` and **`subcaption`**. No new macros are defined. If `\Nlay` or `\nft` are
not defined in the preamble, three uses in `res_param.tex` will fail.

### New labels

`tab:states`, `tab:cycles`, `tab:UA`, `tab:T3e`, `fig:ts`, `fig:composite`,
`fig:composite-stores`, `fig:map-cascade`, `fig:map-geometry`,
`fig:map-field`, `fig:nlay`, `fig:Tm`.

---

## Three things to check before accepting the text

**0. Orientation of the profile figures, which I verified rather than assumed.**
Both half-cycles are recorded in the depth frame (index 0 = top). On charge
the water enters at the top at 160 °C; on discharge it enters at the
*bottom* at 91.1 °C and rises counter-current. The test that settles it is
the PCM temperature against its own melting staircase — the reversed reading
would put the PCM 27 K *above* its melting point at the end of discharge.
Reading the fluid boundary alone is not enough, because the inlet is at
opposite ends in the two half-cycles. Worth a sanity check against your
physical picture of the downcomer.

**1. The v1.0 `T_3e` values in the `\todo` have moved, one of them a lot.**

| store | v1.0 in the paper | v1.3 |
|---|---|---|
| adiabatic | 95.0 °C | **95.03 °C** |
| 1000-well fan | 93.9 °C | **89.90 °C** |
| isolated well | 88.4 °C | **88.04 °C** |
| vertical array, rows | — | **94.46 °C** |

The fan moved because at v1.3 it is computed at the rows packing (57.1 %
loss) rather than the disc (12.8 %). If the sentence in the paper is meant to
quote the disc fan, the number is 93.9 °C and should say *disc*.

**2. A `\chk` worth resolving: is `DT_3C_2C = 55 K` a free choice?**
Within the model a narrower glide improves *both* η_RTE and e′ at every
storage temperature — there is no cell in the 36-point map where 55 K beats
25 K. If the 55 K is set by the PCM selection or the layer grading, that
should be said, because otherwise the reference case is leaving ~12 % of
η_RTE on the table.

**3. The field maps use rows packing throughout, so their 45° column is not
comparable with `tab:formation`'s disc-packing fan rows.** The caption says
so. Worth a glance to make sure that does not read as a contradiction.

---

## Things the paper already has right (no action)

I re-derived two results that turned out to be already in your `results.tex`,
so these are confirmations, not corrections:

- **The rows/disc comparison.** `tab:formation` already carries fan-1000 at
  rows (57.1 %, RTE 0.131) beside the disc (12.8 %, 0.281), and the paragraph
  *"The result is governed by the arrangement of the field"* already draws
  the right conclusion. (The *model report* `boreas_tex/sec/10_results.tex`
  still has the old apples-to-oranges comparison and does need fixing — that
  is a separate document.)
- **The convergence flag.** Your `\chk` in `tab:ref-results` already records
  the 80-cycle limit and the 1e-9 absolute tolerance. For the record, the
  quantitative framing: drift falls from 3.6e5 to 1.22e-9 J/m and plateaus,
  a reduction of 3.2e-15 *relative*. The solve is converged to round-off; the
  tolerance is absolute on a quantity of order 1e5 and so is unreachable by
  construction. Making it relative would fix the flag and cut ~60 of the 80
  cycles. Not changed, because it moves the physics hash.

---

## Regenerating

```
python3 papercalc.py adiabatic isolated fan1000_rows array_rows   # ~2 min
python3 mktex.py          # rebuilds the four tables from papercalc.json
python3 mkmaps.py         # rebuilds the heat maps from the sweep cache
python3 mkfig_store.py    # rebuilds fig:composite-stores
```

`data/boreas_sweep_cache.json` holds all 180 cyclic-steady-state evaluations
keyed by case parameters, so refining a map recomputes only the new cells.
The tables are machine-generated — regenerate rather than hand-edit, or the
paper will drift from the model.

### A trap worth knowing about

`boreas_sweep.py` runs `validate_case` on every point. Without that gate a
`T_m` sweep past 180 °C silently computes on **superheated steam** — water
saturates at 179.9 °C at the shipped 10 bar and CoolProp returns a vapour
enthalpy about ten times the liquid one rather than failing — and the
resulting pinch is meaningless. That artefact produced a false η_ORC optimum
at 185 °C before the gate was added. The check exists in the model; it simply
is not called by `cycle_state_points`.


---

## The animation

`figures/BOREAS_cycle.mp4` (120 frames, 12 fps, ~10 s) and
`figures/BOREAS_cycle.gif` (90 frames, lower resolution, for slides and
anything that will not play video). Three panels with depth on the ordinate
— fluid temperature against the melting staircase, local melted fraction,
overall coefficient — running charge then discharge, with elapsed time and
the mean melted fraction in the header. Axis limits are fixed across the
whole cycle so the two half-cycles can be compared directly and nothing
jumps at the changeover.

One thing to know about the clock: the march records on a roughly
*logarithmic* time grid, which is right for static profiles but would make
constant-rate playback appear to accelerate. The frames are therefore
resampled onto a uniform time grid, so screen time is proportional to
physical time. The consequence is that the sub-minute startup transient goes
past almost instantly — that is the honest depiction, not a dropped frame.

Regenerate with `python3 src/mkanim.py` (needs `rec.pkl`, produced by
`simulate_css_corrected(..., record=True)`; see `src/mkfig_prof.py`).
