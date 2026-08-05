# Handover brief v2 — THUMS journal-paper work

**Supersedes v1.** Paste this into the first message of a Cowork task running **on your
computer**, together with the GitHub token and `THUMS_code_audit_2026-08-05.md`.

---

## Context for the new session

I am upgrading an IHTC conference paper to a journal paper. The system is
power-to-heat-to-power storage: a high-temperature heat pump charges PCM-filled repurposed
oil wells, an ORC discharges them.

Repository: `https://github.com/JaderBarbosaLSU/wellbore-thermal-storage`
All four relevant notebooks are in `thums_core/`:

| File | Role | Status |
|---|---|---|
| `THUMS_Multilayer_BB_param_Jan_12b_converg_improved.ipynb` | **best solver core** | most recently developed (4 Aug 2026) |
| `THUMS_Multilayer_DoE_Feb_22_IHTC.ipynb` | produced the conference-paper numbers | contains a regression — see below |
| `THUMS_Multilayer_DoE_Dec_5.ipynb` | the only notebook with the DoE machinery | older solver core |
| `Pareto_for_THUMS_DoE.ipynb` | post-processing | no actual Pareto front in it |

A prior session produced a full code audit — **read `THUMS_code_audit_2026-08-05.md`
before doing anything else** and do not re-derive it. That audit was written against
`Feb_22`. The section below records what changes now that `Jan_12b` is in the repo.

**Workflow rules** (from `__Colab_GitHub_Token_Workflow.docx`): GitHub is the single source
of truth. You edit and commit; I open the Colab-from-GitHub link and run. I do not edit in
Colab. Version-stamp every change (top text cell, runtime banner, figure captions) and
state the version number with every update. Prefer guards that fail loudly over guards that
substitute an approximation.

**Also note:** GitHub's web editor silently truncated a 2.9 MB notebook to 2 bytes during a
rename earlier in this project. Never use the web editor on these files. Commit through git,
or use Add file → Upload files.

---

## Which notebook is the real solver — resolved

`Jan_12b` is **not** an experiment. It is the most advanced core and it fixes two of the
audit's findings outright. `Feb_22`, which produced the conference-paper numbers, is a
*regression* of it.

### `Jan_12b` fixes, relative to `Feb_22`

- **Audit §2.1 — `__globals__` back-doors: gone.** `D_E_out_HP` and `D_E_in_ORC` are proper
  parameters of `evaluate_Q_ratio_ch` / `evaluate_Q_ratio_dc`. No namespace reach-through
  anywhere in the file.
- **Audit §2.3 — melt-layer conductance: consistent.** `compute_U_i` uses the cylindrical
  form `h_e = k_m/(r_e·ln(1+δ/r_e))`, matching `compute_T_re`.
- New `compute_delta2_fast` and `find_*_fast` solvers, CoolProp state caching
  (`make_cp_state` / `get_props_state` / `compute_h_i_from_state`), a mass-effectiveness
  cache, and guard helpers `is_finite` / `is_pos_finite` / `safe_div`.

### The regression in `Feb_22`

`compute_U_i` in `Feb_22` is a botched edit of the `Jan_12b` version. The cylindrical line
was deleted and the small-δ fallback was left sitting inside the `if` branch:

```python
# Jan_12b — correct
if log_term > 1e-9:
    h_e = k_m / r_e / log_term        # cylindrical shell
else:
    h_e = k_m / delta                 # thin-layer limit

# Feb_22 — the cylindrical line is gone
if log_term > 1e-9:
    h_e = k_m / delta                 # plane slab, used for ALL delta
# (no else — h_e unbound if log_term <= 1e-9)
```

Since `r_e·ln(1+δ/r_e) < δ`, the plane form **underestimates** the melt-layer conductance —
by roughly 34 % at δ ≈ 15.6 mm, more as melting proceeds. The interface temperature driving
the Stefan front was still computed cylindrically, so the two halves of the energy balance
disagreed.

**Consequence for the conference paper.** This error and the `k_w` error (§1.1) push in the
*same* direction: both suppress charging heat transfer. The corrected model should charge
substantially faster and require materially fewer wells. Expect the headline numbers to move
a lot, and in a favourable direction — which is worth knowing before writing any text.

---

## Revised plan

**The solver base is `Jan_12b`, not `Feb_22`.** Do not port fixes into `Feb_22`.

1. **Resolve the `k_w` question** (audit §1.1) — see below. Then fix, re-run the `Jan_12b`
   baseline, and record how far `optimal_N_wells_ch`, `RTE_real` and `m_m_eff` move.
2. **Fix the `RTE` label** (§1.3). Still present in `Jan_12b` at the results dict:
   `'RTE': RTE_real * m_m_eff`. Store both quantities under distinct names.
3. **Silent-zero fallbacks** (§2.2). In `Jan_12b` these live in `compute_delta2_fast`
   (returns `0.0` on several paths) and in the `except` around its call. Note `h_e = 1e9`
   still carries the comment *"Assign a small value, leading to large resistance"* — the
   comment is inverted; 1e9 W/m²·K is essentially zero resistance, i.e. maximum heat
   transfer. Convert these to exceptions and expect breakage.
4. **Delete or validate dead code.** `compute_delta2` (slow) and the non-`_fast` solvers are
   defined but never called in `Jan_12b`. Either delete them or run one validation case
   proving `_fast` agrees with the reference implementation, then delete.
5. **Port the DoE onto the `Jan_12b` core.** It currently exists only in `Dec_5`, on the
   oldest solver. Regenerate the design at Resolution IV first (§1.2) — odd-order generators
   `ABC, ABD, ABE, ABF, ABG, ACD, ACE, ACF, ACG`, verified to give shortest word length 4 in
   the same 128 runs.
6. **Add version stamping and provenance.** `Jan_12b` has none. Stamp version + commit hash
   into the header cell, the runtime banner, every figure caption, and every CSV.
7. **Refactor to `thums_core/` as a package**, then the regression tests: energy closure
   (`∫Q dt` vs `ρ·h_m·V_melt`), `V_melt ≤ V_well` and `m_m_eff ≤ 1`, and zero-loss
   `RTE = 1`. Follow the existing `verification/harness.py` pattern already in this repo.

---

## The `k_w` question — needs the author's answer first

The charging chain has **no steel conductivity anywhere**. In both `Feb_22` and `Jan_12b`,
`evaluate_Q_ratio_ch` declares its conductivity parameter as `k_m_l` and passes it into
`time_profiles_melt`'s `k_w` slot. The discharge twin declares `k_w` and passes it correctly.

So during charging, `R_wall = (r_i/k_w)·ln(r_e/r_i)` uses 0.45 W/m·K instead of 45, and
`fin_m = sqrt(2·h_e/(k_w·fin_t))` gives fin efficiency ≈0.33 instead of ≈0.97.

**What changed between versions:** in `Feb_22` the *plotting* cells passed `k_w` correctly
while the sizing passed `k_m_l` — so figures and table came from different models. In
`Jan_12b` the plotting cells were changed to `k_m_l` as well. Someone harmonised on the PCM
value rather than the steel value.

That harmonisation is why this must be asked rather than assumed. Two readings:

- **A propagated typo** that was later "fixed" by making the inconsistency consistent. Then
  the correct action is to add a `k_w` parameter to the charging chain, mirroring discharge.
- **A deliberate simplification** — e.g. an intent to model a non-metallic liner, or a
  worst-case charging assumption. Then it needs a justifying sentence in the paper, and the
  discharge side needs the same treatment for symmetry.

There is no physical basis for tube-wall conductivity to differ between charge and
discharge, so one of the two paths is wrong either way. Ask the author before editing. If
the answer is "I don't recall," fix it as a typo, run both, and compare — the difference is
the honest measure of how much it mattered.

### The edit, if it is a typo

Mirror the discharge signature. In `Jan_12b`:

- `evaluate_Q_ratio_ch` — add `k_w` before `k_m_l` in the signature; in the
  `time_profiles_melt` call, change the **sixth positional argument only** from `k_m_l` to
  `k_w` (the eleventh stays `k_m_l` — that one really is the PCM).
- `find_N_wells_for_Q_ratio_ch_fast` — same addition, pass through.
- Charging call sites and the standalone charge calls in the plotting cells — add `k_w`.
- Leave the discharge path untouched; it is already correct.

**Verification:** afterwards, searching the file for `T_m_lay, k_m_l, Rf_i_prime` must
return zero matches. At the baseline point, `fin_eff` during charging should move from
≈0.33 to ≈0.97 and `optimal_N_wells_ch` should fall substantially. If nothing moves, the
edit did not take effect — suspect a stale Colab session before concluding the bug was
benign.

---

## Housekeeping

- The "Open in Colab" badge in the first markdown cell of each notebook still points at the
  old repo-root path. Fix in the first commit.
- There is a stale `origin/THUMS` branch at `fc81d83`. Confirm it holds nothing unique,
  then delete it.
- `T_geo` / `geo_grad` exist in `Feb_22` (computed, stored, never used in any balance) and
  do **not** exist at all in `Jan_12b`. Audit §2.4 stands: the formation is adiabatic and
  the PCM is implicitly initialised at its own melting temperature. If the geothermal
  gradient is to remain a DoE factor, it has to actually enter the energy balance.
