# THUMS v0.2 — start here

This archive is a complete, tested Python package that replaces the physics in four
Colab notebooks, plus the documents explaining why.

**Repository:** `https://github.com/JaderBarbosaLSU/wellbore-thermal-storage`
**Built:** 2026-08-05, in a Cowork session that could read the repository but not push to it.

---

## For the assistant picking this up

Unzip at the repository root and commit the lot, then tag `v0.2`:

```bash
unzip thums_v0.2.zip -d /path/to/wellbore-thermal-storage
cd /path/to/wellbore-thermal-storage
git add pyproject.toml thums_core tools notebooks docs verification START_HERE.md README_PACKAGE.md
git commit -m "thums_core v0.2: package skeleton, energy-closure fix, two-constraint sizing"
git tag v0.2 && git push origin main --tags
```

Read in this order before changing anything:

1. `docs/FIX_energy_closure.md` — the substantive finding and the correction
2. `docs/FINDING_compensating_errors.md` — why fixing one defect alone breaks the model
3. `docs/THUMS_code_audit_2026-08-05.md` — the full audit of the original notebooks
4. `docs/THUMS_unified_design_spec.md` — architecture and migration order
5. `docs/THUMS_unification_and_KPI_plan.md` — the KPI registry and post-processing design
6. `docs/Colab_GitHub_Token_Workflow.docx` — the author's own workflow rules; follow them

## Status after v0.2 (done in the follow-up session)

Done:

- package renamed `thums_core/` so the importable tree holds no notebooks; the four
  legacy notebooks stay in `thums/`
- `well_marched.march` is wired into `system.run_cycle`; `Case.front` selects
  `"marched"` (default) or `"legacy"` (reproduces v0.1 exactly, 4e-11)
- `verification/reference/thums_v0_baseline.json` regenerated with the marched front;
  the legacy values it replaced are at tag `v0.2` and in `docs/changes_from_IHTC.md` §6
- `docs/changes_from_IHTC.md` written -- read this before anything else

Still open, in order:

1. make `loss_surplus` an output rather than an assumed 5 % -- it is now measurable
   and it is wrong by 5-22 % once `k_w` is corrected
2. correct the `k_w` argument in the charging chain (largest single mover:
   N_wells 21.7 -> 12.8, binding constraint flips to inventory)
3. set `Case.strict = True` and fix `delta_max`
4. add PCM sensible heat
5. port the DoE onto the unified core using the Resolution IV design in `thums_core/doe.py`

Note `UPLOAD_INSTRUCTIONS.md` describes the manual GitHub drag-and-drop route and is
obsolete now that the package is in git.

---

## What is in here

```
pyproject.toml                  pip install git+https://github.com/... works after this lands
thums_core/
  __init__.py                   version stamp
  errors.py                     typed exceptions; nothing substitutes a value on failure
  config.py                     frozen dataclasses: PCM, Geometry, Cycle, Numerics, Case
  kpi.py                        the KPI registry (label, unit, direction, dependencies)
  doe.py                        Resolution IV fractional factorial + alias enumeration
  results.py                    tidy results schema + run manifest
  stefan.py                     melt front by energy balance  <-- the fix
  well_marched.py               time-marched well, exact closure  <-- the fix
  sizing.py                     N_wells = max(heat transfer, inventory)  <-- new constraint
  system.py                     run_cycle: the strangler boundary over the legacy physics
  _legacy_physics.py            AUTO-GENERATED from the v0.1 notebook; do not edit
  post/effects.py               orthogonal contrasts, Lenth PSE, half-normal plot
  post/pareto.py                real non-dominated sorting
tools/
  extract_legacy.py             regenerate _legacy_physics.py from the notebook
  freeze_baseline.py            step 0: freeze reference results before any change
  diagnose_closure.py           measure the energy-closure error and attribute it
notebooks/
  01_energy_closure.ipynb       thin Colab case study; installs the package from GitHub
verification/reference/
  thums_v0_baseline.json        frozen v0.1 results, defects intact
docs/                           audit, design spec, KPI plan, findings, workflow rules
```

## Verified in the session that built it

- `run_cycle(BASELINE)` reproduces the v0.1 notebook exactly: COP 2.956330409531823,
  eta_ORC 0.23685132027873143, RTE-no-pump 0.4348350503020928
- energy closure: legacy 3–54 % error, marched **1e-16** (machine precision)
- sizing with correct steel conductivity: **13.65 wells**, binding constraint
  **inventory** (heat transfer needs only 8.37) — against ~24 from the legacy model
- Resolution IV design: orthogonal, shortest word length 4, no main effect aliased
  with any two-factor interaction; the legacy design is Resolution III with twelve
  length-3 words
- `notebooks/01_energy_closure.ipynb` runs end to end

Environment: Python 3.11, NumPy 2.4.4, pandas 3.0.2, SciPy 1.17.1, CoolProp 8.0.0.
Note `np.trapezoid` in the legacy physics requires NumPy >= 2.0.

## Known open items

- `system.run_cycle` still uses the legacy front; wiring `march` in is the next step
- the `k_w` argument in the charging chain is still the legacy value by default
  (`Case.charge_uses_wall_conductivity = False`) so the frozen fixture stays valid
- `Case.strict = False` so geometric violations warn rather than raise; flip after the port
- no formation heat loss, no natural convection in the melt, no buoyancy head in the
  1,524 m loop (±7.5 bar, asymmetric between charge and discharge) — see the audit
