# `thums_core/` — first instalment

What is here, what it replaces, and what still has to be built. Companion to
`docs/THUMS_unified_design_spec.md` and `docs/THUMS_unification_and_KPI_plan.md`.

## Status

These are steps 1–3 of the build order: **the layers that change no numbers.**
Nothing here touches the physics. That is deliberate — it is what makes the
later physics fixes attributable to a cause.

| file | replaces | tested |
|---|---|---|
| `thums_core/errors.py` | the silent fallbacks | — |
| `thums_core/kpi.py` | KPI definitions scattered across 3 notebooks | ✔ |
| `thums_core/doe.py` | `build_128run_fractional_16factor` | ✔ Res IV verified |
| `thums_core/results.py` | ad-hoc CSV writing with no provenance | ✔ |
| `thums_core/post/effects.py` | `main_effects`, `compute_main_effects` | ✔ |
| `thums_core/post/pareto.py` | the `Pareto` notebook (which has no front) | ✔ |
| `tools/extract_legacy.py` | — | ✔ |

Still to build: `config.py`, `props.py`, `network.py`, `stefan.py`, `well.py`,
`cycle.py`, `system.py`, `runner.py`, `post/case.py`, `post/fields.py`,
`post/exergy.py`.

## Verified behaviour

Run against the real 16-factor problem, with a synthetic response carrying three
planted main effects and one planted two-factor interaction:

```
Resolution IV (shortest word length 4)
word-length distribution: {4: 43, 6: 96, 8: 207}
columns orthogonal: True
no main effect aliased with a 2FI: True

top 5 effects on eta_rte
  factor    effect  significant      (planted: T_m_C 0.060, h_m 0.036, DT_3C_2C -0.022)
   T_m_C  0.060060         True
     h_m  0.036030         True
DT_3C_2C -0.021937         True

Pareto front over (eta_rte, eps_pcm): 10 non-dominated of 129
dependent-objective guard: ValueError ... ['rho_E'] is derived from E_well
incomplete-design guard:  ValueError ... 1 non-converged run(s); contrasts are not orthogonal
```

The planted effects come back exactly (effect = 2 × coefficient), the planted
`T_m_C × r_e` interaction does **not** contaminate any main effect — which is the
whole point of Resolution IV — and both guards fire.

For comparison, the same module reproduces your legacy design on request:

```
LEGACY: Resolution III (shortest word length 3), 12 words of length 3
T_m_C aliased with: rho_m_s*h_m, rho_m_l*r_e, cp_m_s*fin_t,
                    cp_m_l*fin_L, k_m_s*t_ch, k_m_l*DT_3C_2C
```

`FractionalFactorial` refuses to construct a design below Resolution IV. The
legacy design is reachable only through the explicit
`FractionalFactorial.allow_low_resolution(...)`, so it can be reproduced for the
change log but not used by accident.

## Three design decisions worth knowing

**`weighted_RTE` is gone.** `eta_rte` and `eps_pcm` are separate KPIs. The
efficiency-versus-utilization trade-off is what a Pareto front is for; folding it
into one number required a weighting nobody stated.

**Centre points are not replicates.** `n_centre` defaults to **1**, not 8. The
simulator is deterministic, so repeated centre runs give exactly zero pure error
and cannot support an F test for curvature. Seven of your eight centre runs were
carrying no information. `effects.curvature()` returns the contrast and says so
rather than manufacturing a p-value.

**Effects come from the design, not from group means.** `df[df[f]==high].mean()`
is correct only when the design is balanced and every run converged. Because
failures are written as NaN and `pandas.mean()` skips them silently, one failed
run makes the contrast non-orthogonal with no message. `main_effects` computes
the contrast from the coded matrix and raises if the design is incomplete.

## Environment note

`np.trapezoid` in the legacy code requires NumPy ≥ 2.0. Verified here on
NumPy 2.4.4, pandas 3.0.2, SciPy 1.17.1, CoolProp 8.0.0. On NumPy 1.x it is an
`AttributeError` — the exact Colab-versus-elsewhere trap your workflow note
describes.

## Next

1. `config.py` + `props.py`, then port `network.py` behind them.
2. `system.run_cycle()` driving `_legacy_physics` unchanged — this is the
   strangler step that lets the new schema and post-processing run against the
   old solver, so nothing waits on the physics port.
3. Freeze `verification/reference/thums_v0_baseline.json`, then start
   substituting ported modules one at a time, checking the fixture each time.
