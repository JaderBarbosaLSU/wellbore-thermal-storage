# Which note covers what

The results pack is **cumulative**: `results_pack_v2/` is one directory,
rewritten each round, so the newest zip is the only one needed and every
file in it was generated from the same `pack_runs.json` (V16 enforces
that).

The handoff notes are **not** cumulative. Each answers one review, so
they are kept side by side:

| note | answers | carries |
|---|---|---|
| `HANDOFF_BACK_v2.4.md` | review §1–§7 and §10 | the per-tube rate error and V17; eps_dc relabelled and the below-pinch share; the profile figure fixes; the ring extent macros; the T–s figure at R4; the R3 casing maxima; all 27 section 10 checks, including the three that needed more than a number (the integrator is not explicit; the return-bend double count; the grid is good to a few tenths of a point) |
| `HANDOFF_BACK_v2.8.md` | review v2.8 | the water loop renumbered (1 = field outlet, 4 = field inlet) with zero-default line-loss placeholders, verified as a no-op; the strict lower bound on eta_RTE; the friction discontinuity above Re = 1e6; Table S1's water block. Supersedes the v2.7 note; its content is unchanged and still in the pack |

Read v2.4 first, then v2.8.

Earlier notes (v2.0–v2.3) are in the repository history and are
superseded; nothing in them is still outstanding.
