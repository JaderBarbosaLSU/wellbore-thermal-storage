# Results pack v2.2 — back to the BOREAS paper task

*From the model/notebook task, 2026-10-30. Paste this whole file as the
prompt and attach `BOREAS_results_pack_v22.zip`.*

All of §5 is implemented. **You were right about 5.1, and I checked it
independently before changing anything** — the measurement is below. Every
one of the 111 runs and every derived table has been recomputed.

---

## 5.1 The symmetrization — confirmed, removed

Reproduced here from scratch, at R4, with `spsolve` on both forms:

| L | stations | Δs | ‖G−Gᵀ‖/‖G‖ | symmetrized: negatives | as assembled: negatives |
|---|---|---|---|---|---|
| 1250 | 25 | 50 m | 0.25 | 0 | 0 |
| 1250 | 50 | 25 m | 0.90 | **130** | **0** |
| 1250 | 100 | 12.5 m | 0.89 | **67** | **0** |
| 400 | 8 | 50 m | 0.07 | 0 | 0 |
| 400 | 25 | 16 m | 0.77 | **264** | **0** |
| 400 | 50 | 8 m | 0.76 | **569** | **0** |

Your counts, to the run. The fix is in `Ground.profiles`:

- the system is solved **as assembled**, `A = G + (R_near + G_ii) I`;
- above `_DIRECT_MAX` it uses **BiCGSTAB** with a Jacobi preconditioner.
  CG is gone. It agrees with the direct solve to 1.7e-12 wherever both run;
- the diagonal-dominance comment is corrected. Measured: diagonal
  0.567 m K/W against a largest off-diagonal row sum of **3.73**, so the
  matrix is neither symmetric nor diagonally dominant and CG was doubly
  unjustified.

Two new assertions guard the regression: **V14** (G is asymmetric and is
not symmetrized; zero negatives as assembled) and **V15** (the diagonal is
not dominant). All **16 checks pass**.

### The convergence narrative, withdrawn and replaced

`station_convergence.csv` is redone without the symmetrization, on a
convergence criterion rather than on the absence of negatives. R4,
L = 1250 m:

| Δs | stations | negatives | off-station weight | loss | e′ |
|---|---|---|---|---|---|
| 104 m | 12 | 0 | 25.9 % | 51.75 % | 0.2177 |
| **50 m** | **25** | **0** | **5.2 %** | **52.42 %** | **0.2151** |
| 25 m | 50 | 0 | 25.8 % | 52.42 % | 0.2150 |
| 12.5 m | 100 | 0 | 27.6 % | 52.50 % | 0.2147 |
| 6.25 m | 200 | 0 | 32.7 % | 52.53 % | 0.2146 |

**Zero negatives at every spacing, for R1, R3 and R4, at both lengths.**
"The coarse grid is the well-posed one" is withdrawn, and so is the
suggestion that the ratio rule in `stations_for` was physics — it is now
just the spacing that sets the count.

Δs = 50 m is kept, now justified on convergence: the loss moves 0.001
points and e′ 0.05 % to the next finer grid, 0.12 points across the whole
eightfold refinement. 104 m is too coarse (0.67 points).

**One thing I will not claim.** At L = 400 m the same refinement leaves
about 0.4 points of residual drift (32.55 → 32.95 % from 50 m to 6.25 m),
and the sequence is not monotone. So the discretization is worth **up to a
few tenths of a point of loss**, and I would rather say that than quote a
0.1 pt convergence that only holds at the reference length. It does not
move the cutoff: 400–500 m, unchanged.

---

## 5.2 Reworded — the loss is paid twice

Dropped from the README and from the generated output. `mktables.py` now
prints the comparison itself, so the text cannot drift from the numbers:

| | adiabatic | R4 | |
|---|---|---|---|
| η_ORC | 0.1462 | 0.1401 | **−4 %** |
| η_RTE | 0.3212 | 0.1447 | **−55 %** |
| HTHP input | 3.10 | 6.89 MW per MWe | **2.2×** |
| storage wells per MWe | 18.5 | 37.2 | **2.0×** |

Only the ORC conversion is insensitive. Your figures exactly.

---

## 5.3 CSS cap — flagged

`tab_ref_results` gains a row, **energy convergence, last 5 cycles**: the
relative spread of ΔE_dc over the last five cycles. The adiabatic limit and
R3 hit the 80-cycle cap on a state drift that stalls at ~2e-8 J/m, while
their discharge integral is constant to the last digit. `css_report` in the
notebook prints the same quantity next to the cycle count, so the cap
cannot read as non-convergence anywhere.

---

## 5.4 Tables

- booktabs throughout (`\toprule`/`\midrule`/`\bottomrule`);
- `tab:ref-results`, with the hyphen;
- short column headers (adiabatic / R4 / R3 / R0v), definitions moved into
  the caption, `\footnotesize`, and siunitx units in brackets on the **row**
  labels, since each row has its own unit;
- `tab_T3e` gains three columns, not one: wells per MWe, and also η_RTE and
  the HTHP input, because those are the columns that carry §5.2;
- the 20.5 typo is gone — the table's 20.4 was right, and the note now
  reads from the same generated numbers.

## 5.5 The broken macro

`\RefHTHPInputPerMWe` was nan because `mknumbers.py` read values out of
`intensive.json` by their **LaTeX row labels**, and that key was never
written. Both halves are fixed: the dump now carries plain keys
(`HTHP_in_MW_per_MWe` and the rest), and `mknumbers.py` **raises** if any
macro would expand to nan, inf, an empty string or `--`. The macro is
6.888.

---

## Every macro whose value changed

151 → 161 macros. **19 changed, 10 added, none removed.**

| macro | v2.1 | v2.2 |
|---|---|---|
| `\RefHTHPInputPerMWe` | **nan** | **6.888** |
| `\ConvLossSpread` | 0.29 | 0.12 |
| `\ConvEpSpread` | 0.51 | 0.23 |
| `\ConvNegAtTwoHundred` | 490 | **0** |
| `\RefRgMean` | 0.913 | 0.915 |
| `\RefFourEprime` | 0.2153 | 0.2151 |
| `\RefFourEtaStorage` | 0.4761 | 0.4758 |
| `\RefFourTTwod`, `\RefTtwod` | 116.6 | 116.5 |
| `\RefFourDEch` | 4.352 | 4.351 |
| `\RefFourDEdc` | 2.072 | 2.070 |
| `\RefFourPumping` | 4.85 | 4.86 |
| `\RefFourUAPerMWe` | 3480 | 3481 |
| `\RefFieldMWh` | 97.7 | 97.6 |
| `\LadderTrajectoryStep` | +39.0 | +39.1 |
| `\DutyHPE` | 9.63 | 9.64 |
| `\UAHPEPerMWe` | 1058 | 1059 |
| `\SensThFortyFiveEp` | 0.227 | 0.226 |
| `\SensTopThirtyLoss` | 39.1 | 39.2 |

Apart from the three `Conv*` macros and the broken one, **every change is
in the last digit** — your 0.2 % estimate on the mean R_g,eff was right
(0.913 → 0.915). Nothing in the paper's argument moves.

Added: `\AdiabaticHTHPInputPerMWe`, `\RefChargeElecRatio` (2.2) for §5.2,
and `\{Adiabatic,RefFour,RefThree,RefZeroV\}{CyclesToCSS,EnergyConv}` for
§5.3.

**Unchanged and worth saying so:** the cutoff is 400–500 m for R4, R4_th45
and R4_th70; V7 still regresses to the prototype at 17.39 % and 70.61 %;
V12 and V13 still hold at 1e-16 and 1e-10; no run fails to deliver.

---

## 5.6 Still owed

`BOREAS_parametric.ipynb` at R4 — the maps with L capped at the 5000 ft
shoe, SM Table S1, the water-balance numbers and the glide-shortfall
example. That is the next and last piece of work here.
