# Results pack v2.8 — back to the BOREAS paper task

*From the model/notebook task, 2026-11-03. Attach
`BOREAS_results_pack_v28.zip`. The parametric pack is unchanged from v4.*

All five items are done, and **the renumbering is a verified no-op: zero
macros changed, twelve added, none removed.** 22 checks pass.

---

## The macro diff, computed

| file | changed | added | removed |
|---|---|---|---|
| `boreas_numbers.tex` | **0** | 3 | 0 |
| `boreas_profile_numbers.tex` | **0** | 0 | 0 |
| `boreas_orange_numbers.tex` | **0** | 3 | 0 |
| `boreas_water_numbers.tex` | **0** | 0 | 0 |
| `boreas_ts_numbers.tex` | **0** | 0 | 0 |
| `boreas_video_numbers.tex` | **0** | 6 | 0 |
| `par_numbers.tex` | **0** | 0 | 0 |

Added: the three `…EtaRTEFixedStrict`; `\FrictionCite{Blasius,Nikuradse}`
and `\FrictionReMax`; the six `\LineLoss*`. Nothing renamed in the macro
namespace — the renaming is inside the code, where `DT_M_1D` became
`DT_M_4D`.

## 1. The strict lower bound

| | strict | ≤ averaged | ≤ sliding |
|---|---|---|---|
| **R4** | **0.1333** | 0.1447 | 0.1502 |
| R3 | 0.2565 | 0.2757 | 0.2776 |
| adiabatic | 0.2328 | 0.3212 | 0.3216 |

`\RefFourEtaRTEFixedStrict` and the other two. I took the crude version
you allowed: the below-mean energy is **not converted at all**. So the
bracket at R4 is −7.9 % / +3.8 % about the averaged value, and your
framing is right — it sits between two plants rather than being simply
conservative.

Worth noting the adiabatic case brackets *much* wider (0.2328–0.3216),
because 27 % of its energy arrives below the mean. The bracket is
narrowest where the store is lossiest, which is the opposite of what one
would guess.

## 2. Friction: the discontinuity, and the citations

You are right, and it is a 56 % step: Nikuradse gives 0.01156 at 10⁶
against the constant 0.018. **Extended to 3×10⁶**, its stated range;
above that the run now **raises** rather than extrapolating. The 10⁵ seam
is 0.89 %, as you found.

`\FrictionCiteBlasius` = Blasius 1913, `\FrictionCiteNikuradse` =
Nikuradse 1932, `\FrictionReMax` = 3e6. The forms match what is published
under those names; I have not seen the primary sources, so please check
them against whatever you cite elsewhere before they go in a reference
list.

## 3. Spelling

`mkS1.py` fixed at the generator: "realized", "vapor".

## 4. The renumbering

Implemented as specified. What the code looked like before is worth
recording, because it shaped the work:

- there was **no `T_1d` at all**, and `T_3d` carried the comment *"the
  borehole inlet, state 1d"* — so the old numbering was only ever in a
  comment;
- `T_4c = T_3c = 160.000 °C` already held exactly, so the charge-side
  identity cost nothing;
- the approach was `DT_M_1D`, not `DT_M1D`.

Now:

| | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| charge | 105.00 | 105.00 | 160.00 | 160.00 °C |
| discharge | 116.54 | 116.54 | 91.11 | 91.11 °C |

- `T_1d` is new and is **what the march returns** — the field outlet. The
  evaporator inlet `T_2d` is derived from it as `T_1d − line_dc_dT`. The
  correction pass now hands back `T_1d_realized`, which is what it always
  meant.
- `T_3c = T_4c + line_ch_dT`: the causality runs 4c → 3c, because the
  field inlet is pinned by the top of the cascade and the condenser has to
  work harder to overcome the line.
- `DT_M_1D` → `DT_M_4D`, with the old name accepted, warning, and cleared
  in `__post_init__` so it can never reach a cache key as a second source
  of truth. An explicit `DT_M_4D` always wins.
- Pump heating is written as two explicit identities with a comment, not
  left implicit.
- Line pressure drops are added to the field pumping as a field-level
  duty.
- `tab_states_water.tex` (`tab:states-water`) is new: all eight states
  with T, p, h. `water_states.json` carries the same.
- `\LineLoss{Ch,Dc}` = "not computed", plus the four zero values.
- `layout_fig.py` replaced with your version.

**V19** asserts the four identities hold exactly at zero loss. **V20**
asserts the placeholders do something real at 1 K — and it caught my own
wrong expectation, which is worth passing on: **T_3e does not follow
T_2d one for one.** A 1 K drop at 2d moves T_3e by **0.1826 K**, because
the pinch sits at 15.3 % of the duty and the water there moves by only
that fraction. That is your own evaporator geometry confirming itself
from the other direction, and V20 now checks against 0.153 rather than
against 1.

## 5. The notebooks

Both rebuilt and re-executed end to end.

- Cell 2.1 gains exactly the four line-loss inputs, each 0 and commented
  "placeholders, not yet modeled". Nothing else was added.
- The audit still prints its single line: *REF differs from CASE in 1
  field(s)*.
- A markdown block gives the numbering table and your sentence about the
  coils.
- All eight water temperatures are printed, with the line losses beside
  them, so a student sees the placeholders in the bookkeeping:

```
water loop at the reference case, 1 = field outlet, 4 = field inlet:
  charge     1c  105.00 C  2c  105.00 C  3c  160.00 C  4c  160.00 C
  discharge  1d  116.54 C  2d  116.54 C  3d   91.11 C  4d   91.11 C
  line losses: charge 0.0 K, discharge 0.0 K   (placeholders)
```

- The reference reproduces: loss 0.524152, e′ 0.215101, η_RTE 0.144675,
  T̄_2d 116.5440 °C — pack v2.7 to every published digit.
- The parametric notebook was regenerated, not edited, and
  `par_numbers.tex` is unchanged: 0 changed, 0 added, 0 removed.

### One thing that went wrong, and the check it bought

A builder edit left a dangling string in the model notebook's fetch cell,
and **the notebook was a SyntaxError for a full round**. V18 passed the
whole time — it hashes the *embedded modules*, and says nothing about the
notebook's own cells. The pack was fine, so nothing else noticed.

**V21** now parses every code cell. It is the cheap half of running the
notebook, and it would have caught this immediately. Running it properly
is still the only complete check, and I do that before each ship.

---

## Still open — the authors' items

Unchanged, and all three now have a number attached to the cost of not
deciding:

- **the reference h_m**, 380 kJ/kg at 150 °C. At 260 the field needs 69
  storage wells per MWe instead of 37;
- **the casing annulus data**, worth up to 2.7 points of loss and 5.8 %
  of e′ at R4;
- **the shoe depth per well**, which sets the reference length.
