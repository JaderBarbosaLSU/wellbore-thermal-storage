"""Build notebooks/P2H2P_model.ipynb -- the notebook for someone who will USE
the model rather than audit it.

Design rule, and it is the opposite of the verification notebook's. That one
shows every line of code, extracted from source, because its job is to let a
reader check that the equations and the implementation agree. This one shows
every EQUATION and imports the implementation, because its job is to let a
reader change a parameter and know what will move.

The model arrives as one `%%writefile` cell, so the notebook is self-contained
in Colab without needing access to the repository, and the reader never has to
scroll through it.
"""
import datetime as _dt
import json
import pathlib

STAMP = (_dt.datetime.utcnow() - _dt.timedelta(hours=3)).strftime(
    '%Y-%m-%d %H:%M BRT')
MODEL = pathlib.Path('/tmp/build/thums.py')
OUT = pathlib.Path('/tmp/build/P2H2P_model.ipynb')


def md(t):
    return {'cell_type': 'markdown', 'metadata': {},
            'source': t.strip('\n').splitlines(True)}


def code(t):
    return {'cell_type': 'code', 'metadata': {}, 'execution_count': None,
            'outputs': [], 'source': t.strip('\n').splitlines(True)}


cells = []

# ---------------------------------------------------------------- title
cells.append(md(rf"""
# P2H2P wellbore storage — the model, and how to use it

Latent heat storage in the annulus of a repurposed oil-and-gas well. A
high-temperature heat pump melts a phase-change material during charging; an
organic Rankine cycle recovers the energy during discharging; pressurised water
circulates through finned hairpin tubes in the borehole.

*Model version 0.13 · notebook built {STAMP}*

---

### What this notebook is for

Change a parameter, re-run, see what moves. The equations are here; the code
that evaluates them is imported. If you want to see how a quantity is
*verified*, that is a different notebook — see the last section.

### The calculation, in seven steps

| | step | what is fixed, what is computed |
|---|---|---|
| 1 | thermodynamic cycles | heat pump and ORC evaluated once → COP, $\eta_{{\rm ORC}}$, and the water temperatures at each end of the store |
| 2 | energy budget | from the 1 MWe target and $\eta_{{\rm ORC}}$: how much energy must move each cycle |
| 3 | specify the hardware | $N_{{\rm wells}}$ from latent heat alone; both mass flows from their glides. **Nothing is solved for** |
| 4 | charge | march the exchanger for 10 h with inlet $T_{{4c}}$ |
| 5 | discharge | reverse the flow, march again with inlet $T_{{3d}}$ |
| 6 | repeat | until the store returns to its own initial state (cyclic steady state) |
| 7 | correct once | re-evaluate the ORC at the outlet the store actually produced |

Step 6 is the only iteration. There is no root finding anywhere: the output is
the **deviation** of delivered energy from target, reported rather than driven
to zero. Section 2 explains why it has to be that way.
"""))

# ---------------------------------------------------------------- the module
cells.append(md("""
---
## 0. The model

One cell, written to disk and imported. **You do not need to read it.** The
equations it implements are in §2 and the function that evaluates each one is
named there.

**To change a parameter, go to §3.1 — not here.** Editing a default in this
cell is the one thing that looks like it works and does not: the file on disk
changes, the run does not, and the old answer comes back with no warning.
(The import cell below now reloads, so it *would* work — but §3.1 is one line,
needs no reload, and shows the units.)

Edit this cell only to change the **physics**: a new correlation, a different
resistance network, an extra term. Note that the verification notebook checks
the version in the repository, not your copy.
"""))

cells.append(code('%%writefile thums.py\n' + MODEL.read_text()))

cells.append(md("""
### 0.1 Dependencies

`CoolProp` supplies the water and working-fluid properties and is not in a
stock Colab image, so it is installed here. This is the only cell that needs
the network.
"""))

cells.append(code(r"""
!pip install -q CoolProp
"""))

cells.append(code(r"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import importlib, sys

# RELOAD, do not just import. `%%writefile` above rewrites thums.py on disk,
# but `import` returns the module Python already has in memory: edit the model
# cell, re-run it, re-run this one, and you would silently keep running the OLD
# code. This makes that sequence work.
if 'thums' in sys.modules:
    importlib.reload(sys.modules['thums'])
import thums
from thums import (CASE, validate_case, cycle_state_points, energy_budget,
                   melting_temperatures, T_m_bottom, layer_map,
                   simulate_css, simulate_css_corrected, css_report,
                   march_h, pcm_capacities, pcm_state, unmirror_march,
                   performance_indices, kpi_report, calculate_pressure_drop,
                   exchanger_UA, ua_report,
                   exergy_audit, exergy_report, geothermal_resource,
                   water_h, water_s, stream_exergy_rate)
plt.rcParams.update({'figure.dpi': 110, 'font.size': 9})
pd.set_option('display.width', 200, 'display.max_columns', 30)
print('thums loaded ·', len(open('thums.py').read().splitlines()), 'lines')
"""))

# ---------------------------------------------------------------- assumptions
cells.append(md(r"""
---
## 1. What is modelled, and what is assumed

The store is one hairpin tube per leg-pair in a borehole annulus filled with
PCM, repeated $N_{\rm wells}$ times. Along the tube the PCM is divided into
$N_{\rm lay}$ **layers of different melting temperature** — a cascade — so that
the water meets progressively colder material as it gives up its heat.

Each of the $n_s$ segments is a lumped control volume. The assumptions that
matter:

| | assumption | where it bites |
|---|---|---|
| **H1** | the PCM in a segment is a single lumped state | no radial profile inside the melt |
| **H2** | no axial conduction between segments | the idle periods do nothing at all |
| **H3** | the melt region is a concentric annulus | **saturated**: the fronts of neighbouring tubes touch over 99 % of the well |
| **H5** | one melt-front thickness per segment | see H1 |
| **H6** | the borehole wall is adiabatic | no ground losses, so $\eta_{\rm storage}\equiv1$ at CSS |
| **H7** | the water is quasi-steady | the first ~38 min of each half-cycle is not resolved |
| **H8** | the cascade is laid along the developed length | the two legs of a hairpin carry different $T_m$ at the same depth |
| **H9** | no azimuthal exchange between the legs | they are up to 48.9 K apart at the wellhead |

H3, H6 and H9 are the ones to be careful about; §6 quantifies them.
"""))

# ---------------------------------------------------------------- equations
cells.append(md(r"""
---
## 2. The equations that are actually solved

### 2.1 The state variable

Each segment carries **one scalar**: $E'$, the enthalpy per unit tube length of
the PCM belonging to that segment, measured from a datum of *fully solid
material at its own melting temperature*. Everything else is recovered from it.

With $A_{\rm avail}$ the PCM cross-section per tube,

$$E'_{\rm lat}=\rho\,h_mA_{\rm avail},\qquad C_s=\rho c_{p,s}A_{\rm avail},\qquad C_l=\rho c_{p,l}A_{\rm avail}$$

the constitutive curve is single-valued in $E'$:

$$T_{\rm pcm}=\begin{cases}
T_m+E'/C_s, & E'<0 &\text{(subcooled solid)}\\
T_m, & 0\le E'\le E'_{\rm lat} &\text{(two-phase)}\\
T_m+(E'-E'_{\rm lat})/C_l, & E'>E'_{\rm lat} &\text{(superheated liquid)}
\end{cases}$$

$$A_{\rm melt}=\frac{\min\bigl(\max(E',0),\,E'_{\rm lat}\bigr)}{\rho h_m},
\qquad \varepsilon_{\rm local}=\frac{A_{\rm melt}}{A_{\rm avail}}\in[0,1]$$

The saturation is **the definition, not a guard**: a segment cannot melt PCM it
does not own, and the extra energy goes into superheat. So
$\varepsilon_{\rm local}\in[0,1]$ cannot be violated by any code path.

> *Why not carry the melted area instead?* Because a saturated segment then has
> nowhere to put further heat: the model must either melt material that is not
> there, or discard the heat and break its own energy balance. It is not a
> marginal effect — in this case the store ends every charge fully molten, and
> 9.3 % of the energy passing through it each cycle is carried as superheat or
> subcooling.

`pcm_capacities`, `pcm_state`

### 2.2 The segment equation

Energy balance on the PCM control volume, with the water quasi-steady (H7):

$$\frac{{\rm d}E'_i}{{\rm d}t}=K_i\left(T_{f,i}-T_{{\rm pcm},i}\right),
\qquad
T_{f,i+1}=T_{f,i}-\frac{q'_i\,\Delta z}{\dot m c_p}$$

The conductance closing it is **not** $U_iP_i$. Deriving it from the two
control-volume balances gives the effectiveness form

$$\boxed{\;K_i=\frac{\dot m c_p\left(1-e^{-{\rm NTU}_i}\right)}{\Delta z}\;}
\qquad
{\rm NTU}_i=\frac{2\pi r_i U_i \Delta z}{\dot m c_p}$$

The two agree only as ${\rm NTU}\to0$; here they differ by 4 % to 33 %. Only
the effectiveness form carries the capacity-rate limit, which is what stops
heat flowing from cold to hot as the wall conductance improves.

`compute_U_i`, `march_h`

### 2.3 Integrating it

The equation is nonlinear — $T_{\rm pcm}(E')$ is piecewise — but **affine within
each branch**, so each branch is integrated in closed form:

$$E'(t+\Delta t)=\begin{cases}
E'+K(T_f-T_m)\Delta t & \text{plateau}\\[2pt]
E'_{\rm eq}+\bigl(E'-E'_{\rm eq}\bigr)e^{-K\Delta t/C} & \text{sensible branches}
\end{cases}$$

and a step that would cross a branch boundary is **split at the crossing**. The
exponential cannot overshoot its equilibrium, so the scheme is unconditionally
stable and energy closure is structural.

That matters here: the time grid is logarithmic, and its last step exceeds the
explicit-Euler stability limit by a factor of 7.7. An explicit implementation
diverges.

`advance_segment`

### 2.4 The cascade, and where each half-cycle enters

Layer $l$ sits $\Delta T_{4C,M}$ below the charging water at *its own* leading
face, so the layers are spaced $\Delta T_{3C,2C}/N_{\rm lay}$ apart and the
coldest one is

$$T_m^{\rm bot}=T_m^{\rm top}-\Delta T_{3C,2C}\,
\frac{N_{\rm lay}-1}{N_{\rm lay}}$$

**Only the two inlets may be prescribed.** Both outlets are results of the heat
transfer:

$$T_{4c}=T_m^{\rm top}+\Delta T_{4C,M}\quad(\text{superheat, drives melting}),
\qquad
T_{3d}=T_m^{\rm bot}-\Delta T_{M,1D}\quad(\text{subcooling, drives freezing})$$

The discharge marches from the opposite end, so both the cascade and the state
mirror with the march index; results come back in depth indexing.

`T_m_bottom`, `melting_temperatures`, `layer_map`, `unmirror_march`

### 2.5 Why sizing is a simulation and not a solve

At cyclic steady state the state returns to itself and the store is adiabatic,
so

$$\eta_{\rm storage}({\rm CSS})\equiv1$$

**identically**, for every well count and every flow. The charge and discharge
energy requirements therefore collapse into one equation, which fixes the ratio
of the two flow rates and says nothing about $N_{\rm wells}$. There is no root
to find.

So the hardware is specified and the answer is reported:

$$\frac{\Delta E_{\rm delivered}}{\Delta E_{\rm required}}
=\frac{\text{glide the water actually achieves}}{\Delta T_{3C,2C}}$$

— because the flow was pinned by the *assumed* glide. **The deviation is the
glide shortfall**, and it is the same number on both half-cycles.

One consequence worth internalising before you sweep anything: the deviation
cancels out of the round-trip efficiency exactly. A $+1.2\,\%$ deviation means
$+1.2\,\%$ more electricity out *and* $+1.2\,\%$ more in. It is a statement
about **plant size**, not about performance.

`simulate_css`, `simulate_css_corrected`
"""))

# ---------------------------------------------------------------- couplings
cells.append(md(r"""
---
## 3. The parameters, and how they are coupled

Read this before sweeping anything. Almost every parameter does more than one
job, and the second job is usually the one that surprises you.

### 3.0.1 The storage side

```
ΔT_3C,2C ──────────→ CHARGING mass flow only    ṁ_ch = Q̇_out,HP / cp / ΔT_3C,2C
 charging glide       ...and nothing else, now that ΔT_cascade is its own field

ΔT_3D,2D ──────┬───→ DISCHARGING mass flow      ṁ_dc = Q̇_in,ORC / cp / ΔT_3D,2D
 discharge glide└──→ T_2d (provisional) ───────────────────┐
                                                            │
ΔT_cascade ────┬───→ cascade spacing  (= ΔT_cascade/N_lay)  │
 the ladder    └───→ T_m,bottom ──→ T_3d ──────────────────→┤
                                                            │
N_lay ─────────┬───→ cascade spacing                        │
               └───→ T_m,bottom ──→ T_3d ──────────────────→┤
                                                            │
ΔT_M,1D ───────────→ T_3d  (discharge inlet) ──────────────→┤
                                                            ▼
                                            ┌───────────────────────────┐
                                            │  ORC EVAPORATOR PINCH     │
                                            │  T_3e = min(hot-end rule, │
                                            │             pinch)        │
                                            └───────────┬───────────────┘
ΔT_pinch,ORC ──────────────────────────────────────────→┤
                                                        ▼
                                  η_ORC ──→ energy budget ──→ N_wells
```

### 3.0.15 Where the machine losses act — and why this changed at v0.12

Read this if you are comparing numbers against anything older than v0.12.

**What the model used to do.** Every expansion and compression was
*isentropic*: `s_5e = s_3e`, `s_1h = s_2h`, and so on. The 0.85 efficiencies
were then applied as **multipliers on the work, downstream of the cycle**:

$$\dot Q_{\rm in,ORC}=\frac{\dot W_{\rm el,out}}{\eta_{\rm ORC}\,\eta_{\rm turb}\,\eta_{\rm gen}}
\qquad
\dot W_{\rm el,in}=\frac{\dot Q_{\rm out,HP}}{{\rm COP}\,\eta_{\rm comp}\,\eta_{\rm mot}}$$

So the efficiencies **were** in the model — but the *state points* were not
affected by them.

**Why that was not good enough.** Two reasons, one small and one fatal.

*The small one.* For the ORC efficiency the approximation is good: the
downstream multiplier gives 0.14772 against 0.14806 computed properly — 0.2 %.
For the **heat pump** it is much worse, because the closure runs backwards from
a fixed discharge state: 2.3455 against 2.4731, a **5.4 % under-estimate** of
the COP. An irreversible compressor landing on the same discharge needs a
*colder* suction, and the two-stage machine recovers part of the extra work
through the flash separator — which a single multiplier cannot represent.

*The fatal one.* An isentropic machine is **reversible**, so an exergy balance
across it returns **zero destruction**. The 15 % of work the multiplier removed
appeared nowhere as an irreversibility. A component-wise exergy map built on
those state points would have shown perfect turbines and compressors and
blamed the entire loss on the exchangers.

**What changed.** The efficiencies now act on the state points:

$$h_{5e}=h_{3e}-\eta_{t}\bigl(h_{3e}-h_{5e,s}\bigr)
\qquad
h_{2e}=h_{1e}+\frac{h_{2e,s}-h_{1e}}{\eta_{p}}$$

A real turbine leaves **hotter** than an isentropic one; a real pump absorbs
**more**. For the heat pump, the discharge is still held on the dew line at
$T_{2h}$, so the suction state is now found by bisection — it carries less
superheat than before (14.25 K against 26.53 K), and `epsilon_IHX_2` follows.

**What moved, at the design point:**

| | v0.11a | v0.12 | |
|---|---|---|---|
| $\eta_{\rm ORC}$ reported | 0.17379 | 0.14806 | was the *ideal* cycle |
| COP reported | 2.7594 | 2.4731 | was the *ideal* cycle |
| $N_{\rm wells}$ | 15.659 | 15.625 | −0.2 % |
| **$\eta_{\rm RTE}$** | **0.30033** | **0.31721** | **+5.6 %** |
| $UA$ total | 2277.7 | 2253.3 | −1.1 % |
| CSS deviation | +4.527 % | +4.527 % | *identical* |
| thermal/well, $\rho_E$ | — | — | *identical* |

Two things to take from that table. First, the η_ORC and COP rows are **not**
a degradation: they used to be the reversible cycle's numbers, and are now the
real machine's. The like-for-like comparison is η_RTE. Second, η_RTE went
**up**, because the old treatment over-penalised the heat pump.

**Everything on the storage side is untouched to the last digit**, which is the
check that this change stayed where it belongs.

**Try it yourself.** Set `eta_turb_s = eta_pump_s = eta_comp_s = 1.0` in §3.1
and re-run: you will recover 0.173787 and 2.759421 exactly — the reversible
cycles, to six decimals. That round trip is the regression test for this
change.

### 3.0.2 The plant side

```
T_m ───────────┬───→ the whole cascade ladder
               └───→ T_4c ─┐
ΔT_4C,M ───────────→ T_4c ─┤  (T_3c = T_4c: same water, nothing in between)
                           ├──→ T_2h = T_3c + ΔT_2H,3C ──┐
ΔT_2H,3C ──────────────────┘                              │
                                                          ├──→ COP ──→ budget
ΔT_3A,4A ──────┬───→ T_4a ──→ T_13h = T_4a − ΔT_pinch,HPE─┤
ΔT_pinch,HPE ──┘                                          │
ΔT_sub ───────────────────────────────────────────────────┘

T_sink, ΔT_sink,glide, ΔT_E,sink ──→ T_1e ──→ η_ORC
```

### 3.0.3 Everything else

```
geometry (r_e, fins, ────→ U_i → NTU → K      → rate only, not inventory
 L_well, D_well)     └───→ V_well             → inventory only, via N_wells

h_m, ρ ──────────────────→ inventory (N_wells ∝ 1/ρ h_m)
k_s, k_l ────────────────→ rate, through the melt-layer resistance
c_p,s, c_p,l ────────────→ how much energy the sensible branches carry
```

---

**The old trap is gone — but only because §3.1 says so.** Until v0.9b one field
set the charging flow, the discharging flow, the cascade spacing *and*
`T_m,bottom`. Sweeping the glide moved the whole ladder underneath you. Now
§3.1 writes `DT_3C_2C`, `DT_3D_2D` and `DT_cascade` out separately, so sweeping
`DT_3C_2C` moves **only the charging flow**:

| `DT_3C_2C` | `T_m,bottom` | span | `T_3d` |
|---|---|---|---|
| 40 K | 101.111 °C | 48.889 K | 91.111 °C |
| 55 K | 101.111 °C | 48.889 K | 91.111 °C |
| 70 K | 101.111 °C | 48.889 K | 91.111 °C |

⚠️ **This holds for `case`, not for `CASE`.** `DT_cascade=None` in the raw
dataclass means *follow `DT_3C_2C`*, which restores the old coupling. Always
sweep from `case` (which `sweep()` does by default), never from `CASE`.

**There are no clean knobs.** The previous version of this section claimed
$\Delta T_{4C,M}$ and $\Delta T_{M,1D}$ each moved one inlet and nothing else.
Both claims are false, and the reasons are worth knowing:

| sweep | what you expect | what also moves |
|---|---|---|
| `DT_4C_M` 5→20 K | charge inlet | $T_{3c}=T_{4c}$, so $T_{2h}$ rises 165→180 °C and **COP falls 2.884 → 2.531** |
| `DT_M_1D` 6→18 K | discharge inlet | $T_{3d}$ drops, so $T_{2d}$ and the pinch drop, and **η_ORC falls 0.1815 → 0.1574** |

The second one is the whole point of the pinch work: **anything that moves
$T_{2d}$ or $T_{3d}$ reaches the ORC**, because the evaporating temperature is
set by the composite curves and those start at $T_{3d}$.

So: sweep one parameter, read *both* columns of the KPI table, and when a
result surprises you come back to these three diagrams and follow the arrows.
"""))

cells.append(md(r"""
---
### 3.1 The case

**This is the only cell you need to edit.** `Case` is frozen, so variants are
made with `CASE.with_(field=value)` rather than by assigning to fields.

> **Do not edit the model cell to change a parameter.** It will look as though
> nothing happened. `%%writefile` rewrites `thums.py` on disk, but Python keeps
> the module it already imported, so the run continues with the old values and
> reports them without complaint. The import cell above now reloads, so the
> sequence *edit model cell → re-run it → re-run the import cell* does work —
> but the parameters below are the intended place, and need no reload at all.
>
> **Watch the units.** `h_m` is in **J/kg**, so a PCM with a latent heat of
> 180 kJ/kg is `180_000.0`. Several fields are in SI base units where the
> literature quotes kJ or mm; the comment after each one is authoritative.
"""))

cells.append(code(r"""
# =============================================================================
#  THE DESIGN POINT.  Edit the numbers here -- NOT in the model cell.
#  They are written out rather than left as defaults so that changing one is a
#  one-line edit with the units in front of you.
# =============================================================================
case = CASE.with_(
    # --- PCM -----------------------------------------------------------------
    T_m_C    = 150.0,        # melting temperature of the TOP layer      [C]
    h_m      = 380_000.0,    # latent heat of fusion             [J/kg]  <- J!
    rho_s    = 1550.0,       # solid density                            [kg/m3]
    rho_l    = 1450.0,       # liquid density                           [kg/m3]
    cp_s     = 1280.0,       # solid specific heat                    [J/kg/K]
    cp_l     = 1800.0,       # liquid specific heat                   [J/kg/K]
    k_s      = 0.60,         # solid conductivity                      [W/m/K]
    k_l      = 0.45,         # liquid conductivity                     [W/m/K]
    N_lay    = 9,            # cascade layers along the well               [-]

    # --- the three glides ----------------------------------------------------
    #  These were ONE field until v0.9b.  Two of the equalities they imply are
    #  forced by energy conservation and cannot be changed; the other two are
    #  design choices and are now yours.  See the note below the cell.
    DT_3C_2C = 55.0,         # CHARGING water glide, HTHP condenser         [K]
    DT_3D_2D = 55.0,         # DISCHARGING water glide, ORC evaporator      [K]
    DT_cascade = 55.0,       # cascade grading; SPAN = this * (N_lay-1)/N_lay
    #                          -> span 48.889 K.  MUST EXCEED 45 K, see below

    # --- borehole approaches -------------------------------------------------
    DT_4C_M  = 10.0,         # charge inlet above the hottest layer         [K]
    DT_M_1D  = 10.0,         # discharge inlet below the coldest layer      [K]

    # --- pinch / approach temperatures in the four exchangers ----------------
    DT_pinch_ORC = 5.0,      # ORC evaporator, minimum gap ANYWHERE         [K]
    #                          (interior pinch; this is what sets T_3e)
    DT_2D_3E = 12.0,         # ORC evaporator, hot-end approach             [K]
    #                          T_3e = min(T_2d - this, pinch-feasible).
    #                          INERT at this value: the pinch always binds.
    DT_pinch_HPE = 5.0,      # HTHP evaporator, cold-end approach           [K]
    #                          T_13h = T_source - DT_3A_4A - this
    DT_3A_4A = 10.0,         # how far the SOURCE is cooled                 [K]
    #                          PROMOTED at v0.13. This was a quiet approach
    #                          parameter; it is now the central design
    #                          variable of the source side. It sets BOTH the
    #                          evaporating temperature (so COP falls as it
    #                          rises) and the geothermal flow (which falls
    #                          much faster). Raising it 10 -> 20 K costs 7.8 %
    #                          of COP and halves the number of geothermal
    #                          wells. See the exergy table in section 4.3.
    DT_2H_3C = 10.0,         # HTHP condenser, hot-end approach             [K]
    DT_sub   = 2.0,          # refrigerant subcooling at condenser exit     [K]
    DT_E_sink = 5.0,         # ORC condenser approach, to the sink OUTLET   [K]
    DT_sink_glide = 0.0,     # sink temperature rise; 0 = infinite reservoir[K]
    #                          ABOVE ~5 K the ORC condenser curves CROSS
    T_source_C = 60.0,       # heat-pump source                             [C]
    T_sink_C   = 20.0,       # ORC sink (the ocean, for THUMS)              [C]

    # --- machine efficiencies ------------------------------------------------
    #  ISENTROPIC, and applied INSIDE the cycles from v0.12: they move the
    #  state points, not just the work.  Set all three to 1.0 to recover the
    #  reversible cycles exactly (useful for seeing what they cost).
    eta_turb_s = 0.85,       # ORC turbines, isentropic                     [-]
    eta_pump_s = 0.85,       # ORC pumps, isentropic                        [-]
    eta_comp_s = 0.85,       # heat-pump compressors, isentropic            [-]
    #  ELECTRICAL / MECHANICAL -- correctly OUTSIDE the working fluid
    ElG_eff = 0.95,          # ORC generator                                [-]
    ElH_eff = 0.95,          # compressor motor                             [-]

    # --- the geothermal resource (new at v0.13) ------------------------------
    #  The 60 C source is a FINITE geothermal flow from producing wells in the
    #  cluster, not an unlimited reservoir. The flow itself is NOT an input --
    #  it is derived, exactly as N_wells is:
    #      m_geo = Q_src / (h(T_source) - h(T_source - DT_3A_4A))
    #  and then divided by the per-well yield below to give N_geo. That puts
    #  the uncertainty in a field-measured quantity instead of in a
    #  plant-level flow rate nobody has measured.
    m_dot_geo_well = 12.5,   # yield of one geothermal producer          [kg/s]
    #                          PROVISIONAL. From a 4500 m / 206 C repurposed
    #                          well -- hotter and deeper than ours. Published
    #                          values for repurposed oil and gas wells run
    #                          from about 1 kg/s to this. N_geo is directly
    #                          proportional to it, so treat N_geo as a scaling
    #                          result and not as a prediction.
    T_reinject_min_C = 25.0, # lower bound on the reinjection temperature   [C]
    #                          This is what bounds DT_3A_4A from ABOVE. The
    #                          default is slack on purpose: it is a
    #                          placeholder for a formation constraint we do
    #                          not yet have a number for. WARNS, never raises.
    exergy_convention = 'resource',
    #                          'resource': charge all the exergy the well
    #                              lifts; the reinjected stream is booked as a
    #                              LOSS attributed to no component. Right for
    #                              a dedicated well, which is paid for whether
    #                              or not its exergy is used.
    #                          'stripped': charge only what the evaporator
    #                              removes; reinjection is free.
    #                          They give OPPOSITE guidance on DT_3A_4A --
    #                          'stripped' falls monotonically and pushes the
    #                          design toward the largest possible well count;
    #                          'resource' has an interior optimum near 20 K.

    # --- windows -------------------------------------------------------------
    t_ch     = 10.0,         # charging window                              [h]
    t_dc     = 10.0,         # discharging window                           [h]

    # --- geometry ------------------------------------------------------------
    num_fins = 24,           # fins per leg                                 [-]
    fin_L    = 0.0075,       # fin radial length                            [m]
    fin_t    = 0.0015,       # fin thickness                                [m]

    # --- numerics (not design variables) -------------------------------------
    n_segments = 100,        # segments along the developed tube
    n_times    = 40,         # logarithmic time levels per half-cycle
)

# =============================================================================
#  How many cycles the CSS loop may take before giving up.
#
#  WATCH FOR THE BANNER.  If a run prints
#
#        *** NOT CONVERGED ***
#
#  then the store had NOT returned to its own state when the loop stopped, and
#  every number below that line is meaningless -- the deviation, the melted
#  fractions, all of it.  The run does not fail, it just reports.  Raise this
#  and re-run.
#
#  It is here rather than buried in the calls because it is a number you WILL
#  have to change.  Cycles to convergence grow with the well count and with
#  anything that makes each well work less hard: the design point needs about
#  71, N = 17 needs 142, N = 20 does not converge in 400 at all.
# =============================================================================
N_CYCLES = 250

validate_case(case)
"""))

cells.append(md(r"""
### The three glides — which you may change, and which you may not

Three temperature spans in this plant are nearly the same size. Until v0.9b a
single field set all three, which hid the fact that they are not equally free:

| span | field | status |
|---|---|---|
| HTHP condenser water | `DT_3C_2C` | **specified** |
| borehole, charging | — | **forced** equal to `DT_3C_2C` |
| ORC evaporator water | `DT_3D_2D` | **specified** |
| borehole, discharging | — | **forced** equal to `DT_3D_2D` |
| PCM cascade span | `DT_cascade` | **specified** (span = `DT_cascade`·(N−1)/N) |

The two **forced** rows are not assumptions and you cannot sweep them apart.
`T_3c = T_4c`: the water leaves the heat-pump condenser and enters the borehole
with nothing in between, so its rise in one and its fall in the other are the
same number — energy conservation on a closed loop.

The other two are yours. Setting all three to 55.0 reproduces the old
behaviour exactly.

**`DT_cascade` has a lower bound, and it is tight.** With the span one layer
short of the glide, the approach between the water and the layer it is melting
is exactly `DT_4C_M` at each layer's *leading* face, decaying to 3.889 K at its
trailing face. That sawtooth **is** the cascade. Narrow the span and the
approach at the far end of the well closes:

$$\text{span} \;>\; \text{glide} - \Delta T_{4C,M} \;=\; 45\ \text{K}$$

against a span of 48.889 K — only **3.889 K of margin**. Below it the driving
difference inverts and `validate_case` raises. This is the direction that will
tempt you, because a narrower cascade raises `T_m,bottom` and lets the ORC boil
hotter: η_RTE climbs 0.300 → 0.351 while the residual melt goes 0.012 → 0.196
and the plant misses target by 33 %. Both `DT_cascade` and `DT_3D_2D` have
interior optima and **nobody has looked for them**.

### The pinch temperatures

`DT_pinch_ORC` is the one that matters. The ORC evaporator pinches in its
*interior*, not at an end, because the working fluid takes 76 % of its heat at
one temperature while the water glides — so an approach at the hot end is not a
constraint at all. `DT_2D_3E` is left in the block above because it *is* a real
constraint if you raise it, but at 12 K it never binds.

`DT_sink_glide` is 0 by default, which models the sink as an infinite
reservoir. Give it a value and the ORC condenser acquires an interior pinch at
its desuperheating corner (≈93 % of the duty) which **crosses at about 5 K** —
and the model does not check that one. Leave it at 0 unless you mean it.

> **`N_CYCLES` is the one numerical setting you will actually have to touch.**
> Everything else in the block above is physics. This one is a give-up limit on
> the cyclic-steady-state loop, and if it is too low the run still prints a full
> set of results — they are simply not results, because the store had not
> returned to its own state. The `*** NOT CONVERGED ***` banner is the only
> thing standing between you and a plausible-looking wrong answer, so read it.
>
> If `validate_case` warns that `n_segments` is not a multiple of `N_lay`, the
> layer boundaries do not fall on segment boundaries and one layer is short.
> It is a small effect here — worth about $0.1$ percentage points on the
> deviation — but it is the first thing to rule out if a sweep in `N_lay` looks
> ragged.
"""))

# ---------------------------------------------------------------- plant
cells.append(md(r"""
---
## 4. Running one design

### 4.1 The plant

Steps 1 and 2: the two cycles, then how much energy has to move.
"""))

cells.append(code(r"""
rank, hp, T = cycle_state_points(case)
Eb = energy_budget(case, rank['rank_eff'], hp['hp_cop'], T)
K = 273.15

print('CYCLES')
print(f"  COP (heat pump)        {hp['hp_cop']:8.4f}")
print(f"  eta_ORC                {rank['rank_eff']:8.4f}")
print()
print('WATER AT THE STORE        [C]')
print(f"  charge   in  T_4c      {T['T_4c']-K:8.3f}   = T_m,top + DT_4C_M")
print(f"           out T_2c      {T['T_2c']-K:8.3f}   provisional")
print(f"  discharge in T_3d      {T['T_3d']-K:8.3f}   = T_m,bot - DT_M_1D")
print(f"           out T_2d      {T['T_2d']-K:8.3f}   provisional")
print()
lay = np.array(melting_temperatures(case)[0]) - K
print('CASCADE                   [C]')
print('  ' + '  '.join(f'{v:.1f}' for v in lay))
print(f"  spacing {case.DT_3C_2C/case.N_lay:.3f} K    "
      f"T_m,bottom {T_m_bottom(case)-K:.3f} C")
print()
print('ENERGY PER CYCLE          [kJ, whole field]')
print(f"  required by the ORC    {Eb['D_E_in_ORC']:12.5e}")
print(f"  to be stored           {Eb['D_E_out_HP']:12.5e}")
"""))

cells.append(md(r"""
### 4.2 The store, marched to cyclic steady state

Steps 3–7. This is the whole calculation.
"""))

cells.append(code(r"""
# n_cycles: convergence slowed in v0.7 (the design point needs ~75 cycles,
# and it grows with the well count). The default of 80 is no longer
# comfortable -- see the NOT CONVERGED banner in the report below.
css = simulate_css_corrected(case, record=True, n_cycles=N_CYCLES)
css_report(case, css)
"""))

cells.append(md(r"""
### 4.2.1 Performance indicators

The indicators of the factorial study, evaluated at cyclic steady state.
These are the numbers to collect across a parametric run.
"""))

cells.append(code(r"""
kpi = performance_indices(case, css)
kpi_report(case, css, kpi)
"""))

cells.append(md(r"""
### 4.2.2 What the efficiency costs in hardware

Every efficiency above was bought with exchanger area, and until v0.10 nothing
said how much. That made every approach temperature a **free lunch**: tighten
`DT_pinch_ORC` and $\eta_{\rm ORC}$ rises with nothing to pay.

The conductance each exchanger needs follows from its composite curves —
the same construction as the pinch:

$$UA \;=\; \int_0^{Q}\frac{{\rm d}Q}{\Delta T(Q)},
\qquad
\Delta T_{\rm eff}\;\equiv\;\frac{Q}{UA}$$

$\Delta T_{\rm eff}$ is the single approach that would need the same $UA$. For a
counterflow exchanger with both streams sensible it **is** the log-mean
difference (verified to $4\times10^{-8}$); with a phase change on one side it is
the correct generalisation and the LMTD is not.
"""))

cells.append(code(r"""
ua = ua_report(case, css)
"""))

cells.append(md(r"""
Two things in that table are worth your attention.

**The ORC condenser dominates.** It needs more conductance than the other three
put together, because it rejects 5.9 MW across an approach of about 5 K to a
reservoir. `DT_E_sink` is a parameter nobody has swept, and it is the largest
single lever on total plant area.

**The pinch is no longer free.** Sweeping `DT_pinch_ORC`:

| `DT_pinch_ORC` | $T_{3e}$ | $\eta_{\rm RTE}$ | $UA$ evaporator | $UA$ total |
|---|---|---|---|---|
| 1 K | 100.3 °C | 0.3141 | 721.7 | 2490.4 |
| 2 K | 99.2 °C | 0.3107 | 582.3 | 2373.1 |
| **5 K** | **95.9 °C** | **0.3003** | **415.9** | **2277.7** |
| 10 K | 90.3 °C | 0.2823 | 311.9 | 2308.3 |

Going from 10 K to 1 K buys **11.3 %** of round-trip efficiency and costs
**131 %** more evaporator conductance. Note the total has a shallow minimum
near the default — tightening the evaporator shrinks the plant, which shrinks
the other three. That minimum is an optimisation nobody has run.

**For a DoE, report `UA_total_kW_K` beside `eta_RTE` always.** A design that
wins on efficiency alone may simply be one that specified a bigger exchanger.
"""))

cells.append(md(r"""
### 4.2.3 Where the work is actually lost

Everything above is the **first law**, and the first law cannot see what the
store costs. At cyclic steady state the PCM returns to its own initial state,
so

$$\eta_{\rm storage} \;\equiv\; 1$$

*by construction*. Every joule that goes down the well comes back up. That is
not a result — it is a statement that the energy balance is blind here. The
store's whole cost is a **degradation of temperature**, and only the second
law prices it.

Two things you need before the table makes sense.

**The dead state is the sink, $T_0 = T_{\rm sink} = 20$ °C.** It is the lowest
temperature in the system, so every exergy comes out positive and the heat
dumped at the ORC condenser is genuinely unrecoverable — which is what a dead
state should mean.

**There are *two* reservoirs, not one.** The ORC condenser rejects to the sink.
The heat-pump evaporator does **not** draw from the sink: its source is the
60 °C geothermal water. So the geothermal stream is a second *input*, and the
balance the audit checks is

$$W_{\rm el,in} + \mathcal{E}_{\rm geo}
  \;=\; W_{\rm el,out} + \mathcal{E}_{\rm reinj} + \sum_j I_j$$
"""))

cells.append(code(r"""
ex = exergy_report(case, css)
"""))

cells.append(md(r"""
#### Reading that output

**The three efficiencies at the bottom are the headline.** $\eta_{\rm RTE}$
books the geothermal heat as free. It is not — a dedicated producer is drilled
and paid for whether or not you use its exergy. Price it and 0.317 becomes
$\psi = 0.246$. Report both.

**Three readings of the table, two of them not what we expected.**

*The exchangers beat the machines.* Condenser, evaporators, ORC condenser and
borehole together are **49.4 %** of all destruction; the compressors and
turbines are **19.0 %**. This is a heat-transfer plant, not a turbomachinery
plant — convenient, because heat transfer is what this model resolves and
turbomachinery is what it merely stipulates.

*A valve destroys more than the high-pressure compressor.* The throttle into
the separator is third on the list, above either compressor. That is a design
signal: an expander is worth more here than any plausible gain in compressor
efficiency.

*The borehole is only 6.8 %.* Fourth, and lower than we guessed. The case for
the cascade and for tight approach temperatures rests on **well count and
hardware**, not on lost work. Do not oversell the store.

#### Two things about the check itself — read these

**The GATES are the test; the global sum is not.** Each machine's destruction
is summed from entropy generated *inside* its components and checked against a
boundary balance built from stream states alone. Those are independent, and
they agree to $\sim10^{-15}$.

The *global* line is an **identity**. Because $T_{3c} = T_{4c}$, the borehole
term — taken as the difference of the two water-stream exergies — cancels the
two exchanger terms exactly, and both sides reduce to the same expression. The
first version of this audit printed that identity as a closure check reading
`0.00e+00` and it was reported as a pass. **A residual of exactly zero over
twenty terms of order $10^3$ is not machine precision; it is a tell.** A check
that cannot fail cannot pass.

So the borehole number is correct but **unverified**. Its independent value
needs an entropy curve $s(E')$ alongside the existing $E'(T)$ and an integral
over $z$ and $t$. That is the next piece of work and it is not in v0.13.

**What the gate caught.** It did not close on first construction. The water
mass flow was set as $\dot Q/(c_p\,\Delta T)$ with $c_p$ at the mean
temperature; over a 55 K rise that is 0.096 % inconsistent with the model's own
enthalpies. Invisible to every energy balance — *both sides use the same wrong
flow* — but it fabricated 7.17 kW of exergy. And the formula was written
**twice**: patching one copy changed no reported number at all, which looked
like confirmation and meant nothing, because the live copy was elsewhere.
"""))

cells.append(md(r"""
### 4.2.4 The geothermal side, which we were not counting

The 60 °C source is a **finite** flow from producing wells in the cluster. Until
v0.13 the model drew whatever it needed and no output recorded how much.

The flow is **derived**, not prescribed — the same choice already made for
`N_wells`:

$$\dot m_{\rm geo}
 = \frac{\dot Q_{\rm src}}{h(T_{\rm source}) - h(T_{\rm source}-\Delta T_{3A4A})},
 \qquad
 N_{\rm geo} = \dot m_{\rm geo}\,/\,\dot m_{\rm geo,well}$$

Prescribing a plant-level flow would put the uncertainty in a number nobody has
measured. Deriving it moves the uncertainty into the **per-well yield**, which
is field data. The second law loses nothing by this: the source exergy is
homogeneous of degree one in $\dot m_{\rm geo}$, so it cancels.

`m_dot_geo_well = 12.5` kg/s is **provisional** — read `N_geo` as a scaling
result, not a prediction.

**`DT_3A_4A` is now the central design variable of the source side.** It sets
both the evaporating temperature (so COP falls as it rises) *and* the
geothermal flow (which falls much faster):

| `DT_3A_4A` | COP | $\eta_{\rm RTE}$ | $\dot m_{\rm geo}$ rel. | resource used | $\psi_{\rm stripped}$ | $\psi$ |
|---|---|---|---|---|---|---|
| 5 K | 2.582 | 0.345 | 1.000 | 22.6 % | 0.295 | 0.197 |
| **10 K** | **2.473** | **0.330** | **0.466** | **42.6 %** | **0.288** | **0.245** |
| 15 K | 2.373 | 0.317 | 0.289 | 59.7 % | 0.281 | 0.260 |
| 20 K | 2.280 | 0.305 | 0.202 | 73.9 % | 0.274 | **0.264** |
| 25 K | 2.194 | 0.293 | 0.151 | 85.2 % | 0.267 | 0.263 |
| 30 K | 2.114 | 0.283 | 0.117 | 93.3 % | 0.261 | 0.260 |

**The two exergy conventions give opposite advice.** `'stripped'` falls
monotonically — take as little $\Delta T$ as possible — which is exactly the
configuration needing the *most* geothermal wells. `'resource'` has an interior
optimum near 20 K. The default is `'resource'`, because a dedicated well is
paid for either way.

**The design point at 10 K is badly placed.** Moving to 20 K costs 7.8 % of COP
and cuts the geothermal flow to 43 % of its present value. At plant level that
is not a close call. What stops you is `T_reinject_min_C` — a formation
constraint we have no measured value for yet, which is why it warns rather than
raises.

**One trap.** Because $\dot m_{\rm geo}$ is derived rather than capped, nothing
stops a parametric study buying COP with an arbitrarily large geothermal flow.
**Carry `N_geo` and `wells_per_MWe` in every DoE table.** A study that ranks on
`eta_RTE` alone will find the free-source corner and recommend it.

**Try it yourself.** Set `DT_3A_4A = 20.0` in §3.1 and re-run §4.2–4.2.4. Watch
`eta_RTE` fall, `psi` rise, and `N_geo` halve.
"""))

cells.append(md(r"""
| indicator | meaning | what moves it |
|---|---|---|
| $\eta_T$ | discharge thermal $\to$ net electric | the power block only: $\eta_{\rm ORC}\eta_{\rm turb}\eta_{\rm gen}$. The store cannot change it |
| $\eta_{T,\rm eff}$ | the same, less the discharge pumping | tube diameter, flow rate, well depth |
| $\eta_{\rm RTE}$ | round trip, both parasitics charged | $\mathrm{COP}$, $\eta_{\rm ORC}$ and the parasitics. **Not** the deviation |
| $\varepsilon_{\rm RTE}$ | $\eta_{\rm RTE}\times$ the fraction of the store that cycles | anything that leaves PCM unused |
| $\Delta E_{\rm therm}$ | thermal energy per well per cycle [kWh] | $h_m$, $\rho$, $V_{\rm well}$ |
| $\Delta E_{\rm elec}$ | net electric energy per well [kWh] | the above, times $\eta_{T,\rm eff}$ |

Two warnings about reading these across a factorial table.

**$\eta_{\rm RTE}$ does not move with the deviation, at all.** A field that
delivers 1 % above target also drew 1 % more in, because both flows are pinned
by their glides. If a sweep changes the deviation and leaves $\eta_{\rm RTE}$
alone, it changed the *size* of the plant and nothing else.

**$\varepsilon_{\rm RTE}$ is weighted by the CYCLED fraction**, $\varepsilon$
at end of charge minus $\varepsilon$ at end of discharge — not by the melted
fraction as in the first-cycle study. The melted fraction saturates at $1$ under
this formulation, because a segment that has melted everything puts further
energy into superheat where $\varepsilon$ cannot see it; at this design point it
is exactly $1.0000$, which would make $\varepsilon_{\rm RTE}$ identical to
$\eta_{\rm RTE}$ and useless. The cycled fraction does not saturate: a store
that fills completely and half empties returns $0.5$, which is the statement you
want.
"""))

cells.append(md(r"""
**Reading it.** Three lines carry the result.

- `DEVIATION` — how far the delivered energy sits from the 1 MWe target. Not an
  error: the hardware was specified, not solved for. Positive means the
  latent-only sizing rule overshoots.
- `realised discharge glide / assumed` — the same number. The flow was pinned
  by the assumed glide, so this *is* the deviation.
- `end of discharge  mean eps` — how much of the store fails to refreeze. If it
  is near zero the store is **inventory** limited; if it is well above zero and
  *rises* when you add wells, it is **rate** limited.
"""))

# ---------------------------------------------------------------- plots
cells.append(md(r"""
### 4.3 Profiles along the well

Folded onto depth: a hairpin goes down and comes back, so the developed
coordinate runs to $2L_{\rm well}$ while the depth runs to $L_{\rm well}$ and
back. **Solid = down leg, dotted = return leg.** Because the cascade is laid
along the developed length (H8), the two legs carry different $T_m$ at the same
depth.
"""))

cells.append(code(r"""
ch, dc = css['charge'], css['discharge']
A_avail = case.V_well/(case.num_tubes*case.L_tube)
s = ch['z']; half = len(s)//2
z_dn, z_up = s[:half], case.L_tube - s[half:]
Tm = np.array(css['T_m_seg']) - 273.15

def legs(a):
    a = np.asarray(a); return (z_dn, a[:half]), (z_up, a[half:])
def Tf_seg(a):
    a = np.asarray(a); return 0.5*(a[:-1] + a[1:])

targets = [0.01, 0.1, 0.5, 1.0, 3.0, 10.0]      # hours
def pick(t_s): return [int(np.argmin(np.abs(t_s - th*3600))) for th in targets]
def tlabel(t):
    return f'{t:.0f} s' if t < 60 else (f'{t/60:.0f} min' if t < 3600
                                        else f'{t/3600:.2f} h')

cmap = plt.cm.viridis(np.linspace(0.15, 0.95, len(targets)))
fig, ax = plt.subplots(2, 3, figsize=(14, 8))
for row, (r, ttl) in enumerate(((ch, 'CHARGE'), (dc, 'DISCHARGE'))):
    h = r['history']; ks = pick(r['t_hist'])
    for c_, k in zip(cmap, ks):
        for col, ser in ((0, Tf_seg(h['T_fluid'][k])-273.15),
                         (1, h['A_melt'][k]/A_avail),
                         (2, h['T_pcm'][k]-273.15)):
            (a1, v1), (a2, v2) = legs(ser)
            ax[row][col].plot(v1, a1, color=c_, lw=1.3,
                              label=(tlabel(r['t_hist'][k]) if col == 0 else None))
            ax[row][col].plot(v2, a2, color=c_, lw=1.3, ls=':')
    for col in (0, 2):
        (a1, v1), (a2, v2) = legs(Tm)
        ax[row][col].plot(v1, a1, 'k--', lw=1); ax[row][col].plot(v2, a2, 'k:', lw=1)
    ax[row][0].set_ylabel(f'{ttl}\ndepth [m]')
    ax[row][0].set_xlabel('water [°C]   (dashed: $T_m$)')
    ax[row][1].set_xlabel(r'$\varepsilon_{local}$'); ax[row][1].set_xlim(-0.02, 1.02)
    ax[row][2].set_xlabel(r'$T_{pcm}$ [°C]   (dashed: $T_m$)')
    ax[row][0].legend(fontsize=7, title='time', title_fontsize=7)
    for a in ax[row]: a.invert_yaxis(); a.grid(alpha=.3)
plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
The third column is the one to dwell on. Where the solid curve sits **above**
the dashed $T_m$ the PCM is superheated liquid; below it, subcooled solid. An
area-based model can only ever draw the dashed line.
"""))

cells.append(code(r"""
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
for a, r, ttl in ((ax[0], ch, 'charge'), (ax[1], dc, 'discharge')):
    F = r['history']['A_melt']/A_avail
    im = a.pcolormesh(r['t_hist']/3600, s/2.0, F.T, shading='auto',
                      cmap='inferno', vmin=0, vmax=1)
    a.set_xlabel('time [h]'); a.set_title(ttl + ' (CSS)')
    plt.colorbar(im, ax=a, label=r'$\varepsilon_{local}$')
ax[0].set_ylabel('developed length / 2 [m]'); ax[0].invert_yaxis()
plt.tight_layout(); plt.show()

gap = np.max(np.abs(dc['history']['A_melt'][-1] - ch['history']['A_melt'][0]))
print(f'the cycle closes on itself: max|eps difference| = {gap/A_avail:.2e}')
"""))

cells.append(md(r"""
At cyclic steady state the **right edge of the discharge panel is the left edge
of the charge panel** — the state has returned to itself. That is what CSS
means, and the printed residual is the check.
"""))

# ---------------------------------------------------------------- sweeps
cells.append(md(r"""
---
## 5. Parametric studies

`sweep` runs one parameter and returns a table. The ORC correction is done
**once**, at the design point, and held across the sweep: it is a fairer
comparison — every case is then judged at a common plant operating point — and
it halves the cost. The realised outlet moves by less than a kelvin over most
sweeps, so the two agree anyway.

Start with `record=False` (the default): recording is only needed for profile
plots and roughly doubles the time per point.
""" ))

cells.append(code(r"""
def sweep(field, values, base=None, N=None, correct_once=True,
          n_cycles=None):
    '''Run one parameter over a list of values and collect the results.

    Returns a DataFrame with one row per value. Add whatever else you need to
    the `row` dict below -- everything in the result of simulate_css is
    available.
    '''
    base = base or case
    n_cycles = N_CYCLES if n_cycles is None else n_cycles
    T_2d_fixed = None
    if correct_once:
        T_2d_fixed = simulate_css_corrected(
            base, N=N, n_cycles=n_cycles)['T_2d_realised']
    rows = []
    for v in values:
        c = base.with_(**{field: v})
        bad = [m for lvl, m in validate_case(c, verbose=False) if lvl == 'error']
        if bad:
            print(f'  {field}={v}: SKIPPED -- {bad[0]}');  continue
        r = simulate_css(c, N=N, T_2d=T_2d_fixed, n_cycles=n_cycles)
        kpi_v = performance_indices(c, r)
        if not r['converged']:
            print(f'  {field}={v}: *** NOT CONVERGED *** at n_cycles='
                  f'{n_cycles}, drift {r["history"][-1]["drift"]:.1e}.'
                  f' Raise N_CYCLES; this row is not usable.')
        rows.append({field: v,
                     'converged':    r['converged'],
                     'N_wells':      r['N_wells'],
                     'deviation_%':  100*r['deviation'],
                     'MWe':          r['W_el_out_implied']/1000.0,
                     'glide_dc_K':   r['glide_dc'],
                     'eta_RTE':      kpi_v['eta_RTE'],
                     'UA_tot_kW_K':  kpi_v['UA_total_kW_K'],
                     'UA_orce_kW_K': kpi_v['UA_orce_kW_K'],
                     'pinch_min_K':  kpi_v['pinch_min_K'],
                     'T_2d_real_C':  r['T_2d_realised']-273.15,
                     'eps_end_ch':   r['charge']['eps_local'].mean(),
                     'eps_end_dc':   r['discharge']['eps_local'].mean(),
                     'cycles':       r['cycles']})
        print(f"  {field}={v}: deviation {100*r['deviation']:+.3f} %")
    return pd.DataFrame(rows)


def plot_sweep(df, field, ycols=('deviation_%', 'eps_end_dc')):
    fig, ax = plt.subplots(1, len(ycols), figsize=(5.5*len(ycols), 3.8))
    ax = np.atleast_1d(ax)
    for a, y in zip(ax, ycols):
        a.plot(df[field], df[y], 'o-')
        a.set_xlabel(field); a.set_ylabel(y); a.grid(alpha=.3)
        if y == 'deviation_%': a.axhline(0, color='crimson', ls='--', lw=1)
    plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
### 5.1 Worked example: the discharge-inlet subcooling

$\Delta T_{M,1D}$ is one of the two clean knobs — it moves the discharge inlet
and nothing else — so it is the right place to start.
"""))

cells.append(code(r"""
df = sweep('DT_M_1D', [6.0, 8.0, 10.0, 12.0])
display(df.round(4))
plot_sweep(df, 'DT_M_1D')
"""))

cells.append(md(r"""
Two things to take from this, both of which generalise.

**The residual melt fraction is the diagnostic.** It falls to zero by 15 K: the
store now refreezes completely, so adding driving force stops helping and the
limit has changed from a *rate* limit to an *inventory* limit. Watch
`eps_end_dc` in every sweep — it tells you which lever you are pulling.

**Do not tune to zero deviation.** The curve crosses zero somewhere near
9.3 K. Choosing that value would convert the one honestly reported output of
the procedure back into a satisfied constraint. $\Delta T_{M,1D}=10$ K is used
because it is symmetric with the charging superheat — a principle, not a fit.
"""))

cells.append(md(r"""
### 5.2 A cookbook

| to study… | change | watch | remember |
|---|---|---|---|
| driving force on discharge | `DT_M_1D` | `deviation_%`, `eps_end_dc` | clean knob, moves one inlet |
| driving force on charge | `DT_4C_M` | `eps_end_ch` | clean knob; also moves the whole cascade |
| water flow rate | `DT_3C_2C` | `glide_dc_K`, `N_wells` | **also** moves the cascade spacing and the discharge inlet |
| cascade resolution | `N_lay` | `eps_end_dc` | **also** moves `T_m,bottom`; keep `n_segments` a multiple of it |
| fin geometry | `num_fins`, `fin_L`, `fin_t` | `eps_end_dc`, $U_i$ | rate only; `V_well` barely moves |
| PCM conductivity | `k_s`, `k_l` | `eps_end_dc` | rate only |
| PCM latent heat | `h_m` | `N_wells` | inventory only, $N\propto1/h_m$ |
| PCM sensible capacity | `cp_s`, `cp_l` | the $T_{\rm pcm}$ panel | changes self-levelling, not the energy much |
| window length | `t_ch`, `t_dc` | `deviation_%` | the rate limit is a *rate* limit |
| grid resolution | `n_segments`, `n_times` | everything | convergence check, not a design variable |

A sweep that moves `deviation_%` a lot and `eps_end_dc` barely is changing the
plant's size. One that moves `eps_end_dc` is changing its physics.
"""))

# ---------------------------------------------------------------- limits
cells.append(md(r"""
---
## 6. What this model does not contain

Stated plainly, because a parametric study can wander into the region where one
of these dominates.

1. **No ground heat loss** (H6). The store is adiabatic, so
   $\eta_{\rm storage}\equiv1$ at CSS by construction, and no sweep can produce
   a storage loss. A real field loses to the formation.
2. **H3 is saturated.** The melt fronts of neighbouring tubes are in contact
   over 99 % of the exchanger at the design point, so the annular conduction
   resistance is evaluated at the limit of its validity. Sweeps that *increase*
   melting — more fins, higher $k$, longer charge — go further into that
   region, not out of it.
3. **H9 is unquantified but not small.** The two legs of a hairpin are up to
   48.9 K apart at the wellhead and are assumed perfectly insulated from each
   other. The neglected leak is of order 25 % of the useful duty.
4. **No validation against experiment or a higher-fidelity model.** Energy
   closure is structural — it proves the implementation is consistent, not that
   it is right. This is the largest open item.
5. **Cycle idealisations that remain.** The expansions, compressions and
   pumpings *do* carry isentropic efficiencies as of v0.12 (§3.1), and they act
   on the state points. What is still idealised: no pressure drop anywhere in
   either cycle, no mechanical loss distinct from the isentropic one, and the
   heat-pump discharge is held exactly on the dew line, so the condenser has no
   desuperheating duty.
6. **No volume change on melting**, and the melt-front shape around the fins is
   not modelled.

A deviation of a few per cent sits *inside* the uncertainty of (1)–(3). Report
it as a comparison between cases, not as an absolute.

---

### Where the rest is

| | |
|---|---|
| `P2H2P_verification.ipynb` | fixtures, regression checks, the retired first-cycle sizing path, version history |
| `DESIGN_NOTES.md` | why each superseded approach was superseded, and what the change was worth |
| project report | the full model statement, every equation with the function that evaluates it |
| manuscript | the formulation and its numerics, written for publication |
"""))

# --- guard -------------------------------------------------------------------
# The paper defines macros like \\Aav and \\Elat in its preamble. A notebook has
# no preamble: MathJax renders an unknown macro as a red error, and it is easy
# to paste an equation across from the manuscript without noticing. Fail the
# build rather than ship it.
import re as _re
_MATHJAX_OK = set("""
frac dfrac tfrac begin end cases dcases align aligned text textbf emph mathrm
mathbf mathcal boxed left right bigl bigr Bigl Bigr big Big quad qquad
rho eta varepsilon epsilon delta Delta lambda mu pi sigma tau theta phi omega
alpha beta gamma Gamma Omega Phi Psi Sigma Lambda
psi chi xi nu kappa zeta iota upsilon varphi vartheta varrho varsigma
Theta Xi Upsilon Pi

dot ddot hat bar tilde vec overline underline sqrt sum prod int oint lim
max min inf sup log ln exp sin cos tan
approx equiv propto sim simeq cong neq ne le leq ge geq ll gg
to rightarrow leftarrow Rightarrow leftrightarrow mapsto
in notin subset supset cup cap emptyset forall exists
infty partial nabla cdot cdots ldots dots vdots times pm mp
hline midrule toprule bottomrule
rm it bf sf tt scriptstyle displaystyle limits nolimits
""".split())


def _check_math(cells):
    bad = {}
    for i, c in enumerate(cells):
        if c['cell_type'] != 'markdown':
            continue
        src = ''.join(c['source'])
        for m in _re.findall(r'\\([A-Za-z]+)', src):
            if m not in _MATHJAX_OK:
                bad.setdefault(m, []).append(i)
    if bad:
        for m, where in sorted(bad.items()):
            print(f'  UNKNOWN MACRO \\{m} in cells {where}')
        raise SystemExit('markdown uses macros MathJax will not know; '
                         'write them out or add them to _MATHJAX_OK')


_check_math(cells)

nb = {'cells': cells,
      'metadata': {'kernelspec': {'display_name': 'Python 3',
                                  'language': 'python', 'name': 'python3'},
                   'language_info': {'name': 'python'}},
      'nbformat': 4, 'nbformat_minor': 5}
OUT.write_text(json.dumps(nb, indent=1))
nc = len(cells)
ncode = sum(1 for c in cells if c['cell_type'] == 'code')
print(f'wrote {OUT}')
print(f'  {nc} cells: {ncode} code / {nc-ncode} markdown')
print(f'  model embedded: {len(MODEL.read_text().splitlines())} lines in 1 cell')
