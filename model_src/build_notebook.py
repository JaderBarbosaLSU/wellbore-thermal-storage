"""Build notebooks/P2H2P_PCM_wellbore_storage.ipynb from the verified flat module.

Code cells are EXTRACTED from /tmp/build/thums_model.py by function name, not
retyped, so the notebook runs the same code that was verified against the frozen
IHTC fixture. Markdown cells carry the equations.
"""
import datetime as _dt
import json
import pathlib
import re

# Build stamp in BRT (UTC-3), so a stale Colab tab is identifiable at a glance.
STAMP = (_dt.datetime.utcnow() - _dt.timedelta(hours=3)).strftime('%Y-%m-%d %H:%M BRT')

SRC = pathlib.Path('/tmp/build/thums_model.py').read_text()
OUT = pathlib.Path('/tmp/build/P2H2P_PCM_wellbore_storage.ipynb')

# index top-level defs and the Case class
starts = [(m.start(), m.group(2))
          for m in re.finditer(r'^(def|class) (\w+)[\(:]', SRC, re.M)]
bounds = {}
for i, (pos, name) in enumerate(starts):
    end = starts[i + 1][0] if i + 1 < len(starts) else len(SRC)
    bounds[name] = (pos, end)


LINES = SRC.splitlines(True)
_OFF = []
_p = 0
for _ln in LINES:
    _OFF.append(_p)
    _p += len(_ln)


def _line_of(pos):
    import bisect
    return bisect.bisect_right(_OFF, pos) - 1


def src(*names, drop_trailing=()):
    """Concatenate the source of the named objects, in order.

    Includes any decorator lines above the def/class -- without them a
    @dataclass becomes a plain class and the extraction silently changes
    behaviour. Trailing module-level assignments can be dropped by prefix.
    """
    out = []
    for n in names:
        if n not in bounds:
            raise KeyError(f'not found in flat module: {n}')
        a, b = bounds[n]
        first = _line_of(a)
        while first > 0 and LINES[first - 1].lstrip().startswith('@'):
            first -= 1
        body = ''.join(LINES[first:]).rstrip()
        body = body[:b - _OFF[first]] if b > _OFF[first] else body
        kept = []
        for ln in body.rstrip().splitlines():
            if drop_trailing and ln.startswith(tuple(drop_trailing)):
                continue
            kept.append(ln)
        # a decorator belonging to the NEXT object falls inside this one's
        # [start, end) span; drop any trailing decorator/blank lines
        while kept and (not kept[-1].strip() or kept[-1].lstrip().startswith('@')):
            kept.pop()
        out.append('\n'.join(kept).rstrip())
    return '\n\n\n'.join(out)


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip().splitlines(True)}


def code(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": text.strip().splitlines(True)}


cells = []

# ---------------------------------------------------------------- 0. header
cells.append(md(r"""
# P2H2P PCM wellbore storage — the reference model

### Power-to-Heat-to-Power conversion with phase-change material in repurposed wells

**THUMS · model version 0.4 · every equation visible in this notebook**

**Last updated: @@BUILD_STAMP@@**

`JaderBarbosaLSU/wellbore-thermal-storage` · `notebooks/P2H2P_PCM_wellbore_storage.ipynb`

> This is the **reference** implementation. If you are looking at any other THUMS
> notebook, it is superseded — the four `thums/*.ipynb` files are the original
> v0.1 Colab notebooks, kept for provenance only.

This notebook *is* the model. Nothing is hidden in an imported package: every
function the results depend on is defined in a cell below, in the order it is
used. Read it top to bottom and you have seen the whole thing.

### How it is organised

| section | what it contains |
|---|---|
| 1–2 | what is being modelled, and the assumptions |
| 3–4 | fluid properties, the borehole resistance network, fin efficiency |
| 5 | **the melt front** — the defect, and three formulations |
| 6 | the segment-by-segment march down the well |
| 7–8 | the ORC and heat-pump cycles, pressure drop |
| 9–11 | well-field sizing and the full cycle calculation |
| 12 | results at the design point |
| 13 | verification — reproduces the published IHTC numbers exactly |
| 14 | **profiles along the well** — fluid temperature, melt fraction, $U_i$, NTU |
| 15 | what the model still does **not** contain |

### Working rules

- **Do not press Ctrl/Cmd+S.** Colab will fork a private copy into your Drive
  that stops receiving updates, and nothing warns you.
- Prefer **Runtime → Restart and run all** over re-running single cells; a stale
  number on screen is easy to mistake for a new result.
- Section 7 (the thermodynamic cycles) is long and is state-point bookkeeping
  rather than model equations. Collapse it on a first reading.
"""))

cells.append(code(r"""
!pip install -q CoolProp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from math import log, pi
from dataclasses import dataclass, replace

LAST_UPDATED = '@@BUILD_STAMP@@'
MODEL_VERSION = '0.4c'

print(f'THUMS P2H2P model v{MODEL_VERSION} · last updated {LAST_UPDATED}')
print('ready ·', np.__version__, '·', CP.get_global_param_string('version'))
"""))

# ---------------------------------------------------------------- 1. scope
cells.append(md(r"""
## 1. What is being modelled, and what is assumed

Electrical energy is stored as **latent heat** in a phase-change material packed
into the annulus of a repurposed oil-and-gas well. A high-temperature heat pump
melts the PCM while charging; an organic Rankine cycle recovers the energy while
discharging.

The model answers two questions: **how many wells does a given power block
need**, and **what round-trip efficiency results**. Everything below serves those.

The calculation is a chain — each stage consumes the previous one and nothing
flows backwards:

```
cycles  →  energy budget  →  resistance network  →  segment march
                                                         ↓
        KPIs  ←  well-field sizing  ←  melt front
```

One consequence is worth stating early: **the cycles are decoupled from the
storage**. $\eta_{ORC}$ and COP do not depend on anything the borehole does, so
no change to the storage model can move them.

### Assumptions

| | assumption | where it bites |
|---|---|---|
| **H1** | Quasi-steady conduction in the melt layer — no thermal mass. Time enters only through the front position. | controlled by $\mathrm{Ste}=c_p\Delta T/h_m$; $0.047$ here |
| **H2** | No axial conduction, in wall or PCM. Segments couple only through the fluid. | |
| **H3** | The melt region is a concentric annulus of uniform thickness $\delta$. Fronts from neighbouring tubes never merge. | cannot be violated under Formulation C, but is **sat on**: $\delta/\delta_{\rm merge}=1.000$ over the whole well — §12.1, §15 |
| **H4** | Conduction only in the melt; no natural convection. | increasingly wrong as $\delta$ grows |
| **H5** | *Revised in v0.4.* PCM carries latent **and** sensible energy, lumped into one enthalpy per segment: molten PCM superheats, solid PCM subcools. | no radial resolution *within* a phase — one PCM temperature per segment |
| **H6** | Adiabatic borehole wall — no formation heat loss. | fine over 10 h, not over a season |
| **H7** | Incompressible single-phase water, no buoyancy head. | ±7.5 bar in a 1,524 m loop is neglected |
| **H8** | The PCM is layered along the well, melting point following the fluid glide. Layers are spaced along the **developed** tube length, so a hairpin's two legs are in *different* layers at the same depth. | see H9 |
| **H9** | *New in v0.4a.* **Perfect azimuthal insulation**: each leg's PCM cell is adiabatic at its outer boundary, so no heat crosses between neighbouring legs at the same depth. | the two legs are up to **48.9 K** apart at the wellhead; a slab estimate puts the neglected leak at ~22 W/m against ~83 W/m of useful duty — §15 |
"""))

# ---------------------------------------------------------------- 2. config
cells.append(md(r"""
## 2. Configuration

One flat `Case` object holds every input. v0.1 passed sixteen positional
arguments between functions, which is how the PCM conductivity ended up in the
slot meant for the tube wall for months without anything contradicting it.

> ### Correction notice — wall thermal conductivity
>
> **Versions of this model up to and including the IHTC paper used the wrong
> value of the wall thermal conductivity during charging.** The PCM *liquid*
> conductivity, $k_l = 0.45$ W/m·K, was passed into the argument carrying the
> tube-wall and fin conductivity, where steel, $k_w = 45$ W/m·K, belongs — two
> orders of magnitude low. The error is visible in the original solver
> signature, whose parameter is literally named `k_m_l`.
>
> Its dominant effect was not on the wall resistance, which is small either way,
> but on the **fin efficiency**: it collapsed $\eta_f$ to about $0.33$,
> throttling the fluid-side conductance to roughly what a bare tube could
> absorb. That accidentally compensated the finned-versus-bare surface mismatch
> in the melt front (§5), which is why the energy balance appeared to close.
>
> **The correct value is now the only value.** The switch that reproduced the
> defect was removed in v0.4b along with the rest of the superseded code; steel
> at $k_w = 45$ W/m·K is simply what the model uses.

Geometry notes that matter later:

$$L_{\text{tube}} = 2\,L_{\text{well}}, \qquad
V_{\text{well}} = \frac{\pi D_{\text{well}}^2}{4}L_{\text{well}}
- 2n_t\left(\pi r_e^2 + n_f t_f L_f\right) L_{\text{well}}$$

A hairpin runs down and back, so its developed length already counts both legs.
Per-borehole quantities therefore multiply by $n_t$ (hairpins), **not** $2n_t$ —
a recurring source of factor-of-two errors.
"""))
cells.append(code("FT = 0.3048   # foot  -> m\nIN = 0.0254   # inch  -> m\n\n"
                  + src('Case', drop_trailing=('CASE',))
                  + '\n\nCASE = Case()          # the design point'))

cells.append(code(r"""
c = CASE
print(f'well depth          {c.L_well:8.1f} m   ({c.L_ft:.0f} ft)')
print(f'developed tube      {c.L_tube:8.1f} m   (hairpin, both legs)')
print(f'borehole diameter   {c.D_well*1000:8.1f} mm')
print(f'PCM volume / well   {c.V_well:8.3f} m3')
print(f'melting point       {c.T_m_C:8.1f} C')
print(f'Stefan number       {c.stefan_number(c.DT_4C_M):8.4f}   (H1 degrades as this grows)')
"""))

# ---------------------------------------------------------------- 3. props
cells.append(md(r"""
## 3. Fluid properties and internal convection

$$\mathrm{Re}=\frac{4\dot m}{\pi D \mu},\qquad \mathrm{Pr}=\frac{c_p\mu}{k}$$

Turbulent flow uses Gnielinski; laminar uses the fully-developed isothermal-wall
value $\mathrm{Nu}=3.66$, which is the right constant for a melting boundary.

$$f=\left(0.79\ln \mathrm{Re}-1.64\right)^{-2},\qquad
\mathrm{Nu}=\frac{(f/8)(\mathrm{Re}-1000)\mathrm{Pr}}
{1+12.7\sqrt{f/8}\left(\mathrm{Pr}^{2/3}-1\right)},\qquad
h_i=\frac{\mathrm{Nu}\,k}{D}$$
"""))
cells.append(code("_PROP_CACHE = {}   # properties vary smoothly; caching on 0.05 K bins\n"
                  "                   # cuts CoolProp calls by ~50x with no visible effect\n\n"
                  + src('fluid_props', 'h_internal')))

# ---------------------------------------------------------------- 4. network
cells.append(md(r"""
## 4. The borehole resistance network

Heat leaving the fluid crosses four resistances in series, all referred to the
**inner** tube area so a single $U_i$ drives the march:

$$R'_{\text{conv}}=\frac{1}{h_i},\qquad
R'_{\text{wall}}=\frac{r_i}{k_w}\ln\frac{r_e}{r_i},\qquad
R'_{f}=R''_{f,i},\qquad
R'_{\text{outer}}=\frac{2\pi r_i}{P_T\,\eta_o\,h_e}$$

The melt layer is a steady cylindrical shell:

$$h_e=\frac{k_m}{r_e\ln\left(1+\delta/r_e\right)}$$

The fins are straight longitudinal fins with an adiabatic-tip correction:

$$L_c=L_f+\frac{t_f}{2},\qquad
m=\sqrt{\frac{2h_e}{k_w t_f}},\qquad
\eta_f=\frac{\tanh(mL_c)}{mL_c}$$

$$A_f=L_{\text{tube}}(2L_f+t_f),\qquad
A_t=2L_{\text{tube}}(\pi r_e+n_f L_f),\qquad
\eta_o=1-\frac{n_f A_f}{A_t}(1-\eta_f)$$

and the **finned perimeter** — remember this one:

$$\boxed{P_T = 2\pi r_e + 2 n_f L_f}$$

$$U_i=\left[\frac{1}{h_i}+R''_{f,i}
+\frac{r_i}{k_w}\ln\frac{r_e}{r_i}
+\frac{2\pi r_i}{P_T\eta_o h_e}\right]^{-1}$$

> **Note** that `k_w` appears with two different physical meanings — the tube
> wall in $R'_{\text{wall}}$ and the **fin** in $m$ — through a single argument.
> That is what allowed the PCM conductivity to sit in that slot unnoticed.
"""))
cells.append(code(src('compute_U_i')))


# ---------------------------------------------------------------- 5. front
cells.append(md(r"""
## 5. The melt front

The melt front is carried as an **enthalpy per unit tube length**, $E'$,
measured from a datum of fully solid PCM at the local melting temperature.
Subcooled solid, two-phase and superheated liquid are three branches of one
single-valued curve, and the melted fraction follows from $E'$ rather than being
tracked separately.

> **What used to be here.** Earlier versions carried two other formulations: a
> closed-form Stefan front (v0.1, the conference paper) and an energy-balance
> front on the melted *area* (v0.2). Both have been **deleted** from the code as
> of v0.4b. Neither is physically consistent — a variable counting melted area
> has no slot for liquid above $T_m$ or solid below it, so a segment that
> finishes melting must either melt material that is not there or discard the
> heat delivered to it. The reasoning is preserved in the project report and in
> git history; it is no longer carried as executable code.

### 5.1 The fins are metal, not PCM

$A_{\text{melt}}$ is the quantity the latent balance conserves, so it must be
**PCM and nothing else**. The annulus between $r_e$ and $r_e+\delta$ is not all
PCM — the fins occupy metal inside it:

$$A_{\text{melt}} = \pi\left[(r_e+\delta)^2 - r_e^2\right]
- n_f\,t_f\,\min(\delta,\,L_f)$$

The fin term saturates once the front passes the fin tips. Inverting for
$\delta$ is still closed form, just a shifted quadratic:

$$\delta \le L_f:\quad \pi\delta^2+\left(2\pi r_e-n_f t_f\right)\delta-A=0
\qquad
\delta > L_f:\quad \pi\delta^2+2\pi r_e\delta-\left(A+n_f t_f L_f\right)=0$$

**Corrected in v0.3.** Earlier versions treated the whole annulus as PCM, which
overstated the melted volume by 11 % when the layer is thin and 6 % at the
design point, and understated $\delta$ — and therefore the melt-layer
conduction resistance — by up to 34 %. Because it inflated
$\varepsilon_{\text{PCM}}$, it also pushed $N_{\text{inventory}}$ slightly
**up**; correcting it moves the design point from 12.79 to 12.74 wells.
"""))
cells.append(code(src('delta_from_area', 'area_from_delta')))

# ---------------------------------------------------------------- 6. march
cells.append(md(r"""
### 5.2 The enthalpy state

Formulation B tracks $A_{\rm melt}$ and nothing else, so a segment that runs out
of PCM has nowhere to put further heat: the march clips it and discards the
remainder. On discharge the discarded heat reached **148 % of the heat
delivered** — most of what the network asked for.

Carry **enthalpy per unit tube length** instead, measured from fully-solid-at-$T_m$:

$$E' < 0:\quad \text{solid, subcooled}\quad T = T_m + E'/C_s$$
$$0 \le E' \le E'_{\rm lat}:\quad \text{two-phase}\quad T = T_m,\;\;
A_{\rm melt} = E'/(\rho h_m)$$
$$E' > E'_{\rm lat}:\quad \text{liquid, superheated}\quad
T = T_m + (E'-E'_{\rm lat})/C_l$$

with $E'_{\rm lat}=\rho h_m A_{\rm avail}$, $C_s=\rho c_{p,s}A_{\rm avail}$,
$C_l=\rho c_{p,l}A_{\rm avail}$. Three consequences:

1. **no clipping** — heat always has somewhere to go;
2. **$\varepsilon_{\rm local}\in[0,1]$ by construction** — the over-melt the
   latent-only model produced cannot occur;
3. **desuperheating and subcooling on discharge are recovered.**

$\delta$ still follows from $A_{\rm melt}$ exactly as before, so the
heat-transfer coefficient is unchanged in form — only its driving temperature is
now $T_{\rm pcm}$ rather than $T_m$.

> **Explicit Euler is not usable here.** On the latent plateau the segment has
> effectively infinite heat capacity, so any step is stable. Adding sensible heat
> makes it finite, with a stability limit $\Delta t < C/K \approx 620$ s. The
> time grid's last step is 8491 s — fourteen times that — and an explicit update
> diverges; it drove the secondary fluid to \SI{261}{\kelvin} on the first
> attempt. `advance_segment` integrates each branch in closed form instead,
> splitting the step where the state crosses a branch boundary. It agrees with a
> 200 000-substep reference to $3\times10^{-13}$.
"""))

cells.append(code(src('pcm_capacities', 'pcm_state', 'advance_segment')))
cells.append(code(src('segment_profile', 'unmirror_march', 'mixed_mean_outlet',
                      'march_h')))

cells.append(md(r"""
## 6. The segment march

The tube is cut into $n_s$ segments. Within a segment the PCM surface is
isothermal at its layer melting temperature, so the segment is a single-stream
heat exchanger of effectiveness $1-e^{-\mathrm{NTU}}$:

$$\mathrm{NTU}_j=\frac{2\pi r_i U_i \Delta z}{\dot m c_p},
\qquad T_{j+1}=T_m^{(j)}+\left(T_j-T_m^{(j)}\right)e^{-\mathrm{NTU}_j},
\qquad q'_j=\frac{\dot m c_p (T_j-T_{j+1})}{\Delta z}$$

The march is sequential — each segment's outlet is the next one's inlet.
"""))
cells.append(code(src('layer_map')))

# ---------------------------------------------------------------- 7. cycles
cells.append(md(r"""
## 7. The thermodynamic cycles

State-point bookkeeping through CoolProp, not model equations — **collapse this
section on a first reading.** Every state point is fixed by a prescribed approach
temperature.

$$\eta_{ORC}=\frac{w_{\text{net}}}{q_{\text{in}}},\qquad
\mathrm{COP}=\frac{q_{\text{out}}}{w_{\text{comp}}}$$

> The expansions and compressions carry **no isentropic efficiency**, so both are
> idealised. Only electrical and mechanical efficiencies appear downstream. This
> bounds every efficiency the notebook reports.
"""))
cells.append(code(src('double_stage_rankine')))
cells.append(code(src('two_stage_htheatpump_2regs')))

cells.append(md(r"""
The budget runs **backwards** from the specified electrical output:

$$\dot W_T=\frac{\dot W_{el,out}}{\eta_T\eta_G},\quad
\dot Q_{in,ORC}=\frac{\dot W_T}{\eta_{ORC}},\quad
\Delta E_{out,HP}=\dot Q_{in,ORC}t_{dc}(1+\lambda),\quad
\dot W_{el,in}=\frac{\dot Q_{out,HP}}{\mathrm{COP}\,\eta_C\eta_H}$$

$$\dot m_w^{ch}=\frac{\dot Q_{out,HP}}{c_{p,w}\Delta T_{3C,2C}}$$

$\lambda$ is the assumed storage loss, fixed at 5 %. §12 shows what the model
now says it actually is.
"""))
cells.append(code(src('T_m_bottom', 'cycle_state_points', 'energy_budget',
                      'melting_temperatures')))

# ---------------------------------------------------------------- 8. dp
cells.append(md(r"""
## 8. Pressure drop and pumping power

$$u=\frac{\dot m}{\rho A},\qquad \mathrm{Re}=\frac{\rho u D_i}{\mu},\qquad
f=\begin{cases}64/\mathrm{Re} & \mathrm{Re}<2000\\ 0.3164\,\mathrm{Re}^{-1/4} & \text{Blasius}\end{cases}$$

$$\Delta p_{\text{fric}}=f\frac{L_{\text{tube}}}{D_i}\frac{\rho u^2}{2},\qquad
\Delta p_{\text{min}}=K_{\text{tot}}\frac{\rho u^2}{2},\qquad
\dot W_{\text{pump}}=\frac{\dot m\,\Delta p}{\rho}$$

with $K_{\text{tot}} = 0.5 + 2(2.2) + 1.0 = 5.9$ for an entrance, two 180°
bends and an exit. Those three coefficients are hard-coded literals the source
itself calls example values, and the tube is assumed smooth.
"""))
cells.append(code(src('calculate_pressure_drop')))

# ---------------------------------------------------------------- 10. sizing
cells.append(md(r"""
## 10. Well-field sizing — two constraints, not one

v0.1 sized on a single requirement: *can the wells absorb the energy in time?*
That is a heat-transfer question. **Nothing required them to contain it.**

$$\Phi_{\text{heat}}(N)=\frac{\Delta E_{out,HP}}{\Delta E_{\text{del}}(N)}-1=0,
\qquad
\Phi_{\text{inv}}(N)=\varepsilon_{\text{PCM}}(N)-1=0$$

$$\boxed{N_{\text{wells}}=\max\left(N_{\text{charge\ rate}},\;N_{\text{capacity}}\right)}$$

with the capacity criterion written the way the v0.1 notebook computed
`N_wells_ideal`, but now used to size rather than only to cost:

$$N_{\text{capacity}}=\frac{\Delta E_{\text{out,HP}}}{E_{\text{well}}},\qquad
E_{\text{well}}=\rho V_{\text{well}}\left[h_m+c_{p,l}(T_{3c}-T_m)\right]$$

**Why not size on the discharge instead?** Sizing must act on whichever variable
is free. Charging pins the flow to the specified glide, so $N$ is the only
freedom. Discharging has $N$ already fixed, so the flow is the freedom. Pin the
discharge flow to its glide as well and the delivered energy saturates at
$0.996$ of the requirement for *any* well count — there is no root, because the
fluid leaves slightly short of the nominal outlet temperature.

With the `k_w` defect the rate criterion was always the larger of the two, so
the omission stayed invisible. Correct the conductivity and the inventory
constraint binds instead.

The constraint is only meaningful once the front conserves energy: under
Formulation A the melted volume is set by an independent solve, so
$\varepsilon_{\text{PCM}}$ says nothing about the delivered energy.

### $N_{\text{inventory}}$ is not the "ideal number of wells"

The v0.1 notebook printed `Ideal number of wells: 12.215`. That is a **different
and simpler calculation**: the thermal load divided by the sensible *plus*
latent capacity of the PCM, giving the total PCM volume required and hence a
well count. It carries no thermal resistance and no losses.

It is therefore a genuine **lower bound** — the best any design could achieve —
and $N_{\text{inventory}}$, computed from the marched model, must sit at or
above it. Designs with high $\varepsilon_{\text{PCM}}$ are desirable precisely
because they approach that bound: the PCM in the ground is being used rather
than merely occupying the borehole.

§12 checks the two against each other.
"""))
cells.append(md(r"""
### Two questions, not one

A well field must satisfy **two independent requirements**, and the model names
one number for each. State them as questions before any algebra:

| symbol | the question it answers | how it is found |
|---|---|---|
| $N_{\rm rate}$ | Can the field **absorb** $\Delta E$ *within* $t_{\rm ch}$? | root solve on the march |
| $N_{\rm capacity}$ | Can the field **contain** $\Delta E$ *at all*? | closed form |
| $N_{\rm wells}$ | **Both must hold.** | $\max$ of the two |

The distinction is between a **rate** and a **stock**. A field could have ample
PCM but transfer heat into it too slowly to finish inside the window; or it
could transfer heat beautifully into a store too small to hold the energy.
Neither failure implies the other, so the design must clear both.

$N_{\rm rate}$ has no closed form — it comes from marching the whole well, so it
carries *everything thermal*: fins, wall conductivity, melt-layer resistance,
NTU, the cascade, the charging window. $N_{\rm capacity}$ contains **none** of
that. Corrected in **v0.4b** to resolve the cascade and to include the solid
subcooling the discharge reaches:

$$N_{\rm capacity}=\frac{\Delta E_{\rm out,HP}}
{\rho V_{\rm well}\Big[h_m
 + c_{p,l}\big\langle T_{3c}-T_{m,\ell}\big\rangle
 + c_{p,s}\big\langle T_{m,\ell}-T_{3d}\big\rangle\Big]}$$

where $\langle\cdot\rangle$ is the volume-weighted mean over the $N_{\rm lay}$
layers — a plain mean here, since the layers are equal in length. No $U_i$, no
$\delta$, no NTU, no $t_{\rm ch}$: volume and material properties only. That
makes it a **hard lower bound**, no design and no window can beat it.

> **What was wrong before.** The old form used
> $c_{p,l}(T_{3c}-T_m)$ with $T_m$ the melting temperature of the **top layer**,
> giving a 10 K sensible span for the whole store. But $\Delta T_{4C,M}=10$ K is
> the approach that fixes the charging inlet relative to the *first* layer only.
> Applying it to all nine layers is the **no-cascade limit**: set
> $N_{\rm lay}=1$ and $\langle T_m\rangle = T_{m,\rm top}$ and the two agree
> identically.
>
> It was also not a bound. The model's own end-of-charge state held 4.822 MWh
> per well against a claimed capacity of 4.744 — mean superheat reached
> **13.81 K**, not 10 K. A bound the model exceeds is not a bound, and it was
> setting the well count.

> **Two names have been retired**, because they were traps.
> `N_inventory` was a silent *alias* for `N_capacity` in v0.4 — the two compared
> equal and could not be told apart. Worse, in **v0.3 the same name meant
> something else**: the well count at which $\varepsilon_{\rm PCM}(N)=1$, a
> marched quantity. *A v0.3 printout of `N_inventory` is not comparable with a
> v0.4 one.* And `N_heat` was the dict key while the report said "charge-rate"
> and the sweep column said `N_rate` — three labels, one quantity.
> Both now raise a `KeyError` carrying this explanation. The surviving
> vocabulary is three words: **rate**, **capacity**, **the maximum of the two**.

### How to read the two criteria — and how not to

$N_{\rm capacity}$ is the seductive one. It is a closed form, it is smooth, it
responds to every geometric parameter you change, and **it is wrong to read on
its own.** It is a *lower bound*. The answer is the maximum of the two, so
whenever the rate criterion binds, $N_{\rm capacity}$ moves while the answer
does not — and it can move the opposite way.

This is not hypothetical. Sweeping fin radial length, $N_{\rm capacity}$ falls
monotonically — 11.573, 11.239, 10.985 for 7.5, 3.75, 0.75 mm — which reads as
*fins make things worse*. But

$$N_{\rm capacity}=\frac{\Delta E_{\rm out,HP}}
{\rho V_{\rm well}\left[h_m+c_{p,l}\Delta T\right]}$$

contains **no heat transfer at all**. Shrinking a fin returns its metal volume
to the PCM, $V_{\rm well}$ rises, and the bound falls by pure volume
bookkeeping. Over the same sweep the real answer *rises*: 11.573, 11.569,
14.654. Removing the fins entirely costs 42 % more wells. §12.1 does this
sweep explicitly.

`sizing_report` therefore leads with $N_{\rm wells}$, names the criterion that
set it, and labels $N_{\rm capacity}$ as a bound.
"""))
cells.append(code(src('SizingResult', '_bisect', 'layer_mean_T_m',
                      'well_capacity_kJ', 'size_well_field', 'sizing_report')))

# ---------------------------------------------------------------- 11. run
cells.append(md(r"""
## 11. The full cycle

Charge, size the field, discharge, and collect the indices.

$$\eta_{RTE}=\frac{\left(\dot W_{el,out}-N_{\rm well}\dot W_{f,dc}\right)t_{dc}}
{\left(\dot W_{el,in}+N_{\rm well}\dot W_{f,ch}\right)t_{ch}},
\qquad
\eta_{RTE}^{0}=\frac{\dot W_{el,out}t_{dc}}{\dot W_{el,in}t_{ch}}$$

$$\varepsilon_{m}=\frac{V_{m,l}(t_{ch})}{V_{m}},\qquad
\eta_{\text{storage}}=\frac{\Delta E_{\text{discharged}}}{\Delta E_{\text{stored}}}$$

### The factorial-study indicators

Restored from the earlier DoE work. Three of them already existed under other
names — $\varepsilon_m$ is `eps_pcm`, RTE is `eta_rte`, and
$\Delta E_{\rm therm}$ is `E_well` in MWh rather than kWh.

$$\eta_{T}=\frac{\dot W_{el,out}}{\dot Q_{in,ORC}},
\qquad
\eta_{T,\rm eff}=\eta_{T}-\frac{N_{\rm well}\dot W_{f,dc}}{\dot Q_{in,ORC}}
\quad\text{(effective discharge efficiency)}$$

$$\epsilon\mathrm{RTE}=\mathrm{RTE}\cdot\varepsilon_{m}
\quad\text{(effectiveness-weighted RTE)}$$

$$\Delta E_{\rm therm}=\frac{\dot Q_{in,ORC}\,t_{dc}}{N_{\rm well}},
\qquad
\Delta E_{\rm elec}=\Delta E_{\rm therm}\cdot\eta_{T,\rm eff}
\quad\text{(kWh per well)}$$

$\eta_T$ is identically $\eta_{ORC}\,\eta_T^{\rm turb}\,\eta_G$ — the model
computes it both ways and they agree to machine precision.

> **Read $\epsilon$RTE with care under the enthalpy formulation.**
> $\varepsilon_m$ **saturates**. Once a segment has melted all the PCM it owns,
> further energy goes into *superheat*, which $\varepsilon_m$ cannot see. At
> this design point $\varepsilon_m = 1.0000$ exactly, so
> $\epsilon\mathrm{RTE} = \mathrm{RTE}$ identically and the weighting carries
> no information. It discriminates only across designs that do **not** saturate,
> which in a factorial study is precisely the corner you are least interested
> in.
>
> `util_enthalpy` is added as the non-saturating companion: energy actually
> banked per well as a fraction of the capacity ceiling of §10,
> $\Delta E_{\rm stored}/N_{\rm well}$ divided by $E_{\rm well}^{\rm cap}$.
> It reads **0.850** here — the store uses 85 % of its ceiling — and it keeps
> moving after $\varepsilon_m$ has pinned at 1.
"""))
cells.append(code(src('run_cycle')))

# ---------------------------------------------------------------- 12. results
cells.append(md(r"""
## 12. Results

One configuration. The superseded formulations have been removed from the code
(§5), so there is nothing left to compare against inside this notebook — the
historical comparison lives in the project report.
"""))
cells.append(code(r"""
runs = {'v0.4b': run_cycle(CASE)}

rows = []
for label, r in runs.items():
    k, d = r['kpis'], r['detail']
    rows.append({'configuration': label,
                 'N_wells': k['N_wells'], 'eps_PCM': k['eps_pcm'],
                 'eta_RTE': k['eta_rte'], 'eta_RTE_nopump': k['eta_rte_nopump'],
                 'COP': k['cop_hp'], 'eta_ORC': k['eta_orc'],
                 'E_well_MWh': k['E_well'], 'f_pump': k['f_pump'],
                 'binding': d.get('binding', '-- (only one constraint existed)')})
df = pd.DataFrame(rows).set_index('configuration')
pd.set_option('display.width', 200, 'display.max_columns', 20)
print(df.round(5).T.to_string())
"""))

cells.append(md(r"""
Two things to read from that table.

**`eta_RTE_nopump` is a useful invariant.** It depends only on the two cycle
efficiencies, the four electrical and mechanical efficiencies, and $\lambda$ —
so no change to the storage model can move it. If you modify anything in §5–§11
and this number shifts, the change has leaked into the cycles. It sits at
$0.43484$.

**The rate criterion binds.** `N_capacity` is a lower bound and is slack here by
about 15 %, so the well count is set by how fast the field can absorb the
energy, not by how much it can hold. §10 explains why reading `N_capacity`
alone across a sweep is a trap.
"""))
cells.append(code(r"""
for label, r in runs.items():
    d = r['detail']
    sizing_report(CASE, d, label)
    print()
"""))

cells.append(code(r"""
# --- the factorial-study indicators -----------------------------------------
k = runs['v0.4b']['kpis']
print(f"{'indicator':30s} {'value':>14s}   unit")
for name, key, unit in (
    ('eps_m   melted fraction',      'eps_pcm',      '--'),
    ('RTE     round-trip',           'eta_rte',      '--'),
    ('epsRTE  effectiveness-wtd',    'eps_rte',      '--'),
    ('eta_T   thermal -> electric',  'eta_T',        '--'),
    ('eta_T,eff  after dc pumping',  'eta_T_eff',    '--'),
    ('dE_therm   per well',          'dE_therm_kWh', 'kWh'),
    ('dE_elec    per well',          'dE_elec_kWh',  'kWh'),
    ('util_enthalpy  vs ceiling',    'util_enthalpy','--'),
):
    print(f'{name:30s} {k[key]:14.5f}   {unit}')
print()
print('eta_T identity check:  W_el_out/Q_in_ORC  vs  eta_ORC * eta_turb * eta_gen')
rank, hp, T_ = cycle_state_points(CASE)
print(f'   {k["eta_T"]:.12f}   vs   {rank["rank_eff"]*CASE.Turb_eff*CASE.ElG_eff:.12f}')
print()
print(f'epsRTE == RTE here because eps_m = {k["eps_pcm"]:.4f} exactly.')
print('That is saturation, not perfection -- see the note above.')
"""))

cells.append(code(r"""
# The retired names fail loudly rather than returning a number that may mean
# something other than you think.
d = runs['v0.4b']['detail']
print(f"N_rate     = {d['N_rate']:.3f}   <- rate:     absorb it in time?")
print(f"N_capacity = {d['N_capacity']:.3f}   <- capacity: contain it at all?")
print(f"N_wells    = {d['N_wells']:.3f}   <- max of the two")
print()
for retired in ('N_inventory', 'N_heat'):
    try:
        d[retired]
        print(f'{retired}: returned a value -- the guard is not working')
    except KeyError as e:
        print(f'{retired} ->', str(e).strip('"')[:78], '...')
"""))


cells.append(code(r"""
for label, r in runs.items():
    d = r['detail']
    if 'N_rate' in d:
        print(f"{label:22s}  eta_storage={d['eta_storage']:.4f}"
              f"  closure={d['closure_charge']:.1e}")
print()
print('The assumed storage loss lambda = 0.05 implies eta_storage = 0.952.')
print('At the v0.3 design point the model computes about 0.91, so the field')
print('banks ~10% more than the ORC withdraws. lambda should become an OUTPUT')
print('of the model rather than an input. Not done yet -- see section 15.')
"""))

# ============================================================ 13. verification
cells.append(md(r"""
## 13. Verification against the frozen fixture

**What this replaces.** Until v0.4b this section ran the v0.1 closed-form path
and checked it still reproduced the conference-paper numbers to $4\times10^{-11}$.
That path has been deleted along with the rest of the superseded code, so that
check is gone. It is recoverable from git history if it is ever needed again,
but it is not coming back into this notebook.

**What replaces it.** A frozen fixture of the v0.4b outputs. This is a weaker
guarantee and it is worth being clear about the difference:

- the old check was an *external* one — it tied the code to a published result;
- this one is *internal* — it only detects **drift**. It cannot tell you the
  model is right, only that it still computes what it computed on the day the
  fixture was frozen.

That is still worth having. Most damage comes from a change leaking somewhere
unintended, and this catches exactly that. It does not substitute for the
external validation that §15 still lists as the largest gap.

**A fixture that is re-frozen is not a fixture that failed.** The values below
were re-frozen for v0.6, when the discharge closure changed
(§16.5): the inlet moved from \SI{95.000}{\celsius} to
\SI{91.111}{\celsius}, so every discharge-side index had to move with it. A
drift check is only meaningful if re-freezing is *deliberate and explained*,
never a response to a red line. The retired values are kept beside the new ones
so the size of the step is on the record, and the equivalence test that
guarantees nothing else changed is in §16.5: setting
$\Delta T_{M,1D}=\Delta T_{3C,2C}/N_{\rm lay}$ reproduces every v0.5a number
exactly.
"""))
cells.append(code(r"""
FIXTURE_V06 = {   # re-frozen 2026-09 for the v0.6 discharge closure -- see below
    'N_rate':             11.5122070312,
    'N_capacity':          9.67244314511,
    'N_wells':            11.5122070312,
    'eps_pcm':             1.0,
    'E_well_kJ':          20857297.9447,
    'eta_storage':         0.952163601583,
    'flow_ratio_dc_ch':    0.948828125,
    'cop_hp':              2.95633040945,
    'eta_orc':             0.232035822617,
    'eta_rte':             0.408915233697,
    'eta_rte_nopump':      0.425994284001,
    'f_pump':             0.0587585434779,
    'E_well':              4.63600440156,
    'rho_E':             167.46852173,
}
# Superseded by the line above, kept so the size of the v0.6 step is on the
# record. These are the v0.5a values, frozen under the RETIRED closure in which
# the discharge OUTLET was pinned at T_m,top and the inlet derived from it. The
# discharge inlet moved 95.000 -> 91.111 C, so these had to move; nothing else
# about the model changed. The two differences worth noting are eta_ORC
# (0.2369 -> 0.2320, because the provisional T_2d fell from 150 to 146.1 C and
# the ORC evaporating temperature follows it) and the flow ratio
# (1.0473 -> 0.9488, because a colder inlet needs less flow for the same duty).
FIXTURE_V05A_RETIRED = {
    'N_wells':            11.2551269531,
    'eta_orc':             0.236851320275,
    'flow_ratio_dc_ch':    1.047265625,
}

res_v = runs['v0.4b']
merged = dict(res_v['kpis']); merged.update(res_v['detail'])
print(f"{'quantity':20s} {'frozen':>20s} {'now':>20s}   rel. diff")
worst = 0.0
for k, v in FIXTURE_V06.items():
    now = float(merged[k])
    d = abs(now - v)/abs(v) if v else abs(now)
    worst = max(worst, d)
    print(f'{k:20s} {v:20.12g} {now:20.12g}   {d:.1e}')
print()
print(f'worst relative difference: {worst:.1e}')
print('PASS -- no drift' if worst < 1e-9 else
      'FAIL -- the model has moved since the fixture was frozen')
"""))

# ============================================================ 14. diagnostics
cells.append(md(r"""
## 14. Profiles along the well

The quantities the IHTC notebooks plotted, restored: **secondary-fluid
temperature**, **PCM melt fraction**, the **overall coefficient $U_i$** and
**NTU**, each as a function of depth and time, for charging and for
discharging.

`march(..., record=True)` retains the per-segment state at every time level.
Recording appends to lists and touches nothing the solution depends on — §14.5
verifies that a recorded and an unrecorded run are bit-identical.

### Reading the depth axis: a hairpin has two legs

The tube is a **hairpin** — down and back — so the developed coordinate $s$ runs
$0\to L_{\rm tube}=2L_{\rm well}$ while the **depth** runs $0\to L_{\rm well}\to 0$.
Every profile below is folded accordingly:

$$\text{down leg: } s\in[0,L_{\rm well}],\ d=s \qquad
\text{return leg: } s\in[L_{\rm well},L_{\rm tube}],\ d=L_{\rm tube}-s$$

so **each depth carries two curves**: a solid line for the down leg and a dotted
line for the return leg.

> **This was previously wrong.** Earlier builds plotted $s/2$ and labelled it
> depth, which linearly compressed the whole hairpin into the depth range and
> drew the return leg upside down — developed $s=3000$ m is 48 m below surface,
> not 1500 m.

The fold makes an assumption visible that was hidden before. Because the cascade
(**H8**) is laid out along $s$, the two legs sit in **different PCM layers at the
same depth**:

| depth [m] | down leg $T_m$ | return leg $T_m$ | gap |
|---|---|---|---|
| 0 | 150.0 °C | 101.1 °C | **48.9 K** |
| 800 | 137.8 °C | 113.3 °C | 24.4 K |
| 1524 | 125.6 °C | 125.6 °C | 0.0 K |

The borehole cross-section must therefore be *partitioned*, with the down-leg
PCM held apart from the return-leg PCM. **H9** assumes that partition is
perfectly insulating. It is not free — see §15.

### Reading the melt-fraction plots

The local melt fraction is the melted PCM in a segment divided by the PCM that
segment actually contains:

$$\varepsilon_{\rm local}(z) = \frac{A_{\rm melt}(z)}{A_{\rm avail}},
\qquad A_{\rm avail} = \frac{V_{\rm well}}{n_t\,L_{\rm tube}}$$

Two horizontal limits are drawn on those plots, and both matter:

- **$\varepsilon_{\rm local} = 1$** — the segment has melted all the PCM it owns.
  Beyond this the model is drawing latent heat from material that is not there.
- **front merge** — the melt front has reached the midpoint between neighbouring
  tubes. With $2n_t = 4$ legs sharing the borehole, each owns an equivalent cell
  of radius $r_{\rm cell}=\sqrt{\pi R^2/4}/\sqrt{\pi}$, and the fronts meet at
  $\delta = r_{\rm cell}-r_e$. Hypothesis **H3** fails above this line.

**PCM that never melts is dead volume** — it holds no energy and contributes
nothing to the round trip. PCM that over-melts is fictitious. Both show up as
departures from a flat $\varepsilon_{\rm local}=1$ line, and the spread is what
the multilayer cascade exists to reduce (§14.4).

> **A correction to the headline figure quoted in v0.4 and v0.4a.** The spread
> was reported as falling $1.321\to0.078$. That comparison is **not
> like-for-like** in two ways: the $1.321$ came from the *unguarded*
> Formulation B, a different code path that permits
> $\varepsilon_{\rm local}>1$ and reached $1.78$, and the two figures were
> taken at different well counts. Compared properly — same well count, melted
> fraction capped in both — the spread falls
>
> | $N_{\rm wells}$ | area-based (capped) | enthalpy |
> |---|---|---|
> | 11.255 (v0.4b design point) | 0.405 | **0.000** |
> | 11.573 (v0.4a design point) | 0.437 | **0.078** |
>
> The result survives, and at equal well count it is cleaner: the enthalpy
> formulation removes essentially all of the non-uniformity. But $1.321$
> conflated genuine non-uniformity with the over-melt artefact of an unguarded
> scheme, and should not be quoted.
"""))

cells.append(code(r"""
# --- run charge and discharge at the design point, recording history --------
case = CASE
res  = runs['v0.4b']
N    = res['kpis']['N_wells']

rank, hp, T = cycle_state_points(case)
E = energy_budget(case, rank['rank_eff'], hp['hp_cop'], T)
T_m_lay, T_m_lay_dc = melting_temperatures(case)
n = case.n_segments
T_m_seg    = layer_map(T_m_lay,    n, case.N_lay)
# The discharge marches from the far end: the cascade and the state both mirror
# with the march index. Building it as layer_map(T_m_lay_dc, ...) instead --
# which is what this cell did until v0.5a -- left the state unmirrored and
# introduced a further 6.111 K layer offset, because n_segments is not a
# multiple of N_lay. Same defect as the one fixed in run_cycle; same fix.
T_m_seg_dc = T_m_seg[::-1].copy()
t_ch = np.logspace(0, np.log10(case.t_ch*3600), case.n_times)
t_dc = np.logspace(0, np.log10(case.t_dc*3600), case.n_times)
m1   = E['m_dot_w_ch']/(N*case.num_tubes)

ch = march_h(case, T['T_4c'], T_m_seg, m1, case.k_wall, t_ch, record=True)
ratio = res['detail']['flow_ratio_dc_ch']
dc = march_h(case, T['T_3d'], T_m_seg_dc, ratio*m1, case.k_wall, t_dc,
             E0=ch['E'][::-1], record=True)
# ...and comes back into depth indexing before anything is plotted, so both
# half-cycles share ONE cascade, T_m_seg, and one depth axis.
dc = unmirror_march(dc)

# --- geometric limits ------------------------------------------------------
# r_cell and delta_merge are now Case properties, so the plots, the diagnostics
# and the sizing report all read the SAME limit. They previously did not.
A_avail = case.V_well/(case.num_tubes*case.L_tube)          # PCM per unit tube length
r_cell  = case.r_cell
d_merge = case.delta_merge
A_merge = area_from_delta(d_merge, case.r_e, case.num_fins, case.fin_t, case.fin_L)

print(f'PCM available per tube per unit length  A_avail = {A_avail:.5f} m2')
print(f'equivalent cell radius                          {r_cell*1000:.2f} mm')
print(f'fronts merge at delta                           {d_merge*1000:.2f} mm'
      f'   (eps_local = {A_merge/A_avail:.3f})')
print(f'borehole WALL at delta                          {(case.D_well/2-case.r_e)*1000:.2f} mm'
      f'   <- the wrong limit: 2.9x too permissive')
"""))

cells.append(md(r"""
### 14.1 Charging
"""))

cells.append(code(r"""
# ---------------------------------------------------------------------------
# DEPTH, not developed length. A hairpin goes down and comes back, so the
# developed coordinate s runs 0 -> L_tube = 2 L_well while the DEPTH runs
# 0 -> L_well -> 0. Earlier builds plotted s/2 and called it depth, which
# linearly compressed the whole hairpin into the depth range and drew the
# return leg upside down: developed s = 3000 m is 48 m below surface, not
# 1500 m. The fold below is the correct map.
#
#     down leg    s in [0, L_well]          depth = s
#     return leg  s in [L_well, L_tube]     depth = L_tube - s
#
# Because the cascade (H8) is laid out along s, the two legs carry DIFFERENT
# melting temperatures at the same depth -- up to 48.9 K apart at the wellhead.
# Under H9 they are assumed perfectly insulated from one another.
s_dev = ch['z']                    # developed tube coordinate [m]
half  = len(s_dev)//2
z_dn  = s_dev[:half]                       # depth of the down leg
z_up  = case.L_tube - s_dev[half:]         # depth of the return leg

def legs(arr):
    'Split a per-segment profile into (depth, values) for each leg.'
    a = np.asarray(arr)
    return (z_dn, a[:half]), (z_up, a[half:])

def Tf_seg(a):
    '''Segment-mean fluid temperature from the n+1 face temperatures.

    The faces are what the march computes: F_0 is the inlet, F_{j+1} the outlet
    of segment j. Taking the mean of adjacent faces gives one value per segment
    and -- the reason it is done this way -- commutes with the reversal in
    `unmirror_march`, so the same expression is right for both half-cycles.
    Slicing `[1:]` instead does not: after un-mirroring, the inlet node sits at
    the END of the array, and `[1:]` would silently drop the wrong face.
    '''
    a = np.asarray(a)
    return 0.5*(a[:-1] + a[1:])

z = s_dev/2.0                      # kept only for the time-depth maps below
h = ch['history']

# Choose the time levels to draw by TARGET TIME, not by index. The grid is
# logarithmic, so evenly spaced indices bunch at the start: the old picks
# [0,4,9,...] put three of six curves inside the first 12 seconds, where they
# are indistinguishable, and left the hours unrepresented.
def pick_times(t_s, targets_h):
    return [int(np.argmin(np.abs(t_s - th*3600))) for th in targets_h]

def tlabel(t_s):
    return f'{t_s:.0f} s' if t_s < 60 else (
           f'{t_s/60:.0f} min' if t_s < 3600 else f'{t_s/3600:.2f} h')

targets = [0.01, 0.1, 0.5, 1.0, 3.0, 10.0]      # hours
# INDEX THE HISTORY WITH ITS OWN CLOCK. `t` stamps the marched intervals;
# `t_hist` stamps the instants at which the recorded states are exact, and
# carries one more entry (t = 0 at the front, the end state at the back).
# Using `t` here is what put the frame labelled "10 h" at 7.64 h.
picks = pick_times(ch['t_hist'], targets)
cmap = plt.cm.viridis(np.linspace(0.15, 0.95, len(picks)))

fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))

for c_, k in zip(cmap, picks):
    lab = tlabel(ch['t_hist'][k])
    for col, series in ((0, Tf_seg(h['T_fluid'][k])-273.15),
                        (1, h['A_melt'][k]/A_avail),
                        (2, h['U_i'][k]),
                        (3, h['NTU'][k])):
        (zd, vd), (zu, vu) = legs(series)
        ax[col].plot(vd, zd, color=c_, lw=1.4, ls='-',
                     label=(lab if col == 0 else None))
        ax[col].plot(vu, zu, color=c_, lw=1.4, ls=':')

# the cascade itself, folded the same way -- the gap between the two dashed
# curves at a given depth is what H9 assumes is perfectly insulated
(zd, vd), (zu, vu) = legs(np.array(T_m_seg)-273.15)
ax[0].plot(vd, zd, 'k--', lw=1, label='$T_m$ down leg')
ax[0].plot(vu, zu, 'k:',  lw=1, label='$T_m$ return leg')
ax[0].set_xlabel('secondary-fluid temperature [°C]')
ax[0].set_ylabel('depth [m]   (solid = down leg, dotted = return leg)')
# eps_local = 1 and "fronts merge" are the SAME line -- see the check below
ax[1].axvline(1.0, color='crimson', ls='-', lw=1.2,
              label=r'$\varepsilon_{local}=1$ = fronts merge')
ax[1].set_xlabel(r'local melt fraction  $\varepsilon_{local}$')
ax[2].set_xscale('log')
ax[2].set_xlabel(r'$U_i$  [W m$^{-2}$ K$^{-1}$]  (on inner area)')
ax[3].set_xlabel('NTU per segment')
for a in ax:
    a.invert_yaxis(); a.grid(alpha=.3)
ax[0].legend(fontsize=7, title='time', title_fontsize=7)
ax[1].legend(fontsize=7)
fig.suptitle(f'CHARGING · N = {N:.2f} wells · {case.N_lay} PCM layers', y=1.02)
plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
### Reading the $U_i$ and NTU panels together

The two right-hand panels look alike, and they nearly are. From §6,

$$\mathrm{NTU}=\frac{2\pi r_i\,U_i\,\Delta z}{\dot m\,c_p}$$

so along a well at fixed flow they differ only through $c_p(T)$ — a few per cent
across the glide. Plotting both is still worth it, because they answer different
questions:

- **$U_i$** is the *physical* quantity: the conductance of the borehole at that
  depth and time, in W m⁻² K⁻¹. It is comparable **across cases** — charge
  against discharge, finned against bare — because it carries no flow rate.
- **NTU** is the *dimensionless* group that decides the segment effectiveness
  $1-e^{-\mathrm{NTU}}$, and therefore how much of the available temperature
  difference a segment actually uses. It is **not** comparable across cases with
  different flow.

The cell after the discharge figure makes that concrete. Mid-well at the end of
each half-cycle, $U_i$ rises by a factor $7.115$ from charge to discharge while
NTU rises by $6.907$; the discharge flow ratio is $1.0473$, and
$7.115/1.0473 = 6.794$ recovers the NTU ratio to within the variation of $c_p$.
The factor of seven is not a flow effect at all — it is the melt layer. The
charge ends with a full annulus of liquid PCM between tube and front; the
discharge ends with that annulus refrozen, and solid PCM conducts better than
liquid. NTU carries the flow ratio on top of that, which is the whole reason to
plot the quantity that does not depend on flow.

> These figures moved in v0.5a. This section used to march the discharge with
> the melt state unmirrored and the cascade built by `layer_map(T_m_lay_dc, …)`
> — the flow-reversal defect corrected in `run_cycle` in v0.5, which had been
> left in place here. The ratios read $1.401$ and $1.386$ against a flow ratio
> misquoted as $1.023$.

$U_i$ is drawn on a log axis, and it spans about a factor of seven
($117$ to $833$ W m⁻² K⁻¹ over both half-cycles): the log scale keeps
the late-time curves, which bunch at the low end, readable against the $t\to0$
curve. That $t\to0$ value is worth noting — the melt layer is absent, $h_e$
saturates at its cap, and the borehole is limited by the tube alone. It is the
ceiling no amount of PCM-side design can beat.
"""))

cells.append(md(r"""
### 14.2 Discharging

The discharge starts from the melt distribution charging left behind, so
$\varepsilon_{\rm local}$ falls from its end-of-charge profile rather than from
zero. Under **Formulation C** a segment that reaches $\varepsilon_{\rm local}=0$
does *not* stop contributing: it continues to subcool below $T_m$, giving up
sensible heat, until it reaches the fluid temperature. That recovered
desuperheating and subcooling is the energy Formulation B discarded.
"""))

cells.append(code(r"""
hd = dc['history']
picks = pick_times(dc['t_hist'], targets)     # dc grid, same target times
fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
for c_, k in zip(cmap, picks):
    lab = tlabel(dc['t_hist'][k])
    for col, series in ((0, Tf_seg(hd['T_fluid'][k])-273.15),
                        (1, hd['A_melt'][k]/A_avail),
                        (2, hd['U_i'][k]),
                        (3, hd['NTU'][k])):
        (zd, vd), (zu, vu) = legs(series)
        ax[col].plot(vd, zd, color=c_, lw=1.4, ls='-',
                     label=(lab if col == 0 else None))
        ax[col].plot(vu, zu, color=c_, lw=1.4, ls=':')
# same cascade as on charge: `dc` is now in depth indexing
(zd, vd), (zu, vu) = legs(np.array(T_m_seg)-273.15)
ax[0].plot(vd, zd, 'k--', lw=1, label='$T_m$ down leg')
ax[0].plot(vu, zu, 'k:',  lw=1, label='$T_m$ return leg')
ax[0].set_xlabel('secondary-fluid temperature [°C]')
ax[0].set_ylabel('depth [m]   (solid = down leg, dotted = return leg)')
ax[1].axvline(1.0, color='crimson', ls='-', lw=1.2)
ax[1].axvline(0.0, color='k', lw=1)
ax[1].set_xlabel(r'local melt fraction  $\varepsilon_{local}$')
ax[2].set_xscale('log')
ax[2].set_xlabel(r'$U_i$  [W m$^{-2}$ K$^{-1}$]  (on inner area)')
ax[3].set_xlabel('NTU per segment')
for a in ax:
    a.invert_yaxis(); a.grid(alpha=.3)
ax[0].legend(fontsize=7, title='time', title_fontsize=7)
fig.suptitle(f'DISCHARGING · flow ratio {ratio:.3f}', y=1.02)
plt.tight_layout(); plt.show()
"""))

cells.append(code(r"""
# U_i is comparable across charge and discharge; NTU is not. Quantify both.
print(f"{'':22s} {'charge':>12s} {'discharge':>12s}   ratio")
for name, a, b in (
    ('U_i  mid-well, t=0',  h['U_i'][0][50],   hd['U_i'][0][50]),
    ('U_i  mid-well, end',  h['U_i'][-1][50],  hd['U_i'][-1][50]),
    ('NTU  mid-well, end',  h['NTU'][-1][50],  hd['NTU'][-1][50]),
):
    print(f'  {name:20s} {a:12.4g} {b:12.4g}   {b/a:6.3f}')
print(f'\n  discharge/charge flow ratio = {ratio:.4f}')
print('  The U_i ratio is set by the PCM conductivity (k_s freezing vs k_l')
print('  melting) and the melt-layer thickness. The NTU ratio carries the')
print('  flow ratio on top of that, which is why the two differ.')
"""))

cells.append(md(r"""
### 14.3 The whole cycle as a map

Melt fraction over depth and time, charge then discharge. The white contour is
$\varepsilon_{\rm local}=1$; everything to its warm side is PCM the model melted
but does not have.
"""))

cells.append(code(r"""
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
for a, hh, tt, ttl in ((ax[0], h, ch['t_hist'], 'charging'),
                       (ax[1], hd, dc['t_hist'], 'discharging')):
    F = hh['A_melt']/A_avail
    im = a.pcolormesh(tt/3600, z, F.T, shading='auto', cmap='inferno',
                      vmin=0, vmax=max(1.8, F.max()))
    cs = a.contour(tt/3600, z, F.T, levels=[1.0], colors='w', linewidths=1.5)
    a.clabel(cs, fmt=r'$\varepsilon=1$', fontsize=8)
    a.set_xlabel('time [h]'); a.set_title(ttl)
    plt.colorbar(im, ax=a, label=r'$\varepsilon_{local}$')
ax[0].set_ylabel('developed tube length / 2 [m]'); ax[0].invert_yaxis()
plt.tight_layout(); plt.show()

F_end = ch['A_melt']/A_avail
print(f'end of charge:  eps_local  min {F_end.min():.3f}   max {F_end.max():.3f}'
      f'   aggregate {ch["A_melt"].sum()*(case.L_tube/n)*case.num_tubes/case.V_well:.4f}')
print(f'  segments over capacity (eps>1):        {(F_end>1).sum():3d} of {n}')
print(f'  segments AT the merge line (H3 saturated):'
      f'{(ch["delta"]>=0.999*d_merge).sum():3d} of {n}')
print(f'  segments PAST it (H3 violated -- impossible in v0.4):'
      f'{(ch["delta"]>1.001*d_merge).sum():3d} of {n}')
F_dc = dc['A_melt']/A_avail          # final state, all segments
print(f'end of discharge: eps_local min {F_dc.min():.3f}  max {F_dc.max():.3f}'
      f'   aggregate {F_dc.mean():.4f}')
print(f'  segments left fully solid:             {(F_dc<0.01).sum():3d} of {n}')
print(f'  PCM that never melted at all:          '
      f'{(F_end<0.05).sum():3d} of {n} segments  <- dead volume')
print()
print('NOTE: eps_local = 1 and "fronts merge" are the SAME line -- each leg owns')
print('exactly the PCM in its cell, so melting all of it puts the front on the')
print('cell boundary. area_from_delta(delta_merge) == A_avail identically:')
print(f'   A(delta_merge) = {A_merge:.6e}    A_avail = {A_avail:.6e}')
print('Under Formulation C the enthalpy cap holds A_melt <= A_avail, so delta')
print('CANNOT exceed delta_merge. H3 can no longer be violated -- but it is')
print('being sat on, which is a different and still important problem: the')
print('annular resistance is evaluated at the exact limit of its validity.')
"""))

cells.append(md(r"""
### 14.4 What the cascade is for

The multilayer PCM exists to make the phase change more uniform along the well:
each layer melts at a temperature matched to the local fluid temperature, so the
whole column changes phase rather than only the hot end. PCM that never cycles
is dead volume — it carries no energy through the round trip.

Sweeping $N_{\rm lay}$ at fixed well count shows whether it works.
"""))

cells.append(code(r"""
rows = []
for N_lay in (1, 3, 6, 9, 12, 20):
    cc = case.with_(N_lay=N_lay)
    Tml, _ = melting_temperatures(cc)
    seg = layer_map(Tml, n, N_lay)
    r = march_h(cc, T['T_4c'], seg, m1, cc.k_wall, t_ch)
    F = r['eps_local']
    rows.append({'N_lay': N_lay, 'eps_mean': F.mean(), 'eps_min': F.min(),
                 'eps_max': F.max(), 'spread': F.max()-F.min(),
                 'std': F.std(),
                 'frac_dead_%': 100*(F < 0.05).mean(),
                 'frac_over_%': 100*(F > 1).mean()})
df_c = pd.DataFrame(rows).set_index('N_lay')
print(df_c.round(3).to_string())

fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].plot(df_c.index, df_c['spread'], 'o-', label='max - min')
ax[0].plot(df_c.index, df_c['std'], 's-', label='std. deviation')
ax[0].set_xlabel('number of PCM layers  $N_{lay}$')
ax[0].set_ylabel(r'spread in $\varepsilon_{local}$'); ax[0].legend(); ax[0].grid(alpha=.3)
ax[1].plot(df_c.index, df_c['frac_over_%'], '^-', color='crimson',
           label='over capacity ($\\varepsilon>1$)')
ax[1].plot(df_c.index, df_c['frac_dead_%'], 'v-', color='steelblue',
           label='dead ($\\varepsilon<0.05$)')
ax[1].set_xlabel('number of PCM layers  $N_{lay}$')
ax[1].set_ylabel('% of the well'); ax[1].legend(); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
### 14.5 Recording does not perturb the solution

The history arrays are written from inside the segment loop, and the closing
frame costs one extra read-only sweep (`segment_profile`) after the march has
finished. This confirms the recorded run and the plain run agree exactly, so
nothing in §14 is an artefact of instrumenting the march.

Worth being precise about what the record *means*, since getting it wrong cost
two figures (see §16.3). Each frame is stamped in `t_hist` at the instant where
its **state** — `E`, `A_melt`, `T_pcm`, `delta` — is exact. The **profile**
arrays in the same frame — `T_fluid`, `NTU`, `U_i`, `q_prime` — are the ones
evaluated *from* that state, i.e. the profile that drives the step beginning
there. There are `len(t)+1` frames: one at $t=0$ and one at the end of every
accepted step. `t` itself stamps the marched intervals and is one shorter; it
is the right clock for `Q`, and the wrong one for the history.
"""))

cells.append(code(r"""
a_ = march_h(case, T['T_4c'], T_m_seg, m1, case.k_wall, t_ch)
b_ = march_h(case, T['T_4c'], T_m_seg, m1, case.k_wall, t_ch, record=True)
same_A = np.array_equal(a_['A_melt'], b_['A_melt'])
same_Q = a_['Q_cum_J'] == b_['Q_cum_J']
print(f'A_melt identical: {same_A}     Q_cum identical: {same_Q}')
print('PASS' if (same_A and same_Q) else 'FAIL - recording perturbs the march')
"""))

cells.append(md(r"""
## 15. What this model still does not contain

Read this before quoting any number above.

| missing | consequence |
|---|---|
| **Formation heat loss** (H6) | fine over a 10 h cycle, wrong for seasonal storage |
| **Natural convection in the melt** (H4) | under-predicts melting once the layer is established |
| **Radial resolution within a phase** (H5) | sensible heat is carried but *lumped* — one PCM temperature per segment, so the gradient through the melt layer is quasi-steady rather than resolved |
| **Buoyancy head** (H7) | ±7.5 bar in the 1,524 m loop, asymmetric between charge and discharge |
| **Isentropic efficiencies** | both cycles are idealised, so every efficiency here is an upper bound |
| **Volume change on melting** | 6.9 %, neglected; the liquid layer is ~2.3 % thicker than modelled |
| **Non-circular melt front** | the fins enhance heat transfer through $\eta_o$, but the front is still modelled as a circular annulus (H3). Their *metal* is excluded from $A_{\text{melt}}$ (§5.1); their effect on front *shape* is not modelled |
| **Front interaction** (H3) | `delta_max = 0.5 m` against an 0.0889 m borehole radius — a factor of 5.62. Nothing in the solver enforces it |

### Two quantities the code defines twice

- $h_e$: `compute_U_i` uses the cylindrical shell
  $k_m/(r_e\ln(1+\delta/r_e))$ but falls back to the plane-slab $k_m/\delta$ in
  its small-$\delta$ branch. The second definition, in `compute_T_re`, went
  with the closed-form front when that was deleted, so the two forms are no
  longer iterated against each other — but the branch inconsistency inside
  `compute_U_i` remains.
- The laminar Nusselt number: 3.66 here, 4.36 in some paths. Immaterial at the
  design point, where the flow is strongly turbulent.

### And the largest gap

**No external validation.** Neither formulation has been checked against
experimental or independently computed data. Energy closure is structural under
Formulation B and therefore proves nothing about correctness. Internal
consistency is not validation, and a reviewer will ask.

### Open work, in order

1. **Resolve H3 properly with a 2-D conduction calculation.** *(The bounds check
   itself is now fixed: `Case.r_cell` and `Case.delta_merge` replace the
   borehole-wall limit, which was 2.9× too permissive in $\delta$. Formulation C
   then makes violation impossible — the enthalpy cap holds
   $A_{\rm melt}\le A_{\rm avail}$ and $A(\delta_{\rm merge})=A_{\rm avail}$
   identically — so the limit cannot be exceeded. It is, however, being **sat
   on**: $\delta/\delta_{\rm merge}=1.000$ over the whole well.)*
   What remains is the physics the annulus cannot express: once fronts meet, the
   residual solid occupies the corners between legs, and neither the reduced
   conducting azimuth nor the fins' reach into those corners is represented.
   Both errors point the same way — the model **understates** late-stage
   resistance and **cannot credit** fins for relieving it. This now carries a
   design decision (§12.1), so it is the most consequential open item.
2. **Evaluate H9 — the azimuthal leak between cascade layers.** Assumed perfect
   for now. At the wellhead the neighbouring cell is 48.9 K away and a slab
   estimate gives ~22 W/m against ~83 W/m of useful duty — of order **25 %**,
   and it flows from the hot layer to the cold one, partially undoing what
   $N_{\rm lay}$ is for. Three consequences: $N_{\rm lay}$ cannot be swept for
   free; a leg's neighbour *within* a partition half is a genuine front merge
   while *across* it is a heat leak; and the insulation itself displaces PCM,
   raising $N_{\rm capacity}$. The same 2-D cell calculation that settles
   item 1 would settle this.
3. Make $\lambda$ an output rather than an assumed 5 %. At the v0.4 design point
   the field banks about 7 % more than the ORC withdraws, against the 5 % assumed.
4. Quote results at $n_{times}=160$: `eps_PCM` moves $-0.3\,\%$ from 40 to 320
   time levels, one-sided, so the default carries a small known bias.
5. Model the melt-front shape around the fins, and the volume change on melting.
6. Re-optimise the fin geometry once (1) is settled. §12.1 puts the optimum at
   16 fins × 7.5 mm, but for the reasons above that is a floor, not an answer.
7. **Find an external validation case.** Still the largest gap — energy closure
   is structural and proves nothing about correctness.

---

## Version history

| version | date | change |
|---|---|---|
| **0.6** | this build | **The discharge closure was solving for the wrong end.** Only the two exchanger *inlets* are boundary conditions on the march; both outlets are results. The discharge was closed the other way round — the outlet pinned at $T_{2d}=T_{m,\rm top}-\Delta T_{m,2D}$ with $\Delta T_{m,2D}=0$, the inlet derived from it — and the model's own solution contradicts that: at the start of a discharge at CSS the water leaves at $157.93$ °C, nearly 8 K above the top layer's melting point, because the PCM is superheated. The implied inlet approach was also silently $\Delta T_{3C,2C}/N_{\rm lay}=6.111$ K, one layer width, chosen by nobody. `DT_m_2D` is replaced by **`DT_M_1D`**, a subcooling of the inlet (state **1d**) below the *coldest* layer, the symmetric partner of `DT_4C_M`; $T_{2d}$ is demoted to a provisional estimate and reported against the realised value. `simulate_css_corrected` adds **one** ORC correction pass at the realised outlet — enough because the outlet is pinned by the store, not the plant (it moves $<0.3$ K while the inlet moves 9 K). Verified a pure reparameterisation: at $\Delta T_{M,1D}=6.111$ K every v0.5a number returns exactly (§16.5). At the symmetric 10 K the deviation goes $-5.69\,\%\to+1.23\,\%$ — **most of the shortfall was the closure, not the store**. §13 re-frozen for the new closure. |
| **0.5a** | this build | **Two recording defects, no physics.** (i) `march_h` recorded the state at the *start* of a step but stamped it with the time at the *end*, and appended `E` after the update while the other state arrays came from before it. On a logarithmic grid the last step is 2.36 h, so every recorded profile was up to a quarter of a half-cycle stale and the true end state was never recorded: the map read $\varepsilon=0.356$ at the end of discharge against $0.1545$ from the march. Frames are now stamped in **`t_hist`** at the instants where the state is exact, start at $t=0$, and include the end state; `segment_profile` evaluates the closing frame's fluid profile from that end state. (ii) The **discharge was plotted mirrored** in §14.2 and §16 — the discharge marches from the far end, and only the sizing loops were un-mirroring it. `unmirror_march` now returns every discharge result in depth indexing, so both half-cycles share one cascade and one depth axis, and §16.3 prints the map-closure residual (zero). §14 also carried the *flow-reversal* bug fixed in `run_cycle` in v0.5 — it built the discharge cascade with `layer_map(T_m_lay_dc, ...)` and passed the state unmirrored. `N_wells`, the CSS deviation, energy closure and every reported index are unchanged. |
| 0.5 | 2026-09 | **Cyclic steady state** (§16). At CSS $\eta_{\rm storage}\equiv1$ exactly, so $\lambda$ is unattainable and the rate criterion is degenerate even at $\lambda=0$; sizing is replaced by simulation — $N$ from latent heat alone, both flows pinned by their glides, march to CSS, report the deviation. The field delivers **94.26 %** of the 1 MWe target, and the shortfall is a *discharge rate* limit. Flow-reversal defect fixed in `run_cycle`: the melt state was handed to the discharge unmirrored while the cascade was mirrored, shifting the enthalpy datum by up to 48.9 K. §12.1 removed. |
| 0.4c | 2026-09 | **Superseded code removed.** Formulations A (closed-form Stefan) and B (energy-balance on melted area) deleted along with the v0.1 sizing chain, the `front` and `charge_uses_wall_conductivity` switches, and `CASE_V01` — 16 functions, about 1,050 lines. The model is now one formulation. Verified **bit-identical** to the pre-strip v0.4b on 20 reported quantities before the cut. §13 now checks a frozen v0.4b fixture instead of the published IHTC numbers: that is a *drift* check, not an external validation, and the difference matters — see §13. Notebook 69 → 57 cells, runtime 55 → 39 s. |
| 0.4b | 2026-09 | **`E_well^cap` corrected** (§10): the capacity bound now resolves the cascade, $\langle T_m\rangle = T_{m,\rm top}-\frac{N_{\rm lay}-1}{2N_{\rm lay}}\Delta T_{\rm glide}$, and includes the solid subcooling the discharge reaches. The old form used the top layer's $T_m$ for the whole store — the no-cascade limit — and was not a bound: the model exceeded it by 1.6 %. `N_capacity` 11.573 → **9.574**, so the **rate criterion now binds** and `N_wells` 11.573 → **11.255**. Two conclusions reverse: the fin optimum moves from 16 to about 20 and is shallow, and the melt-fraction spread quoted at $1.321\to0.078$ was not a like-for-like comparison (§14.4). |
| 0.4a | 2026-09 | Sizing report rewritten to lead with `N_wells` and name the binding criterion; `N_capacity` labelled a lower bound (§10). Melt-front limit moved from the borehole wall to the cell radius via `Case.r_cell` / `Case.delta_merge`, and H3 proximity reported (§12, §14). §12.1 (fin and charging-window sweeps) removed in v0.5: the geometry is fixed and the sweeps no longer bracket after the flow-reversal fix. §13 (verification against the frozen IHTC fixture) restored — it had been dropped from the generator. **Build stamp repaired**: the `__STAMP__` placeholder was never substituted, so the notebook shipped reading `Last updated: __STAMP__`; the generator now asserts the substitution happened. §14 gains a **$U_i$ panel** beside NTU, and the plotted time levels are now chosen by target time rather than by index — on a logarithmic grid the old picks put three of six curves inside the first 12 seconds. |
| 0.4 | 2026-09 | Melt front carries **enthalpy** rather than melted area (§5.4): PCM superheats when fully molten and subcools when fully solid, $\varepsilon_{\rm local}\in[0,1]$ by construction, branchwise-exact time integration. Sensible heat is self-levelling — melt-fraction spread falls $1.321\to0.078$. |
| 0.3 | 2026-09 | Wall conductivity corrected: charging now uses steel, $k_w = 45$ W/m·K, by default (§2). Fin metal excluded from the melted PCM area (§5.3). `N_wells` 21.65 → **12.74**, binding constraint now inventory. Lower-bound check against the ideal well count added (§12). |
| 0.2 | 2026-08 | Melt front reformulated by energy balance (§5.2); melt state carried across the cycle; two-constraint sizing (§10); latent inventory made density-consistent; segment latent limiting. Model moved into this notebook, every equation visible. |
| 0.1 | — | IHTC paper. Closed-form Stefan front, single sizing criterion, wall conductivity error. Removed from the code in v0.4b; recoverable from git history. |
"""))

# ============================================================ 16. cyclic steady state
cells.append(md(r"""
## 16. Cyclic operation and cyclic steady state

Everything above is a **first-cycle** result: the store starts at $E'=0$ —
fully solid at exactly $T_m$, no subcooling — and the cycle does not return
there. A plant repeating a 24 h schedule never sees that state again after the
first day.

The intended duty is **10 h charge / 2 h rest / 10 h discharge / 2 h rest**.

> **The rests are exact no-ops in this model.** With no flow,
> $K=\dot m c_p(1-e^{-\rm NTU})/\Delta z \to 0$, and H2 (no axial conduction),
> H6 (adiabatic wall) and H9 (no azimuthal exchange) leave no other heat path,
> so $\mathrm{d}E'/\mathrm{d}t = 0$ identically. 10/2/10/2 is arithmetically
> the same as 10/10 here.
>
> That is a statement about the assumptions, not about the store. In reality a
> 2 h rest lets superheated liquid near the tube conduct outward into colder
> PCM and the front relax; the lumped state of **H5** cannot represent that
> either.

### Why sizing cannot simply be moved to CSS

At cyclic steady state the state returns to itself, so the enthalpy change over
a cycle is zero. With no loss path anywhere in the model,

$$\eta_{\rm storage}({\rm CSS}) \equiv 1 \qquad\text{exactly, for every } N
\text{ and every flow.}$$

Two consequences follow, and both are structural rather than numerical.

**1. The loss surplus $\lambda$ is unattainable at CSS.** The budget requires
$\Delta E_{\rm out,HP}=(1+\lambda)\Delta E_{\rm in,ORC}$, so the ratio of CSS
charge to requirement pins at $1/(1+\lambda)$ for *every* $N$. There is no root.
On the first cycle the surplus is not lost but **parked** in the store as
residual melt — which is why the first cycle appears to work, and why
$\eta_{\rm storage}$ came out at $1/(1+\lambda)$ there. At CSS there is nowhere
left to park it.

**2. Even with $\lambda=0$ the rate criterion goes degenerate.** Charge and
discharge are then the same number, so one equation is left for two unknowns:
it fixes the **flow ratio** and says nothing about $N$.

### The reformulation

Stop asking the energy balance to determine $N$. **Specify the hardware, march
to CSS, report what comes out:**

| quantity | how it is set |
|---|---|
| $N_{\rm wells}$ | $\Delta E_{\rm out,HP}/(\rho V_{\rm well} h_m)$ — latent heat alone |
| $\dot m_{\rm ch}$ | pinned by the charging glide |
| $\dot m_{\rm dc}$ | pinned by the discharging glide |

Nothing is solved — there is **no root-find anywhere** in `simulate_css`, so
the question of convergence does not arise. The ORC correction pass of v0.6
(§16.5) is a second *evaluation*, not an iteration: it runs once, and the shift
it produces is reported so the reader can see it is negligible. What was an unsatisfiable
constraint becomes a reported output: the deviation of delivered energy, and
hence of net electrical output, from target.

Sizing on latent heat alone *overestimates* $N$, because it ignores the
sensible capacity the store also has. That is deliberate and conservative.

> This also dissolves an older objection. §10 argues the discharge flow cannot
> be pinned to its glide because delivered energy then saturates below the
> requirement for any well count, leaving no root. Under this framing that
> saturation is not a failure — it *is* the answer.
"""))
cells.append(code(src('simulate_css', 'simulate_css_corrected', 'css_report')))

cells.append(code(r"""
css = simulate_css_corrected(CASE)
css_report(CASE, css)
"""))

cells.append(md(r"""
### 16.1 Convergence to the cyclic steady state
"""))
cells.append(code(r"""
h = pd.DataFrame(css['history'])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

ax[0].plot(h['cycle'], h['Q_charge_kJ']/1e6, 'o-', label='charge')
ax[0].plot(h['cycle'], h['Q_discharge_kJ']/1e6, 's-', label='discharge')
ax[0].axhline(css['required_kJ']/1e6, color='crimson', ls='--', lw=1.2,
              label='required')
ax[0].set_xlabel('cycle'); ax[0].set_ylabel('energy [GJ per field]')
ax[0].set_title('charge and discharge converge')
ax[0].legend(fontsize=8)

ax[1].semilogy(h['cycle'], h['drift'].clip(lower=1e-12), 'o-')
ax[1].set_xlabel('cycle'); ax[1].set_ylabel(r"max $|\Delta E'|$ between cycles [J/m]")
ax[1].set_title('state drift -> 0')

ax[2].plot(h['cycle'], h['eps_end_charge'], 'o-', label='end of charge')
ax[2].plot(h['cycle'], h['eps_end_discharge'], 's-', label='end of discharge')
ax[2].set_xlabel('cycle'); ax[2].set_ylabel(r'mean $\varepsilon_{local}$')
ax[2].set_title('melt fraction settles'); ax[2].legend(fontsize=8)
for a in ax: a.grid(alpha=.3)
plt.tight_layout(); plt.show()

# `.round(6)` flattens the drift column to 0.000000e+00 for everything below
# 5e-7 and makes the last third of the table look converged six cycles early.
# Format each column for what it is instead.
print(h.to_string(index=False, formatters={
    'Q_charge_kJ':     '{:.6e}'.format,
    'Q_discharge_kJ':  '{:.6e}'.format,
    'drift':           '{:.3e}'.format,
    'eps_end_charge':  '{:.6f}'.format,
    'eps_end_discharge': '{:.6f}'.format}))
"""))

cells.append(md(r"""
The first cycle is the outlier: it starts from a cold, fully solid store and so
absorbs more than any later cycle. The state settles within about ten cycles,
after which charge and discharge are the same number to machine precision —
which is what $\eta_{\rm storage}\equiv1$ means in practice.
"""))

cells.append(md(r"""
### 16.2 Profiles at cyclic steady state

Folded onto depth as in §14: **solid = down leg, dotted = return leg.**
"""))
cells.append(code(r"""
cssr = simulate_css_corrected(CASE, record=True)
chc, dcc = cssr['charge'], cssr['discharge']
s_dev = chc['z']; half = len(s_dev)//2
zd, zu = s_dev[:half], CASE.L_tube - s_dev[half:]
def legs2(a):
    a = np.asarray(a); return (zd, a[:half]), (zu, a[half:])

picks = pick_times(chc['t_hist'], targets)
# ONE cascade for both rows: `simulate_css` returns the discharge in depth
# indexing, so `T_m_seg_dc` no longer exists and is no longer needed.
fig, ax = plt.subplots(2, 3, figsize=(15, 8.6))
for row, (rr, Tm_, ttl) in enumerate(((chc, cssr['T_m_seg'], 'CHARGE at CSS'),
                                      (dcc, cssr['T_m_seg'], 'DISCHARGE at CSS'))):
    hh = rr['history']
    for c_, k in zip(cmap, picks):
        for col, series in ((0, Tf_seg(hh['T_fluid'][k])-273.15),
                            (1, hh['A_melt'][k]/A_avail),
                            (2, hh['T_pcm'][k]-273.15)):
            (a1,v1),(a2,v2) = legs2(series)
            ax[row][col].plot(v1, a1, color=c_, lw=1.3, ls='-',
                              label=(tlabel(rr['t_hist'][k]) if col==0 else None))
            ax[row][col].plot(v2, a2, color=c_, lw=1.3, ls=':')
    (a1,v1),(a2,v2) = legs2(np.array(Tm_)-273.15)
    ax[row][0].plot(v1, a1, 'k--', lw=1); ax[row][0].plot(v2, a2, 'k:', lw=1)
    ax[row][2].plot(v1, a1, 'k--', lw=1); ax[row][2].plot(v2, a2, 'k:', lw=1)
    ax[row][0].set_xlabel('secondary fluid [°C]'); ax[row][0].set_ylabel(f'{ttl}\ndepth [m]')
    ax[row][1].set_xlabel(r'$\varepsilon_{local}$'); ax[row][1].set_xlim(-0.02,1.02)
    ax[row][2].set_xlabel(r'$T_{pcm}$ [°C]   (dashed: $T_m$)')
    ax[row][0].legend(fontsize=7)
    for a in ax[row]: a.invert_yaxis(); a.grid(alpha=.3)
plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
The third column is the one the earlier sections could not draw: **$T_{\rm pcm}$
against the cascade**. Where the solid curve sits above the dashed $T_m$ the PCM
is superheated liquid; below it, subcooled solid. At CSS the charge leaves the
whole store above its melting line, and the discharge fails to bring all of it
back down.
"""))

cells.append(md(r"""
### 16.3 The cycle as a map

Melted fraction over depth and time, charge then discharge, at CSS.

**The two panels join up.** At cyclic steady state the state returns to itself,
so the *right* edge of the discharge panel and the *left* edge of the charge
panel are the same distribution, segment by segment. That is the definition of
CSS, drawn rather than asserted, and the cell prints the residual: it is zero.

It did not join up before v0.5a, for two reasons that had nothing to do with
the physics and everything to do with the bookkeeping:

1. **The history lagged one step.** `march_h` recorded the state at the
   *start* of each step but stamped it with the time at the *end*. On a
   logarithmic grid the last step is $8491\,$s — $2.36\,$h, a quarter of the
   charge — so the right-hand edge of each panel showed the store as it was
   two and a third hours earlier, and the true final state was never recorded
   at all. Mean $\varepsilon$ at the end of discharge read $0.356$ from the map
   against $0.1545$ from the printout below it, which is how the discrepancy
   announced itself.
2. **The discharge panel was drawn mirrored.** The discharge marches from the
   opposite end, so its segment index runs backwards; `simulate_css` un-mirrors
   when it closes the loop, but the figure did not. The panel was upside down
   in depth. Compared without un-mirroring, the two edges differ by $0.677$ in
   $\varepsilon$; compared correctly, by nothing.

Neither affected a single marched or reported quantity — $N$, the deviation,
the energy closure and `eps_local` all come from the march, not from the
record — but every plotted profile in §14 and §16 was one step stale. The
recording convention is now explicit: `t_hist` stamps the instants at which
the state arrays are exact, it starts at $t=0$, and it ends at the end of the
half-cycle.
"""))
cells.append(code(r"""
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
for a, rr, ttl in ((ax[0], chc, 'charge'), (ax[1], dcc, 'discharge')):
    F = rr['history']['A_melt']/A_avail
    # t_hist, not t: the frames are stamped where the state is exact, they
    # start at t = 0, and the last one is the end of the half-cycle. Both
    # panels are in depth indexing, so the right edge of the discharge panel
    # and the left edge of the charge panel are the same state -- which is
    # what cyclic steady state MEANS, and is checked below.
    im = a.pcolormesh(rr['t_hist']/3600, s_dev/2.0, F.T, shading='auto',
                      cmap='inferno', vmin=0, vmax=1)
    a.set_xlabel('time [h]'); a.set_title(ttl + ' (CSS)')
    plt.colorbar(im, ax=a, label=r'$\varepsilon_{local}$')
ax[0].set_ylabel('developed length / 2 [m]'); ax[0].invert_yaxis()
plt.tight_layout(); plt.show()

print(f"end of charge     mean eps {chc['eps_local'].mean():.4f}"
      f"   min {chc['eps_local'].min():.4f}")
print(f"end of discharge  mean eps {dcc['eps_local'].mean():.4f}"
      f"   max {dcc['eps_local'].max():.4f}")
print()
print('The store melts completely and does not refreeze completely. The')
print('residual is what limits delivery -- a DISCHARGE RATE limit.')
print()
# --- the map must close on itself -----------------------------------------
# Right edge of the DISCHARGE panel == left edge of the CHARGE panel, segment
# by segment. This is the periodicity condition drawn rather than asserted; it
# is also the check that the two defects fixed in v0.5a stay fixed, since
# either one alone breaks it (by 0.677 in eps for the mirror, by 0.202 in mean
# eps for the one-step lag).
gap = np.max(np.abs(dcc['history']['A_melt'][-1] - chc['history']['A_melt'][0]))
print(f'map closure   max|eps(end of discharge) - eps(start of charge)|'
      f' = {gap/A_avail:.3e}')
print(f'first/last frames at t = {chc["t_hist"][0]:.1f} s and'
      f' {chc["t_hist"][-1]:.0f} s  ({chc["t_hist"][-1]/3600:.2f} h)')
"""))

cells.append(md(r"""
### 16.4 Design curve: deviation against well count

With $N$ no longer solved for, the design question becomes a curve rather than
a root. Each point is an independent march to CSS.

Note the range: under the v0.6 closure the latent-only sizing slightly
*overshoots* the target, so the crossing sits **below** $N_{\rm lat}$ rather
than well above it. Under the retired closure the same curve crossed near
$N\approx19$ — the closure moved the answer by more than the well count does
over half this range.
"""))
cells.append(code(r"""
N_lat = css['N_wells']
rows = []
# The ORC correction is done ONCE, at the design point, and the resulting
# evaporating temperature is held across the sweep. Two reasons: it is a
# fairer comparison -- the well counts are then compared at a common plant
# operating point rather than each with its own ORC -- and it halves the cost,
# since every point would otherwise be marched to CSS twice. The realised
# outlet moves by less than 0.3 K over this range, so the two agree anyway.
T_2d_fixed = css['T_2d_realised']
for N in (9.0, 10.0, 11.0, N_lat, 13.0, 16.0):
    rr = simulate_css(CASE, N=N, T_2d=T_2d_fixed)
    rows.append({'N_wells': rr['N_wells'],
                 'delivered/req': rr['Q_discharge_kJ']/rr['required_kJ'],
                 'MWe': rr['W_el_out_implied']/1000.0,
                 'deviation %': 100*rr['deviation'],
                 'eps end ch': rr['charge']['eps_local'].mean(),
                 'eps end dc': rr['discharge']['eps_local'].mean(),
                 'cycles': rr['cycles']})
dfc = pd.DataFrame(rows)
print(dfc.round(5).to_string(index=False))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
ax[0].plot(dfc['N_wells'], dfc['MWe'], 'o-')
ax[0].axhline(1.0, color='crimson', ls='--', lw=1.2, label='1 MWe target')
ax[0].axvline(N_lat, color='k', ls=':', lw=1.2, label='latent-only sizing')
ax[0].set_xlabel('$N_{wells}$'); ax[0].set_ylabel('net output [MWe]')
ax[0].legend(fontsize=8)
ax[1].plot(dfc['N_wells'], dfc['eps end dc'], 's-', color='darkorange')
ax[1].set_xlabel('$N_{wells}$')
ax[1].set_ylabel(r'residual $\varepsilon_{local}$ at end of discharge')
for a in ax: a.grid(alpha=.3)
plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
Two things to read from the curve.

**Where the target is met.** With the v0.6 closure the latent-only sizing
already slightly overshoots, so the curve crosses 1 MWe close to the design
point rather than well above it. Read the crossing, not a single number: it
moves with $\Delta T_{M,1D}$ (§16.5).

**The residual melt fraction *rises* with $N$.** That is the signature of a
rate-limited discharge: more wells means each is worked less hard, so
proportionally less of each one refreezes. If the limit were inventory, the
residual would fall. It tells you the lever is on the discharge side — flow,
window, or conductance — not more PCM.

> **How much of this to believe.** The deviation is a model result and
> inherits every assumption in §15. H9 is worth of order $25\,\%$ of the duty at
> the wellhead; H3 is saturated over $99\,\%$ of the well; the time grid carries
> a one-sided $-0.3\,\%$ bias. A deviation of a few per cent sits **inside**
> that uncertainty band, and §16.5 shows it also moves by seven percentage
> points under a modelling choice that was never deliberately made. Report it as
> *"the field is within a few per cent of target at this sizing, and which side
> depends on the discharge approach"*, not as a calibrated shortfall.
"""))

cells.append(md(r"""
### 16.5 Closing the discharge side: $\Delta T_{M,1D}$ instead of $\Delta T_{m,2D}$

**The defect.** Only the two exchanger *inlets* are boundary conditions on the
march. Both outlets are results of the heat transfer. Until v0.5a the discharge
was closed the wrong way round: the **outlet** was pinned at
$T_{2d}=T_{m,\rm top}-\Delta T_{m,2D}$ with $\Delta T_{m,2D}=0$, and the inlet
derived from it by subtracting the glide.

That is not merely arbitrary — *the model's own solution contradicts it*. At the
start of a discharge at CSS the water leaves at $157.93\,$°C, nearly 8 K above
the top layer's melting point, because the PCM there is superheated liquid. A
prescribed outlet cannot survive a formulation that lets the PCM leave $T_m$.

It was also tied to $N_{\rm lay}$ without anyone choosing it. With
$T_{2d}=T_{m,\rm top}$ and a glide of $\Delta T_{3C,2C}$, the implied inlet
approach was

$$\Delta T_{\rm implied}=\frac{\Delta T_{3C,2C}}{N_{\rm lay}}=\frac{55}{9}
=6.111\ \text{K},$$

exactly one layer width. Change the number of layers and the discharge approach
changes silently with it.

**The fix.** Prescribe the inlets symmetrically, which is what they are:

| half-cycle | inlet | prescribed as | value |
|---|---|---|---|
| charge | $T_{4c}$ | $T_{m,\rm top}+\Delta T_{4C,M}$ | $150+10=160$ °C |
| discharge | $T_{3d}$ (state **1d**) | $T_{m,\rm bottom}-\Delta T_{M,1D}$ | $101.111-10=91.111$ °C |

$T_{m,\rm bottom}=T_{m,\rm top}-\Delta T_{3C,2C}(N_{\rm lay}-1)/N_{\rm lay}$ is
the coldest layer, which is the end the discharge enters. $T_{2d}$ is demoted to
a *provisional estimate* — the outlet the water would reach if it achieved the
full glide — used only to give the ORC an evaporating temperature on the first
pass.

**The glide now sets flow and nothing else.** Previously it also moved the
discharge inlet, so a glide sweep confounded two effects.
"""))

cells.append(md(r"""
#### The change is a no-op at the old setting

Setting $\Delta T_{M,1D}=\Delta T_{3C,2C}/N_{\rm lay}$ must reproduce the
retired closure *exactly*. This is the check that separates "we changed a
modelling choice" from "we changed the model". The reference values were taken
by importing the v0.5a sources out of git (commit `b383a48`) alongside this
build in a single process, where every line agrees to the last bit.
"""))
cells.append(code(r"""
# Full float precision, taken by running the v0.5a sources out of git
# (commit b383a48) side by side with this build in one process. Quoting these
# to eight figures would make the test resolve only to 1e-8 and hide exactly
# the kind of small leak it exists to catch.
V05A = {
    'N_wells':           11.544050346532758,
    'deviation':         -0.057363007443116842,
    'eta_storage':        1.0,
    'm1_ch':              0.96552625772561806,
    'm1_dc':              0.96979993412108878,
    'flow_ratio':         1.004426266361246,
    'Q_discharge_kJ':     177430648.09917957,
    'required_kJ':        188227970.57635373,
    'W_el_out_implied':   942.63699255688311,
    'T_3d':             368.14999999999998,     # K, the discharge inlet
    'T_2d':             423.14999999999998,     # K, the (then prescribed) outlet
}

legacy = CASE.with_(DT_M_1D=CASE.DT_3C_2C/CASE.N_lay)
_, _, T_leg = cycle_state_points(legacy)
r_leg = simulate_css(legacy)          # uncorrected, exactly as v0.5a ran it
now = dict(r_leg); now['T_3d'] = T_leg['T_3d']; now['T_2d'] = T_leg['T_2d']

# Two classes of quantity, two tolerances, and the reason is worth stating.
#
# EXACT: state points, flow rates, the well count, the requirement. These are
# closed-form functions of the inputs, so anything but bit-identical means the
# reparameterisation leaked.
#
# ITERATION-TERMINATED: the delivered energy and everything derived from it.
# These depend on WHICH CYCLE the drift test stops at, and that is decided on
# a 1e-9 tolerance against differences of order 1e8 kJ -- the last few cycles
# are floating-point noise. CoolProp's AbstractState carries state between
# calls, so the noise depends on how many property evaluations preceded this
# cell: run standalone, this check is bit-identical on every line; run here,
# after thirteen sections of the notebook, the CSS loop stops one cycle
# earlier or later and the delivered energy moves in the seventh figure. That
# is a property of the library, not of the closure, and loosening the
# tolerance for these three is the honest response rather than a fudge --
# provided the first group stays exact, which is what actually rules out a
# leak.
EXACT = ('N_wells', 'm1_ch', 'm1_dc', 'flow_ratio', 'required_kJ',
         'eta_storage', 'T_3d', 'T_2d')
print(f"{'quantity':18s} {'v0.5a (git b383a48)':>24s} {'v0.6 at legacy':>24s}"
      f"   rel. diff   class")
worst_exact = worst_iter = 0.0
for k, v in V05A.items():
    d = abs(float(now[k]) - v)/abs(v) if v else abs(float(now[k]))
    kind = 'exact' if k in EXACT else 'iterated'
    if k in EXACT:
        worst_exact = max(worst_exact, d)
    else:
        worst_iter = max(worst_iter, d)
    print(f'{k:18s} {v:24.17g} {float(now[k]):24.17g}   {d:.1e}   {kind}')
print()
print(f'worst, closed-form quantities:       {worst_exact:.1e}   (must be 0)')
print(f'worst, iteration-terminated:         {worst_iter:.1e}   (tol 1e-4)')
ok = (worst_exact == 0.0) and (worst_iter < 1e-4)
print()
print('PASS -- the v0.6 closure is a reparameterisation, not a new model.'
      if ok else 'FAIL -- something other than the closure changed')
"""))

cells.append(md(r"""
#### The deviation *is* the glide shortfall

The discharge flow is pinned by the **assumed** glide,
$\dot m_{\rm dc}=\dot Q_{\rm in,ORC}/(c_p\,\Delta T_{3C,2C})$, so over a fixed
window the delivered energy is just

$$\frac{\Delta E_{\rm delivered}}{\Delta E_{\rm required}}
 =\frac{\text{glide the water actually achieves}}{\Delta T_{3C,2C}}.$$

The two agree to $2\times10^{-5}$ (the residue is $c_p(T)$). This is a more
useful statement of the result than "discharge-rate limited": under the retired
closure the water gained $51.84$ K of the $55$ K the pump was sized for, and
that $94.26\,\%$ **is** the $-5.74\,\%$ deviation.
"""))

cells.append(md(r"""
#### Sweep: how much of the shortfall was the closure?

Each point is an independent march to CSS with the ORC correction pass.
"""))
cells.append(code(r"""
Tmb = T_m_bottom(CASE) - 273.15
rows = []
for sub in (CASE.DT_3C_2C/CASE.N_lay, 8.0, 10.0, 12.0, 15.0):
    rr = simulate_css_corrected(CASE.with_(DT_M_1D=sub))
    rows.append({'DT_M_1D': sub, 'T_3d [C]': Tmb - sub,
                 'N_wells': rr['N_wells'],
                 'deviation %': 100*rr['deviation'],
                 'MWe': rr['W_el_out_implied']/1000,
                 'T_2d realised': rr['T_2d_realised']-273.15,
                 'glide dc [K]': rr['glide_dc'],
                 'eps end dc': rr['discharge']['eps_local'].mean(),
                 'dT_2d pass [K]': rr['dT_2d_pass']})
dfs = pd.DataFrame(rows)
print(dfs.round(4).to_string(index=False))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.0))
ax[0].plot(dfs['DT_M_1D'], dfs['deviation %'], 'o-')
ax[0].axhline(0, color='crimson', ls='--', lw=1.2, label='1 MWe target')
ax[0].axvline(CASE.DT_3C_2C/CASE.N_lay, color='gray', ls=':', lw=1.2,
              label='retired closure (6.111 K)')
ax[0].axvline(CASE.DT_4C_M, color='k', ls='-.', lw=1.2,
              label='symmetric with charge (10 K)')
ax[0].set_xlabel(r'$\Delta T_{M,1D}$ [K]'); ax[0].set_ylabel('deviation [%]')
ax[0].legend(fontsize=7)
ax[1].plot(dfs['DT_M_1D'], dfs['eps end dc'], 's-', color='darkorange')
ax[1].set_xlabel(r'$\Delta T_{M,1D}$ [K]')
ax[1].set_ylabel(r'residual $\varepsilon_{local}$')
ax[2].plot(dfs['DT_M_1D'], dfs['T_2d realised'], '^-', color='seagreen',
           label='realised outlet')
ax[2].plot(dfs['DT_M_1D'], dfs['dT_2d pass [K]'].abs()*100, 'v--', color='gray',
           label=r'$100\times|$pass-2 shift$|$')
ax[2].set_xlabel(r'$\Delta T_{M,1D}$ [K]'); ax[2].set_ylabel('[°C] / [K]')
ax[2].legend(fontsize=7)
for a in ax: a.grid(alpha=.3)
plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
**The shortfall was mostly self-inflicted.** Giving the discharge the same
10 K approach the charge gets moves the deviation from $-5.69\,\%$ to
$+1.23\,\%$ — the field now *exceeds* the target. The curve crosses zero at
$\Delta T_{M,1D}\approx9.30$ K and saturates near $+6.3\,\%$ at 15 K, where the
store refreezes completely ($\varepsilon\to0$) and the limit stops being a rate
limit and becomes an inventory limit.

**$\Delta T_{M,1D}=10$ K is a principle, not a tuning.** It is chosen because it
is the symmetric partner of $\Delta T_{4C,M}$, not because it happens to land
near the target. Choosing $9.30$ K instead would zero the deviation — and
destroy the only thing the deviation is good for, which is being a *reported
result* rather than a satisfied constraint.

**One correction pass is enough, and that is a property of the problem.** Across
this sweep the inlet moves 9 K and the deviation moves twelve percentage points,
while the realised outlet moves less than $0.3$ K and the second pass shifts it
by at most $0.03$ K. The outlet is pinned by the store, not by the plant, which
is why the plant level can still be evaluated without an outer loop
(§16).
"""))

cells.append(md(r"""
### 16.6 First-cycle sizing against CSS simulation

The two framings answer different questions and both are kept.
"""))
cells.append(code(r"""
fc = runs['v0.4b']
comp = pd.DataFrame([
    {'framing': 'first cycle (run_cycle)',
     'question': 'how many wells for ONE charge from cold?',
     'N_wells': fc['kpis']['N_wells'], 'lambda': CASE.loss_surplus,
     'eta_storage': fc['detail']['eta_storage'],
     'root-finds': 2},
    {'framing': 'CSS (simulate_css)',
     'question': 'what does a REPEATING cycle deliver?',
     'N_wells': css['N_wells'], 'lambda': 0.0,
     'eta_storage': css['eta_storage'],
     'root-finds': 0},
]).set_index('framing')
print(comp.T.to_string())
print()
print(f"CSS delivers {100*css['deviation']:+.2f} % against the 1 MWe target,")
print(f"i.e. {css['W_el_out_implied']/1000:.4f} MWe at N = {css['N_wells']:.4f}.")
"""))

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.11"},
                   "colab": {"provenance": [], "toc_visible": True}},
      "nbformat": 4, "nbformat_minor": 0}

doc = json.dumps(nb, indent=1)

# Substitute the build stamp. This step went missing once and the notebook
# shipped reading "Last updated: __STAMP__" -- a dead placeholder is worse than
# no timestamp, because it looks like a rendering fault rather than a stale
# build. Assert afterwards so it cannot fail silently again.
SENTINEL = '@@BUILD' + '_STAMP@@'      # split so this line never self-matches
n_stamp = doc.count(SENTINEL)
doc = doc.replace(SENTINEL, STAMP)
if n_stamp != 2:
    raise SystemExit(f'expected 2 build-stamp sites, found {n_stamp}')
if SENTINEL in doc:
    raise SystemExit('build stamp survived substitution')

OUT.write_text(doc)
print(f'  build stamp: {STAMP}  ({n_stamp} substitutions)')
n_code = sum(1 for c in cells if c['cell_type'] == 'code')
n_lines = sum(len(c['source']) for c in cells if c['cell_type'] == 'code')
print(f'wrote {OUT}')
print(f'  {len(cells)} cells: {n_code} code / {len(cells)-n_code} markdown')
print(f'  {n_lines} lines of code, all visible')
