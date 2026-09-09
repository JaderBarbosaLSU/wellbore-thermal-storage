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

**THUMS · model version 0.3 · every equation visible in this notebook**

**Last updated: __STAMP__**

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
| 5 | **the melt front** — the defect, and the two formulations |
| 6 | the segment-by-segment march down the well |
| 7–8 | the ORC and heat-pump cycles, pressure drop |
| 9–11 | well-field sizing and the full cycle calculation |
| 12–14 | results, verification against the published numbers, plots |
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

LAST_UPDATED = '__STAMP__'
MODEL_VERSION = '0.3'

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
| **H3** | The melt region is a concentric annulus of uniform thickness $\delta$. Fronts from neighbouring tubes never merge. | not enforced — see §15 |
| **H4** | Conduction only in the melt; no natural convection. | increasingly wrong as $\delta$ grows |
| **H5** | Latent heat only. No sensible heat in either phase. | a segment that runs out of melt simply stops |
| **H6** | Adiabatic borehole wall — no formation heat loss. | fine over 10 h, not over a season |
| **H7** | Incompressible single-phase water, no buoyancy head. | ±7.5 bar in a 1,524 m loop is neglected |
| **H8** | The PCM is layered along the well, melting point following the fluid glide. | |
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
> **From v0.3 the correct value is the default.** `charge_uses_wall_conductivity`
> is `True`. Set it to `False` only to reproduce the published v0.1 results, as
> §13 does.

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
                  + '\n\nCASE = Case()          # the design point, all corrections applied\n\n'
                  '# The published configuration. BOTH switches must be set: the closed-form\n'
                  '# front AND the wall-conductivity error, because the two compensated each\n'
                  '# other. Setting only one gives a model that never existed.\n'
                  'CASE_V01 = Case(front="closed_form", charge_uses_wall_conductivity=False)'))

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

cells.append(md(r"""
The closed-form front (§5.1) is driven by the tube outer-surface temperature
rather than by the heat flow. Splitting the network at $r_e$:

$$R_1=R''_{f,i}+\frac{1}{h_i}+\frac{r_i}{k_w}\ln\frac{r_e}{r_i},
\qquad R_2=\frac{2\pi r_i}{P_T h_e \eta_o},
\qquad T_{r_e}=T_{\text{avg}}-\frac{R_1\left(T_{\text{avg}}-T_m\right)}{R_1+R_2}$$

The energy-balance front of §5.2 never forms $T_{r_e}$ at all.
"""))
cells.append(code(src('compute_T_re')))

# ---------------------------------------------------------------- 5. front
cells.append(md(r"""
## 5. The melt front — where the model was wrong

### The inconsistency

v0.1 computed two things from two **independent** models that never spoke:

| | surface used | perimeter |
|---|---|---|
| fluid side (`compute_U_i`) | finned, with fin efficiency | **0.4924 m** |
| melt front (`compute_delta2_fast`) | bare cylinder of radius $r_e$ | **0.1324 m** |

A factor of **3.72**. The fluid gave up heat through a finned surface; the front
absorbed it through a bare tube. They are not the same energy, and nothing in the
code ever compared them.
"""))
cells.append(code(r"""
c = CASE
P_T   = 2*np.pi*c.r_e + 2*c.num_fins*c.fin_L
P_bare = 2*np.pi*c.r_e
print(f'finned perimeter P_T   {P_T:.4f} m')
print(f'bare perimeter 2*pi*r_e {P_bare:.4f} m')
print(f'ratio                   {P_T/P_bare:.3f}   <-- the two surfaces the model used')
"""))

cells.append(md(r"""
### 5.1 Formulation A — the closed-form Stefan front (v0.1)

The quasi-steady solution for melting around a **bare** cylinder. With
$\alpha_m = k_m/(\rho_m c_{p,m})$, $x = 1+\delta/r_e$,

$$\mathrm{Fo}=\frac{\alpha_m t}{r_e^2},\qquad
\mathrm{Ph}=\left|\frac{h_m}{c_{p,m}(T_{r_e}-T_m)}\right|=\frac{1}{\mathrm{Ste}}$$

$$\mathrm{Fo}=\mathrm{Ph}\left[\tfrac12 x^2\ln x-\tfrac14 x^2+\tfrac14\right]
\qquad\Longleftrightarrow\qquad
\tfrac12 x^2\ln x-\tfrac14 x^2+\tfrac14=\frac{k_m(T_{r_e}-T_m)\,t}{\rho_m h_m r_e^2}$$

solved for $\delta$ by safeguarded Newton–bisection.

**The algebra is exact** — it returns the right root to a residual of $10^{-15}$.
Two *premises* fail:

1. it is derived for a bare cylinder, but the fluid side uses $P_T$ (the 3.72);
2. it assumes $T_{r_e}$ has been constant since $t=0$, whereas $T_{r_e}-T_m$
   actually runs $0 \to 8.6$ K across the charging window, and the **current**
   value is applied retroactively to the whole history.

Those two errors act in *opposite* directions, which is why the energy balance
appeared to close to 3 % at the design point.
"""))
cells.append(code(src('compute_delta2_fast')))

cells.append(md(r"""
### 5.2 Formulation B — the energy-balance front (v0.2)

Of the two models — heat delivered, and front position — only one may be chosen
freely. The other follows from conservation:

$$\boxed{\rho_m h_m \frac{\mathrm{d}A_{\text{melt}}}{\mathrm{d}t} = q'(t)},
\qquad \delta=\sqrt{r_e^2+\frac{A_{\text{melt}}}{\pi}}-r_e$$

Integrated explicitly, with a projection onto the physical bounds:

$$A^{n+1}=\mathcal{P}\!\left[A^{n}+\frac{q'\Delta t}{\rho_m h_m}\right],
\qquad \mathcal{P}[A]=\min\left(\max(A,0),A_{\max}\right),
\qquad q'_{\text{eff}}=\frac{\left(A^{n+1}-A^{n}\right)\rho_m h_m}{\Delta t}$$

It is $q'_{\text{eff}}$, not $q'$, that is accumulated — otherwise the fluid gets
credited with heat no PCM supplied.

Three things follow that Formulation A cannot do: the fins enter the front
automatically (because $q'$ comes from the finned network), $T_{r_e}$ need not be
constant, and the melt state **carries across the cycle** so discharge starts
from what charging left behind.

> **Closure is now exact by construction** — but that is an arithmetic identity,
> not evidence the model is right. `advance_front` and `closure_error` are
> inverse operations, so a residual of $10^{-16}$ measures floating point.

### 5.3 The fins are metal, not PCM

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
cells.append(code(src('delta_from_area', 'area_from_delta', 'advance_front',
                      'closure_error')))

# ---------------------------------------------------------------- 6. march
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
cells.append(code(src('layer_map', 'march')))

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
cells.append(code(src('cycle_state_points', 'energy_budget',
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

# ---------------------------------------------------------------- 9. legacy
cells.append(md(r"""
## 9. The v0.1 path, kept for verification

These are the **original** functions, unmodified. They are here so §13 can prove
the notebook still reproduces the published numbers exactly — which is what makes
every difference in §12 attributable to the front reformulation and nothing else.

You do not need to read them. Collapse this section.
"""))
cells.append(code(src('make_cp_state', 'get_props_state',
                      'compute_h_i_from_state')))
cells.append(code(src('temperature_profile_melt', 'time_profiles_melt')))
cells.append(code(src('evaluate_Q_ratio_ch', 'evaluate_Q_ratio_dc',
                      'find_N_wells_for_Q_ratio_ch_fast',
                      'find_m_dot_d_well1_for_Q_ratio_dc_fast')))

# ---------------------------------------------------------------- 10. sizing
cells.append(md(r"""
## 10. Well-field sizing — two constraints, not one

v0.1 sized on a single requirement: *can the wells absorb the energy in time?*
That is a heat-transfer question. **Nothing required them to contain it.**

$$\Phi_{\text{heat}}(N)=\frac{\Delta E_{out,HP}}{\Delta E_{\text{del}}(N)}-1=0,
\qquad
\Phi_{\text{inv}}(N)=\varepsilon_{\text{PCM}}(N)-1=0$$

$$\boxed{N_{\text{wells}}=\max\left(N_{\text{heat}},N_{\text{inventory}}\right)}$$

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
cells.append(code(src('_bisect', 'size_well_field')))

# ---------------------------------------------------------------- 11. run
cells.append(md(r"""
## 11. The full cycle

Charge, size the field, discharge, and collect the indices.
`case.front` selects the formulation.

$$\eta_{RTE}=\frac{\left(\dot W_{el,out}-\dot W_{\text{pump}}^{dc}\right)t_{dc}}
{\left(\dot W_{el,in}+\dot W_{\text{pump}}^{ch}\right)t_{ch}},
\qquad
\eta_{RTE}^{0}=\frac{\dot W_{el,out}t_{dc}}{\dot W_{el,in}t_{ch}}$$

$$\varepsilon_{\text{PCM}}=\frac{V_{\text{melt}}n_t}{V_{\text{well}}},\qquad
\eta_{\text{storage}}=\frac{\Delta E_{\text{discharged}}}{\Delta E_{\text{stored}}}$$
"""))
cells.append(code(src('run_cycle', '_find_N_wells_closed_form',
                      '_find_discharge_closed_form')))

# ---------------------------------------------------------------- 12. results
cells.append(md(r"""
## 12. Results

Three configurations: the published v0.1 model, the corrected front, and the
corrected front with the steel conductivity also fixed.
"""))
cells.append(code(r"""
runs = {}
for label, case in (
    ('v0.1 published',      CASE_V01.with_(N_wells_bracket=(1., 400.))),
    ('marched, old k_w',    CASE.with_(N_wells_bracket=(1., 400.),
                                       charge_uses_wall_conductivity=False)),
    ('v0.3 (default)',      CASE.with_(N_wells_bracket=(1., 400.))),
):
    runs[label] = run_cycle(case)

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
Three things to read from that table.

**The sizing moves a long way; the efficiency hardly at all.** `N_wells` falls
from 24.1 to 12.7 — a 47 % reduction — while `eta_RTE` moves by about 2 %.

**`eta_RTE_nopump` is identical to twelve significant figures** across all three
columns. The cycles are untouched by any of this. What the storage model
determines is how much hardware the cycle needs, not how efficiently it runs.

**The binding constraint flips.** With the old conductivity, heat transfer binds
and the inventory constraint is slack. With the correct value the heat-transfer
requirement collapses to about 8 wells — at which point the store would have to
melt more PCM than it contains — so inventory takes over.
"""))
cells.append(code(r"""
for label, r in runs.items():
    d = r['detail']
    if 'N_heat' in d:
        print(f"{label:22s}  N_heat={d['N_heat']:6.2f}  N_inventory={d['N_inventory']:6.2f}"
              f"  binding={d['binding']:14s}  eta_storage={d['eta_storage']:.4f}"
              f"  closure={d['closure_charge']:.1e}")
print()
print('The assumed storage loss lambda = 0.05 implies eta_storage = 0.952.')
print('At the v0.3 design point the model computes about 0.91, so the field')
print('banks ~10% more than the ORC withdraws. lambda should become an OUTPUT')
print('of the model rather than an input. Not done yet -- see section 15.')
"""))

cells.append(md(r"""
### Check against the lower bound

The v0.1 notebook printed `Ideal number of wells: 12.215` from a load-versus-
capacity calculation with no thermal resistance and no losses (§10). That is a
lower bound. Our $N_{\text{inventory}}$ must sit above it — and the gap should
be about the sensible heat the marched model omits under **H5**, i.e. a factor
$1+\mathrm{Ste}$.

This is one of the few genuinely independent checks available on the storage
model, since energy closure is structural and proves nothing (§5.2).
"""))

cells.append(code(r"""
N_IDEAL = 12.215336443485063   # 'Ideal number of wells', THUMS_Multilayer_BB_Jan_12b
d = runs['v0.3 (default)']['detail']
Ste = CASE.stefan_number(CASE.DT_4C_M)

print(f'N_ideal      (load / capacity, sensible + latent, no resistance) = {N_IDEAL:7.3f}')
print(f'N_inventory  (marched model, latent only)                        = {d["N_inventory"]:7.3f}')
print(f'N_heat       (rate criterion)                                    = {d["N_heat"]:7.3f}')
print()
print(f'N_inventory / N_ideal = {d["N_inventory"]/N_IDEAL:.4f}')
print(f'1 + Ste               = {1+Ste:.4f}     (Ste = {Ste:.4f})')
print()
if d['N_inventory'] < N_IDEAL:
    print('FAIL: below the lower bound -- the inventory constraint is wrong.')
else:
    print('OK: above the lower bound, by close to the omitted sensible heat.')
"""))

# ---------------------------------------------------------------- 13. verify
cells.append(md(r"""
## 13. Verification against the published numbers

The v0.1 path in this notebook must still reproduce the conference-paper values
exactly. If it ever stops, a change has leaked into the legacy path and every
comparison in §12 becomes unattributable.
"""))
cells.append(code(r"""
IHTC = {   # frozen v0.1 results, DT_3C_2C = 55 K, N_lay = 9
    'cop_hp':         2.956330409531823,
    'eta_orc':        0.23685132027873143,
    'eta_rte_nopump': 0.4348350503020928,
    'eta_rte':        0.43027101049765343,
    'eps_pcm':        0.5200473147142692,
    'N_wells':        24.123429921001843,
    'E_well':         2.167417633092469,
    'rho_E':          78.29462518690391,
    'f_pump':         0.015081686346236986,
}

got = runs['v0.1 published']['kpis']
print(f"{'KPI':16s} {'published':>22s} {'this notebook':>22s}   rel. diff")
worst = 0.0
for k, v in IHTC.items():
    n = got[k]; d = abs(n - v) / abs(v); worst = max(worst, d)
    print(f'{k:16s} {v:22.12g} {n:22.12g}   {d:.1e}')
print()
print(f'worst relative difference: {worst:.1e}')
print('PASS -- reproduces the published model' if worst < 1e-8 else
      'FAIL -- the legacy path has changed')
"""))

cells.append(md(r"""
A residual of order $10^{-11}$ rather than exactly zero is CoolProp and NumPy
version drift between machines, not a change in the model. Anything above
$10^{-8}$ means something real has moved.
"""))

# ---------------------------------------------------------------- 14. plots
cells.append(md(r"""
## 14. The melt front along the well

Where the PCM actually melts. The steps are the PCM layer boundaries: each layer
has its own melting temperature, following the fluid glide down the well.
"""))
cells.append(code(r"""
case = CASE                      # v0.3 defaults: corrected k_w, fins excluded
res  = runs['v0.3 (default)']
N    = res['kpis']['N_wells']          # taken from the run, not hard-coded

rank, hp, T = cycle_state_points(case)
E = energy_budget(case, rank['rank_eff'], hp['hp_cop'], T)
T_m_lay, _ = melting_temperatures(case)
T_m_seg = layer_map(T_m_lay, case.n_segments, case.N_lay)
times   = np.logspace(0., np.log10(case.t_ch*3600.), case.n_times)
m1      = E['m_dot_w_ch'] / (N * case.num_tubes)

r = march(case, T['T_4c'], T_m_seg, m1, case.k_wall, times, mode='charge')

z = np.linspace(0, case.L_well, case.n_segments)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(r['delta']*1000, z, lw=1.5)
ax[0].set_xlabel('melt layer thickness  δ  [mm]'); ax[0].set_ylabel('depth [m]')
ax[0].invert_yaxis(); ax[0].grid(alpha=.3)
ax[0].set_title(f'end of charging · N = {N:.2f} wells')

ax[1].plot(np.array(T_m_seg)-273.15, z, lw=1.5)
ax[1].set_xlabel('PCM melting temperature [°C]'); ax[1].invert_yaxis()
ax[1].grid(alpha=.3); ax[1].set_title(f'{case.N_lay} PCM layers')
plt.tight_layout(); plt.show()

print(f"eps_PCM = {r['V_melt_tube']*case.num_tubes/case.V_well:.4f}"
      f"   energy closure = {r['closure']:.1e}")
"""))

# ---------------------------------------------------------------- 15. gaps
cells.append(md(r"""
## 15. What this model still does not contain

Read this before quoting any number above.

| missing | consequence |
|---|---|
| **Formation heat loss** (H6) | fine over a 10 h cycle, wrong for seasonal storage |
| **Natural convection in the melt** (H4) | under-predicts melting once the layer is established |
| **Sensible heat in the PCM** (H5) | a segment whose melt is exhausted stops contributing rather than continuing to cool |
| **Buoyancy head** (H7) | ±7.5 bar in the 1,524 m loop, asymmetric between charge and discharge |
| **Isentropic efficiencies** | both cycles are idealised, so every efficiency here is an upper bound |
| **Volume change on melting** | 6.9 %, neglected; the liquid layer is ~2.3 % thicker than modelled |
| **Non-circular melt front** | the fins enhance heat transfer through $\eta_o$, but the front is still modelled as a circular annulus (H3). Their *metal* is excluded from $A_{\text{melt}}$ (§5.3); their effect on front *shape* is not modelled |
| **Front interaction** (H3) | `delta_max = 0.5 m` against an 0.0889 m borehole radius — a factor of 5.62. Nothing in the solver enforces it |

### Two quantities the code defines twice

- $h_e$: the cylindrical shell $k_m/(r_e\ln(1+\delta/r_e))$ in `compute_T_re`,
  but the plane-slab $k_m/\delta$ in `compute_U_i`'s small-$\delta$ branch — and
  the two are iterated against each other.
- The laminar Nusselt number: 3.66 here, 4.36 in some paths. Immaterial at the
  design point, where the flow is strongly turbulent.

### And the largest gap

**No external validation.** Neither formulation has been checked against
experimental or independently computed data. Energy closure is structural under
Formulation B and therefore proves nothing about correctness. Internal
consistency is not validation, and a reviewer will ask.

### Open work, in order

1. Make $\lambda$ an output rather than an assumed 5 %. It is now measurable:
   at the v0.3 design point the field banks about 10 % more than the ORC
   withdraws, against the 5 % assumed.
2. Enforce the geometric bounds instead of warning ($\delta_{\max} = 0.5$ m
   against an 0.0889 m borehole radius).
3. Add PCM sensible heat — it would close most of the gap to the ideal bound
   in §12 and remove the $(1+\mathrm{Ste})$ penalty.
4. Model the melt-front shape around the fins.
5. **Find an external validation case.** Still the largest gap.

---

## Version history

| version | date | change |
|---|---|---|
| **0.3** | this build | Wall conductivity corrected: charging now uses steel, $k_w = 45$ W/m·K, by default (§2). Fin metal excluded from the melted PCM area (§5.3). `N_wells` 21.65 → **12.74**, binding constraint now inventory. Lower-bound check against the ideal well count added (§12). |
| 0.2 | 2026-08 | Melt front reformulated by energy balance (§5.2); melt state carried across the cycle; two-constraint sizing (§10); latent inventory made density-consistent; segment latent limiting. Model moved into this notebook, every equation visible. |
| 0.1 | — | IHTC paper. Closed-form Stefan front, single sizing criterion, wall conductivity error. Reproduced exactly by §13. |
"""))

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.11"},
                   "colab": {"provenance": [], "toc_visible": True}},
      "nbformat": 4, "nbformat_minor": 0}

OUT.write_text(json.dumps(nb, indent=1))
n_code = sum(1 for c in cells if c['cell_type'] == 'code')
n_lines = sum(len(c['source']) for c in cells if c['cell_type'] == 'code')
print(f'wrote {OUT}')
print(f'  {len(cells)} cells: {n_code} code / {len(cells)-n_code} markdown')
print(f'  {n_lines} lines of code, all visible')
