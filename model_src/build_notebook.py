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
| 12 | results — and **12.1, do the fins earn their place?** |
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
MODEL_VERSION = '0.4'

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
### 5.4 Enthalpy state — sensible heat in both phases (v0.4)

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
cells.append(code(src('march_h')))

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
cells.append(code(src('_bisect', 'well_capacity_kJ', 'size_well_field',
                      'sizing_report')))

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
    ('v0.4 (default)',      CASE.with_(N_wells_bracket=(1., 400.))),
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
        sizing_report(CASE, d, label)
        print()
"""))

cells.append(code(r"""
for label, r in runs.items():
    d = r['detail']
    if 'N_heat' in d:
        print(f"{label:22s}  eta_storage={d['eta_storage']:.4f}"
              f"  closure={d['closure_charge']:.1e}")
print()
print('The assumed storage loss lambda = 0.05 implies eta_storage = 0.952.')
print('At the v0.3 design point the model computes about 0.91, so the field')
print('banks ~10% more than the ORC withdraws. lambda should become an OUTPUT')
print('of the model rather than an input. Not done yet -- see section 15.')
"""))

cells.append(md(r"""
### Check against the v0.1 ideal well count

The v0.1 notebook printed `Ideal number of wells: 12.215` from

$$m_{\rm ideal}=\frac{\Delta E_{\rm out,HP}}{h_m+c_{p,l}(T_{3c}-T_m)},
\qquad m_{\rm well}=V_{\rm well}\,\rho_l$$

which is exactly the capacity criterion of §10 — but evaluated with $\rho_l$,
whereas the model uses $\rho_s$ (see `Case.rho_latent`).

> **A correction to what this notebook said at v0.3.** It previously reported
> $N_{\rm inventory}/N_{\rm ideal}=1.0428$ as agreeing with $1+\mathrm{Ste}$, and
> attributed the gap to omitted sensible heat. That compared **inconsistent
> bases** and the agreement was a coincidence: two offsetting differences, a
> factor 0.955 from the sensible credit and 0.936 from the density.
>
> Under v0.4 the question dissolves — `N_capacity` *is* that formula on the
> model's own conventions, so the comparison below is an identity, not a test.
"""))

cells.append(code(r"""
d = runs['v0.4 (default)']['detail']
dT = cycle_state_points(CASE)[2]['T_3c'] - CASE.T_m
E = energy_budget(CASE, *[cycle_state_points(CASE)[k][n] for k, n in
                          ((0,'rank_eff'), (1,'hp_cop'))], cycle_state_points(CASE)[2])
D = E['D_E_out_HP'] * 1000.0          # J

print(f"{'variant':52s} {'N':>9s}")
for lbl, cap, rho in (
    ("v0.1 as coded:  h_m + cp_l*dT,  rho_l",  CASE.h_m + CASE.cp_l*dT, CASE.rho_l),
    ("latent only,                    rho_s",  CASE.h_m,                CASE.rho_s),
    ("h_m + cp_l*dT,  rho_s  <- v0.4 convention", CASE.h_m + CASE.cp_l*dT, CASE.rho_s),
):
    print(f"  {lbl:50s} {D/(cap*rho*CASE.V_well):9.3f}")
print(f"\n  model N_capacity (should equal the last row)  {d['N_capacity']:9.3f}")
print(f"  sensible credit  cp_l*dT/h_m = {CASE.cp_l*dT/CASE.h_m:.4f}  (= Ste)")
print(f"  density factor   rho_l/rho_s = {CASE.rho_l/CASE.rho_s:.4f}")

ident = abs(d['N_capacity'] - D/((CASE.h_m+CASE.cp_l*dT)*CASE.rho_s*CASE.V_well))
print('\nIDENTITY HOLDS' if ident < 1e-6 else f'MISMATCH {ident:.2e}')
"""))

cells.append(md(r"""
## 12.1 Do the fins earn their place?

A worked example of the trap described in §10, and a real design result.

Fins exist because the PCM conducts badly: $k_l = 0.45$ W/m·K against
$45$ W/m·K for the steel they are made of. But fins are metal, and metal in the
borehole is PCM *not* in the borehole. So they trade **capacity** for **rate**,
and the two sizing criteria pull in opposite directions.

Read the `N_wells` column. Reading either criterion alone gives the wrong sign.
"""))

cells.append(code(r"""
def sweep(param_sets, bracket=(3., 400.)):
    out = []
    for label, kw in param_sets:
        c = CASE.with_(N_wells_bracket=bracket, ratio_bracket=(0.05, 8.0), **kw)
        r = run_cycle(c); d = r['detail']
        out.append({'config': label, 'V_well_m3': c.V_well,
                    'N_capacity': d['N_capacity'], 'N_rate': d['N_heat'],
                    'N_wells': d['N_wells'], 'binding': d['binding'],
                    'eps_PCM': d['eps_pcm'],
                    'delta/delta_merge': d['merge_proximity_max']})
    return pd.DataFrame(out).set_index('config')

fins = sweep([(f'{nf:2d} fins x {1000*fL:.2f} mm', dict(num_fins=nf, fin_L=fL))
              for nf, fL in ((24, 0.0075), (24, 0.00375), (24, 0.00075),
                             (16, 0.0075), (12, 0.0075), (8, 0.0075),
                             (0, 1e-9))])
print(fins.round(4).to_string())
print()
print('N_capacity * V_well is CONSTANT -- that column is volume bookkeeping:')
print((fins['N_capacity']*fins['V_well_m3']).round(6).to_string())
"""))

cells.append(md(r"""
The `N_capacity` column falls monotonically as fin metal is removed, and
`N_capacity * V_well` is constant to machine precision — proof that the column
carries no heat transfer whatsoever. The `N_wells` column rises by 42 %.

**The optimum is where the two criteria cross**, because that is where neither
resource is wasted: 16 fins × 7.5 mm gives 11.348, marginally better than the
24 × 7.5 mm design point at 11.573. So the instinct to remove fin metal is
directionally right — worth about 2 %, not the 5 % the bound suggests — and it
reverses sharply past the crossing.
"""))

cells.append(code(r"""
tab = []
for t in (4., 6., 8., 10., 14., 20.):
    row = {'t_ch [h]': t}
    for lbl, kw in (('bare', dict(num_fins=0, fin_L=1e-9)), ('24 fins', {})):
        c = CASE.with_(t_ch=t, N_wells_bracket=(3., 400.),
                       ratio_bracket=(0.05, 8.0), **kw)
        row[lbl] = run_cycle(c)['detail']['N_wells']
    row['fin gain'] = row['bare']/row['24 fins']
    tab.append(row)
print(pd.DataFrame(tab).set_index('t_ch [h]').round(3).to_string())
"""))

cells.append(md(r"""
### Why the answer depends on the charging window

Fins buy **rate**, so their value is set by how much rate you need. Break-even
is around 20 h: give the store long enough and a bare tube melts the same PCM
with none of the volume penalty.

The quasi-steady Stefan solution says where this comes from. Putting
$x = r_{\rm cell}/r_e$ in the closed form of §5.1, a **bare** tube at
$\Delta T = 10$ K reaches the cell boundary in

$$t=\frac{\rho h_m r_e^2}{k_l\Delta T}
\left[\tfrac12 x^2\ln x-\tfrac14 x^2+\tfrac14\right]\approx 6.4\ \text{h}$$

comfortably inside the 10 h window — which is exactly why bare *nearly* works
here, and why the fins look marginal. Halve the window and they are not
marginal at all.

> **Two reasons to distrust this result in the direction of more fins, not
> fewer.** Both are limitations of the model, listed in §15.
>
> 1. **H3 is being sat on.** §14 shows $\delta/\delta_{\rm merge}=1.000$ over
>    the whole well at the finned design point. The annular resistance
>    $\ln(1+\delta/r_e)/2\pi k$ assumes the full azimuth conducts; past merge it
>    does not, and the last solid sits in the corners between legs on a longer,
>    narrower path. Late-stage resistance is understated.
> 2. **Fins enter $U_i$ only as added surface at the tube** — through $\eta_f$
>    and $P_T$. They are *not* represented as radial conduction paths reaching
>    into that corner PCM, which is the single thing a fin does best in this
>    regime. The model cannot credit it.
>
> Neither is a reason to dismiss the sweep. Both are reasons to treat 16 fins as
> a floor rather than an optimum, and to settle it with a 2-D conduction
> calculation before committing to a tube design.
"""))

# ============================================================ 13. verification
cells.append(md(r"""
## 13. Verification against the published numbers

The v0.1 path in this notebook must still reproduce the conference-paper values
exactly. If it ever stops, a change has leaked into the legacy path and every
comparison in §12 becomes unattributable.

> This cell is the reason §9 keeps the original functions verbatim. It is the
> only thing standing between "the model changed because we corrected it" and
> "the model changed and we do not know why."
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
"""))

cells.append(code(r"""
# --- run charge and discharge at the design point, recording history --------
case = CASE
res  = runs['v0.4 (default)']
N    = res['kpis']['N_wells']

rank, hp, T = cycle_state_points(case)
E = energy_budget(case, rank['rank_eff'], hp['hp_cop'], T)
T_m_lay, T_m_lay_dc = melting_temperatures(case)
n = case.n_segments
T_m_seg    = layer_map(T_m_lay,    n, case.N_lay)
T_m_seg_dc = layer_map(T_m_lay_dc, n, case.N_lay)
t_ch = np.logspace(0, np.log10(case.t_ch*3600), case.n_times)
t_dc = np.logspace(0, np.log10(case.t_dc*3600), case.n_times)
m1   = E['m_dot_w_ch']/(N*case.num_tubes)

ch = march_h(case, T['T_4c'], T_m_seg, m1, case.k_wall, t_ch, record=True)
ratio = res['detail']['flow_ratio_dc_ch']
dc = march_h(case, T['T_3d'], T_m_seg_dc, ratio*m1, case.k_wall, t_dc,
             E0=ch['E'], record=True)

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
picks = pick_times(ch['t'], targets)
cmap = plt.cm.viridis(np.linspace(0.15, 0.95, len(picks)))

fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))

for c_, k in zip(cmap, picks):
    lab = tlabel(ch['t'][k])
    for col, series in ((0, h['T_fluid'][k][1:]-273.15),
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

The cell after the discharge figure makes that concrete. Here the discharge flow
ratio is only $1.023$, so the two comparisons very nearly agree — mid-well at
end of run, $U_i$ rises by $1.401$ from charge to discharge while NTU rises by
$1.386$, and $1.401/1.023 = 1.370$ recovers the NTU ratio to within the
variation of $c_p$. The distinction is small at this design point but is not
guaranteed to stay small, which is the reason to plot the quantity that does not
depend on the flow.

$U_i$ is drawn on a log axis, though it spans only about a factor of seven
($117$ to $835$ W m⁻² K⁻¹ over both half-cycles): the log scale keeps
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
picks = pick_times(dc['t'], targets)          # dc grid, same target times
fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
for c_, k in zip(cmap, picks):
    lab = tlabel(dc['t'][k])
    for col, series in ((0, hd['T_fluid'][k][1:]-273.15),
                        (1, hd['A_melt'][k]/A_avail),
                        (2, hd['U_i'][k]),
                        (3, hd['NTU'][k])):
        (zd, vd), (zu, vu) = legs(series)
        ax[col].plot(vd, zd, color=c_, lw=1.4, ls='-',
                     label=(lab if col == 0 else None))
        ax[col].plot(vu, zu, color=c_, lw=1.4, ls=':')
(zd, vd), (zu, vu) = legs(np.array(T_m_seg_dc)-273.15)
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
for a, hh, tt, ttl in ((ax[0], h, ch['t'], 'charging'),
                       (ax[1], hd, dc['t'], 'discharging')):
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

The history arrays are written from inside the segment loop. This confirms the
recorded run and the plain run agree exactly, so nothing in §14 is an artefact
of instrumenting the march.
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
| **0.4a** | this build | Sizing report rewritten to lead with `N_wells` and name the binding criterion; `N_capacity` labelled a lower bound (§10). Melt-front limit moved from the borehole wall to the cell radius via `Case.r_cell` / `Case.delta_merge`, and H3 proximity reported (§12, §14). New §12.1: fin sweep, charging-window sweep, and why reading `N_capacity` alone inverts the answer. §13 (verification against the frozen IHTC fixture) restored — it had been dropped from the generator. **Build stamp repaired**: the `__STAMP__` placeholder was never substituted, so the notebook shipped reading `Last updated: __STAMP__`; the generator now asserts the substitution happened. §14 gains a **$U_i$ panel** beside NTU, and the plotted time levels are now chosen by target time rather than by index — on a logarithmic grid the old picks put three of six curves inside the first 12 seconds. |
| 0.4 | 2026-09 | Melt front carries **enthalpy** rather than melted area (§5.4): PCM superheats when fully molten and subcools when fully solid, $\varepsilon_{\rm local}\in[0,1]$ by construction, branchwise-exact time integration. Sensible heat is self-levelling — melt-fraction spread falls $1.321\to0.078$. |
| 0.3 | 2026-09 | Wall conductivity corrected: charging now uses steel, $k_w = 45$ W/m·K, by default (§2). Fin metal excluded from the melted PCM area (§5.3). `N_wells` 21.65 → **12.74**, binding constraint now inventory. Lower-bound check against the ideal well count added (§12). |
| 0.2 | 2026-08 | Melt front reformulated by energy balance (§5.2); melt state carried across the cycle; two-constraint sizing (§10); latent inventory made density-consistent; segment latent limiting. Model moved into this notebook, every equation visible. |
| 0.1 | — | IHTC paper. Closed-form Stefan front, single sizing criterion, wall conductivity error. Reproduced exactly by §13. |
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
