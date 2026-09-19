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

*Model version 0.6 · notebook built {STAMP}*

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
                   performance_indices, kpi_report, calculate_pressure_drop)
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

Read this before sweeping anything. Several parameters do more than one job,
and the second job is usually the one that surprises you.

```
                 ┌─→ both mass flows          (ṁ = Q̇ / cp / ΔT_3C,2C)
ΔT_3C,2C ────────┼─→ cascade spacing          (= ΔT_3C,2C / N_lay)
  the glide      └─→ T_m,bottom               → discharge inlet → everything

                 ┌─→ cascade spacing
N_lay ───────────┴─→ T_m,bottom               → discharge inlet

ΔT_4C,M ───────────→ charge inlet only
ΔT_M,1D ───────────→ discharge inlet only     ← the two clean knobs

T_m ───────────────→ the whole cascade, AND the HP condensing temperature
                      (T_2h = T_4c + ΔT_2H,3C) → COP

geometry (r_e, fins, ────→ U_i → NTU → K      → rate only, not inventory
 L_well, D_well)     └───→ V_well             → inventory only, via N_wells

h_m, ρ ──────────────────→ inventory (N_wells ∝ 1/ρ h_m)
k_s, k_l ────────────────→ rate, through the melt-layer resistance
c_p,s, c_p,l ────────────→ how much energy the sensible branches carry
```

**The trap.** Change the glide expecting a mass-flow effect and you have also
moved the melting-temperature ladder and the discharge inlet. If you want to
vary *only* the flow, vary the glide and hold $T_m^{\rm bot}-\Delta T_{M,1D}$
fixed by adjusting $\Delta T_{M,1D}$ to compensate.

**The two clean knobs** are $\Delta T_{4C,M}$ and $\Delta T_{M,1D}$: each moves
exactly one inlet and nothing else.
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

    # --- operation -----------------------------------------------------------
    DT_3C_2C = 55.0,         # water glide: sets BOTH flows AND the cascade [K]
    DT_4C_M  = 10.0,         # charge inlet above the hottest layer         [K]
    DT_M_1D  = 10.0,         # discharge inlet below the coldest layer      [K]
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

validate_case(case)
"""))

cells.append(md(r"""
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
css = simulate_css_corrected(case, record=True, n_cycles=250)
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
def sweep(field, values, base=None, N=None, correct_once=True):
    '''Run one parameter over a list of values and collect the results.

    Returns a DataFrame with one row per value. Add whatever else you need to
    the `row` dict below -- everything in the result of simulate_css is
    available.
    '''
    base = base or case
    T_2d_fixed = None
    if correct_once:
        T_2d_fixed = simulate_css_corrected(base, N=N, n_cycles=250)['T_2d_realised']
    rows = []
    for v in values:
        c = base.with_(**{field: v})
        bad = [m for lvl, m in validate_case(c, verbose=False) if lvl == 'error']
        if bad:
            print(f'  {field}={v}: SKIPPED -- {bad[0]}');  continue
        r = simulate_css(c, N=N, T_2d=T_2d_fixed, n_cycles=250)
        rows.append({field: v,
                     'N_wells':      r['N_wells'],
                     'deviation_%':  100*r['deviation'],
                     'MWe':          r['W_el_out_implied']/1000.0,
                     'glide_dc_K':   r['glide_dc'],
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
5. **Idealised cycles.** The expansions and compressions carry no isentropic
   efficiency, so COP and $\eta_{\rm ORC}$ are optimistic.
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
