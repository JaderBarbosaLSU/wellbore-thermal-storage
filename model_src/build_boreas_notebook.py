"""Build notebooks/BOREAS_model.ipynb.

The model arrives as one `%%writefile` cell so the notebook is self-contained
in Colab, exactly as P2H2P_model.ipynb is. That one is NOT replaced and does
not change: BOREAS is a new series, and the whole point of starting one was
that the hairpin model stays runnable.
"""
import datetime as _dt
import json
import pathlib

STAMP = (_dt.datetime.utcnow() - _dt.timedelta(hours=3)).strftime(
    '%Y-%m-%d %H:%M BRT')
MODEL = pathlib.Path('/tmp/build/boreas.py')
# boreas_plots.py is a source file, not a build product -- it is hand-written
# and lives next to this script, unlike boreas.py which is assembled.
PLOTS = pathlib.Path(__file__).parent / 'boreas_plots.py'
OUT = pathlib.Path('/tmp/build/BOREAS_model.ipynb')


def md(t):
    return {'cell_type': 'markdown', 'metadata': {},
            'source': t.strip('\n').splitlines(True)}


def code(t):
    return {'cell_type': 'code', 'metadata': {}, 'execution_count': None,
            'outputs': [], 'source': t.strip('\n').splitlines(True)}


cells = []

cells.append(md(rf"""
# BOREAS — borehole array storage with a conducting formation

Latent-heat storage in repurposed oil-and-gas wells. **A new series, not a new
version of P2H2P/THUMS** — that model is frozen at v0.16 and still runs; its
notebook is untouched.

*v1.4 · built {STAMP}*

---

## What changed, and why the two changes are one change

### 1. The geometry

THUMS put $n$ **hairpins** in the hole — U-tubes, down one leg and back up the
other, both exchanging with the PCM. That turns out to be self-contradictory,
and the contradiction is quantitative.

Charge water enters at $T_m + 10$ K and glides 55 K down, so the **return leg
runs about 45 K below the melting point**. Graded by depth it would refreeze
whatever the down leg had just melted. So THUMS graded the cascade along the
**tube path** instead — which works, but then two different melting points
have to coexist at the same depth, centimeters apart, in contiguous PCM with
no wall between them. Conduction across that gap is

| hairpins | leg-to-leg short circuit | % of a ~120 kW duty |
|---|---|---|
| 1 | 43 kW | 36 % |
| 2 | 62 kW | 52 % |
| 3 | **92 kW** | **77 %** |

a short circuit straight across the cascade the model exists to represent —
and a term THUMS never had.

**BOREAS separates the two jobs.** `n_ft` fin tubes all run *downward* and all
exchange; one **insulated downcomer** returns the whole flow without
exchanging. They meet in a manifold at the bottom of the same hole — no
lateral, no new drilling, the same class of object as the U-bend it replaces.

### 2. The formation

THUMS hypothesis **H6** made the borehole wall adiabatic, so $\eta_{{storage}}$
was identically 1 — not approximately, *exactly*, by construction. BOREAS
conducts to the rock.

### Why the second needed the first

With a hairpin, **every depth is occupied by two tube segments**, so one rock
temperature would have had two PCM segments claiming it. Only when the return
stops exchanging does $z$ become a depth. In BOREAS `L_tube = L_well`, not
twice it, and the map is one-to-one.

Three things follow for free:

- the cascade runs **hot at the bottom**, which *aligns* it with the
  geothermal gradient instead of fighting it
- the PCM becomes **pourable in layers** (path-graded PCM was never
  manufacturable)
- the downcomer's leak goes **into** the store on charge — energy misplaced,
  not energy lost
"""))

# ------------------------------------------------------------------ setup
cells.append(md("""
---
## 1. Setup

`CoolProp` is not in a stock Colab image. The model is written to disk by the
next cell, so this notebook needs nothing from the repository.
"""))

cells.append(code("!pip install -q CoolProp"))
cells.append(code('%%writefile boreas.py\n' + MODEL.read_text()))
GROUND = pathlib.Path(__file__).parent / 'boreas_ground.py'
cells.append(code('%%writefile boreas_ground.py\n' + GROUND.read_text()))
cells.append(md(r"""
The slot coordinates come from the CalGEM WellSTAR export. They are fetched
rather than embedded: the file is the field, and a copy pasted into a
notebook is a copy that goes stale.
"""))
cells.append(code(r"""
import os, urllib.request
CSV = 'THUMS_wells_WellSTAR_clean.csv'
if not os.path.exists(CSV):
    urllib.request.urlretrieve(
        'https://raw.githubusercontent.com/JaderBarbosaLSU/'
        'wellbore-thermal-storage/notebook-model/field_ladder/' + CSV, CSV)
import boreas_ground as G
print('slot data ready')
"""))

cells.append(code(r"""
import importlib, sys
for _m in ('boreas',):
    sys.modules.pop(_m, None)     # %%writefile changed the file on disk
import boreas as B
import numpy as np, matplotlib.pyplot as plt, pandas as pd
print(B.__doc__.strip().splitlines()[0])
c = B.Case()
print(f'\ndefault case: {c.n_ft} fin tubes + 1 downcomer, '
      f'{c.L_ft:.0f} ft, T_m = {c.T_m_C:.0f} C')
"""))

# --------------------------------------------------------------- geometry
### --------------------------------------------- THE DESIGN POINT
cells.append(md(r"""
---
## 2. The reference case (R4) — **edit here, nowhere else**

Every parameter you are likely to change lives in the next cell. Nothing
below it hard-codes a number: each later section starts from `CASE` and
says, in one visible `.with_(...)`, what it is changing and why.

The split is deliberate:

* **`CASE`** is the well, the material and the plant — what you are
  building. It says nothing about how many of them there are or how they
  are arranged on the surface.
* **`FIELD`** is the arrangement — which island, which trajectory, how
  deep the wells kick off. It is kept separate because the *same* well
  gives losses from 3 % to 69 % depending only on this, and because the
  model's own defaults describe an **isolated** well, which is a bound and
  not a design.

`Case` has 83 fields. The ones written out below are those worth touching;
the rest keep their defaults, and the audit cell that follows prints
everything you have changed so nothing hides.
"""))

cells.append(code(r"""
# =============================================================================
#  THE DESIGN POINT.  Edit the numbers HERE -- not in the model cell, and not
#  in the sections below.  Units are in front of every value.
# =============================================================================
CASE = B.Case(
    # --- PCM and cascade -----------------------------------------------------
    T_m_C      = 150.0,      # melting temperature of the HOTTEST layer     [C]
                             #   the cascade is graded with the geotherm, so
                             #   this layer sits at the BOTTOM of the well
    h_m        = 380_000.0,  # latent heat of fusion                     [J/kg]
    k_s        = 0.60,       # solid conductivity                       [W/m/K]
    k_l        = 0.45,       # liquid conductivity                      [W/m/K]
    N_lay      = 9,          # cascade layers along the well                [-]
    DT_3C_2C   = 55.0,       # cascade grading; SPAN = this*(N_lay-1)/N_lay [K]
                             #   -> 48.889 K.  DT_cascade overrides it if set;
                             #   setting DT_cascade EQUAL to this is a no-op.

    # --- the well ------------------------------------------------------------
    #   NOTE  the storage length is NOT here. It is L_md in FIELD below,
    #   in meters of MEASURED depth, because that is what the trajectory
    #   and the casing shoe are quoted in. Setting L_ft here as well is an
    #   error and the guard below says so rather than silently winning.
    n_ft       = 3,          # fin tubes per well, all downward             [-]
    num_fins   = 24,         # fins per tube                                [-]
    fin_L      = 0.0035,     # fin length                                   [m]
    fin_t      = 0.0015,     # fin thickness                                [m]
    t_ins      = 0.008,      # DOWNCOMER insulation (the casing is bare)    [m]
    k_ins      = 0.007,      # downcomer insulation, VIT                [W/m/K]

    # --- casing: THUMS upper zone, 8-5/8 in 32# (Truex & Hunter, SPE-2562)
    #   D_in is the casing INSIDE diameter and is what bounds the PCM. The
    #   201.2 mm ID carries 46 % more material per tube than the 177.8 mm
    #   legacy string, which is worth +36 % on e' -- a larger effect than
    #   the choice of ground model. Set D_in = 7.0, r_casing_o = 0.098 and
    #   r_bore = 0.110 to recover the legacy reference.
    D_in       = 201.2 / 25.4,  # 8-5/8 in 32# inside diameter             [in]
    r_casing_o = 0.21910 / 2,   # 219.1 mm outside diameter                 [m]
    r_bore     = 0.31120 / 2,   # 12-1/4 in hole                            [m]
    k_cement   = 0.7,           # annulus fill: mud above the TOC       [W/m/K]
                                #   cement (1.0) instead moves the loss by
                                #   0.09 pts and e' by 0.1 %, so the top of
                                #   cement is not worth chasing

    # --- formation -----------------------------------------------------------
    T_surface_C   = 18.0,    # surface temperature                          [C]
    grad_K_per_km = 45.0,    # geothermal gradient                       [K/km]
    k_rock        = 2.5,     # rock conductivity                        [W/m/K]
    rho_cp_rock   = 2.34e6,  # rock volumetric heat capacity          [J/m3/K]
    t_operation_yr = 10.0,   # field age at which R_g is evaluated          [yr]

    # --- plant ---------------------------------------------------------------
    W_dot_el_out = 1000.0,   # net electrical output                      [kWe]
    t_ch         = 10.0,     # charge duration                               [h]
    t_dc         = 10.0,     # discharge duration                            [h]
    T_source_C   = 60.0,     # geothermal source                             [C]
    T_sink_C     = 20.0,     # heat sink                                     [C]
    P            = 10.0e5,   # secondary-loop pressure                      [Pa]
                             #   water saturates at 179.9 C here, which CAPS
                             #   T_m near 165 C; see section 10.5
    DT_pinch_ORC = 5.0,      # ORC evaporator pinch                          [K]
    DT_2H_3C     = 10.0,     # HTHP condenser approach                       [K]

    # --- numerics ------------------------------------------------------------
    n_segments = 99,         # a MULTIPLE of N_lay, or a layer boundary falls
                             #   inside a segment and validate_case warns
)

# =============================================================================
#  THE FIELD.  Measured slot positions and real trajectories, not an
#  idealized arrangement. The old packing / n_islands / fan_angle_deg
#  description is retired: it could not represent where the wells actually
#  are, and at this spacing the difference is not small.
#
#    island      : slot coordinates from the WellSTAR export
#    trajectory  : 'vertical' or 's_curve' (build-and-hold, SPE-2562)
#    kop_rule    : 'staggered' outer rows kick off shallow, inner deep
#                  'fixed'     every well at kop_ft
#    L_md        : storage length in MEASURED depth, bounded by the shoe
# =============================================================================
FIELD = dict(
    island     = 'Grissom',   # 363 slots after dedup and stray removal   [-]
    trajectory = 's_curve',   # 'vertical' | 's_curve'
    kop_rule   = 'staggered', # 'staggered' | 'fixed'
    theta_max  = 55.0,        # maximum inclination                     [deg]
    kop_ft     = 300.0,       # used only when kop_rule = 'fixed'        [ft]
    L_md       = 1250.0,      # storage length in measured depth          [m]
    ds_station = 50.0,        # superposition station spacing             [m]
                              #   ONE spacing for every run of every sweep,
                              #   not a count and not a per-run search.
                              #   50 m is the shipped value: zero negative
                              #   strengths at every length from 100 to
                              #   2000 m, for R1, R3 and R4. Refining to
                              #   25 m introduces them while moving the
                              #   loss by 0.04 points, so the finer grid is
                              #   worse posed and not more accurate.
                              #   See section 5.2 and station_convergence.csv.
)

# The storage length reaches Case through L_ft, and ONLY from here.
REF = CASE.with_(L_ft=FIELD['L_md'] / 0.3048)
G.DS_STATION = FIELD['ds_station']

# ---- no silent overrides ----------------------------------------------
# If L_ft was also set in CASE above and disagrees with L_md, that is an
# error, not a warning: the student would otherwise edit L_ft, see nothing
# change, and have no way to find out why.
_L_ft_default = B.Case().L_ft
if abs(CASE.L_ft - _L_ft_default) > 1e-9 and \
        abs(CASE.L_ft * 0.3048 - FIELD['L_md']) > 1e-6:
    raise ValueError(
        f"two storage lengths: CASE.L_ft = {CASE.L_ft:.1f} ft "
        f"({CASE.L_ft * 0.3048:.1f} m) against FIELD['L_md'] = "
        f"{FIELD['L_md']:.1f} m. Set the length in FIELD['L_md'] only.")

# The superposition profile for that field. Built once here and reused by
# every later cell through `field_ground`, so no section can quietly fall
# back to the idealized geometry.
#
# Cached by specification, including the station spacing. The file for the
# default FIELD is fetched from the repository; editing FIELD above forces
# a fresh solve, which the cell says it is doing.
MD, RG, TVD, FINFO = G.cached_reference_field(
    B, REF, CSV, island=FIELD['island'], trajectory=FIELD['trajectory'],
    kop_rule=FIELD['kop_rule'], theta_max=FIELD['theta_max'],
    kop_ft=FIELD['kop_ft'], L_md=FIELD['L_md'], ds=FIELD['ds_station'])
SLOTS = G.grissom_positions(CSV)        # plan view, for the figures below
print('reference case (R4) set.  REF = CASE at the field length; '
      'MD/RG/TVD describe the field.')
"""))

cells.append(md(r"""
### 2.1 What you actually changed

`Case` carries 83 fields and the cell above writes out 30. This prints
every field of `REF` that differs from the model default, so nothing you
set — or forgot to set — is invisible.
"""))

cells.append(code(r"""
import dataclasses
base = B.Case()
diff = [(f.name, getattr(base, f.name), getattr(REF, f.name))
        for f in dataclasses.fields(base)
        if getattr(base, f.name) != getattr(REF, f.name)]
print(f'{len(diff)} field(s) differ from the model defaults:\n')
print(f"  {'field':18s}{'default':>14}  ->  {'yours':<14}")
for n, d, v in diff:
    print(f'  {n:18s}{str(d):>14}  ->  {str(v):<14}')

# No later cell may correct an input. Every .with_(...) below section 2
# is a VISIBLE sensitivity, never a fix; this prints what section 2 itself
# changed between CASE and REF so an override cannot hide.
_ref_vs_case = [(f.name, getattr(CASE, f.name), getattr(REF, f.name))
                for f in dataclasses.fields(CASE)
                if getattr(CASE, f.name) != getattr(REF, f.name)]
print(f'\nREF differs from CASE in {len(_ref_vs_case)} field(s) '
      f'(should be L_ft alone, set from FIELD["L_md"]):')
for n, a, b in _ref_vs_case:
    print(f'  {n:18s}{str(a):>14}  ->  {str(b):<14}')
if [n for n, _, _ in _ref_vs_case] not in ([], ['L_ft']):
    raise ValueError('section 2 changed something other than the storage '
                     'length between CASE and REF')

print(f'\nstation spacing {FIELD["ds_station"]:.0f} m -> '
      f'{G.Ground.stations_for(FIELD["L_md"], FIELD["ds_station"])} stations '
      f'at L = {FIELD["L_md"]:.0f} m')

print('\nvalidate_case on REF:')
msgs = B.validate_case(REF, verbose=False)
if not msgs:
    print('  (clean)')
for lvl, m in msgs:
    print(f'  [{lvl}] {m[:100]}')
"""))



cells.append(md(r"""
---
## 3. The hole, in cross-section

`n_ft` finned tubes and one insulated downcomer, sized so its velocity matches
theirs. That sizing is the whole reason this geometry costs nothing in
pumping: the downcomer carries the flow of all `n_ft` tubes through `n_ft`
times the area, so the velocity is the same everywhere and the total path is
still $2L_{well}$ — exactly what the hairpin had.

The downcomer runs **off-center**, as in the schematic. A concentric
arrangement does not close for three 1¼″ tubes in a 7″ hole; an offset one
packs comfortably. `validate_case` says so rather than leaving it to a comment.
"""))

cells.append(code(r"""
def pack(case):
    # Place the downcomer tangent to the wall, then fit the fin tubes round it.
    # SOLVED, not eyeballed: a hand-placed sketch had a fin tube overlapping
    # the downcomer by 30 mm, which is the kind of thing a figure will happily
    # show you forever. The constraint is only that no two circles intersect
    # and all sit inside the hole -- and it is checked below the plot, because
    # a schematic nobody verifies is a drawing, not a result.
    R = case.D_well / 2
    a, b = case.r_e + case.fin_L, case.r_id_outer
    d_id, d_ft = R - b, R - a                     # both tangent to the wall
    # smallest angle off the downcomer that clears it
    cosmin = (d_id**2 + d_ft**2 - (a + b)**2) / (2 * d_id * d_ft)
    th_min = np.arccos(np.clip(cosmin, -1, 1))
    span = 2 * np.pi - 2 * th_min
    n = case.n_ft
    if n == 1:
        th = np.array([np.pi])
    else:
        th = th_min + np.linspace(0, span, n)
    return (d_id, 0.0), [(d_ft * np.cos(t), d_ft * np.sin(t)) for t in th], a, b, R


(xid, yid), fts, a, b, R = pack(c)
fig, ax = plt.subplots(figsize=(6, 6))
ax.add_patch(plt.Circle((0, 0), R * 1000, fc='#2ecc71', ec='k', lw=2, zorder=1))
ax.add_patch(plt.Circle((xid * 1000, yid * 1000), b * 1000, fc='w',
                        ec='#1a5276', lw=6, zorder=2))
ax.text(xid * 1000, yid * 1000, 'ID', ha='center', va='center', fontsize=12, zorder=3)
for (x, y) in fts:
    ax.add_patch(plt.Circle((x * 1000, y * 1000), a * 1000, fc='w',
                            ec='#1a5276', lw=5, zorder=2))
    ax.text(x * 1000, y * 1000, 'FT', ha='center', va='center', fontsize=10, zorder=3)
L = R * 1100
ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect('equal'); ax.axis('off')
ax.set_title(f'{c.n_ft} fin tubes + 1 insulated downcomer\n'
             f'{2000*R:.0f} mm hole, PCM = {100*c.V_well/c.V_borehole:.0f} % of it',
             fontsize=12)
plt.show()

# the figure is only worth looking at if it is a real packing
gaps = [np.hypot(x - xid, y - yid) - a - b for (x, y) in fts]
print(f'  smallest downcomer-to-tube gap  {1000*min(gaps):+6.1f} mm')
if len(fts) > 1:
    pp = [np.hypot(fts[i][0]-fts[j][0], fts[i][1]-fts[j][1]) - 2*a
          for i in range(len(fts)) for j in range(i+1, len(fts))]
    print(f'  smallest tube-to-tube gap       {1000*min(pp):+6.1f} mm')
print(f'  fin tube      {2000*c.r_e:.1f} mm OD, {1000*c.fin_L:.1f} mm fins, '
      f'{c.num_fins} of them')
print(f'  downcomer     {2000*c.r_id_bore:.1f} mm bore, {2000*c.r_id_outer:.0f} mm over '
      f'{1000*c.t_ins:.0f} mm insulation (k = {c.k_ins} W/m/K)')
print(f'  r_cell        {1000*c.r_cell:.1f} mm   (melt fronts meet here)')
print(f'  L_tube        {c.L_tube:.0f} m = L_well  <-- not 2x, and that is the point')
"""))

# -------------------------------------------------------------- formation
cells.append(md(r"""
---
## 4. The formation

The thermal circuit from the PCM outward is a series chain, and the terms are
wildly unequal:

$$R_{tot} = \underbrace{R_{PCM}}_{\sim 25\%} + \underbrace{R_{casing}}_{0.1\%}
+ \underbrace{R_{cement}}_{5\%} + \underbrace{R_{ground}(t)}_{\sim 70\%}$$

The ground resistance uses the **exact** large-time Carslaw--Jaeger prefactor
$2e^{-\gamma/2} = 1.49861$, not the $2.5$ carried by v1.0 on an attribution
that could not be substantiated. The difference was a *constant* additive
offset of \SI{0.0326}{} m·K/W that never decayed with time; under the
`max()` rule it moved only the isolated-well bound and the sparse fan.

The casing is thermally invisible and the cement is a rounding error. **The
ground resistance is the model**, and unlike the others it *grows with time*
as the thermal front spreads — which is why `t_operation_yr` is a parameter
and why a quoted loss figure is meaningless without the age beside it.

### The gradient is independent, and only checked

`T_source_C` is a **lower bound** on the formation temperature — the water
cools coming up the production string, and the producing wells may be deeper
than the retrofits. An earlier draft derived the gradient from the source
temperature, which writes a bound as an equality and throws away exactly the
physics the deficit represents. Here both are specified and `validate_case`
enforces only the inequality.
"""))

cells.append(code(r"""
print(f'  R_PCM->casing   {np.log((c.D_well/2)/c.r_e)/(2*np.pi*c.n_ft*0.5*(c.k_s+c.k_l)):8.4f} m K/W')
print(f'  R_cement        {np.log(c.r_bore/c.r_casing_o)/(2*np.pi*c.k_cement):8.4f}')
print(f'  R_ground        {B.ground_resistance(c):8.4f}   <-- at t = {c.t_operation_yr:.0f} yr')
print(f'  R_total         {B.formation_resistance(c):8.4f}')

T_rock, depth = B.formation_profile(c)
T_m_lay, _ = B.melting_temperatures(c)
T_m_seg = B.layer_map(T_m_lay, c.n_segments, c.N_lay)
fig, ax = plt.subplots(figsize=(6.5, 6))
ax.plot(T_m_seg - 273.15, depth, lw=2.5, label='PCM cascade')
ax.plot(T_rock - 273.15, depth, lw=2.5, ls='--', label='undisturbed rock')
ax.fill_betweenx(depth, T_rock - 273.15, T_m_seg - 273.15,
                 where=(T_m_seg > T_rock), alpha=.18, color='#c0392b', label='loses')
ax.fill_betweenx(depth, T_rock - 273.15, T_m_seg - 273.15,
                 where=(T_m_seg <= T_rock), alpha=.3, color='#2471a3', label='GAINS')
ax.invert_yaxis(); ax.set_xlabel('temperature, C'); ax.set_ylabel('depth, m')
ax.set_title(f'cascade hot at the BOTTOM, aligned with the gradient\n'
             f'{c.grad_K_per_km:.0f} K/km from {c.T_surface_C:.0f} C', fontsize=11)
ax.legend(); ax.grid(alpha=.3); plt.show()
"""))

cells.append(md(r"""
### Why hot-at-the-bottom matters

The rock is coldest at the surface and hottest at depth. THUMS ran the cascade
**hot at the top**, so the hottest PCM faced the coldest rock and the coldest
PCM faced rock that — at 8000 ft and a Wilmington gradient — is *hotter than
it is*. Those deep layers never freeze, and that capacity is dead.

The net loss is identical either way (the mean of $T_{PCM}$ is the same), but
the **distribution** is not:

| at 8000 ft, 56.5 K/km | ΔT at top | ΔT at bottom | cannot freeze |
|---|---|---|---|
| hot-top (THUMS) | +132 K | −55 K | **29 %** |
| hot-bottom (BOREAS) | +83 K | −6 K | **7 %** |

This is not a free choice: the downcomer *makes* it hot-at-bottom, because
charge water reaches the bottom before it has given up any heat.
"""))

# ------------------------------------------------------- the field, measured
cells.append(md(r"""
---
## 5. The field — where the wells actually are

Until v1.2 the field was an abstraction: a well count, a wellhead spacing
and a cone half-angle, from which a perimeter and a shielding factor
followed. That description cannot represent **where the wells are**, and at
1.65 m spacing the difference is not a refinement — it is most of the
answer.

So the geometry is now read from the **CalGEM WellSTAR export**: 363 Grissom
slots, deduplicated and stripped of strays. Three facts about them, all
from the data rather than from a description:

* the measured median nearest-neighbor spacing is **1.65 m**, 10 % tighter
  than the nominal 6 ft. The idealized rungs below are built at the measured
  value, so that the step from them to the real slots isolates the *layout*;
* Grissom is a **thick closed ring**, not a filled disc and not a straight
  double row: no slot lies inside $0.3\,R_\mathrm{max}$ and 72 % lie beyond
  $0.7\,R_\mathrm{max}$, against 0.09 and 0.51 for a disc;
* the wells are **directionally drilled**, so measured depth is not true
  vertical depth. At $L = 1250$ m MD the S-curve reaches only about 890 m
  TVD, where the rock is 58 °C rather than the 73 °C a vertical well would
  see. The store is colder *and* better shielded than the vertical
  idealisation says.

### The ground model

Every well in the field is a line source. At a given MD, well $i$ sees its
own wall rise plus the rise its neighbors impose,

$$\Delta T_i = R_\mathrm{near}\,q_i + \sum_k G_{ik}\,q_k, \qquad
G_{ik} = \frac{1}{4\pi k}\,E_1\!\left(\frac{d_{ik}^{2}}{4\alpha t}\right)$$

with $d_{ik}$ the shortest distance from the target to well $k$'s
trajectory. Two ways to close this:

| model | assumption | what it gets wrong |
|---|---|---|
| `flux` | every well emits the **same** $q$ | edge wells are in fact hotter; it averages a temperature that is not uniform |
| `temp` | every well sits at the **same** $\Delta T$ | nothing; the strengths are solved for and edge wells come out stronger |

`temp` is the one used throughout, because the PCM in every well is at the
same melting temperature — equal temperature is the physical boundary
condition, not equal flux. Solving $(R_\mathrm{near}\mathbf{I} +
\mathbf{G})\,\mathbf{q} = \mathbf{1}$ gives
$R_{g,\mathrm{eff}} = 1/\overline{q} - R_\mathrm{near}$.

**The specification degenerates correctly.** With a single well, $G_{ii}$
must reduce to the isolated-well $R_g$ already in the model. It does, to
$8\times10^{-7}$ — which is the check worth having, because it ties the new
ground model to the old one rather than to a reference value.
"""))

cells.append(code(r"""
# the degeneracy check: one well, no neighbors
_c = REF
_g = G.Ground(k=_c.k_rock, rhocp=_c.rho_cp_rock, t_yr=_c.t_operation_yr,
              rb=_c.r_bore)
print(f'  G_ii (new superposition kernel)  {_g.kern(_c.r_bore):.9f} m K/W')
print(f'  R_g  (shipped isolated-well)     {B.ground_resistance(_c):.9f} m K/W')
print(f'  difference                       '
      f'{abs(_g.kern(_c.r_bore) - B.ground_resistance(_c)):.2e}')
"""))

cells.append(md(r"""
### 5.2 One station spacing, fixed, for every run

Stations are the measured depths at which the interaction matrix is
assembled. The method treats each as an infinite line source in the plane
normal to the trajectory, so **a station must represent a segment long
compared with the distance to the neighboring wells**.

Fix the *count* and that condition changes with length: 50 stations is 25 m
apart in a 1250 m hole and 2 m apart in a 100 m one, against a well spacing
of 1.65 m. The short case then returns *negative* strengths — wells
absorbing heat from rock they are hotter than — and $R_{g,\mathrm{eff}}$
goes negative with them.

So the input is a **spacing**, `FIELD['ds_station'] = 50 m`, and the count
follows as $\max(3, \mathrm{round}(L/\Delta s))$. One value for the
ladder, the $L$ sweep and every sensitivity. Three things follow:

* the validity condition is the same at every length;
* the stations sit at the same measured depths whatever $L$ is, so
  $R_{g,\mathrm{eff}}$ at a given MD can be *compared* across lengths —
  that is check V10, and it agrees to 0.00 % between 200 m and 2000 m;
* every run of a sweep is discretized identically, so a difference between
  two runs is the thing being varied and not the grid.

**V11 is an assertion, not a search.** If the fixed spacing produces a
negative strength the run fails and says so. It does not quietly try
another count — that is what made two runs of the same sweep incomparable
in v2.0.

### What the grid study actually says

Refining does *not* converge to a better answer here, and that is worth
seeing rather than assuming. At the reference case:

| stations | $\Delta s$ | negatives | off-station weight | loss | $e'$ |
|---|---|---|---|---|---|
| 25 | 50.0 m | **0** | 5.2 % | 52.39 % | 0.2153 |
| 50 | 25.0 m | 130 | 25.8 % | 52.43 % | 0.2151 |
| 100 | 12.5 m | 67 | 27.6 % | 52.51 % | 0.2147 |
| 200 | 6.25 m | 490 | 32.7 % | 52.68 % | 0.2142 |

The answer barely moves — 0.29 points of loss and 0.51 % of $e'$ across an
eightfold refinement — while the negative strengths appear from 25 m
spacing on. **The coarse grid is the well-posed one**; the fine ones are
cheaper to compute than to justify. 50 m gives zero negatives at every
length from 100 to 2000 m, for R1, R3 and R4.

### Two structural results that make most of this moot

Both measured, in the cell below, not assumed.

**Vertical fields are exactly block diagonal in station.** For parallel
vertical wells the nearest point on a neighbor is at the same depth, so
every coupling stays within its own station — the off-station kernel weight
is 0.00 % — and every block is the *same* matrix, because horizontal
distances do not change with depth. $R_{g,\mathrm{eff}}$ is therefore
independent of both the station count and $L$: R1 gives 15.783 m K/W and R3
gives 3.230, at every length and every count tested. The solver exploits
this directly, which is what makes the disc rung computable at all — at
1.65 m spacing all 363 wells sit inside 18 m, so every target has about
90 000 candidate samples inside the cutoff.

**The S-curve is block tridiagonal.** The wells converge and diverge, so a
neighbor's nearest point lands at most one station away, carrying 5.2 % of
the kernel weight at the reference. This is the only case in which the
count does anything at all.
"""))

cells.append(code(r"""
# the structure claim, measured rather than asserted
import numpy as np
_c, _P = REF, SLOTS
_kopm, _ = G.staggered_kop(_P)
_gr = G.Ground(k=_c.k_rock, rhocp=_c.rho_cp_rock, t_yr=_c.t_operation_yr,
               rb=_c.r_bore)
_Rn = float(np.atleast_1d(B.formation_resistance(_c, np.array([0.0])))[0]
            - B.ground_resistance(_c))
_nst = G.Ground.stations_for(FIELD['L_md'], FIELD['ds_station'])
print(f"ds = {FIELD['ds_station']:.0f} m -> {_nst} stations "
      f"at L = {FIELD['L_md']:.0f} m MD\n")
for _t, _k in (('vertical', None), ('s_curve', _kopm)):
    _md, _X = G.trajectories(_P, FIELD['L_md'], kind=_t, kop=_k,
                             theta_max=FIELD['theta_max'])
    _pr = _gr.profiles(_md, _X, _Rn, nst=_nst)
    print(f"  {_t:9s} off-station kernel weight "
          f"{100*_gr.off_station_weight:5.2f} %   "
          f"negative strengths {_gr.n_negative:5d}   "
          f"R_g,eff mean {_pr['temp'].mean():7.3f} m K/W")
print('\n0.00 % is the block-diagonal case: the count cannot matter there.')
"""))

cells.append(code(r"""
import numpy as np, matplotlib.pyplot as plt
P = SLOTS
kop_m, is_out = G.staggered_kop(P)

fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.0))
for grp, col, lab in ((is_out, '#b03028', 'outer row'),
                      (~is_out, '#23356b', 'inner row')):
    ax[0].plot(P[grp, 0], P[grp, 1], 'o', ms=2.2, color=col, label=lab)
ax[0].set_aspect('equal'); ax[0].legend(frameon=False, markerscale=3)
ax[0].set_xlabel('east / m'); ax[0].set_ylabel('north / m')
ax[0].set_title(f'(a) {len(P)} Grissom slots')

_md2, _X2 = G.trajectories(P, 2000.0, kind='s_curve', kop=kop_m,
                           theta_max=FIELD['theta_max'])
_rad = np.hypot(_X2[:, :, 0], _X2[:, :, 1])
for i in range(0, len(P), 7):
    ax[1].plot(_rad[i], _X2[i, :, 2], lw=0.4, color='#6a6a6a', alpha=.7)
ax[1].axhline(FINFO['tvd_bottom'], color='#b03028', lw=1.2, ls='--')
ax[1].annotate(f"L = {FIELD['L_md']:.0f} m MD", (0, FINFO['tvd_bottom']),
               color='#b03028', fontsize=9, va='bottom')
ax[1].invert_yaxis()
ax[1].set_xlabel('radius from centroid / m'); ax[1].set_ylabel('TVD / m')
ax[1].set_title('(b) S-curve trajectories')

from scipy.spatial import cKDTree
for kind, col in (('vertical', '#23356b'), ('s_curve', '#b03028')):
    m2, X2 = G.trajectories(P, 2000.0, kind=kind, kop=kop_m,
                            theta_max=FIELD['theta_max'])
    nn = [np.median(cKDTree(X2[:, j, :2]).query(X2[:, j, :2], k=2)[0][:, 1])
          for j in range(0, len(m2), 8)]
    ax[2].plot(m2[::8], nn, color=col, lw=1.8, label=kind)
ax[2].set_yscale('log'); ax[2].legend(frameon=False)
ax[2].set_xlabel('MD / m'); ax[2].set_ylabel('median spacing / m')
ax[2].set_title('(c) nearest-neighbor spacing')
plt.tight_layout(); plt.show()

print(f"measured median spacing at the wellhead: "
      f"{np.median(cKDTree(P).query(P, k=2)[0][:, 1]):.3f} m "
      f"(nominal 6 ft = {G.B0_NOMINAL:.2f} m)")
"""))

cells.append(code(r"""
# R_g,eff along the well: vertical against the real trajectory
_MDv, _RGv, _TVDv, _Iv = G.cached_reference_field(
    B, REF, CSV, trajectory='vertical', L_md=FIELD['L_md'], verbose=False)

fig, a = plt.subplots(figsize=(5.6, 4.4))
a.plot(_RGv, _TVDv, color='#23356b', lw=2.0, label='vertical')
a.plot(RG, TVD, color='#b03028', lw=2.0, label=f"S-curve, "
       rf"$\theta_\mathrm{{max}}$ = {FIELD['theta_max']:.0f}$^\circ$")
a.axvline(B.ground_resistance(REF), color='#6a6a6a', lw=1.2, ls='--')
a.annotate('isolated well', (B.ground_resistance(REF), 0), color='#6a6a6a',
           fontsize=9, rotation=90, ha='right', va='top')
a.invert_yaxis(); a.set_xscale('log')
a.set_xlabel(r'$R_{g,\mathrm{eff}}$ / m K W$^{-1}$'); a.set_ylabel('TVD / m')
a.legend(frameon=False)
a.set_title('Shielding along the well, equal-temperature superposition',
            fontsize=10)
plt.tight_layout(); plt.show()

print(f"vertical : mean {_Iv['Rg_mean']:6.2f} m K/W, "
      f"TVD at the shoe {_Iv['tvd_bottom']:.0f} m")
print(f"S-curve  : mean {FINFO['Rg_mean']:6.2f} m K/W, "
      f"TVD at the shoe {FINFO['tvd_bottom']:.0f} m")
print(f"isolated : {B.ground_resistance(REF):6.2f} m K/W")
print(f"the field is {_Iv['Rg_mean']/B.ground_resistance(REF):.0f}x more "
      f"resistive than a lone well when the wells are vertical, and only "
      f"{FINFO['Rg_mean']/B.ground_resistance(REF):.0f}x once they deviate: "
      f"the trajectory gives back most of the shielding.")
"""))

# ----------------------------------------------------------------- ladder
cells.append(md(r"""
---
## 6. Run it — the field ladder

One run does not tell you what a field is worth; the *difference* between
two runs does. So the result is reported as a **ladder**, in which each rung
changes exactly one thing and nothing else moves:

| rung | geometry | what the step from the rung above isolates |
|---|---|---|
| **R0v** | isolated well, vertical | the bound. No field at all |
| **R0i** | isolated well, S-curve | the **trajectory alone**: same MD, colder rock, no neighbors |
| **R1** | disc, vertical, 1.65 m | the **field**: 363 neighbors in the most favorable packing |
| **R2** | straight double row, vertical | the **layout**: same count and spacing, strung out instead of compact |
| **R3** | Grissom slots, vertical | the **real positions**: a ring, not an idealized row |
| **R4** | Grissom slots, S-curve | the **real trajectories**: the **reference case** |

Two rules make the comparison honest, and both were violated in an earlier
draft:

1. **Every rung is at the measured 1.65 m.** Building the idealized rungs at
   the nominal 6 ft mixes a layout change with a spacing change, and that
   understated the layout penalty by 40 % — R2 → R3 is $+4.2$ points, not
   $+3.0$.
2. **Every rung carries the same casing**, the THUMS upper-zone 8-5/8 in
   string of section 2. Its 201.2 mm ID holds 46 % more PCM per tube than the
   legacy 177.8 mm, worth **+36 % on $e'$** — a larger effect than the choice
   of ground model, and a reason not to compare across casings by accident.

The 111 runs behind this section take about two hours, so they are read
from a cached JSON. Set `RECOMPUTE = True` to rebuild any of them locally;
the cache is keyed by the full specification, so nothing is recomputed twice.
"""))

cells.append(code(r"""
import os, json, urllib.request
RECOMPUTE = False            # True -> run the ladder here instead of reading it

PACK = 'pack_runs.json'
if not os.path.exists(PACK):
    urllib.request.urlretrieve(
        'https://raw.githubusercontent.com/JaderBarbosaLSU/'
        'wellbore-thermal-storage/notebook-model/field_ladder/' + PACK, PACK)
RUNS = [r for r in json.load(open(PACK)).values() if 'error' not in r]
print(f'{len(RUNS)} runs in the cache')

def pick(step, L=None, ground='superposition', **kw):
    c = [r for r in RUNS if r['step'] == step and r['ground'] == ground
         and all(abs(r[k] - v) < 1e-9 if isinstance(v, float) else r[k] == v
                 for k, v in kw.items())]
    if L is not None:
        c = [r for r in c if abs(r['L_md'] - L) < 1e-6]
    return sorted(c, key=lambda r: r['L_md'])
"""))

cells.append(code(r"""
import pandas as pd
LADDER = ['R0v', 'R0i', 'R1', 'R2', 'R3', 'R4']
WHAT = {'R0v': 'isolated, vertical', 'R0i': 'isolated, S-curve',
        'R1': 'disc, vertical', 'R2': 'double row, vertical',
        'R3': 'Grissom slots, vertical (R3)', 'R4': 'Grissom slots, S-curve'}

tab, prev = {}, None
for s in LADDER:
    r = pick(s, FIELD['L_md'])
    if not r:
        continue
    r = r[0]
    tab[s] = {'geometry': WHAT[s], 'loss %': 100 * r['loss'],
              'step': '' if prev is None else f"{100*(r['loss']-prev):+.1f}",
              'eta_storage': r['eta_storage'], 'eta_RTE': r['RTE'],
              "e' [kWh/m]": r['e_prime'],
              'T_2d [C]': r['T2d_C'], 'TVD [m]': r['tvd_bottom'],
              'delivers': r['delivers']}
    prev = r['loss']
display(pd.DataFrame(tab).T)
print(f"field at the reference case (R4): "
      f"{pick('R4', FIELD['L_md'])[0]['field_MWh_per_cycle']:.1f} MWh "
      f"per cycle from {pick('R4', FIELD['L_md'])[0]['N_field']:.0f} wells")
"""))

cells.append(code(r"""
fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.0))
steps = [s for s in LADDER if pick(s, FIELD['L_md'])]
loss = [100 * pick(s, FIELD['L_md'])[0]['loss'] for s in steps]
ep = [pick(s, FIELD['L_md'])[0]['e_prime'] for s in steps]
x = np.arange(len(steps))
lab = [WHAT[s].replace(', ', ',\n') for s in steps]
ax[0].bar(x, loss, color=['#b03028' if v > 50 else '#23356b' for v in loss])
ax[1].bar(x, ep, color='#2f6f4f')
for a, v, f in ((ax[0], loss, '{:.1f}'), (ax[1], ep, '{:.3f}')):
    a.set_xticks(x); a.set_xticklabels(lab, fontsize=8)
    for i, vi in enumerate(v):
        a.text(i, vi + 0.02 * max(v), f.format(vi), ha='center', fontsize=8)
ax[0].set_ylabel('formation loss / %')
ax[1].set_ylabel(r"$e'$ / kWh per m of well")
fig.suptitle(f"The field ladder at L = {FIELD['L_md']:.0f} m MD, "
             f"Grissom, N = 363", y=1.0, fontsize=11)
plt.tight_layout(); plt.show()
"""))

cells.append(md(r"""
### Reading the ladder

**The field is worth an order of magnitude, and the trajectory takes a
third of it back.** R0 → R1 is the big step: 363 neighbors cut the loss
from 64 % to 3 %, because each well only has to push heat out through its
share of the cluster boundary. R3 → R4 then puts 39 points back on,
entirely through the trajectory: below the kickoff the wells separate, the
shielding weakens, and the store spends most of its length in rock it does
not share with anybody.

**Layout matters more than compactness suggests.** R1 → R2 is $+7$ points
for the same 363 wells at the same spacing — stringing them into a line
roughly triples the perimeter per well. R2 → R3 is a further $+3$: the real
ring is *better* than a straight double row of the same length, because a
closed loop has no end-caps, but worse than the idealized row's uniformity
allows for.

**Nothing here fails to deliver.** Earlier versions reported short isolated
wells that absorbed heat on discharge, and v1.3 and earlier hid that behind
an `abs()`. With the signed fix and the THUMS casing the regime is gone: an
isolated 100 m well goes from absorbing (81.9 °C, loss 115 %) to delivering
(99.2 °C, loss 86 %). Both changes were needed — the sign fix reports the
legacy case correctly, the casing removes the case.
"""))

cells.append(md(r"""
---
## 7. The gates, watched to fail

A check that has never been seen to fail is not evidence. Both of these were
corrupted deliberately, and the second exists *because* the first was silent.

**The closure cannot see a flipped loss sign.** Negating `q_loss` leaves the
residual at 5e-15 — `q_prime` and `Q_loss_cum` both carry the error and it
cancels algebraically. The closure reads only its own inputs.

So the loss term carries a **Clausius gate**: accumulate
$q_{loss}\,(T_{PCM}-T_{rock})$, which is $K(\Delta T)^2$ under the correct
convention and $-K(\Delta T)^2$ under a flipped one.
"""))

cells.append(code(r"""
import dataclasses, gc
src = open('boreas.py').read()
r_ok = B.simulate_css_corrected(REF)
print(f'{"(uncorrupted)":<28} closure {r_ok["charge"]["closure"]:.2e}   '
      f'clausius {r_ok["charge"]["loss_clausius"]:+.3e}')
for lab, bad in [
        ('drop q_loss from q_prime',
         src.replace('q_prime[i] = q_tot + q_loss', 'q_prime[i] = q_tot')),
        ('flip the loss sign',
         src.replace('q_loss = K_l * (T_eff - T_rock_arr[i])',
                     'q_loss = K_l * (T_rock_arr[i] - T_eff)'))]:
    assert bad != src, lab
    ns = {}
    exec(compile(bad, 'corrupt', 'exec'), ns)
    rb = ns['simulate_css_corrected'](ns['Case'](**dataclasses.asdict(REF)))
    cl, cz = rb['charge']['closure'], rb['charge']['loss_clausius']
    print(f'{lab:<28} closure {cl:.2e} {"FIRES" if cl > 1e-6 else "silent":>6}   '
          f'clausius {cz:+.3e} {"FIRES" if cz < 0 else "silent"}')
    # each corrupted copy is a second full module plus its recorded march.
    # Keeping both alive alongside the real one is about 1 GB, enough to be
    # killed on a small runtime, and nothing below needs them.
    del ns, rb
    gc.collect()
"""))

cells.append(md(r"""
### The adiabatic limit

Driving `k_rock` to zero must recover THUMS hypothesis H6 exactly: no loss,
$\eta_{storage} = 1$.

This gate found a real bug on its first run. `ground_resistance` returned
`0.0` when the line source fell outside its validity, meaning *"the front has
not cleared the hole"*. But zero **resistance** is infinite **conductance** —
it claims the rock is a perfect heat sink. The test returned
$\eta_{storage}=0.20$ and a loss **larger** than the physical case.
"""))

cells.append(code(r"""
r0 = B.simulate_css_corrected(REF.with_(k_rock=1e-9))
print(f'  k_rock -> 0 : Q_loss {r0["Q_loss_kJ"]:.3e} kJ, '
      f'eta_storage {r0["eta_storage"]:.9f}')
print(f'  {"PASS" if r0["Q_loss_kJ"] == 0.0 else "FAIL"}'
      f'   (residual is the CSS drift tolerance, not a loss)')
"""))

# --------------------------------------------------------- the depth question
cells.append(md(r"""
---
## 8. How deep should the storage go?

With a real trajectory the question changes shape. A vertical well gains
from every extra meter: the rock is hotter and the shielding does not
weaken, so $e'$ rises monotonically and there is no optimum — depth is a
pure capital lever. **A deviated well does not.** Past the kickoff each
extra meter of measured depth buys less true vertical depth and sits in
rock the neighbors no longer share, so the marginal return falls away.

The cutoff is defined on the **marginal** round-trip efficiency, not the
average: $\Delta E_\mathrm{out}/\Delta E_\mathrm{in}$ across each interval
in $L$. The storage stops paying where that marginal value first drops
below half its value **over the first 100 m**.

Which interval supplies that reference is not a detail. Because $E(0) = 0$,
the marginal over $0 \to 100$ m *is* $\eta_\mathrm{RTE}(100)$, and it is
larger than the marginal over $100 \to 200$ m. Using the $100 \to 200$
interval instead lowers the threshold from 0.131 to 0.122 and moves the
cutoff a full step, from 400--500 m to 500--600 m. Both are defensible
definitions; they are not the same definition, and the paper reports the
first.
"""))

cells.append(code(r"""
fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.0))
COL = dict(zip(LADDER, ('#6a6a6a', '#9c27a0', '#23356b', '#0a9396',
                        '#2f6f4f', '#b03028')))
cut = {}
for s in LADDER:
    rr = pick(s)
    if len(rr) < 2:
        continue
    L = np.array([r['L_md'] for r in rr])
    ax[0].plot(L, [r['e_prime'] for r in rr], 'o-', ms=3.5, lw=1.5,
               color=COL[s], label=f'{s}: {WHAT[s]}')
    ax[1].plot(L, [r['RTE'] for r in rr], 'o-', ms=3.5, lw=1.5, color=COL[s])
    Eo = np.array([r['E_out_kWh'] for r in rr])
    Ei = np.array([r['E_in_kWh'] for r in rr])
    # Marginals on the interior intervals, plus the first one, 0 -> L[0].
    # E(0) = 0, so the marginal over that first interval is just RTE(L[0]) --
    # and it is the reference the criterion is defined against. Taking
    # marg[0] (the 100 -> 200 m interval) instead is a different and smaller
    # number, and it moves the cutoff a full step, to 500-600 m.
    marg = np.diff(Eo) / np.where(np.diff(Ei) == 0, np.nan, np.diff(Ei))
    Lm = 0.5 * (L[1:] + L[:-1])
    ax[2].plot(Lm, marg, 'o-', ms=3.5, lw=1.5, color=COL[s])
    first = Eo[0] / Ei[0] if abs(L[0] - 100.0) < 1e-6 else marg[0]
    ax[2].plot([0.5 * L[0]], [first], 'o', ms=3.5, color=COL[s])
    thr = 0.5 * first
    below = np.where(marg < thr)[0]
    if len(below):
        cut[s] = (L[below[0]], L[below[0] + 1], first, thr)
        ax[2].plot(Lm[below[0]], marg[below[0]], '*', ms=15, color=COL[s])

for a, lab in zip(ax, (r"$e'$ / kWh m$^{-1}$", r'$\eta_\mathrm{RTE}$',
                       r'marginal $\Delta E_\mathrm{out}/\Delta E_\mathrm{in}$')):
    a.set_xlabel('storage length $L$ / m MD'); a.set_ylabel(lab)
    a.axvspan(FIELD['L_md'], 2000.0, color='#cccccc', alpha=.35, lw=0)
ax[0].legend(frameon=False, fontsize=7.5, ncol=1)
ax[2].set_ylim(bottom=0)
ax[1].annotate('beyond the provisional shoe',
               (FIELD['L_md'], ax[1].get_ylim()[1]), fontsize=8, color='#6a6a6a',
               rotation=90, ha='right', va='top')
plt.tight_layout(); plt.show()

print('cutoff -- marginal RTE first below half its value over the first')
print('100 m, which is RTE(100) because E(0) = 0:')
for s in LADDER:
    if s in cut:
        a_, b_, m0, thr = cut[s]
        print(f'  {s:4s} {WHAT[s]:26s} {a_:5.0f}-{b_:5.0f} m MD   '
              f'(first 100 m {m0:.3f}, threshold {thr:.3f})')
    elif len(pick(s)) > 1:
        print(f'  {s:4s} {WHAT[s]:26s} none within the sweep')
"""))

cells.append(md(r"""
**Vertical wells have no cutoff; deviated ones stop paying between 400 and
600 m MD**, and the reference case's own crossing falls in 400--500. That is
not a statement about heat transfer in a tube — it is a statement about the
kickoff. The S-curve reaches about 890 m TVD at 1250 m
MD, and most of the shortfall is accumulated in the build section, where
the wells are simultaneously diverging from each other and falling behind
in depth.

### How sharp is that number?

The marginal curve is flat where it crosses the threshold — 0.123 against
0.131 over 400--500 m — so the crossing is resolved to one interval of the
$L$ grid and no finer. Within that limit it is now **stable**: the
45$^\circ$ and 70$^\circ$ variants put it in the same 400--500 m interval
as the reference case.

They did not, in v2.0. There the station count was chosen per run, as the
first count that solved with no negative strengths, so two runs of the same
sweep were discretized differently — at $L = 500$ m the reference case
landed on 25 stations and the 45$^\circ$ variant on 6 — and that moved $e'$
by 2.8 % against the 0.4 % the inclination was worth. The cutoff came out
one interval later for both inclination variants, in an order that was not
monotonic in angle. **That was the grid, not physics**, and section 5.2
describes the fixed spacing that removed it.

Two things are solid underneath. Below the kickoff the three inclinations
are now *identical to four decimals*, as they must be — there is nothing to
deviate yet, and in v2.0 they were not. And a deep kickoff removes the
cutoff altogether: `R4_kop300` has none within the sweep.

The grey band marks measured depths beyond the provisional casing shoe.
**The shoe depth is the one genuinely blocking unknown in the geometry**:
the top of cement moves the loss by 0.09 points across its whole plausible
range, the mud conductivity by 0.17, the casing grade by 2.2 % of $e'$ —
all sub-2 %, because the near-well resistance they move is itself only
about 7 % of the total. The shoe is different in kind, because it is a
constraint on $L$ rather than a resistance.
"""))

# ------------------------------------------------------------ sensitivities
cells.append(md(r"""
---
## 8b. What the answer is sensitive to

Four one-at-a-time variations on the design rung R4, each re-solving the
superposition from scratch rather than rescaling anything.
"""))

cells.append(code(r"""
import pandas as pd
SENS = [('R4', 'reference case: staggered KOP, 55 deg, 10 yr'),
        ('R4_kop300', 'every well kicks off at 300 ft'),
        ('R4_th45', 'theta_max = 45 deg'),
        ('R4_th70', 'theta_max = 70 deg'),
        ('R4_top1', 'field age 1 yr'),
        ('R4_top30', 'field age 30 yr')]
base = pick('R4', FIELD['L_md'])[0]
rows, missing = {}, []
for step, what in SENS:
    r = pick(step, FIELD['L_md'])
    if not r:
        missing.append(f'{step} ({what})')
        continue
    r = r[0]
    rows[step] = {'what': what, 'loss %': 100 * r['loss'],
                  'd(loss) [pts]': 100 * (r['loss'] - base['loss']),
                  "e' [kWh/m]": r['e_prime'],
                  "d(e') %": 100 * (r['e_prime'] / base['e_prime'] - 1),
                  'eta_RTE': r['RTE'], 'TVD [m]': r['tvd_bottom']}
display(pd.DataFrame(rows).T.round(3))
if missing:
    print('not in this cache yet, so not shown rather than shown as blank:')
    for m in missing:
        print('  ' + m)
"""))

cells.append(md(r"""
**The kickoff is the largest single lever, and the field age is the widest
spread.** Moving every well to a 300 ft kickoff takes the loss from 52 % to
35 % — more than either inclination change — which is the mechanism behind
the cutoff: a deep kickoff keeps the store inside the tightly spaced
section. Inclination at 45° and 70° gives 50 % and 55 %, a narrow band
either side.

Field age runs 72 % at one year to 39 % at thirty, the widest of anything
tested. That is not a design variable, and it cuts the right way: the
thermal field matures around the store and the loss falls with it. It also
means **any loss quoted without an age is meaningless** — everything here
is at ten years.

### Two ground models, one geometry

The `flux` rows in the cache are the equal-strength alternative at the same
geometry. They are kept because the difference between them is the size of
the modeling choice, and that belongs in the limitations rather than in a
footnote.
"""))

cells.append(code(r"""
rows = {}
for s in ('R1', 'R2', 'R3', 'R4'):
    t = pick(s, FIELD['L_md'], ground='superposition')
    f = pick(s, FIELD['L_md'], ground='flux')
    if t and f:
        rows[s] = {'geometry': WHAT[s],
                   'temp: loss %': 100 * t[0]['loss'],
                   'flux: loss %': 100 * f[0]['loss'],
                   'difference [pts]': 100 * (f[0]['loss'] - t[0]['loss']),
                   "temp: e'": t[0]['e_prime'], "flux: e'": f[0]['e_prime']}
display(pd.DataFrame(rows).T.round(3))
print('equal flux is the optimistic reading: it averages a wall temperature')
print('that is not uniform, and the edge wells it under-weights are the ones')
print('doing the losing.')

# V1/V7: the two legacy-casing rows exist only to regress against the
# prototype, which was computed at the 177.8 mm string. They are NOT a
# result -- they are the check that this implementation reproduces the one
# it replaced once the casing is set back.
print()
for s_, want in (('R3_legacy', 17.4), ('R4_legacy', 70.7)):
    r = pick(s_, FIELD['L_md'])
    if r:
        print(f'{s_:10s} {100*r[0]["loss"]:5.1f} % against the prototype\'s '
              f'{want:4.1f} %')
"""))

cells.append(md(r"""
---
## 9. Change something

Every parameter is a field of `Case`. The ones BOREAS adds:

| | |
|---|---|
| `n_ft` | fin tubes per hole (the downcomer is always one more) |
| `t_ins`, `k_ins` | downcomer insulation |
| `T_surface_C`, `grad_K_per_km` | the geothermal gradient |
| `k_rock`, `rho_cp_rock`, `k_cement` | formation and annulus |
| `t_operation_yr` | how old the field is when you evaluate it |
| `D_in`, `r_casing_o`, `r_bore` | the casing string; `D_in` is what bounds the PCM |

**Nothing below section 2 may correct an input.** Every `.with_(...)` in
the sections above is a *sensitivity* --- a visible, commented variation
on the reference case --- and never a fix for something set in section 2.
The audit cell enforces the one case that used to break this rule: the
storage length, which section 2 set as `L_ft = 5000` and then silently
replaced with 1250 m. There is now one length input, `FIELD['L_md']`, and
setting `L_ft` as well raises.

The field is **not** a `Case` field any more. `n_wells_field`,
`B_surface_m`, `fan_angle_deg` and `packing` still exist for backward
compatibility with v1.2 results, but nothing in this notebook uses them:
the arrangement lives in `FIELD` and reaches the model through
`boreas_ground.field_ground`, which patches `formation_profile` and
`formation_resistance` for the duration of a run. To change the field,
edit `FIELD` in section 2 and re-run — not a `Case` field.

Always run `validate_case` first — it is fast, and it catches the gradient
contradicting the source temperature, a cascade that has been flipped
hot-to-top, a geometry that will not pack, and a well whose deep layers cannot
freeze.
"""))

### ------------------------------------------------- 8. cycle diagrams
cells.append(md(r"""
---
## 10. The surface cycles: $T$--$s$ diagrams and state points

The borehole model and the plant models meet at the water loop, and the only
quantity the store sets is $T_\mathrm{2d}$, the temperature of the water
arriving at the ORC evaporator. Everything else in the two cycles below is
frozen from the hairpin study. So these diagrams are *context*, not a result:
they say what machine the store is feeding, and they make visible how much of
the store's temperature the ORC actually uses.

Three things in these figures are deliberate and worth stating, because each
is a place where a $T$--$s$ plot is usually drawn wrongly:

1. **Isobars are sampled in enthalpy, not temperature.** Inside the dome $T$
   and $p$ are not independent, so a state cannot be recovered from $(T, p)$
   there and a temperature-sampled isobar breaks silently at the saturation
   line --- which is exactly where every heat-exchanger leg lives. Enthalpy is
   monotonic along an isobar through subcooled liquid, dome and superheat
   alike.
2. **Throttles are drawn as isenthalps**, sampled in pressure, not as the
   straight lines a sketch usually shows. The $7\mathrm{h} \to 8\mathrm{h}$
   expansion is strongly curved in $T$--$s$ and the curvature *is* the
   irreversibility.
3. **The connectivity is read from the model's enthalpies, not from the
   numbering.** The ORC is a *reheat* cycle ($3\mathrm{e} \to 5\mathrm{e}$,
   reheat at $p_\mathrm{int}$ to $6\mathrm{e}$, then $\to 4\mathrm{e}$); the
   HTHP is a two-stage *flash-economizer* cycle. Neither is apparent from the
   state list alone.

### A note on fonts

Times New Roman is proprietary and is **not** installed on Colab. Asking
`matplotlib` for it by name there does not fail --- it silently falls back to
DejaVu Sans and the figure comes out in the wrong face. The style helper asks
for Times first and then for its metric-compatible clones (Nimbus Roman,
Liberation Serif), which are present on virtually every Linux image and are
indistinguishable from Times at figure scale. `active_serif()` reports which
one actually resolved, so the fallback is visible rather than silent. On a Mac
with Times installed, you get the real thing.
"""))

cells.append(code('%%writefile boreas_plots.py\n' + PLOTS.read_text()))

cells.append(code(r"""
import importlib, boreas_plots as BP
importlib.reload(BP)
BP.use_paper_style()
print('serif resolved to:', BP.active_serif()[0])

ref = REF

# THE CYCLES ARE DRAWN AT THE OUTLET THE GRISSOM FIELD ACTUALLY DELIVERS.
# cycle_state_points' own T_2d is provisional -- the value the outlet WOULD
# take if the top layer stayed fully molten -- and at the reference case the
# field misses it by nearly 30 K. Drawing the T-s diagram at the provisional
# value while the tables report the realized one is how an earlier version
# ended up with a figure and a state-point table that were not the same case.
import pickle, os
if os.path.exists('rec_field.pkl'):
    REC = pickle.load(open('rec_field.pkl', 'rb'))['rec']
else:
    with G.field_ground(B, MD, RG, TVD):          # ~1 min
        REC = B.simulate_css_corrected(REF, record=True)
    pickle.dump(dict(rec=REC), open('rec_field.pkl', 'wb'))

_prov = B.cycle_state_points(ref)[2]['T_2d'] - 273.15
rank, hp, Tdict = B.cycle_state_points(ref, T_2d=REC['T_2d_realized'])
print(f'T_2d provisional (top layer fully molten)   {_prov:7.2f} C')
print(f"T_2d realized on the Grissom field          "
      f"{REC['T_2d_realized'] - 273.15:7.2f} C"
      f"   ({REC['T_2d_realized'] - 273.15 - _prov:+.1f} K)")
print('every diagram and state-point table below is at the realized value.')
"""))

cells.append(md(r"""
### 10.1 The two cycles side by side
"""))

cells.append(code(r"""
fig, _ = BP.ts_both(ref, rank, hp)
fig.savefig('BOREAS_Ts_cycles.pdf')
fig.show()
"""))

cells.append(md(r"""
### 10.2 Both on one axes, with the water levels

This is the coupling diagram. The two cycles never exchange working fluid ---
they exchange *water* --- so the horizontal dotted lines are the only thing
that makes the pair a single machine rather than two unrelated plots. Four
levels are shown: the sink, the geothermal source feeding the HTHP evaporator,
the HTHP condenser inlet $T_\mathrm{2c}$ (the charge loop), and
$T_\mathrm{2d}$, the store outlet --- the one line the borehole model sets.
"""))

cells.append(code(r"""
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8.0, 5.8))
BP.ts_combined(ref, rank, hp, Tdict, ax=ax)
fig.tight_layout()
fig.savefig('BOREAS_Ts_combined.pdf')
fig.show()
"""))

cells.append(md(r"""
### 10.2b The water on the $T$--$s$ diagram, with every pinch marked

The horizontal lines above are *inlet temperatures*. They say where the water
starts, not what it does, and a pinch cannot be read off them. The water
actually glides --- $T_\mathrm{2d} \to T_\mathrm{3d}$ through the ORC
evaporator, $T_\mathrm{2c} \to T_\mathrm{3c}$ through the HTHP condenser ---
and it is the glide against a *flat* boiling or condensing plateau that
creates the binding approach.

Water cannot be plotted on these axes directly: its entropy is on a scale with
nothing to do with cyclopentane's. What can be plotted is the water
**temperature at each point of the exchanger**, placed at the entropy the
working fluid has at that same point. The two streams are tied together by
cumulative duty, so at a given $Q$ the working fluid sits at $(s, T_\mathrm{wf})$
and the water at $T_\mathrm{w}$; plotting $(s, T_\mathrm{w})$ puts the water
line directly above (evaporator) or below (condenser) the working-fluid path,
and the **vertical gap at every abscissa is the real local approach**.

This is a duty mapping, not a thermodynamic path. The water line is not a
process in cyclopentane's coordinates and its abscissa means nothing for the
water itself. Stated plainly it is an honest and very readable construction;
left unstated it would be a category error, so the axis label says so.
"""))

cells.append(code(r"""
fig, ax = plt.subplots(figsize=(8.6, 6.2))
BP.ts_combined(ref, rank, hp, Tdict, ax=ax, water='glide', model=B)
fig.tight_layout()
fig.savefig('BOREAS_Ts_combined_water.pdf')
fig.show()
"""))

cells.append(md(r"""
Two things are now visible that the level lines hid.

The ORC-evaporator water (black) descends across the diagram and **touches the
working fluid at the bubble-point corner** --- the pinch is there, not at
either end. And the charge water (teal, $T_\mathrm{2c} \to T_\mathrm{3c}$) runs
nearly parallel to it, offset by the storage approach: that offset is the
round trip through the PCM, and it is the only part of this figure the
borehole model can move.

### 10.2c One panel per exchanger

Four exchangers in four corners of the diagram; overlaying them all makes a
figure in which no individual pinch can be read. Each panel below carries the
dome, the working-fluid path, the duty-mapped water line and the pinch arrow.
The recovered pinches should reproduce the specification
($\Delta T_\mathrm{pinch,ORC}$, $\Delta T_\mathrm{E,sink}$,
$\Delta T_\mathrm{2H,3C}$, $\Delta T_\mathrm{pinch,HPE}$) --- if one of them
does not, the overlay and the model disagree and the overlay is wrong.
"""))

cells.append(code(r"""
fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.8))
_, ov = BP.ts_water(B, ref, rank, hp, Tdict, axes=axes)
fig.tight_layout()
fig.savefig('BOREAS_Ts_water.pdf')

print(f"{'exchanger':20s}{'pinch, K':>10}{'specified':>12}")
spec = {'ORC evaporator': ref.DT_pinch_ORC,
        'ORC condenser': ref.DT_E_sink,
        'HTHP condenser': ref.DT_2H_3C,
        'HTHP evaporator': ref.DT_pinch_HPE}
for k, v in ov.items():
    print(f"{k:20s}{v['pinch_K']:>10.3f}{spec[k]:>12.3f}")
fig.show()
"""))

cells.append(md(r"""
### 10.3 ORC evaporator composite curves: the pinch is *interior*

This is the figure that explains the gap. The water glides down from
$T_\mathrm{2d}$ to $T_\mathrm{3d}$ as a nearly straight sensible line; the
working fluid preheats, then goes **flat** at $T_\mathrm{3e}$ for the rest of
the exchanger. The smallest approach therefore occurs neither at the hot end
nor at the cold end but at the corner where boiling starts, and panel (b)
makes the V unmistakable.

The grey dashed curve is the cold composite that the hot-end rule *alone*
would permit, $T_\mathrm{3e} = T_\mathrm{2d} - \Delta T_\mathrm{2D,3E}$. It
does not merely violate the 5 K pinch --- it sits **above** the water over
most of the duty, so heat would have to flow from the cold stream to the hot
one. That is not a slightly optimistic design, it is a second-law violation,
and the efficiency computed from it is unreachable at any exchanger size. The
printed `hot_end_reversed_pct` is the fraction of the load that runs
backwards.

The curves here are the model's own: `orc_composite` calls
`_orc_cold_composite` and `_water_curve` rather than rebuilding them, so what
is plotted is by construction what `orc_pinch` minimizes over. The returned
`pinch_K` should agree with `rank['pinch']` to rounding --- that agreement is
the figure's self-check, and it is worth reading rather than trusting.
"""))

cells.append(code(r"""
fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.7))
_, info = BP.orc_composite(B, ref, rank, Tdict, axes=axes)
fig.tight_layout()
fig.savefig('BOREAS_ORC_composite.pdf')

for k, v in info.items():
    print(f'  {k:24s} {v:10.4f}')
print(f"\n  model's own pinch        {rank['pinch']:10.4f}   <- must agree")
"""))

cells.append(md(r"""
Note where the onset of boiling lands on the duty axis. It is **not**
$(h_\mathrm{bub} - h_\mathrm{10e})/Q_\mathrm{tot}$: that is the main stream's
preheat alone and lands about five points of duty too low, because the
composite also carries the reheat stream ($y$ of the flow, $5\mathrm{e} \to
6\mathrm{e}$), whose entire duty lies *below* the bubble temperature and is
therefore already banked by the time the main stream starts to boil. The break
is read off the assembled composite instead of re-derived --- re-deriving it
is how the marker ends up a few points away from the pinch it is supposed to
coincide with.

### 10.3b Look at the gap between $T_\mathrm{2d}$ and $3\mathrm{e}$

In the reference case the store delivers water at $T_\mathrm{2d}$ but the ORC
evaporates well below it. Which constraint is binding is reported by the model
itself, and it is not always the one you would assume: if the hot-end rule
$T_\mathrm{3e} = T_\mathrm{2d} - \Delta T_\mathrm{2D,3E}$ bound, the approach
would be 12 K. Run the cell and read `T_3e_binding`.
"""))

cells.append(code(r"""
T_2d = Tdict['T_2d'] - 273.15
print(f"water from the store, T_2d        {T_2d:8.2f} C")
print(f"hot-end rule would give T_3e      {rank['T_3e_hot_end'] - 273.15:8.2f} C"
      f"   (approach {T_2d - (rank['T_3e_hot_end'] - 273.15):.1f} K)")
print(f"realized T_3e                     {rank['T_3e'] - 273.15:8.2f} C"
      f"   (approach {T_2d - (rank['T_3e'] - 273.15):.1f} K)")
print(f"binding constraint                {rank['T_3e_binding']}")
print()
T_carnot = 1.0 - (rank['T_1e'] / rank['T_3e'])
print(f"Carnot between T_1e and T_3e      {T_carnot:8.4f}")
print(f"eta_ORC                           {rank['rank_eff']:8.4f}"
      f"   ({rank['rank_eff'] / T_carnot * 100:.0f} % of Carnot)")
"""))

cells.append(md(r"""
### 10.4 State points and the scalar summary
"""))

cells.append(code(r"""
st = BP.state_table(rank, hp)
display(st)
st.to_csv('BOREAS_state_points.csv')

# h is kJ/kg in the model's dictionaries while s is J/kg/K. The mixed units
# are the model's; state_table converts s and leaves h alone. Do not "fix"
# this by dividing h -- that was a real error earlier in this work.
# to_latex routes through the Styler in pandas 2.x and wants
# jinja2 >= 3.1.2. Degrade instead of crashing the notebook on a stack that
# has an older one.
try:
    print(st.to_latex(float_format='%.3f'))
except ImportError as e:
    print(f'(to_latex unavailable: {e})\n')
    print(st.to_string())
"""))

cells.append(code(r"""
display(BP.cycle_summary(rank, hp))
"""))

cells.append(md(r"""
### 10.5 What actually moves the cycles

Worth running before assuming otherwise. The geothermal source temperature
moves the **COP** --- it sets the heat-pump lift --- and *nothing on the ORC
side*. $T_\mathrm{2d}$ is fixed by the cascade specification
($\Delta T_\mathrm{3C,2C}$ and the melting temperature), not by the source, so
$T_\mathrm{3e}$ and $\eta_\mathrm{ORC}$ do not budge. The two cycles are
coupled to the store but decoupled from each other, which is why the DOE found
$\eta_\mathrm{RTE}$ below ground to be almost entirely pumping.
"""))

cells.append(code(r"""
alt = REF.with_(grad_K_per_km=56.5, T_source_C=70.0)
rank_a, hp_a, T_a = B.cycle_state_points(alt)
print(f"T_2d   {Tdict['T_2d'] - 273.15:7.2f} ->{T_a['T_2d'] - 273.15:7.2f} C")
print(f"T_3e   {rank['T_3e'] - 273.15:7.2f} ->{rank_a['T_3e'] - 273.15:7.2f} C")
print(f"etaORC {rank['rank_eff']:7.4f} ->{rank_a['rank_eff']:7.4f}")
print(f"COP    {hp['hp_cop']:7.4f} ->{hp_a['hp_cop']:7.4f}")
"""))

### ------------------------------------------- 9. well profiles
cells.append(md(r"""
---
## 11. Thermal behavior along the well

The figures above are about the surface machines. These are about the store
itself, and they need a RECORDED march:
`simulate_css_corrected(..., record=True)` keeps the fluid and PCM
temperature, the melted area, $U_i$ and $q'$ at every time level, at the
cost of one extra cycle.

### The coordinate, which is the thing to get right

**The solver's `z` is not depth.** Quoting `formation_profile`: on charge
the water enters the fin tubes at the *bottom* --- it arrives there down the
insulated downcomer --- so `z = 0` is the bottom of the well and depth
decreases with `z`. Plotting against the raw `z` and calling it depth draws
the well upside down, which is exactly what an earlier version of these
figures did. `depth_axis()` takes the depth array from the model itself so
there is one definition.

Both half-cycle histories come back in the **charge** frame, index 0 at the
bottom, even though the discharge march is run mirrored. Reading the fluid
boundary condition does not settle that, because the inlet sits at opposite
ends in the two half-cycles. The test that does: PCM cannot be above its own
melting point while it is freezing, so at the end of discharge
`max(T_pcm - T_m)` must be $\le 0$ when the history is read against the
cascade **forward**, and the cell prints both readings rather than asserting
one.

Read the two numbers for what they are. The forward reading is the one that
has to pass, and it does, at `+0.00` K. The reversed reading is *evidence*
only when it fails — on the idealized vertical geometry it came out at
`+26.9` K, a flagrant violation. On the Grissom S-curve it comes out
negative too, because the cascade and the realized PCM temperatures both
fall towards the top and the mirrored comparison happens not to break. **A
reversed reading that passes is not a second confirmation**; it means this
particular test has stopped discriminating for this geometry, and the
forward reading is carrying the check alone.
"""))

cells.append(code(r"""
# THE SAME RECORDED MARCH as section 10, run inside field_ground: every
# profile below is on the measured Grissom field, with the superposition
# resistance R_g,eff(MD) and the rock temperature on TRUE VERTICAL DEPTH.
# Run it outside and the rock at the bottom of the store comes out 15 K too
# hot, because the S-curve reaches only 891 m TVD at 1250 m MD.
rec = REC
ref_well = REF
A_av = BP.pcm_capacities_A(ref_well)
# Depth from the TRAJECTORY, not from the solver coordinate. On an S-curve
# the two differ by a third: 1250 m of hole reaches 891 m of depth.
import numpy as _np
_xi = (_np.arange(99) + 0.5) / 99
_mdv = REF.L_well * (1.0 - _xi)                  # charge frame: 0 = bottom
dseg = _np.interp(_mdv, _np.r_[0.0, MD], _np.r_[0.0, TVD])
dnode = _np.interp(REF.L_well * (1.0 - _np.linspace(0, 1, 100)),
                   _np.r_[0.0, MD], _np.r_[0.0, TVD])
print(f'depth runs {dseg[0]:.0f} m (index 0, bottom) -> {dseg[-1]:.0f} m (top)')

# the orientation check, printed rather than trusted
Tm = rec['T_m_seg'] - 273.15
Tp = rec['discharge']['history']['T_pcm'][-1] - 273.15
print(f'end of discharge, max(T_pcm - T_m): forward {(Tp - Tm).max():+.2f} K, '
      f'reversed {(Tp - Tm[::-1]).max():+.2f} K  <- forward must be <= 0')
"""))

cells.append(code(r"""
fig, ax = plt.subplots(1, 2, figsize=(8.8, 5.8), sharex=True)
BP.well_profiles(rec, ref_well, axes=ax, vertical=True, staircase=True,
                 depth_seg=dseg, depth_node=dnode, offset=-273.15, n_times=6,
                 field='T_fluid', overlay='T_pcm',
                 overlay_label=r'$T_\mathrm{pcm}$',
                 ylabel=r'temperature / $^\circ$C')
fig.tight_layout(); fig.show()
"""))

cells.append(code(r"""
fig, ax = plt.subplots(1, 2, figsize=(8.8, 5.8), sharex=True)
BP.well_profiles(rec, ref_well, axes=ax, vertical=True, field='A_melt',
                 depth_seg=dseg, depth_node=dnode, scale=1.0 / A_av,
                 ylabel=r'melted fraction, $A_\mathrm{melt}/A_\mathrm{av}$')
fig.tight_layout(); fig.show()

fig, ax = plt.subplots(1, 2, figsize=(8.8, 5.8), sharex=True)
BP.well_profiles(rec, ref_well, axes=ax, vertical=True, field='U_i',
                 depth_seg=dseg, depth_node=dnode,
                 ylabel=r'$U_i$ / W m$^{-2}$ K$^{-1}$')
fig.tight_layout(); fig.show()

fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.2))
BP.well_integrals(rec, ref_well, axes=ax)
fig.tight_layout(); fig.show()
"""))

cells.append(md(r"""
### 11.1 Animation of the full cycle

Charge then discharge, three panels, with the frames resampled onto a
**uniform** time grid --- the march records on a roughly logarithmic one,
which would make constant-rate playback appear to accelerate.
"""))

cells.append(code(r"""
from IPython.display import Video
BP.animate_well(rec, ref_well, 'BOREAS_cycle.mp4', fps=12, frames_per_half=60)
Video('BOREAS_cycle.mp4', embed=True, width=760)
"""))


cells.append(code(r"""
mine = REF.with_(
    n_ft=5,                    # more, smaller tubes -> more PCM, more pumping
    grad_K_per_km=56.5,        # central Wilmington, 1949 survey
)
for lvl, m in B.validate_case(mine, verbose=False):
    print(f'  [{lvl}] {m}')

r = B.simulate_css_corrected(mine)
B.css_report(mine, r)
B.kpi_report(mine, r)
"""))

cells.append(md(r"""
---
### What this model will not tell you

- **Anything about the surface machines.** They are frozen here. Round-trip
  efficiency is set by the heat-pump *lift* — $(T_h - T_0)/(T_h - T_{src})$ —
  and so by the geothermal source temperature, not by anything below ground.
- **A multi-year loss curve.** `t_operation_yr` evaluates the ground
  resistance at one age. The decay from year 1 to year 30 comes from running
  it several times, not from an outer loop.
- **An optimum.** This is a model, not a search.
"""))

nb = {'cells': cells,
      'metadata': {'colab': {'provenance': [], 'toc_visible': True},
                   'kernelspec': {'name': 'python3', 'display_name': 'Python 3'},
                   'language_info': {'name': 'python'}},
      'nbformat': 4, 'nbformat_minor': 0}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(nb, indent=1))
n_code = sum(1 for c in cells if c['cell_type'] == 'code')
print(f'{OUT.name}')
print(f'  {len(cells)} cells: {n_code} code / {len(cells) - n_code} markdown')
