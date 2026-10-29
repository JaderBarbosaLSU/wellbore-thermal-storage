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

*v1.3 · built {STAMP}*

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
have to coexist at the same depth, centimetres apart, in contiguous PCM with
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
## 2. The design point — **edit here, nowhere else**

Every parameter you are likely to change lives in the next cell. Nothing
below it hard-codes a number: each later section starts from `CASE` and
says, in one visible `.with_(...)`, what it is changing and why.

The split is deliberate:

* **`CASE`** is the well, the material and the plant — what you are
  building. It says nothing about how many of them there are or how they
  are arranged on the surface.
* **`FIELD`** is the arrangement — count, trajectory, packing. It is kept
  separate because the *same* well gives losses from 0.9 % to 79 %
  depending only on this, and because the model's defaults here describe an
  **isolated** well, which is a bound and not a design.

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
    L_ft       = 5000.0,     # useful storage length                       [ft]
    n_ft       = 3,          # fin tubes per well, all downward             [-]
    num_fins   = 24,         # fins per tube                                [-]
    fin_L      = 0.0035,     # fin length                                   [m]
    fin_t      = 0.0015,     # fin thickness                                [m]
    t_ins      = 0.008,      # DOWNCOMER insulation (the casing is bare)    [m]
    k_ins      = 0.007,      # downcomer insulation, VIT                [W/m/K]
    r_casing_o = 0.098,      # casing outer radius                          [m]
    r_bore     = 0.110,      # drilled hole radius, 7 in                    [m]

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
#  THE FIELD.  How many of that well, and how they are arranged.
#    trajectory : fan_angle_deg = 0 is a vertical array; > 0 deviates below
#                 the kickoff, so the spacing grows with depth
#    packing    : 'disc' a filled circle, 'rows' two straight lines side by
#                 side per island, 'ring' that double line closed into a loop
# =============================================================================
FIELD = dict(
    n_wells_field = 1000.0,  # wells sharing the formation                  [-]
    fan_angle_deg = 0.0,     # 0 = vertical array                         [deg]
    packing       = 'rows',  # 'disc' | 'rows' | 'ring'
    n_islands     = 4,       # separate clusters                            [-]
    B_surface_m   = 1.83,    # wellhead spacing, 6 ft                       [m]
    row_gap_m     = 3.66,    # between the two rows, 12 ft                  [m]
)

REF = CASE.with_(**FIELD)        # the reference case used throughout
print('design point set.  REF =  CASE  +  FIELD')
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

The downcomer runs **off-centre**, as in the schematic. A concentric
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

# ----------------------------------------------------------------- B(z)
cells.append(md(r"""
---
## 5. $B(z)$ — the field is a cone, not an array

THUMS wellheads sit on **6 ft centres**, 12 ft between double rows, and over
1200 wells are drilled *directionally* from four islands with drift angles to
84° across 6500 acres. So the spacing is not a number — it is a function of
depth.

The governing quantity is how long neighbours take to merge thermally,
$r^2/\alpha$:

| spacing | merge time |
|---|---|
| 1.83 m (6 ft) | **36 days** |
| 8 m | 1.9 yr |
| 162 m (reservoir) | 779 yr |

The shallow section is **one thermal body within a month**; the deep section is
a field of isolated wells. In the same well.

**And the two effects oppose each other.** Shallow rock is cold, so the driving
difference is largest exactly where shielding is strongest; deep rock is hot,
so it is smallest exactly where the wells stand alone. That partial
cancellation is specific to directionally-drilled island fields — it does not
arise in purpose-drilled BTES, where spacing is uniform by design.

A shielded well's loss is the cluster perimeter flux shared among $N$ wells,

$$R_{perimeter}(z) = \frac{N\sqrt{\pi\alpha t}}{2\pi k\,R_{field}(z)}$$

which carries no driving temperature and so is a true resistance. The two
regimes combine by taking the **larger** resistance — not a blend: a well
cannot leak more than if it stood alone, a field cannot leak more than its
perimeter allows, and whichever binds governs.
"""))

cells.append(code(r"""
cf = CASE.with_(n_wells_field=1000.0)
z = np.linspace(0, cf.L_well, 200)
Bz = B.well_spacing(cf, z)
Rf = B.formation_resistance(cf, z)
R_iso = B.formation_resistance(B.Case())

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8))
a1.plot(Bz, z, lw=2.5, color='#1a5276'); a1.invert_yaxis()
a1.set_xlabel('spacing B(z), m'); a1.set_ylabel('depth, m'); a1.grid(alpha=.3)
a1.set_title(f'{cf.B_surface_m:.2f} m at the wellhead, '
             f'{Bz[-1]:.0f} m at {cf.L_well:.0f} m', fontsize=11)
a2.plot(Rf, z, lw=2.5, label='field, B(z)')
a2.axvline(R_iso, color='#c0392b', ls='--', lw=2, label='isolated well')
a2.invert_yaxis(); a2.set_xscale('log')
a2.set_xlabel('formation resistance, m K/W'); a2.set_ylabel('depth, m')
a2.legend(); a2.grid(alpha=.3)
a2.set_title('shielding weakens with depth but never stops binding', fontsize=11)
plt.tight_layout(); plt.show()
print(f'  shielding factor: {Rf[0]/R_iso:6.0f}x at the top, '
      f'{Rf[-1]/R_iso:5.1f}x at the bottom')
"""))


# ----------------------------------------------------------------- run
cells.append(md(r"""
---
## 6. Run it

`n_wells_field=None` is the **isolated well** --- the conservative bound,
and the v1.3 default. It is not a pessimistic estimate of the real loss; it
is a different configuration, and one nobody would build.

**What the four words mean.** Two of them describe the *trajectories*, two
the *plan view* of the wellheads:

| term | meaning |
|---|---|
| **array** | the wells stay vertical, $\theta = 0$, so the spacing $B(z)$ is the wellhead spacing at every depth |
| **fan** | the wells are vertical to the kickoff and then deviate at half-angle $\theta$, so $B(z)$ grows with depth |
| **disc** | a circle *filled* with wellheads on $B_0$ centres; $P_0 = 2\pi R_0$ with $R_0 = B_0\sqrt{N/\pi}$ |
| **rows** | two straight lines of wells side by side, $w_\mathrm{row}$ apart, per island; $P_0 = n_\mathrm{isl}\,2(L_\mathrm{row} + w_\mathrm{row})$ |
| **ring** | the same double line closed into a loop, like stands round a stadium; it loses from **both** faces, inner and outer, so $P_0 = n_\mathrm{isl}\,2\,L_\mathrm{row}$ |

So *array* and *fan* are alternatives, and *disc*, *rows* and *ring* are
alternatives; every field case is one of each.

**Packing must be named explicitly, and that is not decoration.** The model
default is `packing='disc'`, which *minimises* perimeter for a given well
count and is therefore the most favourable arrangement that exists. A row
that omits `packing` silently reports the superseded v1.2 geometry: the
vertical array comes out at 0.9 % loss instead of 7.9 %. Comparing a fan
computed at disc packing against an array computed at rows is how an
earlier draft concluded the two were "the same machine".

Fourteen cases, about a minute each.
"""))

cells.append(code(r"""
import numpy as np, pandas as pd
# the geometry, before any thermal calculation
geo = {}
for N in (100., 1000.):
    for pk, ni in (('disc', 1), ('rows', 4), ('ring', 4)):
        cc = CASE.with_(n_wells_field=N, fan_angle_deg=0., packing=pk, n_islands=ni)
        P0 = float(np.atleast_1d(B.field_perimeter(cc, 0.))[0])
        geo[(int(N), pk)] = {'islands': ni, 'P0 [m]': P0,
                             'per well [m]': P0 / N}
g = pd.DataFrame(geo).T.round(3)
g.index.names = ['N', 'packing']
display(g)
print('ring differs from rows only by the two end-caps of each row,')
print('2*w_row per island: 1.6 % of the perimeter at N=1000, 14 % at N=100.')
"""))

cells.append(code(r"""
import pandas as pd
DISC = dict(packing='disc', n_islands=1)
ROWS4 = dict(packing='rows', n_islands=4)
RING4 = dict(packing='ring', n_islands=4)

CONFIGS = [
    ('adiabatic, k_r -> 0',      dict(k_rock=1e-9, **DISC)),
    ('isolated well',            dict(**DISC)),
    ('fan, 100 wells, disc',     dict(n_wells_field=100., **DISC)),
    ('fan, 100 wells, rows',     dict(n_wells_field=100., **ROWS4)),
    ('fan, 1000 wells, disc',    dict(n_wells_field=1000., **DISC)),
    ('fan, 1000 wells, rows',    dict(n_wells_field=1000., **ROWS4)),
    ('fan, 1000, 70 deg, disc',  dict(n_wells_field=1000., fan_angle_deg=70., **DISC)),
    ('fan, 1000, 70 deg, rows',  dict(n_wells_field=1000., fan_angle_deg=70., **ROWS4)),
    ('vertical array, disc',     dict(n_wells_field=1000., fan_angle_deg=0., **DISC)),
    ('vertical array, rows',     dict(n_wells_field=1000., fan_angle_deg=0., **ROWS4)),
    ('fan, 100 wells, ring',     dict(n_wells_field=100., **RING4)),
    ('fan, 1000 wells, ring',    dict(n_wells_field=1000., **RING4)),
    ('fan, 1000, 70 deg, ring',  dict(n_wells_field=1000., fan_angle_deg=70., **RING4)),
    ('vertical array, ring',     dict(n_wells_field=1000., fan_angle_deg=0., **RING4)),
]

rows = {}
for lab, kw in CONFIGS:
    cc = CASE.with_(**kw)        # CASE, not a fresh Case(): one design point
    r = B.simulate_css_corrected(cc)
    k = B.performance_indices(cc, r)
    rows[lab] = {'loss %': 100 * k['f_formation_loss'],
                 'eta_storage': r['eta_storage'],
                 'eta_RTE': k['eta_RTE'],
                 'eps_cycled': k['eps_cycled'],
                 'wells/MWe': k['wells_per_MWe'],
                 'T_2d [C]': r['T_2d_realised'] - 273.15}
    print(f'  {lab:26s} loss {rows[lab]["loss %"]:5.1f} %', flush=True)
display(pd.DataFrame(rows).T.round(4))
"""))

cells.append(md(r"""
Expected, with the default 45 deg fan:

| configuration | loss | $\eta_{storage}$ | $\eta_{RTE}$ |
|---|---|---|---|
| isolated well | 75.6 % | 0.244 | 0.072 |
| fan, 100 wells | 54.9 % | 0.451 | 0.138 |
| **fan, 1000 wells** | **12.8 %** | **0.872** | **0.281** |
| array at 1.83 m | 0.93 % | 0.991 | 0.323 |

Two orders of magnitude between the bounds, from one geometric parameter. The
realistic case loses **12.8 %** — about 15 % relative on round-trip efficiency
against the adiabatic 0.323. Material, and nothing like fatal.

The fan *angle* matters nearly as much as the well count: steepening 45 deg to
70 deg more than doubles the loss, because the wells separate faster with
depth. In a real field that angle is set by where the reservoir is, which makes
it a site parameter rather than a design one.

**This also settles the insulation question.** For an isolated well there is a
shallow optimum near 5 mm of casing insulation, worth about +19 % on delivered
energy before the PCM-volume penalty turns the curve over. In a field the loss
is already under 1 %, so insulation trades roughly 11 % of the PCM for
nothing.

That an optimum exists *at all* is a symptom that the isolated-well scope is
unphysical: a model whose advice is "insulate" is reporting on a configuration
that does not exist.
"""))

# ----------------------------------------------------------------- gates
cells.append(md(r"""
**Closing the rows into a ring buys almost nothing.** The ring loses from
both faces, so its perimeter is the rim length twice over, and against the
double row it saves only the two end-caps --- 1.6 % of the perimeter at
1000 wells. The loss moves from 7.9 % to 7.8 % for the vertical array and
from 57.1 % to 56.9 % for the fan. At 100 wells the saving is larger (14 %
of the perimeter, 65.5 % to 64.9 % loss) because the end-caps are a bigger
share of a shorter row. The arrangement that matters is whether the wells
are strung out at all, not how the line is closed.

**What this table does not contain.** The paper also reports a
*superposition* variant of the two 1000-well rows, in which
`formation_resistance` is replaced by a line-source superposition over the
actual well positions. That is **not** what the shipped code computes, so
those rows cannot be reproduced from this notebook; they came from a
separate calculation and are labelled as such in the paper.
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
import dataclasses
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

# ------------------------------------------------------------- playground
cells.append(md(r"""
---
## 7b. Packing is not spacing

`field_radius` in v1.0--v1.2 turned a well count and a spacing into a
footprint the obvious way, $R_0 = B_0\sqrt{N/\pi}$ --- a compact disc. The
spacing was right (6 ft, from the Oxy deck). The **packing** was not.

A disc *minimises* perimeter for a given count, and shielding goes as
$1/\text{perimeter}$, so it is the most favourable arrangement that exists.
THUMS wellheads sit in **linear double rows along the island rim**, split
across **four islands**, each with its own boundary.

| arrangement | perimeter per well |
|---|---|
| disc | 0.205 m |
| four islands, double rows | **1.859 m** |

**9.06x more perimeter, so 9x the loss.** `packing` is now a `Case` field.
"""))

cells.append(code(r"""
for pk, kw in [('disc', {}), ('rows', dict(n_islands=4))]:
    cc = CASE.with_(**{**FIELD, 'packing': pk, **kw})
    P = float(np.atleast_1d(B.field_perimeter(cc, 0.))[0])
    R = float(np.atleast_1d(B.formation_resistance(cc, np.array([0.])))[0])
    print(f'{pk:>6}: perimeter {P:7.0f} m -> {P/1000:.3f} m per well, '
          f'R_form {R:7.3f} m K/W')
"""))

# ------------------------------------------------------- depth sweep
cells.append(md(r"""
---
## 8. How deep should the storage go?

The 2014 Oxy drilling deck settles the geometry this model was abstracting.
THUMS wellheads sit on **6 ft centres** in double rows, and the wells are
**directionally drilled**: 8200 ft of horizontal displacement, **4175 ft
TVD**, 80° maximum inclination.

Two things follow, and both matter.

**The wells do not fan uniformly — the bundle *sheds* members.** Each well
kicks off shallow and leaves at its own depth; the ones that remain stay
tightly packed. The cone abstraction in section 4 keeps all $N$ wells and
grows their spacing. Reality keeps a tight core and loses wells from it.

**Storage can only occupy the upper, near-vertical, large-diameter cased
section** — which is exactly where the bundle is still tight. So the realistic
spacing over the storage interval is close to the *island* value, not the
fanned one.

That makes depth a genuine open question, and not the one you might expect.
"""))

cells.append(code(r"""
import pandas as pd
rows = {}
for L in (1000., 2000., 3000., 4000.):
    for tag, kw in [('island', dict(n_wells_field=1000., fan_angle_deg=0.)),
                    ('fan 45', dict(n_wells_field=1000.))]:
        cc = CASE.with_(L_ft=L, **kw)
        r = B.simulate_css_corrected(cc)
        k = B.performance_indices(cc, r)
        rows[(tag, L)] = {'loss %': 100*k['f_formation_loss'],
                          'eta_RTE': k['eta_RTE'],
                          'wells/MWe': k['wells_per_MWe'],
                          'kW/well': k['kW_per_well']}
    print(f'{L:.0f} ft done', flush=True)
df = pd.DataFrame(rows).T.round(4)
df.index.names = ['spacing', 'L_ft']
display(df)
"""))

cells.append(md(r"""
### The answer, and it is not a compromise

| $L$, ft | island loss | island $\eta_{RTE}$ | wells/MWe | fan-45 loss | fan-45 $\eta_{RTE}$ |
|---|---|---|---|---|---|
| 1000 | 1.29 % | 0.3253 | 106.0 | 1.29 % | 0.3253 |
| 1500 | 1.24 % | 0.3256 | 73.1 | 2.39 % | 0.3214 |
| 2000 | 1.20 % | 0.3256 | 56.7 | 4.25 % | 0.3146 |
| 3000 | 1.11 % | 0.3253 | 40.4 | 7.96 % | 0.3008 |
| 4000 | 1.02 % | 0.3245 | 32.3 | 10.86 % | 0.2894 |

**At island spacing the loss is flat at about 1 % and $\eta_{RTE}$ is pinned
near 0.325 at every depth.** Shorter wells are even marginally *better*,
because shorter tubes pump more easily.

So there is **no thermodynamic optimum in depth**. Depth is purely a capital
lever: 106 wells per MWe at 1000 ft against 32 at 4000 ft, a factor of 3.3.

> **Go as deep as the near-vertical section allows.** The binding constraint
> is well geometry, not heat transfer.

Note also that the two columns agree exactly at 1000 ft: below a 300 m kickoff
the fan has not opened yet, so there is nothing to shield against.

**A caution on the default.** `L_ft = 5000` exceeds the 4175 ft TVD of the
deepest THUMS well. Anything run at the default depth is a well this field
does not have.

**And on scale.** At 32 wells per MWe, the roughly 1000 active THUMS wells
support of order **30 MWe** — if every one of them were converted.
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
| `n_wells_field` | `None` = isolated; a number = shielded |
| `B_surface_m`, `z_kickoff_m`, `fan_angle_deg` | the directional fan |

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
rank, hp, Tdict = B.cycle_state_points(ref)
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
is plotted is by construction what `orc_pinch` minimises over. The returned
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
print(f"realised T_3e                     {rank['T_3e'] - 273.15:8.2f} C"
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
## 11. Thermal behaviour along the well

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
melting point while it is freezing, and at the end of discharge
`max(T_pcm - T_m)` is `+0.00` K read forward against `+26.9` K reversed.
"""))

cells.append(code(r"""
import pickle, os
if os.path.exists('rec.pkl'):
    rec = pickle.load(open('rec.pkl', 'rb'))
else:
    rec = B.simulate_css_corrected(REF, record=True)   # ~1 min
    pickle.dump(rec, open('rec.pkl', 'wb'))

ref_well = REF
A_av = BP.pcm_capacities_A(ref_well)
dseg, dnode = BP.depth_axis(B, ref_well, 99)
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
