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

*v1.0 · built {STAMP}*

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
cells.append(md(r"""
---
## 2. The hole, in cross-section

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
## 3. The formation

The thermal circuit from the PCM outward is a series chain, and the terms are
wildly unequal:

$$R_{tot} = \underbrace{R_{PCM}}_{\sim 25\%} + \underbrace{R_{casing}}_{0.1\%}
+ \underbrace{R_{cement}}_{5\%} + \underbrace{R_{ground}(t)}_{\sim 70\%}$$

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
## 4. $B(z)$ — the field is a cone, not an array

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
cf = B.Case(n_wells_field=1000.0)
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
## 5. Run it

`n_wells_field=None` is the **isolated well** — the conservative bound, and
the v1.0 default. It is not a pessimistic estimate of the real loss; it is a
different configuration, and one nobody would build.

About a minute per case.
"""))

cells.append(code(r"""
rows = {}
for lab, kw in [('isolated well', {}),
                ('fan, 100 wells', dict(n_wells_field=100.)),
                ('fan, 1000 wells', dict(n_wells_field=1000.)),
                ('array, 1.83 m', dict(n_wells_field=1000., fan_angle_deg=0.))]:
    cc = B.Case(n_segments=99, **kw)
    r = B.simulate_css_corrected(cc)
    k = B.performance_indices(cc, r)
    rows[lab] = {'loss %': 100 * k['f_formation_loss'],
                 'eta_storage': r['eta_storage'],
                 'eta_RTE': k['eta_RTE'],
                 'eps_cycled': k['eps_cycled'],
                 'wells/MWe': k['wells_per_MWe']}
    print(f'  {lab} done')
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
---
## 6. The gates, watched to fail

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
src = open('boreas.py').read()
r_ok = B.simulate_css_corrected(B.Case(n_segments=99))
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
    rb = ns['simulate_css_corrected'](ns['Case'](n_segments=99))
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
r0 = B.simulate_css_corrected(B.Case(n_segments=99, k_rock=1e-9))
print(f'  k_rock -> 0 : Q_loss {r0["Q_loss_kJ"]:.3e} kJ, '
      f'eta_storage {r0["eta_storage"]:.9f}')
print(f'  {"PASS" if r0["Q_loss_kJ"] == 0.0 else "FAIL"}'
      f'   (residual is the CSS drift tolerance, not a loss)')
"""))

# ------------------------------------------------------------- playground
cells.append(md(r"""
---
## 7. Change something

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

cells.append(code(r"""
mine = B.Case(
    n_segments=99,
    n_ft=5,                    # more, smaller tubes -> more PCM, more pumping
    grad_K_per_km=56.5,        # central Wilmington, 1949 survey
    t_operation_yr=10.0,
    n_wells_field=1000.0,
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
