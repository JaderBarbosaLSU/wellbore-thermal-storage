"""Handoff section 10: the orange checks only the code can close.

Resumable: every value is written to orange.json as it is computed, and
a block whose values are already there is skipped.
"""
import json
import os
import sys

import numpy as np

import boreas as B
import boreas_ground as G

CSV = 'THUMS_wells_WellSTAR_clean.csv'
L_MD = 1250.0
THUMS = dict(D_in=201.2 / 25.4, r_casing_o=0.21910 / 2,
             r_bore=0.31120 / 2, k_cement=0.7)
CASE = B.Case(n_segments=99, L_ft=L_MD / 0.3048, t_operation_yr=10.0,
              **THUMS)
STORE = 'orange.json'
M = json.load(open(STORE)) if os.path.exists(STORE) else {}
NOTE = {}


def add(n, v, fmt='{:.4g}', note=''):
    M[n] = v if isinstance(v, str) else fmt.format(v)
    NOTE[n] = note
    json.dump(M, open(STORE, 'w'), indent=1)


def have(*names):
    return all(n in M for n in names)


def field(traj='s_curve', L=L_MD):
    return G.cached_reference_field(
        B, CASE, CSV, trajectory=traj, kop_rule='staggered', theta_max=55.0,
        kop_ft=300.0, L_md=L, ds=G.DS_STATION, verbose=False)[:3]


md, Rg, tvd = field()

# --- 2.2  formation temperature at the ends of the store ---------------
if not have('RefTfTop'):
    with G.field_ground(B, md, Rg, tvd):
        fp = B.formation_profile(CASE)
    Tf = (np.asarray(fp[0]) if isinstance(fp, tuple) else np.asarray(fp))
    Tf = Tf - 273.15
    add('RefTfBottom', float(Tf[0]), '{:.1f}', 'formation T, bottom, C')
    add('RefTfTop', float(Tf[-1]), '{:.1f}', 'formation T, top, C')
    add('RefTvdBottom', float(tvd[-1]), '{:.0f}', 'TVD at the bottom, m')
    print(f'formation {Tf[-1]:.1f} (top) to {Tf[0]:.1f} C (bottom)',
          flush=True)

# --- 3.4  the ORC correction pass --------------------------------------
if not have('RefFourCorrPassDT', 'RefThreeCorrPassDT',
            'AdiabaticCorrPassDT'):
    for key, traj in (('RefFour', 's_curve'), ('RefThree', 'vertical')):
        m_, R_, t_ = field(traj)
        with G.field_ground(B, m_, R_, t_):
            r = B.simulate_css_corrected(CASE)
        add(f'{key}CorrPassDT', abs(r['dT_2d_pass']), '{:.1e}',
            'change in T_2d between the two ORC passes, K')
    r = B.simulate_css_corrected(CASE.with_(k_rock=1e-9))
    add('AdiabaticCorrPassDT', abs(r['dT_2d_pass']), '{:.1e}', '')
    print('correction pass done', flush=True)

# --- 3.1  the step-to-time-constant ratio ------------------------------
if not have('TauRatioLast'):
    with G.field_ground(B, md, Rg, tvd):
        rec = B.simulate_css_corrected(CASE, record=True)
    th = np.asarray(rec['charge']['t_hist'], float)
    dt_last = float(th[-1] - th[-2])
    H = rec['charge']['history']
    U = float(np.asarray(H['U_i'], float)[-1, 0])
    E = np.asarray(H['E'], float)[:, 0]
    Tp = np.asarray(H['T_pcm'], float)[:, 0]
    dE, dT = np.abs(np.diff(E)), np.abs(np.diff(Tp))
    ok = (dT > 1e-9) & (dE > 0)
    Cseg = float(np.median(dE[ok] / dT[ok])) if ok.any() else float('nan')
    Kseg = U * 2.0 * np.pi * CASE.r_i * (CASE.L_well / CASE.n_segments)
    tau = Cseg / Kseg
    add('TauLast', tau, '{:.0f}', 'segment time constant C/K, s')
    add('DtLast', dt_last, '{:.0f}', 'last step on the record grid, s')
    add('TauRatioLast', dt_last / tau, '{:.1f}',
        'dt/tau at the last step, first segment, R4 -- NOT a stability '
        'number, see IntegratorScheme')
    add('IntegratorScheme', 'closed form per enthalpy branch', None,
        'advance_segment solves dE/dt = K(T0 - T_pcm(E)) exactly inside '
        'each branch, so there is no explicit-stability limit')
    print(f'dt/tau {dt_last / tau:.1f} (dt {dt_last:.0f} s, '
          f'tau {tau:.0f} s)', flush=True)

# --- 3.1  the integrator check -----------------------------------------
# Not a CFL check. advance_segment is exact inside a branch, so the only
# error is at BRANCH CROSSINGS. Take one segment over the longest step in
# the grid, straddling the solidus, and compare the closed-form answer
# with the same integrator subdivided 2e5 times.
if not have('IntegratorCheck'):
    E_lat = CASE.h_m * 1.0
    C_s, C_l = 2.0e6, 2.2e6
    T_m, T0, K = 150.0 + 273.15, 165.0 + 273.15, 50.0
    worst = 0.0
    for E0 in (-0.35 * E_lat, 0.0, 0.5 * E_lat, 0.999 * E_lat):
        for dt in (8491.0, 3600.0):
            one, _ = B.advance_segment(E0, T_m, T0, K, dt, E_lat,
                                       C_s, C_l, True)
            n = 200000
            E = E0
            h = dt / n
            for _ in range(n):
                E, _ = B.advance_segment(E, T_m, T0, K, h, E_lat,
                                         C_s, C_l, True)
            scale = max(abs(one), abs(E), 1.0)
            worst = max(worst, abs(one - E) / scale)
    add('IntegratorCheck', worst, '{:.1e}',
        'closed form against 2e5 substeps of itself, worst relative '
        'enthalpy error over four branch states and two step lengths')
    print(f'integrator check {worst:.2e}', flush=True)

# --- 3.5  grid convergence ---------------------------------------------
for tag, kw in (('NsBase', dict(n_segments=99, n_times=40)),
                ('NsFine', dict(n_segments=198, n_times=40)),
                ('NtFine', dict(n_segments=99, n_times=80))):
    if have(f'Conv{tag}Loss'):
        continue
    c = CASE.with_(**kw)
    with G.field_ground(B, md, Rg, tvd):
        rr = B.simulate_css_corrected(c)
        kk = B.performance_indices(c, rr)
    add(f'Conv{tag}Loss', 100 * kk['f_formation_loss'], '{:.3f}', 'loss, %')
    add(f'Conv{tag}Ep', kk['kWh_per_m_well'], '{:.4f}', "e'")
    print(f"  n_s={kw['n_segments']} n_t={kw['n_times']}  "
          f"loss {100 * kk['f_formation_loss']:.3f} %  "
          f"e' {kk['kWh_per_m_well']:.4f}", flush=True)

# --- 4.  the ring, measured --------------------------------------------
if not have('GrissomExtentLong'):
    P = G.grissom_positions(CSV)
    Q = P - P.mean(0)
    _, _, vt = np.linalg.svd(Q, full_matrices=False)
    proj = Q @ vt.T
    ext = proj.max(0) - proj.min(0)
    rad = np.hypot(Q[:, 0], Q[:, 1])
    add('GrissomExtentLong', float(max(ext)), '{:.0f}',
        'extent along the long principal axis of the 363 slots, m')
    add('GrissomExtentShort', float(min(ext)), '{:.0f}', 'the short one, m')
    add('GrissomRminOverRmax', float(rad.min() / rad.max()), '{:.2f}',
        'nearest slot to the centroid over the farthest')
    add('GrissomNslots', float(len(P)), '{:.0f}', '')
    print(f'ring {max(ext):.0f} x {min(ext):.0f} m, '
          f'r_min/r_max {rad.min() / rad.max():.2f}', flush=True)

# --- 7.  the casing unknowns at R3 alone -------------------------------
if not have('SensCasingMaxLossRthree'):
    rows = [v for v in json.load(open('extra_runs.json')).values()
            if v.get('kind') == 'casing' and v.get('step') == 'R3'
            and 'error' not in v]
    ref3 = [v for v in rows if v['tag'] == 'ref'][0]
    add('SensCasingMaxLossRthree',
        max(abs(100 * (v['loss'] - ref3['loss'])) for v in rows), '{:.2f}',
        'largest |d loss| over the three unknowns at R3, points')
    add('SensCasingMaxEpRthree',
        max(abs(100 * (v['e_prime'] / ref3['e_prime'] - 1)) for v in rows),
        '{:.2f}', 'and for e-prime, %')

# --- 2.5  what the pumping model does ----------------------------------
add('PumpEta', 'none', None,
    'hydraulic power only; no pump or motor efficiency is applied')
add('FrictionForm', 'Haaland', None,
    'turbulent friction factor; laminar below the transition')
add('BendAllowance', 'two per well', None,
    'one 180 deg turn at the bottom manifold, counted on the fin-tube '
    'leg only; the downcomer leg is called with n_bends = 0')

# --- 2.2  the surface casing above the storage -------------------------
add('SurfaceCasing', 'not represented', None,
    'the model carries ONE string over the whole storage length, with '
    'the annulus fill k_cement; a second, shallower surface string is '
    'not represented')

with open('orange_notes.json', 'w') as f:
    json.dump(NOTE, f, indent=1)
print(f'\n{len(M)} values in {STORE}')
