"""Run the parametric grids on the Grissom field. Resumable, cached.

One superposition solve per storage length, reused by every point that
shares it (see par_grids). Each point is then a full CSS march inside
`field_ground`, so the rock temperature is on TRUE VERTICAL DEPTH and the
resistance is the measured R_g,eff(MD) -- not the idealized field the old
parametric notebook used.
"""
import json, os, sys, time
import numpy as np
import boreas as B
import boreas_ground as G
import par_grids as A

CSV = 'THUMS_wells_WellSTAR_clean.csv'
CACHE = 'par_runs.json'
THUMS = dict(D_in=201.2/25.4, r_casing_o=0.21910/2, r_bore=0.31120/2,
             k_cement=0.7)
BASE = dict(T_m_C=150.0, h_m=380_000.0, k_s=0.60, k_l=0.45, N_lay=9,
            DT_3C_2C=55.0, n_ft=3, num_fins=24, fin_L=0.0035, fin_t=0.0015,
            t_ins=0.008, k_ins=0.007, T_surface_C=18.0, grad_K_per_km=45.0,
            k_rock=2.5, rho_cp_rock=2.34e6, t_operation_yr=10.0,
            W_dot_el_out=1000.0, t_ch=10.0, t_dc=10.0, T_source_C=60.0,
            T_sink_C=20.0, P=10.0e5, DT_pinch_ORC=5.0, DT_2H_3C=10.0,
            n_segments=99, **THUMS)

_FIELD = {}


def field_for(L):
    """R_g,eff(MD) and TVD(MD) for the R4 geometry at this length."""
    if L not in _FIELD:
        c = B.Case(**{**BASE, 'L_ft': L / A.FT})
        _FIELD[L] = G.cached_reference_field(
            B, c, CSV, trajectory='s_curve', kop_rule='staggered',
            theta_max=55.0, kop_ft=300.0, L_md=L, ds=G.DS_STATION,
            verbose=False)[:3]
    return _FIELD[L]


def auto_P(case):
    """Lowest secondary pressure that keeps the water liquid with margin.

    Sweeping T_m at a fixed 10 bar silently runs into saturation: water
    boils at 179.9 C there, so the top of the cascade has nowhere to go
    and the point is infeasible for a reason that has nothing to do with
    the storage. auto_P raises the pressure instead, so the T_m axis
    measures the PCM and not the loop pressure.
    """
    from CoolProp.CoolProp import PropsSI
    T_hot = case.T_m_C + 273.15 + 20.0
    p = PropsSI('P', 'T', min(T_hot, 640.0), 'Q', 0, 'Water') * 1.10
    return max(case.P, float(p))


def run_point(pt):
    kw = dict(BASE)
    L = pt['L_md']
    kw['L_ft'] = L / A.FT
    use_auto = kw.pop('auto_P', None) or pt.get('auto_P', False)
    for k, v in pt.items():
        if k in ('L_md', 'auto_P'):
            continue
        kw[k] = v
    case = B.Case(**kw)
    if use_auto:
        case = case.with_(P=auto_P(case))
    md, Rg, tvd = field_for(L)
    with G.field_ground(B, md, Rg, tvd):
        r = B.simulate_css_corrected(case)
        k = B.performance_indices(case, r)
    out = dict(pt)
    out.update(
        P_bar=case.P / 1e5,
        loss=k['f_formation_loss'], eta_storage=r['eta_storage'],
        RTE=k['eta_RTE'], e_prime=k['kWh_per_m_well'],
        eps_RTE=k['eps_RTE'], eta_T=k['eta_T'], eta_T_eff=k['eta_T_eff'],
        T2d_C=r['T_2d_realized'] - 273.15,
        T3e_C=r['rank']['T_3e'] - 273.15, eta_ORC=r['rank']['rank_eff'],
        wells_per_MWe=k['wells_per_MWe'], N_geo_per_MWe=k['N_geo_per_MWe'],
        rho_E=k['rho_E_kWh_m3'], f_pump=k['f_pump'], f_bypass=k['f_bypass'],
        dE_therm_kWh=k['dE_therm_kWh'], delivers=bool(k['delivers']),
        tvd_bottom=float(tvd[-1]), Rg_mean=float(np.mean(Rg)),
        cycles=int(r['cycles']),
        span_K=float(case.DT_cascade * (case.N_lay - 1) / case.N_lay)
        if getattr(case, 'DT_cascade', None) else
        float(case.DT_3C_2C * (case.N_lay - 1) / case.N_lay))
    return out


def key_of(name, pt):
    return json.dumps([name, {k: pt[k] for k in sorted(pt)}], sort_keys=True)


def main(budget=150.0):
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    todo = [(n, p) for n, g in A.grids().items() for p in g
            if key_of(n, p) not in cache]
    print(f'{len(cache)} cached, {len(todo)} to run', flush=True)
    t0 = time.time()
    for name, pt in todo:
        if time.time() - t0 > budget:
            print('BUDGET', flush=True)
            break
        try:
            r = run_point(pt)
            print(f"  {name:12s} {str({k: v for k, v in pt.items()}):62s} "
                  f"loss={100*r['loss']:6.2f}% e'={r['e_prime']:.4f}",
                  flush=True)
        except Exception as e:
            r = dict(pt, error=f'{type(e).__name__}: {e}'[:200])
            print(f"  {name:12s} ERROR {r['error'][:90]}", flush=True)
        cache[key_of(name, pt)] = r
        json.dump(cache, open(CACHE, 'w'), indent=0, default=float)
    left = sum(1 for n, g in A.grids().items() for p in g
               if key_of(n, p) not in cache)
    print(f'done; {len(cache)} cached, {left} remaining', flush=True)
    return left


if __name__ == '__main__':
    sys.exit(0 if main(float(sys.argv[1]) if len(sys.argv) > 1 else 150.0)
             == 0 else 1)
