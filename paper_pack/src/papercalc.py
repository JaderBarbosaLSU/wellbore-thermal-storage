"""Numbers for the enthalpy paper's section 4.1, at v1.3."""
import json, sys
import boreas as B

ROWS = dict(packing='rows', n_islands=4)
CONFIGS = {
    'adiabatic':      dict(n_segments=99, k_rock=1e-9),
    'isolated':       dict(n_segments=99),
    'fan1000_rows':   dict(n_segments=99, n_wells_field=1000., fan_angle_deg=45., **ROWS),
    'array_rows':     dict(n_segments=99, n_wells_field=1000., fan_angle_deg=0.,  **ROWS),
}

def one(name):
    c = B.Case(**CONFIGS[name])
    r = B.simulate_css_corrected(c)
    k = B.performance_indices(c, r)
    rank, T, Eb = r['rank'], r['T'], r['budget']
    hp = r.get('hp')
    out = dict(
        name=name,
        T_2d_realised_C=r['T_2d_realised'] - 273.15,
        T_3e_C=rank['T_3e'] - 273.15,
        T_3e_hot_end_C=rank['T_3e_hot_end'] - 273.15,
        T_3e_binding=rank['T_3e_binding'],
        pinch_ORC_K=rank['pinch'],
        eta_ORC=rank['rank_eff'], y_frac=rank['y_frac'],
        eta_RTE=k['eta_RTE'], loss_frac=r['Q_loss_kJ'] / r['Q_charge_kJ'],
        eta_storage=r['eta_storage'],
        N_wells=r['N_wells'],
        m1_ch=r['m1_ch'], m1_dc=r['m1_dc'],
        Q_dot_in_ORC_kW=Eb['Q_dot_in_ORC'], Q_dot_out_HP_kW=Eb['Q_dot_out_HP'],
        W_dot_el_in_kW=Eb['W_dot_el_in'], m_dot_w_ch=Eb['m_dot_w_ch'],
        rank={k: (v if isinstance(v, str) else float(v))
              for k, v in rank.items()},
        Tdict={k: float(v) for k, v in T.items()},
        hpd={k: (v if isinstance(v, str) else float(v))
             for k, v in r['hp'].items()},
    )
    try:
        ua = B.exchanger_UA(c, rank, r['hp'], T, Eb)
        out['UA'] = {kk: {a: float(b) for a, b in vv.items()}
                     for kk, vv in ua.items()}
    except Exception as e:
        out['UA_error'] = f'{type(e).__name__}: {e}'[:200]
    return out

if __name__ == '__main__':
    res = {}
    try:
        res = json.load(open('papercalc.json'))
    except Exception:
        pass
    for n in sys.argv[1:]:
        if n not in res:
            res[n] = one(n)
            json.dump(res, open('papercalc.json', 'w'), indent=1, default=float)
            print('ok', n, flush=True)
