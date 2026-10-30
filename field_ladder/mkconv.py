"""station_convergence.csv from the measured grid."""
import csv, json
rows = json.load(open('station_convergence.json'))
rows.sort(key=lambda r: (r['step'], -r['L_md'], r['nst']))
for r in rows:
    r.setdefault('ds_target', r['ds_m'])
cols = ['step','L_md','ds_target','ds_m','nst','ds_over_B0',
        'n_negative','off_station_weight','Rg_mean','loss','e_prime',
        'RTE','secs']
for r in rows:
    r.setdefault('ds_over_B0', r['ds_m']/1.65)
with open('results_pack_v2/station_convergence.csv','w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
    w.writeheader(); w.writerows(rows)
print(f'wrote results_pack_v2/station_convergence.csv ({len(rows)} rows)')
