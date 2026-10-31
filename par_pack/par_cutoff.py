"""Handoff section 10, last row: the reference PCM at L = 500 m, the cutoff.

Two of the paper's checks wait on this. The L grid of the maps does not
contain 500 m, so it is run here rather than interpolated.
"""
import json, os
import par_runner as R, par_grids as A

CACHE = 'par_cutoff.json'
out = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
for L in (500.0, A.L_REF):
    k = f'{L:.0f}'
    if k in out:
        continue
    out[k] = R.run_point(dict(L_md=L, N_lay=9, n_ft=3))
    json.dump(out, open(CACHE, 'w'), indent=1)
    r = out[k]
    print(f"L = {L:6.0f} m   e' {r['e_prime']:.4f}   RTE {r['RTE']:.4f}   "
          f"wells/MWe {r['wells_per_MWe']:6.2f}   loss {100*r['loss']:5.2f} %",
          flush=True)
if len(out) == 2:
    a = out['500']
    print(f"\nisland output at L = 500 m: "
          f"{363.0 / a['wells_per_MWe']:.1f} MWe from the 363 slots")
