import json, os, time
import pack_runner as P
CACHE = P.CACHE
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
specs = P.build_specs()
t0 = time.time()
BUDGET = float(__import__('sys').argv[1]) if len(__import__('sys').argv)>1 else 1e9
for i, s in enumerate(specs):
    k = P.key_of(s)
    if k in cache:
        continue
    if time.time() - t0 > BUDGET:
        print('BUDGET', flush=True); break
    try:
        r = P.run_one(s)
        print(f"  {i+1:3d}/{len(specs)} {s['step']:10s} L={s['L_md']:6.0f} "
              f"{s['ground'][:4]:5s} nst={r['n_stations']:4d} "
              f"loss={100*r['loss']:6.2f}% e'={r['e_prime']:.4f}", flush=True)
    except Exception as e:
        r = dict(s, error=f'{type(e).__name__}: {e}'[:300])
        print(f"  {i+1:3d}/{len(specs)} {s['step']:10s} L={s['L_md']:6.0f} "
              f"ERROR {r['error'][:90]}", flush=True)
    cache[k] = r
    json.dump(cache, open(CACHE, 'w'), indent=0, default=float)
print(f'done in {time.time()-t0:.0f} s; {len(cache)} runs', flush=True)
