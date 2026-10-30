import json, time, sys, resource
import pack_runner as P
specs = [s for s in P.build_specs()
         if P.key_of(s) not in json.load(open(P.CACHE))]
print('remaining:', [(s['step'], s['L_md']) for s in specs], flush=True)
for s in specs:
    t0 = time.time()
    try:
        r = P.run_one(s)
    except Exception as e:
        print(s['step'], s['L_md'], 'ERROR', type(e).__name__, e, flush=True)
        continue
    c = json.load(open(P.CACHE)); c[P.key_of(s)] = r
    json.dump(c, open(P.CACHE, 'w'), indent=0, default=float)
    print(f"{s['step']:12s} L={s['L_md']:6.0f} nst={r['n_stations']:3d} "
          f"loss={r['loss']:.4f}  {time.time()-t0:.0f} s  peakRSS="
          f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f} MB",
          flush=True)
print('ALLDONE')
