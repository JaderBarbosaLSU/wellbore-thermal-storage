"""Resumable parameter sweeps for BOREAS.

WHY A CACHE AND A TIME BUDGET. One cyclic-steady-state run costs about 17 s
and a background process does not survive between tool calls here, so a 7x7
map is ~15 minutes of foreground work and cannot be launched and left. Every
evaluation is therefore written to a JSON cache keyed by the case parameters,
and `run_grid` stops cleanly when its budget is spent. Re-running the same
call picks up exactly where it left off, so a map is accumulated over as many
calls as it needs and no evaluation is ever repeated.

The cache key is the SORTED parameter dict, not a positional tuple, so adding
a swept variable invalidates only the entries that actually change.
"""
import json
import pathlib
import time

CACHE = pathlib.Path("boreas_sweep_cache.json")

# Bump when the metrics or the feasibility gate change. Entries carry the
# version in their key, so a bump invalidates stale results instead of
# silently serving them. v2 added the validate_case gate, without which a
# T_m sweep past 180 C returned numbers computed on SUPERHEATED STEAM: at
# P = 10 bar water saturates at 179.9 C, CoolProp happily returns a vapour
# enthalpy ~10x the liquid one, and the resulting pinch is meaningless. The
# model already had that check; it lives in validate_case, which
# cycle_state_points never calls.
VERSION = "v2"


def _key(kind, kw):
    return f"{VERSION}|{kind}|" + json.dumps(
        {k: kw[k] for k in sorted(kw)}, sort_keys=True, default=str)


def min_water_pressure(model, case, margin_K=10.0):
    """Lowest loop pressure that keeps the water liquid, plus a margin.

    The hottest water in the plant is the charging inlet T_4c, not T_2d, so
    that is what sets the pressure. Returns Pa, rounded up to the next bar.
    """
    import CoolProp.CoolProp as CP
    import math
    T_hot = case.T_m + case.DT_4C_M + margin_K
    return math.ceil(CP.PropsSI("P", "T", T_hot, "Q", 0, "Water") / 1e5) * 1e5


def feasibility(model, case):
    """(ok, first blocking message). Blocking == validate_case 'error'."""
    for lvl, msg in model.validate_case(case, verbose=False):
        if lvl.lower() in ("error", "fatal"):
            return False, msg
    return True, ""


def load(path=CACHE):
    p = pathlib.Path(path)
    return json.loads(p.read_text()) if p.exists() else {}


def save(cache, path=CACHE):
    pathlib.Path(path).write_text(json.dumps(cache, indent=0, default=float))


# --------------------------------------------------------------- metrics
def _prepare(model, kw):
    """Build the Case, optionally auto-pressurising, and gate it.

    `auto_P=True` raises the water-loop pressure to whatever keeps the loop
    liquid at this T_m. That is a real design cost, not a free fix, so it is
    opt-in and the chosen pressure is reported alongside every result.
    """
    kw = dict(kw)
    auto_P = kw.pop("auto_P", False)
    c = model.Case(**kw)
    if auto_P:
        c = model.Case(**{**kw, "P": max(c.P, min_water_pressure(model, c))})
    ok, why = feasibility(model, c)
    return c, ok, why


def cycle_metrics(model, **kw):
    """Cheap: surface cycles only, no borehole march (~6 s)."""
    c, ok, why = _prepare(model, kw)
    if not ok:
        return dict(feasible=False, why=why[:160], P_bar=c.P / 1e5)
    rk, hp, T = model.cycle_state_points(c)
    return dict(feasible=True, P_bar=c.P / 1e5,
                eta_ORC=rk["rank_eff"], COP=hp["hp_cop"],
                product=rk["rank_eff"] * hp["hp_cop"],
                T_2d_C=T["T_2d"] - 273.15, T_3e_C=rk["T_3e"] - 273.15,
                T_3e_binding=rk["T_3e_binding"],
                y_frac=rk["y_frac"], x_8h=hp["x_8h"],
                p_evap_ORC_bar=rk["p_evape"] / 1e5,
                pinch_K=rk["pinch"])


def css_metrics(model, **kw):
    """Full cyclic steady state plus performance indices (~17 s)."""
    c, ok, why = _prepare(model, kw)
    if not ok:
        return dict(feasible=False, why=why[:160], P_bar=c.P / 1e5)
    r = model.simulate_css_corrected(c)
    k = model.performance_indices(c, r)
    Qch = r["Q_charge_kJ"]
    out = dict(feasible=True, P_bar=c.P / 1e5,
               eta_RTE=k["eta_RTE"],
               eps_RTE=k.get("eps_RTE"),
               kWh_per_m_well=k["kWh_per_m_well"],
               f_pump=k.get("f_pump"),
               eta_storage=r["eta_storage"],
               loss_frac=(r["Q_loss_kJ"] / Qch) if Qch else float("nan"),
               N_wells=r["N_wells"],
               T_2d_realised_C=r["T_2d_realised"] - 273.15,
               converged=bool(r["converged"]))
    return {a: b for a, b in out.items() if b is not None}


KINDS = {"cycle": cycle_metrics, "css": css_metrics}


# ----------------------------------------------------------------- driver
def run_grid(model, kind, grid, fixed=None, budget_s=480.0, cache=None,
             path=CACHE, verbose=True):
    """Evaluate `grid` (list of kwargs dicts), resuming and respecting budget.

    Returns (cache, n_done_this_call, n_remaining). A point that raises is
    stored as {'error': ...} rather than retried forever -- an infeasible
    corner is a result, and silently dropping it would leave a hole in the
    map that looks like a bug in the plot.
    """
    fixed = fixed or {}
    cache = load(path) if cache is None else cache
    fn = KINDS[kind]
    t0 = time.time()
    done = 0
    todo = [kw for kw in grid if _key(kind, {**fixed, **kw}) not in cache]
    for kw in todo:
        if time.time() - t0 > budget_s:
            break
        full = {**fixed, **kw}
        k = _key(kind, full)
        try:
            cache[k] = {**kw, **fn(model, **full)}
        except Exception as e:
            cache[k] = {**kw, "error": f"{type(e).__name__}: {e}"[:200]}
        done += 1
        # Save after EVERY point, not at the end. The call can be killed by a
        # hard tool timeout with no chance to clean up, and a sweep that only
        # persists on completion loses the whole call's work when that
        # happens -- which it did, once, for 120 s of runs.
        save(cache, path)
        if verbose:
            print(f"  [{done}/{len(todo)}] {kw}  "
                  f"({time.time() - t0:.0f} s)", flush=True)
    left = len(todo) - done
    if verbose:
        print(f"done {done}, remaining {left}, elapsed "
              f"{time.time() - t0:.0f} s")
    return cache, done, left


def _worker(args):
    """Evaluate one point in a subprocess. Imports the model by NAME.

    The model module cannot be passed through a process boundary, so workers
    import it themselves; `model_name` is whatever `import` string reaches
    the assembled model ('boreas').
    """
    model_name, kind, kw = args
    import importlib
    model = importlib.import_module(model_name)
    try:
        return kw, KINDS[kind](model, **kw)
    except Exception as e:
        return kw, {"error": f"{type(e).__name__}: {e}"[:200]}


def run_grid_parallel(model_name, kind, grid, fixed=None, budget_s=95.0,
                      workers=3, cache=None, path=CACHE, verbose=True):
    """Same contract as `run_grid`, across processes.

    Four cores are available and a CSS point costs ~30 s, so three workers
    turn a 36-point map from nine tool calls into three. Three and not four:
    memory is ~3 GB and each worker carries its own CoolProp tables.

    Results are written as they land, so a hard tool timeout loses at most
    the points still in flight.
    """
    # multiprocessing.Pool, NOT ProcessPoolExecutor: the executor's shutdown
    # joins futures that have already started, so a budget overrun blocks
    # until the slowest in-flight point finishes and the whole tool call is
    # killed instead. Pool.terminate() stops the workers outright, and
    # because every result is written as it lands, nothing is lost.
    import multiprocessing as mp
    fixed = fixed or {}
    cache = load(path) if cache is None else cache
    todo = [kw for kw in grid if _key(kind, {**fixed, **kw}) not in cache]
    t0 = time.time()
    done = 0
    if not todo:
        return cache, 0, 0
    pool = mp.Pool(processes=workers)
    try:
        it = pool.imap_unordered(
            _worker, [(model_name, kind, {**fixed, **kw}) for kw in todo])
        while True:
            left_s = budget_s - (time.time() - t0)
            if left_s <= 0:
                break
            try:
                kw_full, res = it.next(timeout=left_s)
            except mp.TimeoutError:
                break
            except StopIteration:
                break
            kw = {k: v for k, v in kw_full.items() if k not in fixed}
            cache[_key(kind, kw_full)] = {
                **kw, **{k: v for k, v in res.items() if k not in fixed}}
            done += 1
            save(cache, path)
            if verbose:
                print(f"  [{done}/{len(todo)}] {kw} "
                      f"({time.time() - t0:.0f} s)", flush=True)
    finally:
        pool.terminate()
        pool.join()
    save(cache, path)
    left = len(todo) - done
    if verbose:
        print(f"done {done}, remaining {left}, "
              f"elapsed {time.time() - t0:.0f} s")
    return cache, done, left


def collect(kind, grid, fixed=None, path=CACHE):
    """Pull a finished grid out of the cache as a list of row dicts."""
    fixed = fixed or {}
    cache = load(path)
    rows = []
    for kw in grid:
        v = cache.get(_key(kind, {**fixed, **kw}))
        if v is not None:
            rows.append(v)
    return rows


def grid2d(xname, xs, yname, ys):
    return [{xname: float(x), yname: (float(y) if not isinstance(y, int)
                                      else y)} for y in ys for x in xs]
