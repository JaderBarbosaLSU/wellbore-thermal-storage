"""boreas_field.py -- well-field geometry for BOREAS, without editing boreas.py.

Builds the ground resistance seen by a representative storage well from the
actual positions of all N storage wells, at every measured depth (MD), and
runs BOREAS against it.

Geometry
    disc_positions      N wells filling a circle at one well per B0^2
    rows_positions      one straight double row, B0 along a row, w between rows
    grissom_positions   Active+Idle Grissom slots on the ring (WellSTAR)
    trajectories        vertical, or literature S-curves radiating from the
                        island centre (Clark & Barth 1983; Truex & Hunter 1969;
                        Robertson et al. 1987)

Ground models (all evaluated at the age of the field, t_op)
    'per'   perimeter rule, Eq. (Rper) of the paper: max(R_g, N sqrt(pi a t)/(k P))
    'flux'  line sources of EQUAL STRENGTH over the actual positions; the
            field-mean wall temperature rise per unit heat rate
    'temp'  line sources of EQUAL TEMPERATURE: every well at a given MD has the
            same PCM-to-rock difference, the strengths are solved for (edge
            wells lose more), and the field-mean conductance is returned.
            The near-well series resistance (PCM, casing, cement) is inside
            the solve, so the result is consistent with BOREAS' series chain.

Line-source interaction is local and two-dimensional: a target point on well i
interacts with well k through the infinite line source E1(d^2/4at)/(4 pi k),
where d is the SHORTEST distance from the target to well k's trajectory. For
parallel vertical wells this is the classical superposition; for inclined
wells it is the perpendicular distance to the locally straight neighbour.
The ground surface and the well ends are not represented (as in Eq. (Rg)).
"""
import numpy as np
import pandas as pd
from scipy.special import exp1
from scipy.spatial import cKDTree
import scipy.sparse as sp
import scipy.sparse.linalg as spla

FT = 0.3048
# Nominal wellhead spacing, 6 ft and 12 ft, from the Oxy Long Beach
# drilling deck. B0_MEASURED is what the WellSTAR slot coordinates actually
# show for Grissom: median nearest-neighbour 1.65 m over 363 slots, 10 %
# tighter than nominal. The idealised rungs (disc, double row) are built at
# the MEASURED value so that the step from them to the real slot positions
# isolates the LAYOUT and is not contaminated by a spacing change.
B0_NOMINAL, W_DEF = 1.83, 3.66
B0_MEASURED = 1.65
B0_DEF = B0_MEASURED


# ----------------------------------------------------------------- slots --
def disc_positions(N, B0=B0_DEF):
    i = np.arange(N) + 0.5
    R0 = B0 * np.sqrt(N / np.pi)
    r = np.sqrt(i / N) * R0
    th = np.pi * (3 - np.sqrt(5)) * i
    return np.c_[r * np.cos(th), r * np.sin(th)]


def rows_positions(N, B0=B0_DEF, w=W_DEF):
    n1 = int(np.ceil(N / 2)); n2 = N - n1
    x1 = np.arange(n1) * B0; x2 = np.arange(n2) * B0
    P = np.r_[np.c_[x1, np.zeros(n1)], np.c_[x2, np.full(n2, w)]]
    return P - P.mean(0)


def grissom_positions(csv, statuses=("Active", "Idle"), dedupe=0.6, stray=5.0):
    df = pd.read_csv(csv, dtype={"API": str})
    g = df[(df.Island == "Grissom") & df.WellStatus.isin(statuses)]
    lat0 = g.Latitude.mean(); R = 6371000.0
    P = np.c_[np.radians(g.Longitude) * R * np.cos(np.radians(lat0)),
              np.radians(g.Latitude) * R]
    t = cKDTree(P); keep = np.ones(len(P), bool)
    for i, j in sorted(t.query_pairs(dedupe)):
        if keep[i] and keep[j]:
            keep[j] = False
    P = P[keep]
    d, _ = cKDTree(P).query(P, 2)
    P = P[d[:, 1] < stray]                       # drop isolated stray slots
    return P - P.mean(0)


def perimeter_buffer(P, r):
    from shapely.geometry import Point
    from shapely.ops import unary_union
    return unary_union([Point(p).buffer(r, 24) for p in P]).length


# ---------------------------------------------------------- trajectories --
def s_curve(md, kop, build, theta_max, hold_to=None, drop=2.5, theta_end=45.0):
    """Inclination [deg] along MD [m]: vertical to kop, build [deg/100 ft] to
    theta_max, hold; optional drop [deg/100 ft] from hold_to [m] to theta_end."""
    b = build / (100 * FT); dr = drop / (100 * FT)
    inc = np.clip((md - kop) * b, 0.0, theta_max)
    if hold_to is not None and theta_end < theta_max:
        inc = np.where(md > hold_to,
                       np.maximum(theta_max - (md - hold_to) * dr, theta_end), inc)
    return inc


def trajectories(P, L, kind="vertical", ds=5.0, kop=None, build=4.0,
                 theta_max=55.0, hold_to=3000 * FT, center=None):
    """Polylines for all wells: array (N, n, 3) of (x, y, TVD), MD grid md."""
    md = np.arange(0.0, L + 1e-9, ds)
    N = len(P)
    if kind == "vertical":
        X = np.repeat(P[:, None, :], len(md), 1)
        return md, np.concatenate([X, np.broadcast_to(md, (N, len(md)))[..., None]], 2)
    c = P.mean(0) if center is None else center
    az = np.arctan2(P[:, 1] - c[1], P[:, 0] - c[0])
    kop = np.broadcast_to(np.asarray(kop, float), (N,))
    out = np.empty((N, len(md), 3))
    for i in range(N):
        inc = np.radians(s_curve(md, kop[i], build, theta_max, hold_to))
        incm = 0.5 * (inc[1:] + inc[:-1])
        dh = np.r_[0, np.cumsum(np.sin(incm) * ds)]
        dz = np.r_[0, np.cumsum(np.cos(incm) * ds)]
        out[i, :, 0] = P[i, 0] + dh * np.cos(az[i])
        out[i, :, 1] = P[i, 1] + dh * np.sin(az[i])
        out[i, :, 2] = dz
    return md, out


def staggered_kop(P, center=None, outer=(100, 175, 250), inner=(300, 400, 500)):
    """Outer-row slots kick off shallow (long drift goes 'over the top'),
    inner-row slots deep; values cycle along the ring (Clark & Barth 1983)."""
    c = P.mean(0) if center is None else center
    r = np.hypot(*(P - c).T); ang = np.arctan2(*(P - c).T[::-1])
    t = cKDTree(P); _, nb = t.query(P, 12)
    is_out = r >= np.median(r[nb], 1)
    order = np.argsort(ang); kop = np.empty(len(P))
    io = ii = 0
    for j in order:
        if is_out[j]:
            kop[j] = outer[io % len(outer)]; io += 1
        else:
            kop[j] = inner[ii % len(inner)]; ii += 1
    return kop * FT, is_out


# ------------------------------------------------------ ground resistance --
class Ground:
    def __init__(self, k=2.5, rhocp=2.34e6, t_yr=10.0, rb=0.110, cutoff=80.0):
        self.k, self.a = k, k / rhocp
        self.t = t_yr * 365.25 * 86400.0
        self.rb, self.cutoff = rb, cutoff
        self.Rg = np.log(2 * np.exp(-np.euler_gamma / 2) *
                         np.sqrt(self.a * self.t) / rb) / (2 * np.pi * k)

    def kern(self, d):
        return exp1(d ** 2 / (4 * self.a * self.t)) / (4 * np.pi * self.k)

    def pairs(self, md, X, nst=50):
        """For targets at nst MD stations on every well: shortest distance to
        every other well within the cutoff, and the neighbour's station."""
        N, n, _ = X.shape
        st = (np.arange(nst) + 0.5) * md[-1] / nst
        ist = np.clip(np.searchsorted(md, st), 1, n - 1)
        T = X[:, ist, :].reshape(-1, 3)                      # targets
        owner_t = np.repeat(np.arange(N), nst)
        pts = X.reshape(-1, 3)
        owner = np.repeat(np.arange(N), n); idx = np.tile(np.arange(n), N)
        tree = cKDTree(pts)
        rows, cols, dist = [], [], []
        for t0 in range(0, len(T), 400):
            lists = tree.query_ball_point(T[t0:t0 + 400], self.cutoff + 6.0)
            for tt, L in enumerate(lists):
                t_id = t0 + tt
                L = np.asarray(L)
                L = L[owner[L] != owner_t[t_id]]
                if len(L) == 0:
                    continue
                o = owner[L]; j = idx[L]
                p = T[t_id]
                # exact distance to the two polyline segments around each sample
                best = np.full(len(L), np.inf)
                for dj in (-1, 0):
                    a_ = np.clip(j + dj, 0, n - 2)
                    A = X[o, a_]; Bq = X[o, a_ + 1]
                    v = Bq - A; wv = p - A
                    tpar = np.clip((wv * v).sum(1) / np.maximum((v * v).sum(1), 1e-12), 0, 1)
                    dd = np.linalg.norm(wv - tpar[:, None] * v, axis=1)
                    best = np.minimum(best, dd)
                srt = np.lexsort((best, o))
                o_s, b_s, j_s = o[srt], best[srt], j[srt]
                first = np.r_[True, o_s[1:] != o_s[:-1]]
                o_u, d_u, j_u = o_s[first], b_s[first], j_s[first]
                keep = d_u < self.cutoff
                o_u, d_u, j_u = o_u[keep], d_u[keep], j_u[keep]
                s_nb = np.clip(np.searchsorted(st, md[j_u]), 0, nst - 1)
                rows.append(np.full(len(o_u), t_id)); cols.append(o_u * nst + s_nb)
                dist.append(np.maximum(d_u, self.rb))
        return st, np.concatenate(rows), np.concatenate(cols), np.concatenate(dist), N

    @staticmethod
    def stations_for(L, B0=1.65, ratio=10.0, lo=3, hi=50):
        """How many MD stations the superposition can actually support.

        THE METHOD NEEDS THE STATION SPACING TO BE LARGE COMPARED WITH THE
        WELL SPACING. Each station is treated as an infinite line source in
        the plane normal to the trajectory, and that is only meaningful if
        the station represents a segment much longer than the distance to
        the neighbouring wells. Violate it and the solve returns NEGATIVE
        strengths -- wells absorbing heat while hotter than the rock, which
        cannot happen -- and R_g,eff comes out nonsensical or negative.

        Measured: at L = 100 m with 50 stations the spacing is 2.0 m
        against a well spacing of 1.65 m, 8565 of 18150 strengths go
        negative and R_g,eff reads -0.18 m K/W. Cutting to 10 stations
        (11.1 m, ratio 6.7) gives no negatives and 3.16; 5 and 3 stations
        give 3.17 and 3.16, so the answer is converged once the ratio is
        respected. The old fixed nst=50 was silently wrong for every
        trajectory shorter than about 800 m.
        """
        import math
        return int(min(hi, max(lo, math.floor(L / (ratio * B0)) + 1)))

    def profiles(self, md, X, R_near, nst=50, perimeter=None):
        """Field-mean R_g,eff(MD) for 'flux' and 'temp' (and 'per' if given)."""
        st, r, c, d, N = self.pairs(md, X, nst)
        g = self.kern(d)
        M = N * nst
        G = sp.csr_matrix((g, (r, c)), shape=(M, M))
        self_t = self.kern(np.array([self.rb]))[0]
        # equal flux: mean wall rise per unit q
        flux = (self_t + np.asarray(G.sum(1)).ravel()).reshape(N, nst).mean(0)
        # equal temperature: (R_near + self + G) q = 1
        A = (G + G.T) * 0.5 + sp.identity(M) * (R_near + self_t)
        q = spla.spsolve(A.tocsc(), np.ones(M)).reshape(N, nst)
        # V11: every solved strength must be positive. All wells sit at the
        # same temperature and all are hotter than the far field, so every
        # one of them rejects heat; a negative q is the discretisation
        # failing, not a result. Recorded so the caller can assert on it.
        self.n_negative = int((q < 0).sum())
        temp = 1.0 / q.mean(0) - R_near
        out = dict(md=st, flux=np.maximum(flux, self.Rg), temp=np.maximum(temp, self.Rg),
                   flux_raw=flux, temp_raw=temp, q_spread=q.std(0) / q.mean(0),
                   tvd=X[:, np.clip(np.searchsorted(md, st), 0, len(md) - 1), 2].mean(0))
        if perimeter is not None:
            out["per"] = np.maximum(self.Rg, N * np.sqrt(np.pi * self.a * self.t)
                                    / (self.k * perimeter(st)))
        return out


# ------------------------------------------------------------ BOREAS hook --
def run_boreas(B, L, md_prof, Rg_prof, tvd_prof=None, n_segments=99, **case_kw):
    """Run simulate_css_corrected with R_g,eff(MD) and T_f(TVD(MD)).

    Patches B.formation_profile / B.formation_resistance for the duration of
    the call only. The 'depth' passed between them is MD from the wellhead.
    """
    case = B.Case(n_segments=n_segments, L_ft=L / FT, **case_kw)
    fp0, fr0 = B.formation_profile, B.formation_resistance
    R_near = fr0(case) - B.ground_resistance(case)
    md_p = np.asarray(md_prof, float); Rg_p = np.asarray(Rg_prof, float)
    tvd_p = md_p if tvd_prof is None else np.asarray(tvd_prof, float)

    def tvd_of(mdv):
        if tvd_prof is None:
            return np.asarray(mdv, float)
        x = np.r_[0.0, md_p]; y = np.r_[0.0, tvd_p]
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
        mdv = np.asarray(mdv, float)
        return np.where(mdv <= x[-1], np.interp(mdv, x, y), y[-1] + slope * (mdv - x[-1]))

    def fprof(c, n_segments=None, discharge=False):
        n = n_segments or c.n_segments
        xi = (np.arange(n) + 0.5) / n
        mdv = c.L_well * (1.0 - xi)                 # charge frame: 0 = bottom
        T = c.T_surface_C + 273.15 + c.grad_K_per_km / 1000.0 * tvd_of(mdv)
        if discharge:
            return T[::-1].copy(), mdv[::-1].copy()
        return T, mdv

    def fres(c, depth=None):
        if depth is None:
            return fr0(c)
        return R_near + np.interp(np.asarray(depth, float), md_p, Rg_p)

    B.formation_profile, B.formation_resistance = fprof, fres
    try:
        r = B.simulate_css_corrected(case)
        k = B.performance_indices(case, r)
    finally:
        B.formation_profile, B.formation_resistance = fp0, fr0
    Nw = r["N_wells"]
    E_out = (k["W_el_out_kW"] - k["pumping_dc_kW"]) * case.t_dc / Nw   # net, per well
    E_in = (k["W_el_in_kW"] + k["pumping_ch_kW"]) * case.t_ch / Nw
    return dict(L=L, loss=k["f_formation_loss"], RTE=k["eta_RTE"],
                # v1.4 onwards the package reports whether the store
                # actually gives heat back; surface it rather than let a
                # non-delivering rung look like a small positive.
                delivers=bool(k.get("delivers", r.get("delivers", True))),
                eta_storage=k["eta_storage"], e_prime=k["kWh_per_m_well"],
                E_out_kWh=E_out, E_in_kWh=E_in, T2d_C=r["T_2d_realised"] - 273.15,
                N_wells_1MW=Nw, R_near=R_near)


# ======================================================================
#  Making the measured field the DEFAULT, not an add-on.
#
#  A notebook declares the field once and every later cell -- profiles,
#  animation, parametric studies -- inherits it. Without this the
#  idealised R_field(z) leaks back in wherever someone calls
#  simulate_css_corrected directly, which is most places.
# ======================================================================
from contextlib import contextmanager as _ctx


def reference_field(B, case, csv, island='Grissom', trajectory='s_curve',
                    kop_rule='staggered', theta_max=55.0, kop_ft=300.0,
                    L_md=None, verbose=True):
    """Positions, trajectories and R_g,eff(MD) for the design point.

    Returns (md, Rg, tvd, info). Hand the first three to `field_ground`.

    Station counts stride GEOMETRICALLY: stepping down one at a time from
    50 costs ~26 sparse solves at L = 1250 m, slower than the march it
    feeds, and the clean counts are scattered anyway.
    """
    import numpy as _np
    L = float(L_md if L_md is not None else case.L_well)
    if island != 'Grissom':
        raise ValueError(f'no slot data for island {island!r}')
    P = grissom_positions(csv)
    kop = None
    if trajectory == 's_curve':
        kop = (staggered_kop(P)[0] if kop_rule == 'staggered'
               else _np.full(len(P), kop_ft * FT))
    md, X = trajectories(P, L, kind=trajectory, kop=kop, theta_max=theta_max)
    R_near = float(_np.atleast_1d(
        B.formation_resistance(case, _np.array([0.0])))[0]
        - B.ground_resistance(case))
    gr = Ground(k=case.k_rock, rhocp=case.rho_cp_rock,
                t_yr=case.t_operation_yr, rb=case.r_bore)
    nst0 = Ground.stations_for(L)
    cands = sorted({max(3, int(round(nst0 * f)))
                    for f in (1.0, 0.85, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)},
                   reverse=True)
    pr = nst = None
    for cand in cands:
        pr = gr.profiles(md, X, R_near, nst=cand)
        if getattr(gr, 'n_negative', 0) == 0:
            nst = cand
            break
    if nst is None:
        raise RuntimeError(f'no clean station count among {cands} at L={L}')
    info = dict(N=len(P), L_md=L, n_stations=nst, R_near=R_near,
                Rg_mean=float(_np.mean(pr['temp'])),
                Rg_top=float(pr['temp'][0]), Rg_bottom=float(pr['temp'][-1]),
                tvd_bottom=float(pr['tvd'][-1]), island=island,
                trajectory=trajectory, kop_rule=kop_rule,
                theta_max=theta_max, positions=P, kop=kop)
    if verbose:
        print(f"field: {info['N']} {island} slots, {trajectory}, "
              f"L = {L:.0f} m MD -> TVD {info['tvd_bottom']:.0f} m")
        print(f"  R_g,eff {info['Rg_top']:.2f} (top) to "
              f"{info['Rg_bottom']:.2f} (bottom), mean {info['Rg_mean']:.2f} "
              f"m K/W over {nst} stations")
    return pr['md'], pr['temp'], pr['tvd'], info


@_ctx
def field_ground(B, md, Rg, tvd):
    """Run the model against a measured field, not the idealised R_field(z).

    Patches formation_profile and formation_resistance for the duration, so
    ANY call inside -- including simulate_css_corrected(record=True), which
    is what the profile and animation figures need -- sees the
    superposition resistance and the true vertical depth. Restored on exit
    even if the body raises.
    """
    import numpy as _np
    fp0, fr0 = B.formation_profile, B.formation_resistance
    md_a = _np.asarray(md, float)
    Rg_a = _np.asarray(Rg, float)
    tvd_a = _np.asarray(tvd, float)
    x = _np.r_[0.0, md_a]
    y = _np.r_[0.0, tvd_a]
    slope = (y[-1] - y[-2]) / (x[-1] - x[-2])

    def _tvd(mdv):
        mdv = _np.asarray(mdv, float)
        return _np.where(mdv <= x[-1], _np.interp(mdv, x, y),
                         y[-1] + slope * (mdv - x[-1]))

    def fprof(c, n_segments=None, discharge=False):
        n = n_segments or c.n_segments
        xi = (_np.arange(n) + 0.5) / n
        mdv = c.L_well * (1.0 - xi)          # charge frame: z = 0 is bottom
        T = c.T_surface_C + 273.15 + c.grad_K_per_km / 1000.0 * _tvd(mdv)
        if discharge:
            return T[::-1].copy(), mdv[::-1].copy()
        return T, mdv

    def fres(c, depth=None):
        if depth is None:
            return fr0(c)
        rn = float(_np.atleast_1d(fr0(c, _np.array([0.0])))[0]
                   - B.ground_resistance(c))
        return rn + _np.interp(_np.asarray(depth, float), md_a, Rg_a)

    B.formation_profile, B.formation_resistance = fprof, fres
    try:
        yield
    finally:
        B.formation_profile, B.formation_resistance = fp0, fr0
