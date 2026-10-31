"""boreas_field.py -- well-field geometry for BOREAS, without editing boreas.py.

Builds the ground resistance seen by a representative storage well from the
actual positions of all N storage wells, at every measured depth (MD), and
runs BOREAS against it.

Geometry
    disc_positions      N wells filling a circle at one well per B0^2
    rows_positions      one straight double row, B0 along a row, w between rows
    grissom_positions   Active+Idle Grissom slots on the ring (WellSTAR)
    trajectories        vertical, or literature S-curves radiating from the
                        island center (Clark & Barth 1983; Truex & Hunter 1969;
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
wells it is the perpendicular distance to the locally straight neighbor.
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
# show for Grissom: median nearest-neighbor 1.65 m over 363 slots, 10 %
# tighter than nominal. The idealized rungs (disc, double row) are built at
# the MEASURED value so that the step from them to the real slot positions
# isolates the LAYOUT and is not contaminated by a spacing change.
# Above this many unknowns the direct sparse solve is replaced by CG; see
# Ground.profiles. 6000 is below the smallest size at which the direct
# factorisation has been seen to exhaust memory, and the two agree to 1e-12
# wherever both run.
_DIRECT_MAX = 6000

# THE station spacing, in meters of measured depth. One number, used by
# every run of the ladder, the L sweep and the sensitivities. See
# Ground.stations_for for why this is a spacing and not a count.
DS_STATION = 50.0

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
        every other well within the cutoff, and the neighbor's station."""
        N, n, _ = X.shape
        st = (np.arange(nst) + 0.5) * md[-1] / nst
        ist = np.clip(np.searchsorted(md, st), 1, n - 1)
        T = X[:, ist, :].reshape(-1, 3)                      # targets
        owner_t = np.repeat(np.arange(N), nst)
        pts = X.reshape(-1, 3)
        owner = np.repeat(np.arange(N), n); idx = np.tile(np.arange(n), N)
        tree = cKDTree(pts)
        rows, cols, dist = [], [], []
        # Chunked, and the chunk is small on purpose. query_ball_point
        # returns a Python list per target, and near the surface an S-curve
        # puts every one of the 363 wells' ~250 samples inside the cutoff --
        # tens of thousands of indices per target. At a chunk of 400 that
        # intermediate is hundreds of megabytes, and four of these running
        # in parallel is what exhausts the machine. Smaller chunks cost a
        # few per cent in Python overhead and bound the peak.
        for t0 in range(0, len(T), 50):
            lists = tree.query_ball_point(T[t0:t0 + 50], self.cutoff + 6.0)
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
                rows.append(np.full(len(o_u), t_id, dtype=np.int32))
                cols.append((o_u * nst + s_nb).astype(np.int32))
                dist.append(np.maximum(d_u, self.rb))
        if not rows:
            # A single well, or a field so spread that nothing is inside
            # the cutoff: no interactions at all. The caller must still get
            # arrays, and the solve then returns the isolated R_g, which is
            # what V3 checks.
            z32 = np.zeros(0, np.int32)
            return st, z32, z32, np.zeros(0), N
        return st, np.concatenate(rows), np.concatenate(cols), np.concatenate(dist), N

    @staticmethod
    def stations_for(L, ds=None):
        """How many MD stations to assemble the superposition on.

        ONE FIXED STATION SPACING, ds, for every run. Not a count, and not
        a per-run search.

        The method treats each station as an infinite line source in the
        plane normal to the trajectory, so a station must represent a
        segment long compared with the distance to the neighboring wells.
        Fix the COUNT instead and that condition changes with L: 50
        stations is 25 m apart in a 1250 m hole and 2 m apart in a 100 m
        one, against a well spacing of 1.65 m. The short case then returns
        negative strengths -- wells absorbing heat from rock they are
        hotter than -- and R_g,eff goes negative with them.

        Fix the SPACING and the condition is the same at every length, the
        stations sit at the same measured depths whatever L is (which is
        what makes V10 meaningful), and every run of a sweep is discretized
        identically.

        ds = 50 m is the shipped value. Measured over the whole L sweep,
        100 to 2000 m: zero negative strengths at every length for R1, R3
        and R4, with the off-station kernel weight rising smoothly from
        0 to 6.4 %. Refining to 25 m introduces negatives at L = 1250 m
        while moving the loss by 0.04 points, so the fine grid is worse
        posed and not more accurate. See station_convergence.csv.
        """
        return int(max(3, round(L / float(ds if ds else DS_STATION))))

    @staticmethod
    def _is_vertical(X, tol=1e-9):
        """True when every trajectory keeps its surface (x, y)."""
        return bool(np.abs(X[:, :, :2] - X[:, :1, :2]).max() <= tol)

    def profiles(self, md, X, R_near, nst=50, perimeter=None):
        """Field-mean R_g,eff(MD) for 'flux' and 'temp' (and 'per' if given).

        VERTICAL FIELDS ARE SOLVED ONCE, NOT nst TIMES. For parallel
        vertical wells the nearest point on a neighbor is at the same
        depth, so the interaction matrix is exactly block diagonal in
        station (measured: off-station kernel weight 0.00 %) and every
        block is the same matrix -- the horizontal distances do not change
        with depth. R_g,eff is then a single number, independent of both
        nst and L, which is what the measured sweep shows (R1: 15.783 m K/W
        at every length and every count; R3: 3.230).

        This is not an approximation, and it is what makes the disc rung
        tractable: at 1.65 m spacing all 363 wells sit inside 18 m, so
        every target has ~90 000 candidate samples inside the 80 m cutoff
        and doing that 25 times over was not finishing.
        """
        if self._is_vertical(X) and nst > 1:
            # One station, on a two-sample stub of the same trajectories.
            # The kernel sees only the horizontal distance, which is the
            # same at every depth for vertical wells, so the stub gives the
            # identical matrix at a fraction of the cost. Its own tvd is
            # meaningless and is discarded; the real one is taken below.
            one = self.profiles(md[:2], X[:, :2, :], R_near, nst=1)
            st = (np.arange(nst) + 0.5) * md[-1] / nst
            rep = lambda v: np.repeat(np.asarray(v).reshape(-1)[:1], nst)
            out = dict(md=st, flux=rep(one['flux']), temp=rep(one['temp']),
                       flux_raw=rep(one['flux_raw']),
                       temp_raw=rep(one['temp_raw']),
                       q_spread=rep(one['q_spread']),
                       tvd=X[:, np.clip(np.searchsorted(md, st), 0,
                                        len(md) - 1), 2].mean(0))
            self.n_negative = one and self.n_negative
            self.off_station_weight = 0.0
            if perimeter is not None:
                out['per'] = np.maximum(
                    self.Rg, N_ := X.shape[0] * np.sqrt(np.pi * self.a * self.t)
                    / (self.k * perimeter(st)))
            return out
        st, r, c, d, N = self.pairs(md, X, nst)
        g = self.kern(d)
        M = N * nst
        self_t = self.kern(np.array([self.rb]))[0]
        # equal flux: mean wall rise per unit q. Row sums only, so this
        # needs no matrix at all.
        # Structural diagnostic, recorded while the entries are in hand:
        # how much kernel weight couples a target to a neighbor at a
        # DIFFERENT station. Zero means the system is exactly block
        # diagonal in station and the answer cannot depend on nst; a large
        # value means the stations have been refined past the point where
        # a station still represents a line source, which is the same
        # failure that produces negative strengths.
        self.off_station_weight = float(
            g[(r % nst) != (c % nst)].sum() / g.sum()) if len(g) else 0.0
        rowsum = np.bincount(r, weights=g, minlength=M)
        flux = (self_t + rowsum).reshape(N, nst).mean(0)
        del rowsum
        # equal temperature: (R_near + G_ii) q + G q = 1.
        #
        # ASSEMBLED AS IT STANDS. Earlier versions solved
        # (G + G^T)/2 instead, and that was wrong. G is not symmetric for
        # inclined wells and has no reason to be: row (i, s) says "the
        # temperature rise at MY location, on well i at station s, due to
        # well k", and the nearest point on k from there is at some other
        # station s'. The reverse entry is about a different target point
        # and a different distance. Averaging the two mixes rows that are
        # not about the same place.
        #
        # Measured here: ||G - G^T|| / ||G|| is 0.25 at 50 m station
        # spacing and 0.90 at 25 m. The symmetrized system returns 130
        # negative strengths at 25 m and 569 at 8 m; the system as
        # assembled returns NONE at any spacing tested, and its mean
        # R_g,eff moves under 0.8 % across a fourfold refinement. The
        # negative strengths were an artifact of the symmetrization, not
        # of the grid.
        #
        # For parallel vertical wells the nearest point IS at the same
        # depth and the distances are mutual, so G is symmetric there and
        # the symmetrization was a no-op -- which is why this went unseen.
        A = sp.coo_matrix(
            (np.concatenate([g, np.full(M, R_near + self_t)]),
             (np.concatenate([r, np.arange(M, dtype=r.dtype)]),
              np.concatenate([c, np.arange(M, dtype=r.dtype)]))),
            shape=(M, M))
        del g, r, c, d
        b = np.ones(M)
        if M <= _DIRECT_MAX:
            q = spla.spsolve(A.tocsc(), b)
        else:
            # NOT CG. CG needs a symmetric positive-definite matrix, and
            # this one is neither symmetric nor diagonally dominant: the
            # diagonal is R_near + G_ii = 0.567 m K/W against a largest
            # off-diagonal row sum of 3.73, because near the surface all
            # 363 wells are inside the 80 m cutoff. BiCGSTAB with a
            # Jacobi preconditioner, checked against the direct solve
            # wherever both run.
            Ac = A.tocsr()
            del A
            Mp = spla.LinearOperator(
                (M, M), matvec=lambda x, dg=Ac.diagonal(): x / dg)
            q, info = spla.bicgstab(Ac, b, rtol=1e-12, atol=0.0,
                                    maxiter=5000, M=Mp)
            if info != 0:
                raise RuntimeError(
                    f'BiCGSTAB did not converge (info={info}) at M={M}')
        q = q.reshape(N, nst)
        # V11: every solved strength must be positive. All wells sit at the
        # same temperature and all are hotter than the far field, so every
        # one of them rejects heat; a negative q is the discretization
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
                E_out_kWh=E_out, E_in_kWh=E_in, T2d_C=r["T_2d_realized"] - 273.15,
                N_wells_1MW=Nw, R_near=R_near)


# ======================================================================
#  Making the measured field the DEFAULT, not an add-on.
#
#  A notebook declares the field once and every later cell -- profiles,
#  animation, parametric studies -- inherits it. Without this the
#  idealized R_field(z) leaks back in wherever someone calls
#  simulate_css_corrected directly, which is most places.
# ======================================================================
from contextlib import contextmanager as _ctx


def reference_field(B, case, csv, island='Grissom', trajectory='s_curve',
                    kop_rule='staggered', theta_max=55.0, kop_ft=300.0,
                    L_md=None, ds=None, verbose=True):
    """Positions, trajectories and R_g,eff(MD) for the design point.

    Returns (md, Rg, tvd, info). Hand the first three to `field_ground`.

    The station count follows from ONE fixed spacing, `ds` (default
    DS_STATION = 50 m). There is no search: V11 is asserted instead.
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
    nst = Ground.stations_for(L, ds)
    pr = gr.profiles(md, X, R_near, nst=nst)
    # V11, as an assertion and not a search. Every well sits at the same
    # temperature and every one of them is hotter than the far field, so
    # every solved strength must be positive. A negative one is the
    # discretization failing, and the response is to report it, not to
    # quietly try a different station count -- that is what made two runs
    # of the same sweep incomparable in v2.0.
    if gr.n_negative:
        raise RuntimeError(
            f'V11 FAILED: {gr.n_negative} of {len(P) * nst} strengths are '
            f'negative at L = {L:.0f} m with ds = {ds or DS_STATION:.0f} m '
            f'({nst} stations). This is a geometry or kernel problem; do '
            f'not work around it by changing the station count.')
    # tvd_bottom is the TVD at the BOTTOM OF THE STORE, md = L, taken
    # from the trajectory. It used to be pr['tvd'][-1], the TVD at the
    # last STATION CENTRE -- at ds = 50 m that is md = 1225 m, so the
    # reference case was reported 18 m too shallow, 891 against 909.
    info = dict(N=len(P), L_md=L, n_stations=nst,
                ds_station=float(ds or DS_STATION),
                tvd_last_station=float(pr['tvd'][-1]),
                off_station_weight=gr.off_station_weight, R_near=R_near,
                Rg_mean=float(_np.mean(pr['temp'])),
                Rg_top=float(pr['temp'][0]), Rg_bottom=float(pr['temp'][-1]),
                tvd_bottom=float(X[:, -1, 2].mean()), island=island,
                trajectory=trajectory, kop_rule=kop_rule,
                theta_max=theta_max, positions=P, kop=kop)
    if verbose:
        print(f"field: {info['N']} {island} slots, {trajectory}, "
              f"L = {L:.0f} m MD -> TVD {info['tvd_bottom']:.0f} m")
        print(f"  R_g,eff {info['Rg_top']:.2f} (top) to "
              f"{info['Rg_bottom']:.2f} (bottom), mean {info['Rg_mean']:.2f} "
              f"m K/W over {nst} stations")
    return pr['md'], pr['temp'], pr['tvd'], info


# The superposition solve is the expensive step in the whole model: 363
# wells x 50 stations is an 18150-unknown sparse system, about 85 s, and
# reference_field tries several station counts before one solves cleanly.
# That is fine in a batch runner and absurd in the first cell of a
# notebook, so the result is cached by the SPECIFICATION it came from.
_REF_URL = ('https://raw.githubusercontent.com/JaderBarbosaLSU/'
            'wellbore-thermal-storage/notebook-model/field_ladder/')


def _ref_key(case, island, trajectory, kop_rule, theta_max, kop_ft, L,
             ds=None):
    """Everything the profile depends on, and nothing else.

    The case enters only through r_bore (the kernel's inner radius), the
    rock properties and the field age; the PCM, the plant and the fin
    geometry do not touch the ground solve, so changing them must not
    invalidate the cache.
    """
    bits = [island, trajectory, kop_rule, f'{theta_max:.4g}',
            f'{kop_ft:.4g}', f'{L:.6g}', f'{case.r_bore:.6g}',
            f'{case.k_rock:.6g}', f'{case.rho_cp_rock:.6g}',
            f'{case.t_operation_yr:.6g}', f'ds{float(ds or DS_STATION):.4g}']
    return 'fieldref_' + '_'.join(bits).replace('.', 'p').replace('+', '')


def cached_reference_field(B, case, csv, cache_dir='.', fetch=True,
                           verbose=True, **kw):
    """`reference_field`, read from disk or from the repo when possible.

    Returns the same (md, Rg, tvd, info). `info` carries 'cached', which is
    the provenance: True means nothing was solved here. A cache miss is not
    an error -- it recomputes and writes the file -- but it costs minutes,
    so it says so first.
    """
    import os
    import urllib.request
    import numpy as _np
    L = float(kw.get('L_md') or case.L_well)
    key = _ref_key(case, kw.get('island', 'Grissom'),
                   kw.get('trajectory', 's_curve'),
                   kw.get('kop_rule', 'staggered'),
                   float(kw.get('theta_max', 55.0)),
                   float(kw.get('kop_ft', 300.0)), L,
                   kw.get('ds'))
    path = os.path.join(cache_dir, key + '.npz')
    if not os.path.exists(path) and fetch:
        try:
            urllib.request.urlretrieve(_REF_URL + key + '.npz', path)
        except Exception:
            pass
    if os.path.exists(path):
        try:
            d = _np.load(path, allow_pickle=True)
            md, Rg, tvd = d['md'], d['Rg'], d['tvd']
            info = dict(d['info'].item())
            info['cached'] = True
            if verbose:
                print(f"field: {info['N']} {info['island']} slots, "
                      f"{info['trajectory']}, L = {L:.0f} m MD -> TVD "
                      f"{info['tvd_bottom']:.0f} m   [cached]")
                print(f"  R_g,eff {info['Rg_top']:.2f} (top) to "
                      f"{info['Rg_bottom']:.2f} (bottom), mean "
                      f"{info['Rg_mean']:.2f} m K/W over "
                      f"{info['n_stations']} stations")
            return md, Rg, tvd, info
        except Exception:
            pass
    if verbose:
        print('no cached field for this specification -- solving the '
              f'superposition.\n  Written to {path} and reused from then '
              'on.', flush=True)
    md, Rg, tvd, info = reference_field(B, case, csv, verbose=verbose, **kw)
    info['cached'] = False
    try:
        store = {k: v for k, v in info.items() if k not in ('positions', 'kop')}
        _np.savez(path, md=md, Rg=Rg, tvd=tvd, info=_np.array(store, dtype=object))
    except Exception:
        pass
    return md, Rg, tvd, info


@_ctx
def field_ground(B, md, Rg, tvd):
    """Run the model against a measured field, not the idealized R_field(z).

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

    # A TVD mapping for the PLOTTING layer. formation_profile's second
    # return value has to stay MEASURED depth, because fres interpolates
    # R_g,eff against it -- but a figure drawn on that axis shows a
    # 1250 m well where the S-curve reaches 891 m. depth_axis picks this
    # up when present, so every figure and the animation get TVD without
    # the solver's own lookup changing at all.
    B._boreas_tvd_of = _tvd

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
        if hasattr(B, '_boreas_tvd_of'):
            del B._boreas_tvd_of
