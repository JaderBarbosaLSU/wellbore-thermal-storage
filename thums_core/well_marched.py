"""Time-marched single-well model with exact energy closure.

The ONLY difference from `_legacy_physics.temperature_profile_melt` is where the
melt-front position comes from:

    legacy : delta = compute_delta2_fast(...)      independent Stefan solve, bare cylinder
    here   : rho*h_m*dA/dt = q'                    integrated from the delivered heat

Everything else -- internal convection, wall and fouling resistance, fin efficiency,
the melt-layer conduction resistance, the NTU segment march -- is the legacy code,
called unchanged. That isolates the change so the difference in results is
attributable to the front formulation and to nothing else.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import _legacy_physics as L
from . import stefan
from .errors import ThumsPropertyError

CP = L.CP
PI = np.pi


@dataclass
class WellState:
    """Melt cross-sectional area per segment, per tube. Carried across the cycle."""
    A_melt: np.ndarray
    t: float = 0.0

    @classmethod
    def solid(cls, n_segments: int):
        return cls(np.zeros(n_segments))

    def volume(self, dz: float, n_tubes: int) -> float:
        return float(np.sum(self.A_melt * dz)) * n_tubes

    def copy(self):
        return WellState(self.A_melt.copy(), self.t)


@dataclass
class MarchResult:
    state: WellState
    Q_cum_J: float = 0.0
    t: np.ndarray = field(default_factory=lambda: np.array([]))
    Q: np.ndarray = field(default_factory=lambda: np.array([]))
    T_out: np.ndarray = field(default_factory=lambda: np.array([]))
    closure: float = float("nan")
    Q_rejected_J: float = 0.0
    warnings_energy: str = ""
    warnings: list = field(default_factory=list)


_PROP_CACHE: dict = {}


def _fluid_props(fluid, T, P):
    # Properties vary smoothly; caching on 0.05 K bins cuts CoolProp calls by
    # ~50x in a segment march with no visible effect on the result.
    key = (fluid, round(T * 20.0), round(P))
    hit = _PROP_CACHE.get(key)
    if hit is not None:
        return hit
    try:
        rho = CP.PropsSI("D", "T", T, "P", P, fluid)
        mu = CP.PropsSI("V", "T", T, "P", P, fluid)
        k = CP.PropsSI("L", "T", T, "P", P, fluid)
        cp = CP.PropsSI("C", "T", T, "P", P, fluid)
    except ValueError as e:
        raise ThumsPropertyError(f"{fluid} at T={T:.2f} K, P={P:.0f} Pa: {e}") from e
    _PROP_CACHE[key] = (rho, mu, k, cp)
    return rho, mu, k, cp


def _h_i(fluid, T, P, r_i, m_dot):
    rho, mu, k, cp = _fluid_props(fluid, T, P)
    D = 2 * r_i
    Re = 4.0 * m_dot / (PI * D * mu)
    Pr = cp * mu / k
    if Re < 2300:
        Nu = 4.36
    else:
        f = (0.79 * np.log(Re) - 1.64) ** -2
        den = 1 + 12.7 * np.sqrt(f / 8.0) * (Pr ** (2.0 / 3.0) - 1)
        Nu = (f / 8.0) * (Re - 1000) * Pr / den
    return Nu * k / D, cp, Re, Pr


def march(geom, pcm, *, T_inlet, T_m_seg, m_dot, P, fluid, k_wall, Rf_i,
          times, n_segments=100, state=None, mode="charge", strict_bounds=False,
          A_max=None):
    """March one tube of one well over `times`, starting from `state`.

    T_m_seg : melting temperature per segment [K] (the multilayer PCM profile)
    mode    : 'charge' (fluid hotter than PCM) or 'discharge'
    """
    L_tube = geom.L_tube
    dz = L_tube / n_segments
    r_i, r_e = geom.r_i, geom.r_e
    # k_m IS directional: the layer around the tube is liquid while melting and
    # solid while freezing. rho is NOT -- see PCM.rho_latent. Making the latent
    # inventory directional too would hand back rho_s/rho_l = 1.069 times the
    # energy that went in.
    k_m = pcm.k_l if mode == "charge" else pcm.k_s
    rho_m = pcm.rho_latent

    st = WellState.solid(n_segments) if state is None else state.copy()
    # Closure is on the CHANGE in melted area over this march, not on the
    # absolute area -- otherwise a discharge starting from a charged state
    # compares its own extracted heat against the inventory it inherited.
    A0 = st.A_melt.copy()
    res = MarchResult(state=st)
    t_prev = 0.0
    Q_cum = 0.0
    Q_rejected = 0.0
    ts, Qs, Touts = [], [], []

    for t in np.asarray(times, dtype=float):
        dt = t - t_prev
        if dt <= 0:
            t_prev = t
            continue
        delta = stefan.delta_from_area(st.A_melt, r_e)
        T0 = T_inlet
        q_prime = np.zeros(n_segments)
        q_demand = np.zeros(n_segments)
        for i in range(n_segments):
            h_i, cp_d, _, _ = _h_i(fluid, T0, P, r_i, m_dot)
            U_i = L.compute_U_i(h_i, r_i, r_e, k_wall, Rf_i, k_m, float(delta[i]),
                                L_tube, geom.fin_t, geom.fin_L, geom.num_fins)
            NTU = (2 * PI * r_i * U_i * dz) / (m_dot * cp_d)
            NTU = float(np.clip(NTU, -50.0, 50.0))
            T1 = T_m_seg[i] + (T0 - T_m_seg[i]) * np.exp(-NTU)
            q_seg = m_dot * cp_d * (T0 - T1)          # W, positive when melting
            q_demand[i] = q_seg / dz

            # Limit the segment to the latent heat it can actually supply or
            # absorb THIS STEP, and take the fluid temperature from the limited
            # value. Clamping only the front afterwards leaves the fluid warmed
            # (on discharge) by heat no PCM gave up, so downstream segments see
            # an inlet temperature that the energy balance does not support.
            # Only latent heat is modelled: a segment whose melt is exhausted
            # stops contributing rather than cooling sensibly, so `Q_rejected`
            # below is a lower bound on what a sensible-heat model would recover.
            q_avail_melt = st.A_melt[i] * rho_m * pcm.h_m / dt      # W/m, freezing
            if q_seg / dz < -q_avail_melt:
                q_prime[i] = -q_avail_melt
            elif A_max is not None and q_seg / dz > (A_max - st.A_melt[i]) * rho_m * pcm.h_m / dt:
                q_prime[i] = (A_max - st.A_melt[i]) * rho_m * pcm.h_m / dt
            else:
                q_prime[i] = q_seg / dz
            T1 = T0 - q_prime[i] * dz / (m_dot * cp_d)
            T0 = T1

        st.A_melt, q_eff = stefan.advance(st.A_melt, q_prime, dt, rho_m,
                                          pcm.h_m, A_max=A_max)
        st.t = t
        # Accumulate the heat actually taken up as latent heat. Where the segment
        # hit a limit (fully solid on discharge, or fully melted when A_max is
        # set) q_eff < q_prime, and the difference is recorded rather than
        # quietly booked as delivered energy.
        Q_cum += float(np.sum(q_eff) * dz) * dt
        Q_rejected += float(np.sum(q_demand - q_eff) * dz) * dt
        ts.append(t)
        Qs.append(float(np.sum(q_eff) * dz))
        Touts.append(T0)
        t_prev = t

    res.Q_cum_J = Q_cum
    res.Q_rejected_J = Q_rejected
    res.t = np.array(ts)
    res.Q = np.array(Qs)
    res.T_out = np.array(Touts)
    res.closure = stefan.closure_error(Q_cum, st.A_melt - A0, dz, rho_m, pcm.h_m)
    if abs(Q_rejected) > 1e-9 * max(abs(Q_cum), 1.0):
        res.warnings_energy = (
            f"{Q_rejected/1e9:.4f} GJ of the network's heat could not be taken up "
            f"({abs(Q_rejected)/max(abs(Q_cum), 1.0)*100:.2f} % of delivered): "
            "segments reached a physical limit (fully solid, or A_max)")
    # NOTE: L_tube = 2*L_well already covers both legs of ONE hairpin, so the
    # per-borehole multiplier is num_tubes (hairpins), not 2*num_tubes.
    res.warnings = stefan.check_bounds(st.A_melt, r_e, dz, geom.num_tubes,
                                       geom.V_well, geom.D_well, strict=strict_bounds)
    return res
