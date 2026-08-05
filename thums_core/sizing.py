"""Well-field sizing under BOTH binding constraints.

v0.1 sizes on one criterion only: solve `Q_ratio_ch = 1`, i.e. *can the wells
absorb the energy within the charging window?* That is a heat-transfer question.
Nothing anywhere requires the wells to **contain** enough PCM to hold that energy.

The inventory number exists in the legacy code -- the baseline cell prints
`N_wells_ideal = 12.2` as "Ideal number of wells" -- and is then never used again.

With the legacy `k_w` defect the heat-transfer requirement happened to be the
larger of the two, so the missing constraint never bound and the omission stayed
invisible. Correct the conductivity and the heat-transfer criterion collapses
towards N ~ 1.5, where the model melts more PCM than exists.

    N_wells = max(N_heat_transfer, N_inventory)

`N_inventory` is the smallest N for which the melted fraction eps(N) <= 1. It is
well defined only once the melt front conserves energy -- see `thums_core.stefan`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import well_marched as W
from .errors import ThumsConvergenceError


@dataclass
class SizingResult:
    N_wells: float
    N_heat: float
    N_inventory: float
    binding: str          # 'heat_transfer' | 'inventory'
    eps_pcm: float
    Q_ratio: float
    converged: bool = True
    iterations: int = 0
    warnings: tuple = ()


def _evaluate(N, *, geom, pcm, cy, num, T_inlet, T_m_seg, m_dot_total, times,
              k_wall, n_segments, mode="charge"):
    m1 = m_dot_total / (N * geom.num_tubes)
    r = W.march(geom, pcm, T_inlet=T_inlet, T_m_seg=T_m_seg, m_dot=m1, P=cy.P,
                fluid=cy.fluid2, k_wall=k_wall, Rf_i=geom.Rf_i, times=times,
                n_segments=n_segments, mode=mode)
    dz = geom.L_tube / n_segments
    eps = r.state.volume(dz, geom.num_tubes) / geom.V_well
    Q_total = r.Q_cum_J * geom.num_tubes * N / 1000.0        # kJ, whole field
    return eps, Q_total, r


def _bisect(f, lo, hi, tol, max_iter, what):
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        # monotone and same sign: the root is outside the bracket
        raise ThumsConvergenceError(what, iterations=0,
                                    last=(lo if abs(flo) < abs(fhi) else hi))
    for it in range(1, max_iter + 1):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol or (hi - lo) < 1e-3:
            return mid, it
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    raise ThumsConvergenceError(what, iterations=max_iter, last=0.5 * (lo + hi))


def size_well_field(*, geom, pcm, cy, num, T_inlet, T_m_seg, m_dot_total,
                    D_E_required_kJ, times, k_wall, n_segments=60):
    """Return the well count that satisfies heat transfer AND inventory."""
    kw = dict(geom=geom, pcm=pcm, cy=cy, num=num, T_inlet=T_inlet,
              T_m_seg=T_m_seg, m_dot_total=m_dot_total, times=times,
              k_wall=k_wall, n_segments=n_segments)
    lo, hi = num.N_wells_bracket

    # --- inventory: smallest N with eps(N) <= 1 -----------------------------
    # eps falls monotonically with N (same total energy spread over more PCM).
    def f_inv(N):
        eps, _, _ = _evaluate(N, **kw)
        return eps - 1.0

    try:
        N_inv, it_inv = _bisect(f_inv, lo, hi, 1e-3, num.max_iterations,
                                "inventory constraint eps(N) = 1")
    except ThumsConvergenceError:
        eps_lo, _, _ = _evaluate(lo, **kw)
        if eps_lo <= 1.0:
            N_inv, it_inv = lo, 0          # inventory never binds
        else:
            raise
    # --- heat transfer: Q_delivered(N) >= required --------------------------
    def f_heat(N):
        _, Q, _ = _evaluate(N, **kw)
        return D_E_required_kJ / Q - 1.0

    try:
        N_heat, it_heat = _bisect(f_heat, lo, hi, num.tol_Q_ratio,
                                  num.max_iterations, "heat-transfer sizing")
    except ThumsConvergenceError:
        _, Q_lo, _ = _evaluate(lo, **kw)
        if D_E_required_kJ / Q_lo <= 1.0:
            N_heat, it_heat = lo, 0        # heat transfer never binds
        else:
            raise

    N = max(N_heat, N_inv)
    eps, Q_total, r = _evaluate(N, **kw)
    return SizingResult(N_wells=N, N_heat=N_heat, N_inventory=N_inv,
                        binding="inventory" if N_inv >= N_heat else "heat_transfer",
                        eps_pcm=eps, Q_ratio=D_E_required_kJ / Q_total,
                        iterations=it_inv + it_heat, warnings=tuple(r.warnings))
