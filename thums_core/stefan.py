"""Melt-front advance by energy balance, not by an independent closed-form solve.

THE DEFECT THIS REPLACES
------------------------
v0.1 computes the heat leaving the fluid and the position of the melt front from
two *independent* models:

* fluid side  -- `compute_U_i` / `compute_T_re`, using the FINNED outer surface
  `perim_T = 2*pi*r_e + 2*n_fin*fin_L` (0.4924 m at the design geometry) with a
  fin efficiency;
* melt front  -- `compute_delta2_fast`, the quasi-steady Stefan solution for a
  BARE concentric cylinder of radius r_e (perimeter 0.1324 m), driven by
  `T_re - T_m` and assuming `T_re` has been constant since t = 0.

The two surfaces differ by a factor of 3.72, so the heat the fluid gives up and
the latent heat the front absorbs are not the same energy. Measured at the design
point, `Q_in / (rho h_m V_melt)` is 2.02-2.16 with the correct steel conductivity.

THE FORMULATION HERE
--------------------
Only one of those is a free model. The other follows from conservation:

    rho_m * h_m * dA_melt/dt = q'(t)          [W/m of tube]

where q' is whatever the resistance network actually delivers. Integrating this
makes closure exact by construction rather than something to check afterwards,
and the melt-layer thickness follows from the area:

    delta = sqrt(r_e**2 + A_melt/pi) - r_e

Three further defects disappear as a side effect:

1. `T_re` no longer has to be assumed constant from t = 0 -- the state is marched,
   so a varying inlet temperature is handled correctly.
2. The state carries across the cycle. Discharge starts from the melt distribution
   charging left behind, instead of restarting from a bare tube in fully solid PCM.
   The legacy formulation cannot do this at all, which is why round-trip efficiency
   is not bounded by the storage inventory.
3. The fins enter the front advance automatically, because q' comes from the
   finned network.

NOTE ON `h_e = 1e9`
-------------------
The legacy branch `if delta <= 0: h_e = 1e9` is CORRECT as a limit -- a zero-
thickness melt layer has zero conduction resistance. Its comment ("Assign a small
value, leading to large resistance") is inverted, and my earlier audit called the
value wrong; it is not. The real defect is that error paths *route into* this
branch, so a failed solve is reported as maximum heat transfer.
"""

from __future__ import annotations

import numpy as np

from .errors import ThumsGeometryError


def delta_from_area(A_melt, r_e: float):
    """Melt-layer thickness from melted cross-sectional area (concentric annulus)."""
    A = np.maximum(A_melt, 0.0)
    return np.sqrt(r_e ** 2 + A / np.pi) - r_e


def area_from_delta(delta, r_e: float):
    d = np.maximum(delta, 0.0)
    return np.pi * ((r_e + d) ** 2 - r_e ** 2)


def advance(A_melt, q_prime, dt: float, rho_m: float, h_m: float):
    """One explicit step of  rho*h_m*dA/dt = q'.

    q_prime is heat per unit tube length [W/m]; positive melts, negative freezes.
    """
    return np.maximum(A_melt + q_prime * dt / (rho_m * h_m), 0.0)


def check_bounds(A_melt, r_e: float, dz: float, n_tubes_in_hole: int,
                 V_well: float, D_well: float, strict: bool = True):
    """Melt fronts must stay inside the borehole and must not merge.

    Returns a list of messages; raises ThumsGeometryError when strict.
    """
    msgs = []
    delta = delta_from_area(A_melt, r_e)
    r_out = r_e + delta
    if np.any(r_out > D_well / 2.0):
        n = int((r_out > D_well / 2.0).sum())
        msgs.append(f"melt front leaves the borehole in {n} segment(s): "
                    f"max r = {r_out.max()*1000:.1f} mm > {D_well/2*1000:.1f} mm")
    V = float(np.sum(A_melt * dz)) * n_tubes_in_hole
    if V > V_well:
        msgs.append(f"melted volume {V:.3f} m3 exceeds PCM volume {V_well:.3f} m3 "
                    f"(eps = {V/V_well:.3f})")
    if strict and msgs:
        raise ThumsGeometryError("; ".join(msgs))
    return msgs


def closure_error(Q_in_J: float, A_melt, dz: float, rho_m: float, h_m: float) -> float:
    """|Q_in - rho*h_m*V| / Q_in. Zero by construction if `advance` was used."""
    V = float(np.sum(A_melt * dz))
    return abs(Q_in_J - rho_m * h_m * V) / abs(Q_in_J) if Q_in_J else np.nan
