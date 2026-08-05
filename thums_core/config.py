"""Configuration objects.

Frozen dataclasses. Nothing in the package takes a positional list of sixteen
arguments -- that signature is what allowed `k_m_l` to be passed where `k_w`
belonged for months without contradiction.

`Geometry.k_wall` is read by BOTH charge and discharge. During the strangler
phase the legacy call chain still receives whatever the notebook passed, so the
defect is faithfully reproduced; `Case.legacy_charge_k_wall` selects which value
goes into the legacy charging path, so the two variants can be run and compared
rather than argued about.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

FT = 0.3048
IN = 0.0254


@dataclass(frozen=True)
class PCM:
    T_m_C: float = 150.0
    rho_s: float = 1550.0
    rho_l: float = 1450.0
    cp_s: float = 1280.0
    cp_l: float = 1800.0
    k_s: float = 0.60
    k_l: float = 0.45
    h_m: float = 380_000.0
    cost_per_kWh: float = 20.0
    name: str = "Saher 2024 trimodal"

    @property
    def T_m(self) -> float:
        return self.T_m_C + 273.15

    def stefan(self, dT: float) -> float:
        """Ste = cp*dT/h_m -- the quasi-steady melt-layer assumption degrades as this grows."""
        return self.cp_l * dT / self.h_m


@dataclass(frozen=True)
class Geometry:
    L_ft: float = 5000.0
    D_in: float = 7.0
    r_i: float = 0.01623
    r_e: float = 0.02108
    fin_t: float = 0.0015
    fin_L: float = 0.0075
    num_fins: int = 24
    num_tubes: int = 2          # hairpins, so 2*num_tubes legs in the borehole
    k_wall: float = 45.0        # steel -- used in BOTH directions
    Rf_i: float = 1e-3
    geo_grad: float = 0.05
    T_geo_surface_C: float = 18.0

    @property
    def L_well(self) -> float:
        return self.L_ft * FT

    @property
    def L_tube(self) -> float:
        return 2.0 * self.L_well

    @property
    def D_well(self) -> float:
        return self.D_in * IN

    @property
    def V_borehole(self) -> float:
        return np.pi * self.D_well ** 2 / 4.0 * self.L_well

    @property
    def A_cs_fins(self) -> float:
        return self.num_fins * self.fin_t * self.fin_L

    @property
    def V_well(self) -> float:
        """PCM volume: borehole minus tubes and fins."""
        n_legs = 2 * self.num_tubes
        A_tube = np.pi * self.r_e ** 2
        return self.V_borehole - n_legs * (A_tube + self.A_cs_fins) * self.L_well

    def legacy_vector(self) -> list:
        """The geom_par_vector the legacy physics expects.

        NOTE: v0.1 (`Jan_12b`) unpacks ELEVEN elements. The `Feb_22` branch
        unpacks thirteen, having appended `geo_grad` and `T_geo_baseline_K` --
        which it then computes, stores in the results frame, and never uses in
        any energy balance (audit 2.4). So the two branches are not
        interchangeable at this interface, and the geothermal fields below exist
        for the ported model, not for the legacy call.
        """
        return [self.L_well, self.L_tube, self.D_well, self.r_i, self.r_e,
                2 * self.r_i, 2 * self.r_e, self.fin_t, self.fin_L,
                self.num_fins, self.num_tubes]

    def validate(self, numerics: "Numerics", strict: bool = False) -> list[str]:
        """Check geometric admissibility.

        Returns a list of warnings. With strict=True these become exceptions.

        strict is False during the strangler phase so that the legacy behaviour
        is reproduced exactly and the frozen fixture stays valid. It flips to
        True once the port is complete -- at which point `delta_max = 0.5 m`
        against an 0.089 m borehole radius stops being permitted.
        """
        from .errors import ThumsGeometryError
        issues = []
        r_tip = self.r_e + self.fin_L
        if r_tip >= self.D_well / 2:
            issues.append(f"fin tip radius {r_tip:.4f} m exceeds borehole radius "
                          f"{self.D_well/2:.4f} m")
        if numerics.delta_max > self.D_well / 2:
            issues.append(
                f"delta_max = {numerics.delta_max} m exceeds the borehole radius "
                f"{self.D_well/2:.4f} m, so the melt front is not bounded by the "
                "well (audit 2.5)")
        if strict and issues:
            raise ThumsGeometryError("; ".join(issues))
        return issues


@dataclass(frozen=True)
class Cycle:
    fluid: str = "cyclopentane"
    fluid2: str = "Water"
    refrig: str = "cyclopentane"
    P: float = 1e6
    W_dot_el_out: float = 1000.0
    Turb_eff: float = 0.85
    ElG_eff: float = 0.95
    Comp_eff: float = 0.85
    ElH_eff: float = 0.95
    # discharge side
    T_sink_C: float = 20.0
    DT_E_sink: float = 5.0
    DT_2D_3E: float = 12.0
    DT_m_2D: float = 0.0
    # charge side
    T_source_C: float = 60.0
    DT_4C_M: float = 10.0
    DT_3A_4A: float = 10.0
    DT_3A_13H: float = 10.0
    DT_2H_3C: float = 10.0
    DT_sub: float = 2.0
    DT_3C_2C: float = 55.0        # secondary-fluid glide; DT_2D_3D is tied to it
    t_ch: float = 10.0
    t_dc: float = 10.0
    orc_stages: int = 2
    hthp_stages: int = 2
    loss_surplus: float = 0.05    # legacy flat storage loss

    @property
    def DT_2D_3D(self) -> float:
        """Tied to DT_3C_2C by construction.

        The legacy code carried these as two independent variables with a comment
        saying they must be equal -- which is a comment, not a constraint.
        """
        return self.DT_3C_2C


@dataclass(frozen=True)
class Numerics:
    n_segments: int = 100
    n_times: int = 40
    delta_max: float = 0.5
    delta_tol: float = 1e-5
    delta_maxiter: int = 20
    N_wells_bracket: tuple = (3.0, 120.0)
    ratio_bracket: tuple = (0.2, 2.0)
    tol_Q_ratio: float = 1e-3
    max_iterations: int = 40


@dataclass(frozen=True)
class Case:
    pcm: PCM = field(default_factory=PCM)
    geometry: Geometry = field(default_factory=Geometry)
    cycle: Cycle = field(default_factory=Cycle)
    numerics: Numerics = field(default_factory=Numerics)
    N_lay: int = 9
    # Strangler-phase switch. False reproduces the legacy defect (PCM liquid
    # conductivity used as the tube-wall/fin conductivity during charging);
    # True uses Geometry.k_wall on both sides, which is the physical answer.
    charge_uses_wall_conductivity: bool = False
    # strict=False reproduces legacy behaviour (warnings, not errors) so the
    # frozen fixture stays valid during the port. Flip to True afterwards.
    strict: bool = False

    def with_(self, **kw) -> "Case":
        """Return a copy with overrides, e.g. case.with_(N_lay=12)."""
        top = {k: v for k, v in kw.items()
               if k in ("N_lay", "charge_uses_wall_conductivity", "strict")}
        rest = {k: v for k, v in kw.items() if k not in top}
        obj = replace(self, **top) if top else self
        for group in ("pcm", "geometry", "cycle", "numerics"):
            g = getattr(obj, group)
            fields = {k: v for k, v in rest.items() if hasattr(g, k)}
            if fields:
                obj = replace(obj, **{group: replace(g, **fields)})
                rest = {k: v for k, v in rest.items() if k not in fields}
        if rest:
            raise KeyError(f"unknown configuration fields: {sorted(rest)}")
        return obj

    def validate(self) -> list[str]:
        return self.geometry.validate(self.numerics, strict=self.strict)

    def key(self) -> str:
        """Stable hash for caching. Changes whenever any input changes."""
        import hashlib
        import json
        from dataclasses import asdict
        blob = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(blob.encode()).hexdigest()[:16]


BASELINE = Case()   # the v0.1 reference point
