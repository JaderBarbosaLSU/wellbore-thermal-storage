

# ==========================================================================
# config
# ==========================================================================
# One flat Case object. v0.1 passed sixteen positional arguments, which is how
# the PCM conductivity ended up in the slot meant for the tube wall for months
# without anything contradicting it. Named fields make that mistake visible.

from dataclasses import dataclass, replace

FT = 0.3048
IN = 0.0254


@dataclass(frozen=True)
class Case:
    # ---- PCM (Saher 2024 trimodal) ----
    T_m_C: float = 150.0          # melting temperature            [C]
    rho_s: float = 1550.0         # solid density                  [kg/m3]
    rho_l: float = 1450.0         # liquid density                 [kg/m3]
    cp_s: float = 1280.0          # solid specific heat            [J/kg/K]
    cp_l: float = 1800.0          # liquid specific heat           [J/kg/K]
    k_s: float = 0.60             # solid conductivity             [W/m/K]
    k_l: float = 0.45             # liquid conductivity            [W/m/K]
    h_m: float = 380_000.0        # latent heat of fusion          [J/kg]
    cost_per_kWh: float = 20.0

    # ---- borehole and tube geometry ----
    L_ft: float = 5000.0          # well depth                     [ft]
    D_in: float = 7.0             # borehole diameter              [in]
    r_i: float = 0.01623          # tube inner radius              [m]
    r_e: float = 0.02108          # tube outer radius              [m]
    fin_t: float = 0.0015         # fin thickness                  [m]
    fin_L: float = 0.0075         # fin radial length              [m]
    num_fins: int = 24            # fins per leg
    num_tubes: int = 2            # HAIRPINS per borehole (2 legs each)
    k_wall: float = 45.0          # steel conductivity             [W/m/K]
    Rf_i: float = 1e-3            # internal fouling resistance    [m2K/W]

    # ---- power block and operating conditions ----
    fluid: str = "cyclopentane"   # ORC working fluid
    fluid2: str = "Water"         # borehole secondary fluid
    refrig: str = "cyclopentane"  # heat-pump refrigerant
    P: float = 1e6                # secondary-fluid pressure       [Pa]
    W_dot_el_out: float = 1000.0  # ORC net electrical output      [kW]
    Turb_eff: float = 0.85
    ElG_eff: float = 0.95
    Comp_eff: float = 0.85
    ElH_eff: float = 0.95
    T_sink_C: float = 20.0
    DT_E_sink: float = 5.0
    DT_2D_3E: float = 12.0
    DT_m_2D: float = 0.0
    T_source_C: float = 60.0
    DT_4C_M: float = 10.0
    DT_3A_4A: float = 10.0
    DT_3A_13H: float = 10.0
    DT_2H_3C: float = 10.0
    DT_sub: float = 2.0
    DT_3C_2C: float = 55.0        # secondary-fluid glide          [K]
    t_ch: float = 10.0            # charging duration              [h]
    t_dc: float = 10.0            # discharging duration           [h]
    loss_surplus: float = 0.05    # assumed storage loss, lambda
    N_lay: int = 9                # PCM layers along the well

    # ---- numerics ----
    n_segments: int = 100
    n_times: int = 40
    delta_max: float = 0.5
    delta_tol: float = 1e-5
    delta_maxiter: int = 20
    N_wells_bracket: tuple = (3.0, 120.0)
    ratio_bracket: tuple = (0.2, 2.0)
    tol_Q_ratio: float = 1e-3
    max_iterations: int = 40

    # ---- switches ----
    front: str = "energy_balance"        # or "closed_form" (= v0.1)
    # CORRECTED 2026-09. Versions of this model up to and including the IHTC
    # paper passed the PCM LIQUID conductivity (0.45 W/m/K) into the slot that
    # carries the tube-wall and fin conductivity during CHARGING, where steel
    # (45 W/m/K) belongs -- two orders of magnitude low. The error is visible in
    # the original solver signature, whose parameter is literally named `k_m_l`.
    #
    # Its main effect is not on the wall resistance, which is small either way,
    # but on the fin efficiency: it collapsed eta_f to about 0.33, throttling the
    # fluid-side conductance to roughly what a bare tube could absorb. That
    # accidentally compensated the finned-vs-bare surface mismatch in the melt
    # front, which is why the energy balance appeared to close.
    #
    # Set False only to reproduce the published v0.1 results (section 13).
    charge_uses_wall_conductivity: bool = True

    # ---- derived ----
    @property
    def T_m(self):
        return self.T_m_C + 273.15

    @property
    def L_well(self):
        return self.L_ft * FT

    @property
    def L_tube(self):
        """Developed length of ONE hairpin: down and back."""
        return 2.0 * self.L_well

    @property
    def D_well(self):
        return self.D_in * IN

    @property
    def V_borehole(self):
        return np.pi * self.D_well ** 2 / 4.0 * self.L_well

    @property
    def V_well(self):
        """PCM volume in one borehole: the hole, less tubes and fins."""
        n_legs = 2 * self.num_tubes
        A_tube = np.pi * self.r_e ** 2
        A_fins = self.num_fins * self.fin_t * self.fin_L
        return self.V_borehole - n_legs * (A_tube + A_fins) * self.L_well

    @property
    def rho_latent(self):
        """Density for the latent inventory, the SAME in both directions.

        The front is carried as a solid-equivalent area, so V_melt/V_well is a
        melted MASS fraction and the energy banked on charging equals the energy
        available on discharging. Using rho_l to melt and rho_s to freeze -- the
        directional choice that is correct for the conductivity k_m -- would
        return rho_s/rho_l = 1.069 times more energy than went in.
        """
        return self.rho_s

    def geom_vector(self):
        """The eleven positional geometry arguments the v0.1 functions expect."""
        return [self.L_well, self.L_tube, self.D_well, self.r_i, self.r_e,
                2 * self.r_i, 2 * self.r_e, self.fin_t, self.fin_L,
                self.num_fins, self.num_tubes]

    def stefan_number(self, dT):
        """Ste = cp dT / h_m. The quasi-steady melt layer degrades as this grows."""
        return self.cp_l * dT / self.h_m

    def with_(self, **kw):
        return replace(self, **kw)


CASE = Case()          # the design point, all corrections applied

# The published configuration. BOTH switches must be set: the closed-form front
# AND the wall-conductivity error, because the two compensated each other.
CASE_V01 = Case(front="closed_form", charge_uses_wall_conductivity=False)


# ==========================================================================
# front_energy_balance   (the v0.2 correction)
# ==========================================================================

def area_from_delta(delta, r_e, n_f=0, t_f=0.0, L_f=0.0):
    """Melted PCM cross-section for a melt layer of thickness `delta`.

    The annulus between r_e and r_e+delta is not all PCM: the fins occupy metal
    inside it. Their cross-section within the front is

        n_f * t_f * min(delta, L_f)

    which saturates once the front passes the fin tips. Excluding it matters
    because A_melt is the quantity the latent-heat balance conserves, so it must
    be PCM and nothing else. At the design point the fins are 6 % of the melted
    area, rising to 11 % when the layer is thin.
    """
    d = np.maximum(delta, 0.0)
    A_annulus = np.pi * ((r_e + d) ** 2 - r_e ** 2)
    A_fin = n_f * t_f * np.minimum(d, L_f)
    return np.maximum(A_annulus - A_fin, 0.0)


def delta_from_area(A_melt, r_e, n_f=0, t_f=0.0, L_f=0.0):
    """Invert `area_from_delta`. Closed form; the fins only shift the quadratic.

    For delta <= L_f the fin metal grows with the front:

        pi d^2 + (2 pi r_e - n_f t_f) d - A = 0

    and for delta > L_f it is constant at n_f t_f L_f:

        pi d^2 + 2 pi r_e d - (A + n_f t_f L_f) = 0

    With n_f = 0 both collapse to  d = sqrt(r_e^2 + A/pi) - r_e,  the bare
    concentric annulus used before this correction.
    """
    A = np.maximum(np.asarray(A_melt, dtype=float), 0.0)
    # branch 1: front still inside the fins
    b = 2.0 * np.pi * r_e - n_f * t_f
    d_in = (-b + np.sqrt(b * b + 4.0 * np.pi * A)) / (2.0 * np.pi)
    # branch 2: front beyond the fin tips
    d_out = np.sqrt(r_e ** 2 + (A + n_f * t_f * L_f) / np.pi) - r_e
    return np.where(d_in <= L_f, d_in, d_out)


def advance_front(A_melt, q_prime, dt, rho_m, h_m, A_max=None):
    """One explicit step of   rho * h_m * dA/dt = q'.

    THIS IS THE WHOLE CORRECTION. v0.1 solved for the front independently of the
    heat the network delivered, using a bare-cylinder Stefan solution, while the
    fluid gave up its heat through a FINNED surface 3.72 times larger. Here the
    front simply follows the delivered heat, so the two can no longer disagree.

    Returns the new area AND the heat actually taken up. They differ wherever a
    segment hits a physical limit: you cannot freeze PCM that is already solid,
    and you cannot melt more than the segment contains. Accumulating q_prime
    instead of q_eff would credit the fluid with heat no PCM supplied.
    """
    A_raw = A_melt + q_prime * dt / (rho_m * h_m)
    A_new = np.maximum(A_raw, 0.0)
    if A_max is not None:
        A_new = np.minimum(A_new, A_max)
    q_eff = (A_new - A_melt) * (rho_m * h_m) / dt
    return A_new, q_eff


def closure_error(Q_in_J, dA_melt, dz, rho_m, h_m):
    """|Q_in - rho h_m dV| / Q_in.

    Zero by construction when advance_front was used -- the two are inverse
    operations. This is an arithmetic check, NOT evidence the model is right.
    """
    V = float(np.sum(dA_melt * dz))
    return abs(Q_in_J - rho_m * h_m * V) / abs(Q_in_J) if Q_in_J else np.nan


# ==========================================================================
# march   (segment-by-segment integration down the tube)
# ==========================================================================

_PROP_CACHE = {}


def fluid_props(fluid, T, P):
    """Density, viscosity, conductivity, cp. Cached on 0.05 K bins."""
    key = (fluid, round(T * 20.0), round(P))
    hit = _PROP_CACHE.get(key)
    if hit is not None:
        return hit
    rho = CP.PropsSI("D", "T", T, "P", P, fluid)
    mu = CP.PropsSI("V", "T", T, "P", P, fluid)
    k = CP.PropsSI("L", "T", T, "P", P, fluid)
    cp = CP.PropsSI("C", "T", T, "P", P, fluid)
    _PROP_CACHE[key] = (rho, mu, k, cp)
    return rho, mu, k, cp


def h_internal(fluid, T, P, r_i, m_dot):
    """Internal convective coefficient. Gnielinski if turbulent, else laminar.

        Re = 4 m_dot / (pi D mu)
        f  = (0.79 ln Re - 1.64)^-2
        Nu = (f/8)(Re-1000) Pr / [1 + 12.7 sqrt(f/8) (Pr^(2/3) - 1)]
    """
    rho, mu, k, cp = fluid_props(fluid, T, P)
    D = 2 * r_i
    Re = 4.0 * m_dot / (np.pi * D * mu)
    Pr = cp * mu / k
    if Re < 2300:
        Nu = 3.66                      # fully developed, isothermal wall
    else:
        f = (0.79 * np.log(Re) - 1.64) ** -2
        Nu = (f / 8.0) * (Re - 1000) * Pr / (
            1 + 12.7 * np.sqrt(f / 8.0) * (Pr ** (2.0 / 3.0) - 1))
    return Nu * k / D, cp, Re, Pr


def layer_map(T_m_lay, n_segments, N_lay):
    """Melting temperature per segment, on the v0.1 node convention."""
    z = np.linspace(0.0, 1.0, n_segments + 1)
    z_lay = np.linspace(0.0, 1.0, N_lay + 1)
    T_m = np.empty(n_segments + 1)
    for j in range(N_lay):
        T_m[(z >= z_lay[j]) & (z < z_lay[j + 1])] = T_m_lay[j]
    T_m[-1] = T_m_lay[-1]
    return T_m[1:]


def march(case, T_inlet, T_m_seg, m_dot, k_wall, times,
          n_segments=None, A_melt0=None, mode="charge", A_max=None,
          record=False):
    """March one tube over `times`, starting from the melt state `A_melt0`.

    Returns a dict. The melt state is returned so the discharge can START from
    what the charge left behind -- v0.1 could not do this, and restarted every
    discharge from a bare tube in fully solid PCM, which is why its recovered
    energy was not bounded by its stored energy.
    """
    n_segments = n_segments or case.n_segments
    dz = case.L_tube / n_segments
    r_i, r_e = case.r_i, case.r_e
    # k_m IS directional: liquid layer while melting, solid while freezing.
    # rho is NOT -- see Case.rho_latent.
    k_m = case.k_l if mode == "charge" else case.k_s
    rho_m = case.rho_latent

    A = np.zeros(n_segments) if A_melt0 is None else np.array(A_melt0, float)
    A0 = A.copy()
    t_prev, Q_cum, Q_rejected = 0.0, 0.0, 0.0
    ts, Qs, T_outs = [], [], []

    # Optional per-segment history, for the diagnostic plots. Recording only
    # appends to lists; it touches nothing the solution depends on, which the
    # notebook verifies by comparing a recorded and an unrecorded run.
    hist = {k: [] for k in ("T_fluid", "A_melt", "delta", "NTU", "U_i",
                            "q_prime", "q_demand")} if record else None

    for t in np.asarray(times, float):
        dt = t - t_prev
        if dt <= 0:
            t_prev = t
            continue
        delta = delta_from_area(A, r_e, case.num_fins, case.fin_t,
                                case.fin_L)
        T0 = T_inlet
        q_prime = np.zeros(n_segments)
        q_demand = np.zeros(n_segments)
        if record:
            T_prof = np.empty(n_segments + 1); T_prof[0] = T_inlet
            NTU_prof = np.empty(n_segments); U_prof = np.empty(n_segments)

        for i in range(n_segments):
            h_i, cp_d, _, _ = h_internal(case.fluid2, T0, case.P, r_i, m_dot)
            U_i = compute_U_i(h_i, r_i, r_e, k_wall, case.Rf_i, k_m,
                              float(delta[i]), case.L_tube, case.fin_t,
                              case.fin_L, case.num_fins)
            NTU = float(np.clip((2 * np.pi * r_i * U_i * dz) / (m_dot * cp_d),
                                -50.0, 50.0))
            T1 = T_m_seg[i] + (T0 - T_m_seg[i]) * np.exp(-NTU)
            q_seg = m_dot * cp_d * (T0 - T1)          # W, +ve when melting
            q_demand[i] = q_seg / dz

            # Limit the segment to the latent heat it can supply or absorb THIS
            # step, and take the fluid temperature from the limited value.
            # Clamping only the front afterwards leaves the fluid warmed by heat
            # no PCM gave up, which starves every segment downstream.
            q_avail = A[i] * rho_m * case.h_m / dt
            if q_demand[i] < -q_avail:
                q_prime[i] = -q_avail
            elif A_max is not None and q_demand[i] > (A_max - A[i]) * rho_m * case.h_m / dt:
                q_prime[i] = (A_max - A[i]) * rho_m * case.h_m / dt
            else:
                q_prime[i] = q_demand[i]
            T0 = T0 - q_prime[i] * dz / (m_dot * cp_d)
            if record:
                T_prof[i + 1] = T0
                NTU_prof[i] = NTU
                U_prof[i] = U_i

        if record:
            hist["T_fluid"].append(T_prof.copy())
            hist["delta"].append(delta.copy())
            hist["NTU"].append(NTU_prof.copy())
            hist["U_i"].append(U_prof.copy())
            hist["q_prime"].append(q_prime.copy())
            hist["q_demand"].append(q_demand.copy())

        A, q_eff = advance_front(A, q_prime, dt, rho_m, case.h_m, A_max)
        if record:
            hist["A_melt"].append(A.copy())     # state AFTER the step
        Q_cum += float(np.sum(q_eff) * dz) * dt
        Q_rejected += float(np.sum(q_demand - q_eff) * dz) * dt
        ts.append(t)
        Qs.append(float(np.sum(q_eff) * dz))
        T_outs.append(T0)
        t_prev = t

    return {
        "A_melt": A,
        "delta": delta_from_area(A, r_e, case.num_fins, case.fin_t, case.fin_L),
        "Q_cum_J": Q_cum,
        "Q_rejected_J": Q_rejected,
        "t": np.array(ts),
        "Q": np.array(Qs),
        "T_out": np.array(T_outs),
        "V_melt_tube": float(np.sum(A * dz)),
        "closure": closure_error(Q_cum, A - A0, dz, rho_m, case.h_m),
        "z": np.linspace(0.0, case.L_tube, n_segments),
        "history": {k: np.array(v) for k, v in hist.items()} if record else None,
    }
