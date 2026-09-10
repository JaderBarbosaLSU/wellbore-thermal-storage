

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
    # v0.4. True carries PCM enthalpy (latent + sensible in both phases), so a
    # fully melted segment superheats and a fully solid one subcools instead of
    # the march discarding the heat. False pins the PCM at T_m (v0.3 behaviour).
    sensible_heat: bool = True
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

    @property
    def r_cell(self):
        """Radius of the equivalent PCM cell owned by ONE tube leg.

        The 2*num_tubes legs share the borehole cross-section equally, so each
        owns pi R^2 / (2 n_t) of it; writing that as pi r_cell^2 gives

            r_cell = R / sqrt(2 n_t),      R = D_well / 2.

        This is the radius at which a leg's melt front meets its neighbour's.
        It is NOT the borehole wall: with four legs in a 7 in hole r_cell is
        44.45 mm while the wall is at 88.90 mm, so a limit imposed at the wall
        allows delta = 67.82 mm where the real limit is 23.37 mm -- a factor
        2.9 too permissive, and therefore silent at every design point.
        """
        return (self.D_well / 2.0) / np.sqrt(2.0 * self.num_tubes)

    @property
    def delta_merge(self):
        """Melt-layer thickness at which neighbouring fronts touch.

        Hypothesis H3 (concentric annulus, no interaction between legs) is a
        statement about delta < delta_merge. At and beyond this thickness the
        remaining solid sits in the corners between legs, reached through a
        constricted path that the annular resistance ln(1+delta/r_e)/(2 pi k)
        does not represent.

        Under Formulation C this is not a bound that can be VIOLATED --
        area_from_delta(delta_merge) equals A_avail identically, and the
        enthalpy cap holds A_melt <= A_avail -- so delta <= delta_merge always.
        It is a bound the solution can SIT ON, which is what it does here.
        """
        return self.r_cell - self.r_e

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


# ==========================================================================
# enthalpy state   (v0.4: sensible heat in both phases)
# ==========================================================================
# Formulation B tracks A_melt and nothing else, so a segment that runs out of
# PCM has nowhere to put further heat: the march clips it and discards the
# remainder. At the design point the discarded heat on discharge was 148% of
# the heat delivered -- i.e. most of what the network asked for.
#
# The fix is to carry ENTHALPY per unit tube length instead of melted area, with
# the melted area recovered from it. One scalar per segment, measured from the
# fully-solid-at-T_m datum:
#
#     E' < 0            solid, SUBCOOLED     T = T_m + E'/C_s
#     0 <= E' <= E'_lat two-phase            T = T_m,  A_melt = E'/(rho h_m)
#     E' > E'_lat       liquid, SUPERHEATED  T = T_m + (E'-E'_lat)/C_l
#
# Three things follow, and the third is the reason to do it:
#
#   1. no clipping anywhere -- heat always has somewhere to go;
#   2. the melt fraction is bounded in [0,1] BY CONSTRUCTION, so the
#      eps_local > 1 that the latent-only model produced cannot occur;
#   3. desuperheating and subcooling during discharge are recovered, which is
#      the energy the latent-only model threw away.
#
# delta still follows from A_melt exactly as before, so the heat-transfer
# coefficient is unchanged in form -- only its driving temperature is now
# T_pcm rather than T_m.


def pcm_capacities(case):
    """Per unit tube length: latent capacity and the two sensible capacities."""
    A_avail = case.V_well / (case.num_tubes * case.L_tube)   # PCM area per tube [m2]
    E_lat = case.rho_latent * case.h_m * A_avail             # J/m to melt it all
    C_s = case.rho_latent * case.cp_s * A_avail              # J/m/K, solid
    C_l = case.rho_latent * case.cp_l * A_avail              # J/m/K, liquid
    return A_avail, E_lat, C_s, C_l


def pcm_state(E, T_m_seg, case):
    """Enthalpy per unit tube length -> (melted area, PCM temperature).

    The inverse of the enthalpy curve above. `A_melt` saturates at the segment's
    own PCM content, which is what makes eps_local <= 1 structural.
    """
    A_avail, E_lat, C_s, C_l = pcm_capacities(case)
    T_m_seg = np.asarray(T_m_seg, dtype=float)
    A = np.clip(E, 0.0, E_lat) / (case.rho_latent * case.h_m)
    if not case.sensible_heat:                    # v0.3 behaviour: PCM pinned at T_m
        return A, T_m_seg.copy()
    T = np.where(E < 0.0, T_m_seg + E / C_s,
                 np.where(E > E_lat, T_m_seg + (E - E_lat) / C_l, T_m_seg))
    return A, T



def advance_segment(E, T_m, T0, K, dt, E_lat, C_s, C_l, sensible):
    """Integrate ONE segment's enthalpy over dt, exactly, branch by branch.

    The segment obeys   dE'/dt = K (T0 - T_pcm(E'))   with T0 held over the step.
    That is linear inside each branch of the enthalpy curve, so it can be solved
    in closed form rather than stepped:

      plateau  T_pcm = T_m constant           -> dE'/dt constant, E' linear in t
      solid    T_pcm = T_m + E'/C_s           -> exponential relaxation to
      liquid   T_pcm = T_m + (E'-E_lat)/C_l      E'_eq = C (T0 - T_m) [+E_lat]

    Explicit Euler is NOT usable here. On the latent plateau the segment has
    effectively infinite heat capacity, so any step is stable; the moment
    sensible heat is added the capacity becomes finite and the stability limit
    is dt < C/K, about 620 s at the design point. The time grid's last step is
    8491 s, fourteen times that, and an explicit update diverges -- it drove the
    secondary fluid to 261 K on the first attempt.

    Returns (E_new, q_avg) with q_avg = (E_new - E)/dt, so the heat taken from
    the fluid is exactly the enthalpy the PCM gained: closure stays structural.
    """
    if K <= 0.0 or dt <= 0.0:
        return E, 0.0
    E0, t_left = E, dt
    for _ in range(6):   # at most a few branch crossings
        if t_left <= 0.0:
            break
        if not sensible:                     # v0.3: PCM pinned at T_m
            E = E + K * (T0 - T_m) * t_left
            break
        if 0.0 <= E <= E_lat:                # ---- latent plateau: linear ----
            rate = K * (T0 - T_m)
            if rate == 0.0:
                break
            t_edge = ((E_lat - E) / rate) if rate > 0 else (E / -rate)
            if t_edge >= t_left:
                E = E + rate * t_left
                break
            # Land on the boundary and step just PAST it. Without the nudge the
            # next iteration re-enters the plateau at zero rate and the segment
            # sticks at E_lat for ever -- which showed up as every segment
            # reporting a melt fraction of exactly 1.000 and no superheat.
            eps = 1e-9 * max(E_lat, 1.0)
            E = (E_lat + eps) if rate > 0 else -eps
            t_left -= t_edge
        else:                                # ---- sensible: exponential ----
            if E < 0.0:
                C, E_ref, lo, hi = C_s, 0.0, -np.inf, 0.0
            else:
                C, E_ref, lo, hi = C_l, E_lat, E_lat, np.inf
            E_eq = E_ref + C * (T0 - T_m)
            E_new = E_eq + (E - E_eq) * np.exp(-K * t_left / C)
            if lo <= E_new <= hi:
                E = E_new
                break
            edge = hi if E_new > hi else lo   # crosses back onto the plateau
            num, den = E_eq - edge, E_eq - E
            if den == 0.0 or num / den <= 0.0:
                E = edge
                break
            t_edge = -C / K * np.log(num / den)
            if t_edge >= t_left or t_edge <= 0.0:
                E = E_new
                break
            E = edge
            t_left -= t_edge
    return E, (E - E0) / dt


def march_h(case, T_inlet, T_m_seg, m_dot, k_wall, times,
            n_segments=None, E0=None, record=False):
    """Segment march on the ENTHALPY state. Replaces `march` from v0.4.

    Identical to `march` in every respect except what is integrated: the state
    is E' [J per m of tube] rather than A_melt [m2], so no heat is ever clipped
    and the PCM temperature is free to leave T_m.
    """
    n_segments = n_segments or case.n_segments
    dz = case.L_tube / n_segments
    r_i, r_e = case.r_i, case.r_e
    A_avail, E_lat, C_s, C_l = pcm_capacities(case)

    E = np.zeros(n_segments) if E0 is None else np.array(E0, float)
    E_start = E.copy()
    t_prev, Q_cum = 0.0, 0.0
    ts, Qs, T_outs = [], [], []
    hist = {k: [] for k in ("T_fluid", "T_pcm", "A_melt", "delta", "NTU",
                            "U_i", "q_prime", "E")} if record else None

    for t in np.asarray(times, float):
        dt = t - t_prev
        if dt <= 0:
            t_prev = t
            continue

        A, T_pcm = pcm_state(E, T_m_seg, case)
        delta = delta_from_area(A, r_e, case.num_fins, case.fin_t, case.fin_L)
        T0 = T_inlet
        q_prime = np.zeros(n_segments)
        if record:
            T_prof = np.empty(n_segments + 1); T_prof[0] = T_inlet
            NTU_prof = np.empty(n_segments); U_prof = np.empty(n_segments)

        for i in range(n_segments):
            # the layer next to the tube is liquid wherever any melt exists
            k_m = case.k_l if E[i] > 0.0 else case.k_s
            h_i, cp_d, _, _ = h_internal(case.fluid2, T0, case.P, r_i, m_dot)
            U_i = compute_U_i(h_i, r_i, r_e, k_wall, case.Rf_i, k_m,
                              float(delta[i]), case.L_tube, case.fin_t,
                              case.fin_L, case.num_fins)
            NTU = float(np.clip((2 * np.pi * r_i * U_i * dz) / (m_dot * cp_d),
                                -50.0, 50.0))
            # conductance of the segment to the fluid, per unit tube length:
            # q' = K (T0 - T_pcm), from the NTU effectiveness
            K = m_dot * cp_d * (1.0 - np.exp(-NTU)) / dz
            E[i], q_prime[i] = advance_segment(E[i], T_m_seg[i], T0, K, dt,
                                               E_lat, C_s, C_l,
                                               case.sensible_heat)
            # the fluid gives up exactly what the PCM took
            T0 = T0 - q_prime[i] * dz / (m_dot * cp_d)
            if record:
                T_prof[i + 1] = T0; NTU_prof[i] = NTU; U_prof[i] = U_i

        if record:
            hist["T_fluid"].append(T_prof.copy()); hist["T_pcm"].append(T_pcm.copy())
            hist["A_melt"].append(A.copy());       hist["delta"].append(delta.copy())
            hist["NTU"].append(NTU_prof.copy());   hist["U_i"].append(U_prof.copy())
            hist["q_prime"].append(q_prime.copy()); hist["E"].append(E.copy())

        Q_cum += float(np.sum(q_prime) * dz) * dt
        ts.append(t); Qs.append(float(np.sum(q_prime) * dz)); T_outs.append(T0)
        t_prev = t

    A, T_pcm = pcm_state(E, T_m_seg, case)
    dE = float(np.sum(E - E_start) * dz)
    delta_end = delta_from_area(A, r_e, case.num_fins, case.fin_t, case.fin_L)

    # ---- H3 proximity ----------------------------------------------------
    # delta can never EXCEED delta_merge here (see Case.delta_merge), so the
    # useful diagnostic is not a violation flag but how close to the limit the
    # solution runs, and over how much of the well. The annular conduction
    # resistance is least trustworthy where this approaches 1.
    d_merge = case.delta_merge
    prox = delta_end / d_merge if d_merge > 0 else np.full_like(delta_end, np.nan)

    return {
        "E": E, "A_melt": A, "T_pcm": T_pcm,
        "delta": delta_end,
        "Q_cum_J": Q_cum,
        "closure": abs(Q_cum - dE) / abs(Q_cum) if Q_cum else np.nan,
        "t": np.array(ts), "Q": np.array(Qs), "T_out": np.array(T_outs),
        "V_melt_tube": float(np.sum(A * dz)),
        "eps_local": A / A_avail,
        "delta_merge": d_merge,
        "merge_proximity": prox,
        "merge_proximity_max": float(np.max(prox)),
        "f_near_merge": float(np.mean(prox > 0.90)),
        "E_sensible_J": float(np.sum(np.where(E > E_lat, E - E_lat,
                                              np.where(E < 0, E, 0.0))) * dz),
        "z": np.linspace(0.0, case.L_tube, n_segments),
        "history": {k: np.array(v) for k, v in hist.items()} if record else None,
    }


def well_capacity_kJ(case, T_charge_in, T_discharge_in):
    """Energy ONE well can hold, latent plus sensible -- the capacity criterion.

    This is the quantity the v0.1 notebook computed as `N_wells_ideal` and then
    used only for costing. Written here as an energy per well so it can size the
    field directly:

        E_well = rho V_well [ h_m + cp_l (T_charge_in - T_m) ]

    The sensible term is the zero-resistance limit: with no thermal resistance
    the liquid would reach the charging fluid inlet temperature. `cycle=True`
    adds the solid subcooling the discharge can also reach, which is the full
    swing between the charged and discharged states.
    """
    m_well = case.rho_latent * case.V_well
    E = case.h_m
    if case.sensible_heat:
        E += case.cp_l * max(T_charge_in - case.T_m, 0.0)
    return m_well * E / 1000.0                      # kJ
