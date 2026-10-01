

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
    fin_L: float = 0.0035         # fin radial length              [m]
    num_fins: int = 24            # fins per fin tube
    # ---- THE DOWNCOMER GEOMETRY. This is what BOREAS is for. ------------
    # THUMS put n hairpins in the hole: each a U-tube, down one leg and back
    # up the other, both legs exchanging with the PCM. That geometry is
    # self-contradictory and the contradiction is quantitative.
    #
    # Charge water enters at T_m + DT_4C_M and glides DT_3C_2C down, so the
    # RETURN leg runs ~45 K BELOW the melting point. Graded by depth it
    # refreezes what the down leg just melted. So THUMS graded the cascade
    # along the TUBE PATH instead -- which works, but requires two different
    # melting points to coexist at the same depth, centimetres apart, with
    # contiguous PCM and no wall between them. Conduction between those legs
    # is 43 kW at one hairpin and 92 kW at three, against a ~120 kW duty, and
    # it is a short circuit straight across the cascade the model exists to
    # represent. THUMS never had that term.
    #
    # BOREAS separates the two jobs. n_ft FIN TUBES all run DOWNWARD and all
    # exchange; ONE INSULATED DOWNCOMER returns the whole flow without
    # exchanging. They meet in a manifold at the bottom of the same hole --
    # no lateral, no new drilling, the same class of object as the U-bend it
    # replaces. Consequences:
    #
    #   z becomes DEPTH, one-to-one (L_tube = L_well, not twice it)
    #   the legs are PARALLEL, not series
    #   the cascade grades with depth, so the PCM is POURABLE in layers
    #   the cascade runs HOT AT THE BOTTOM, which aligns it with the
    #     geothermal gradient instead of fighting it -- see cascade_hot_bottom
    #
    # Sizing the downcomer for equal velocity keeps the pressure drop
    # identical to the hairpin's: same velocity everywhere, same total path
    # 2*L_well. The nested-coaxial alternative costs 3.3x on f_pump, which is
    # the one quantity Study 2 showed governs efficiency below ground.
    n_ft: int = 3                 # FIN TUBES per borehole (all downward)
    t_ins: float = 0.008          # downcomer insulation thickness [m]
    k_ins: float = 0.007          # downcomer insulation, VIT       [W/m/K]
    k_wall: float = 45.0          # steel conductivity             [W/m/K]
    Rf_i: float = 1e-3            # internal fouling resistance    [m2K/W]

    # ---- the formation. THUMS hypothesis H6 said the wall was adiabatic,
    # which made eta_storage identically 1 -- a number no reviewer of a
    # borehole storage paper will accept. These are what break it.
    T_surface_C: float = 18.0     # ground surface / mean annual air [C]
    grad_K_per_km: float = 45.0   # geothermal gradient              [K/km]
    k_rock: float = 2.5           # formation conductivity        [W/m/K]
    rho_cp_rock: float = 2.34e6   # formation volumetric heat cap [J/m3/K]
    k_cement: float = 1.0         # annular cement                [W/m/K]
    r_casing_o: float = 0.098     # casing outer radius             [m]
    r_bore: float = 0.110         # drilled hole radius             [m]
    t_operation_yr: float = 10.0  # age of the field when evaluated  [yr]
    # The gradient and the source temperature are specified INDEPENDENTLY and
    # checked against each other, never derived one from the other. 60 C at
    # the wellhead is a LOWER BOUND on the formation temperature -- the water
    # cools coming up -- and the producing wells may be deeper than these.
    # Writing that bound as an equality is what an earlier draft did, and it
    # threw away exactly the physics the deficit represents.
    L_producing_ft: float = 8000.0   # depth the geothermal water comes from

    # ---- THE WELL FIELD IS A CONE, NOT AN ARRAY -------------------------
    # THUMS wellheads sit on 6 ft centres with 12 ft between double rows,
    # and over 1200 wells are drilled DIRECTIONALLY from four islands with
    # drift angles to 84 deg, reaching 6500 acres. So the spacing is not a
    # number, it is a function of depth: metres apart at the wellhead, of
    # order 150 m at reservoir depth.
    #
    # That matters more than any other single parameter here. Neighbours
    # merge thermally in r^2/alpha, which at 6 ft is 36 DAYS and at 160 m is
    # eight centuries. The shallow section is one thermal body within a month
    # and barely loses at all; the deep section is a field of isolated wells.
    #
    # And the two effects run opposite ways with depth. Shallow rock is cold,
    # so the driving difference is largest exactly where the shielding is
    # strongest; deep rock is hot, so the driving difference is smallest
    # exactly where the wells stand alone. Neither a single spacing nor a
    # single resistance can represent that, which is why B is B(z).
    #
    # n_wells_field = None recovers the isolated well, which is the
    # conservative bound and the v1.0 default behaviour.
    n_wells_field: float = None   # wells sharing the formation; None = alone
    B_surface_m: float = 1.83     # wellhead spacing, 6 ft                [m]
    z_kickoff_m: float = 300.0    # vertical above the kickoff point      [m]
    fan_angle_deg: float = 45.0   # effective half-angle of the fan     [deg]

    # ---- power block and operating conditions ----
    fluid: str = "cyclopentane"   # ORC working fluid
    fluid2: str = "Water"         # borehole secondary fluid
    refrig: str = "cyclopentane"  # heat-pump refrigerant
    P: float = 1e6                # secondary-fluid pressure       [Pa]
    W_dot_el_out: float = 1000.0  # ORC net electrical output      [kW]
    # ---- machine efficiencies -------------------------------------------
    # ISENTROPIC efficiencies, applied INSIDE the cycles from v0.12: the
    # expansions, compressions and pumpings move the state points themselves.
    # Before v0.12 the state points were isentropic and 0.85 was applied as a
    # multiplier on the work downstream. That is accurate for the efficiency
    # -- 0.6 % at the design point -- but it leaves every state point ideal,
    # so the turbine exit, the regenerator split, the composite curves and
    # the evaporator pinch were all the reversible cycle's. It also makes a
    # component-wise exergy balance impossible: an isentropic machine
    # destroys nothing, so the lost work appears nowhere. See DN-21.
    #
    # Setting all three to 1.0 recovers the pre-v0.12 cycles exactly.
    eta_turb_s: float = 0.85      # ORC turbines, isentropic
    eta_pump_s: float = 0.85      # ORC pumps, isentropic
    eta_comp_s: float = 0.85      # heat-pump compressors, isentropic
    # ELECTRICAL / MECHANICAL, and correctly applied outside the working
    # fluid: these are not thermodynamic irreversibilities of the cycle.
    ElG_eff: float = 0.95         # ORC generator
    ElH_eff: float = 0.95         # compressor motor
    # RETIRED at v0.12. These were the isentropic efficiencies applied as
    # downstream multipliers; keeping them would double-count. validate_case
    # raises if a Case still names either.
    Turb_eff: float = None
    Comp_eff: float = None
    T_sink_C: float = 20.0
    # Approach at the ORC condenser, measured to the sink OUTLET. With
    # DT_sink_glide = 0 the sink is an infinite reservoir and outlet = inlet,
    # which is the model as built. Referencing it to the inlet -- as this did
    # before v0.9b -- is the same wrong-end pattern that made DT_3A_13H
    # unsafe; it was merely masked by the reservoir assumption. See DN-18.
    DT_E_sink: float = 5.0
    # Temperature rise of the sink stream through the ORC condenser. Zero is
    # the reservoir idealisation. Give it a value and the condenser acquires
    # an INTERIOR pinch at the desuperheat corner, ~93 % of the duty, which
    # crosses at about 5.4 K -- the ORC condenser is not structurally safe
    # the way the HTHP condenser is.
    DT_sink_glide: float = 0.0
    DT_2D_3E: float = 12.0
    # Minimum water-to-working-fluid gap ANYWHERE in the ORC evaporator, not
    # just at its hot end. DT_2D_3E alone cannot keep the two composite curves
    # apart, because the ORC takes most of its heat at one temperature while
    # the water glides; see `orc_pinch`. Raise this and the ORC boils lower
    # and yields less; lower it and the exchanger grows. See DN-15.
    DT_pinch_ORC: float = 5.0
    # Discharge-inlet subcooling below the COLDEST cascade layer, state 1d in
    # the plant diagram; the symmetric partner of DT_4C_M. See DN-8.
    DT_M_1D: float = 10.0
    T_source_C: float = 60.0
    DT_4C_M: float = 10.0
    # How far the SOURCE stream is cooled in the heat-pump evaporator. This is
    # the specification: it sets the source flow needed for a given duty.
    DT_3A_4A: float = 10.0
    # Minimum approach in the heat-pump evaporator, at its COLD end -- where
    # the source leaves and the refrigerant enters. The evaporating
    # temperature is DERIVED from it,
    #     T_13h = T_source_C - DT_3A_4A - DT_pinch_HPE,
    # so the approach can no longer be negative by arithmetic. Unlike the ORC
    # evaporator the cold stream here is isothermal with NO preheat kink, so
    # the gap is monotonic in Q and an end approach is genuinely sufficient:
    # no composite-curve search is needed. See DN-17.
    DT_pinch_HPE: float = 5.0
    # RETIRED at v0.9. The evaporating temperature used to be set by
    # T_13h = T_source_C - DT_3A_13H, independently of DT_3A_4A, so the
    # approach was the DIFFERENCE of two free fields and nothing checked its
    # sign. Left here only so that an old Case naming it fails loudly in
    # validate_case rather than being silently ignored.
    DT_3A_13H: float = None
    DT_2H_3C: float = 10.0
    DT_sub: float = 2.0
    # The CHARGING water glide, across the HTHP condenser. Because
    # T_3c = T_4c -- the water leaves the condenser and enters the borehole
    # with no state change between -- this is identically the borehole
    # charging glide. That equality is energy conservation on a closed loop,
    # not an assumption, and cannot be relaxed.
    DT_3C_2C: float = 55.0        # charging water glide           [K]
    # The DISCHARGING water glide, across the ORC evaporator, and likewise
    # identically the borehole discharging glide. Until v0.9b this was not a
    # field at all: both mass flows were divided by DT_3C_2C, which silently
    # asserted that the two half-cycles glide by the same amount. They need
    # not. None keeps the old behaviour exactly. See DN-18.
    DT_3D_2D: float = None        # discharging water glide        [K]
    # The cascade grading parameter. Layer spacing is DT_cascade / N_lay and
    # the top-to-bottom SPAN is DT_cascade * (N_lay - 1) / N_lay, one layer
    # short of DT_cascade itself.
    #
    # Setting this equal to the charging glide -- which None does -- is not
    # arbitrary: it is the unique choice that makes the approach between the
    # water and the layer it is melting equal to DT_4C_M at the LEADING FACE
    # of every layer, decaying to DT_4C_M - DT_cascade/N_lay at each trailing
    # face. That sawtooth is the cascade.
    #
    # Prescribing it independently is legitimate, but it is bounded BELOW,
    # not above: the approach at the far end of the well is
    #     charging : span - (glide_ch - DT_4C_M)
    #     discharge: span - (glide_dc - DT_M_1D)
    # so a span SMALLER than glide minus the approach inverts the driving
    # difference at that end. At the default the margin is only 3.889 K.
    DT_cascade: float = None      # cascade grading parameter      [K]
    t_ch: float = 10.0            # charging duration              [h]
    t_dc: float = 10.0            # discharging duration           [h]

    # ---- the geothermal resource (new at v0.13) --------------------------
    # Until v0.13 the 60 C source was effectively unlimited: the model drew
    # whatever flow the duty needed and no output recorded how much that was.
    # The exergy audit makes it a priced resource, so the flow is now reported
    # and converted into a well count.
    #
    # m_geo is NOT an input. It is DERIVED, exactly as N_wells is:
    #     m_geo = Q_src / (h(T_source) - h(T_source - DT_3A_4A))
    # and the uncertainty is carried by the per-well yield below, which is a
    # field-measured quantity, rather than by a plant-level flow that nobody
    # has ever measured. See DN-22.
    #
    # PROVISIONAL DEFAULT. 12.5 kg/s is the producer in the 4500 m / 206 C
    # repurposed-well case study, which is a hotter and deeper well than ours;
    # the published range for repurposed oil and gas wells runs from about
    # 1 kg/s up to that figure, and co-produced water from depleted fields can
    # exceed it. Treat N_geo as proportional to this number, not as a result.
    m_dot_geo_well: float = 12.5  # per producing well             [kg/s]
    # Lower bound on the reinjection temperature, T_4a = T_source - DT_3A_4A.
    # This is what bounds DT_3A_4A from ABOVE, and the second law wants to run
    # against that bound: raising DT_3A_4A cuts m_geo far faster than it cuts
    # COP. The default is deliberately slack -- it is a placeholder for a
    # formation constraint (silica scaling, injectivity) that the project does
    # not yet have a number for. validate_case warns, it does not raise.
    T_reinject_min_C: float = 25.0
    # Which exergy is charged for the geothermal stream.
    #   'resource'  everything the well lifts, relative to T_sink, with the
    #               reinjected stream booked as a named LOSS attributed to no
    #               component. The right convention for a dedicated well, which
    #               is paid for whether or not its exergy is used.
    #   'stripped'  only what the evaporator removes; reinjection is free.
    # The two give OPPOSITE guidance on DT_3A_4A: 'stripped' falls monotonically
    # and pushes the design toward the largest possible well count, 'resource'
    # has an interior optimum near 20 K. See DN-23.
    exergy_convention: str = "resource"
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
    # The ONLY remaining switch. True carries PCM enthalpy (latent + sensible in
    # both phases), so a fully melted segment superheats and a fully solid one
    # subcools. False pins the PCM at T_m, which is retained not as a model but
    # as the control case for the self-levelling comparison: it isolates what
    # the sensible branches do, at a fixed well count.
    sensible_heat: bool = True

    # Which material lies between the tube and the phase front, and therefore
    # which shell the heat has to cross. See `conduction_shell` and DN-13.
    #   'directional'  melting -> liquid shell (k_l); freezing -> frozen
    #                  shell (k_s). The physical arrangement in both directions.
    #   'melt_side'    the melted thickness always, with k_l whenever any melt
    #                  exists. Correct while melting, INVERTED while freezing.
    #                  Retained only to reproduce results published before the
    #                  correction; it is not a model.
    front_geometry: str = "directional"

    # ---- derived ----
    @property
    def glide_dc_spec(self):
        """Specified discharging glide; falls back to the charging glide."""
        return self.DT_3C_2C if self.DT_3D_2D is None else self.DT_3D_2D

    @property
    def cascade_param(self):
        """Cascade grading parameter; falls back to the charging glide."""
        return self.DT_3C_2C if self.DT_cascade is None else self.DT_cascade

    @property
    def cascade_span(self):
        """Top-to-bottom melting range, one layer short of cascade_param."""
        return self.cascade_param * (self.N_lay - 1) / self.N_lay

    @property
    def T_m(self):
        return self.T_m_C + 273.15

    @property
    def L_well(self):
        return self.L_ft * FT

    @property
    def L_tube(self):
        """Developed length of one fin tube.

        EQUAL TO L_well, where THUMS had twice it. The fin tubes only go
        down; the return is the downcomer and does not exchange. This single
        line is what makes z a depth rather than a path length, and so what
        lets a depth-dependent formation temperature be attached at all.
        """
        return self.L_well

    @property
    def r_id_bore(self):
        """Downcomer bore, sized so its velocity matches the fin tubes'.

        The downcomer carries the WHOLE flow that the n_ft tubes carry
        between them, so equal velocity means equal total area. That is the
        condition under which this geometry has the same pressure drop as the
        hairpin it replaces -- same velocity everywhere, same total path.
        """
        return self.r_i * np.sqrt(float(self.n_ft))

    @property
    def r_id_outer(self):
        """Downcomer outside radius: bore + wall + insulation."""
        wall = self.r_e - self.r_i
        return self.r_id_bore + wall + self.t_ins

    @property
    def D_well(self):
        return self.D_in * IN

    @property
    def V_borehole(self):
        return np.pi * self.D_well ** 2 / 4.0 * self.L_well

    @property
    def V_well(self):
        """PCM volume in one borehole: the hole, less tubes and fins."""
        A_tube = np.pi * self.r_e ** 2
        A_fins = self.num_fins * self.fin_t * self.fin_L
        A_id = np.pi * self.r_id_outer ** 2
        occupied = self.n_ft * (A_tube + A_fins) + A_id
        return self.V_borehole - occupied * self.L_well

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

        The n_ft fin tubes share the PCM-bearing cross-section equally. The
        downcomer occupies area but exchanges no heat with the PCM it
        displaces, so it is SUBTRACTED before the division rather than
        counted as a participant -- a fin tube's cell is bounded by its
        neighbours and by the downcomer's outer wall, not by a share of it.

            pi r_cell^2 = (pi R^2 - pi r_id_outer^2) / n_ft

        It is NOT the borehole wall. With three legs in a 7 in hole r_cell is
        about 40 mm while the wall is at 88.90 mm, so a melt-layer limit
        imposed at the wall would be a factor 2 too permissive and therefore
        silent at every design point.
        """
        R = self.D_well / 2.0
        A_pcm_annulus = np.pi * (R ** 2 - self.r_id_outer ** 2)
        return np.sqrt(A_pcm_annulus / (np.pi * self.n_ft))

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
        """The eleven positional geometry arguments the geometry helpers expect."""
        return [self.L_well, self.L_tube, self.D_well, self.r_i, self.r_e,
                2 * self.r_i, 2 * self.r_e, self.fin_t, self.fin_L,
                self.num_fins, self.n_ft]

    def stefan_number(self, dT):
        """Ste = cp dT / h_m. The quasi-steady melt layer degrades as this grows."""
        return self.cp_l * dT / self.h_m

    def with_(self, **kw):
        return replace(self, **kw)


CASE = Case()          # the design point



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
    A_avail = case.V_well / (case.n_ft * case.L_tube)   # PCM area per tube [m2]
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


def segment_profile(case, T_inlet, T_m_seg, m_dot, k_wall, E, n_segments=None):
    """Instantaneous fluid and conductance profile for a frozen PCM state.

    The inner loop of `march_h` without the time advance: it evaluates the
    resistance network, the segment conductance and the fluid temperature
    profile for the state `E` as it stands, and returns the instantaneous
    q' = K (T_f - T_pcm).

    It exists so that the last recorded frame of a march can carry a fluid
    profile belonging to the *end* state rather than to the state one step
    earlier. Everything it computes is already computed inside `march_h`;
    this is one extra sweep per recorded march, not one per time level.
    """
    n_segments = n_segments or len(E)
    dz = case.L_tube / n_segments
    A, T_pcm = pcm_state(E, T_m_seg, case)
    delta = delta_from_area(A, case.r_e, case.num_fins, case.fin_t, case.fin_L)
    A_avail, E_lat = pcm_capacities(case)[:2]
    d_melt, k_melt = conduction_shell(case, A, A_avail, True, E, E_lat)
    d_froz, k_froz = conduction_shell(case, A, A_avail, False, E, E_lat)

    T_prof = np.empty(n_segments + 1)
    T_prof[0] = T_inlet
    NTU_prof = np.empty(n_segments)
    U_prof = np.empty(n_segments)
    q_prime = np.empty(n_segments)

    T0 = T_inlet
    for i in range(n_segments):
        melting = T0 > T_pcm[i]
        d_path = float(d_melt[i] if melting else d_froz[i])
        k_m = float(k_melt[i] if melting else k_froz[i])
        h_i, cp_d, _, _ = h_internal(case.fluid2, T0, case.P, case.r_i, m_dot)
        U_i = compute_U_i(h_i, case.r_i, case.r_e, k_wall, case.Rf_i, k_m,
                          d_path, case.L_tube, case.fin_t,
                          case.fin_L, case.num_fins)
        NTU = float(np.clip((2 * np.pi * case.r_i * U_i * dz) / (m_dot * cp_d),
                            -50.0, 50.0))
        K = m_dot * cp_d * (1.0 - np.exp(-NTU)) / dz
        q_prime[i] = K * (T0 - T_pcm[i])
        T0 = T0 - q_prime[i] * dz / (m_dot * cp_d)
        T_prof[i + 1] = T0
        NTU_prof[i] = NTU
        U_prof[i] = U_i

    return dict(T_fluid=T_prof, NTU=NTU_prof, U_i=U_prof, q_prime=q_prime,
                A_melt=A, T_pcm=T_pcm, delta=delta, E=np.array(E, float))


def bulk_shape_factor(case):
    """Dimensionless resistance of a SINGLE-PHASE cell, mean to tube surface.

    A cell that is entirely one phase has no front. Heat entering at the tube is
    stored throughout its volume, so the temperature varies with radius and the
    lumped state T_pcm is the VOLUME MEAN. The resistance that belongs with that
    mean is not the surface-to-surface annulus value.

    Take the cell as an annulus r_e <= r <= r_o, adiabatic at r_o (its neighbour
    is identical, so the boundary is a symmetry plane), storing sensible heat at
    a uniform volumetric rate s = rho c dT/dt. Quasi-steady, the heat crossing
    radius r is what the material beyond it stores,

        -k 2 pi r dT/dr = s pi (r_o^2 - r^2),

    which integrates to

        T(r_e) - T(r) = (s/2k) [ r_o^2 ln(r/r_e) - (r^2 - r_e^2)/2 ].

    Averaging over the volume and dividing by the total Q' = s pi (r_o^2-r_e^2)
    gives R'_bulk = S/(2 pi k) with, writing beta = r_o/r_e,

        S = [ beta^4 ln(beta) - beta^4/2 + beta^2/2 - (beta^2-1)^2/4 ]
            / (beta^2 - 1)^2 .

    S depends on GEOMETRY ALONE -- k cancels -- so one shape factor serves both
    phases and only the conductivity changes between them. For the default cell,
    beta = 2.1086 and S = 0.34672, against ln(beta) = 0.74604 for the
    surface-to-surface annulus: the mean-to-surface resistance is smaller by a
    factor 2.152, the cylindrical analogue of the familiar 1/3 for a slab with
    uniform generation.

    Verified against numerical quadrature of the same integral to six decimals.
    """
    beta = case.r_cell / case.r_e
    b2 = beta * beta
    return ((b2 * b2 * np.log(beta) - b2 * b2 / 2.0 + b2 / 2.0
             - (b2 - 1.0) ** 2 / 4.0) / (b2 - 1.0) ** 2)


def bulk_equivalent_delta(case):
    """The shell thickness that reproduces the bulk resistance.

    The network in `compute_U_i` takes a thickness, not a resistance, and wraps
    it in the fin efficiency and the finned-area ratio. Feeding it an equivalent
    thickness therefore keeps the single-phase branches inside exactly the same
    machinery as the two-phase ones, rather than bolting a second path onto it:

        ln(1 + delta_eq/r_e) = S    =>    delta_eq = r_e (e^S - 1).

    Like S it is a property of the geometry alone -- the same delta_eq serves
    subcooled solid and superheated liquid, and only k changes. For the default
    cell it is 8.74 mm, against 0 and 23.37 mm, which are the two values the
    model used before this was worked out.
    """
    return case.r_e * (np.exp(bulk_shape_factor(case)) - 1.0)


def conduction_shell(case, A_melt, A_avail, melting, E=None, E_lat=None):
    """Thickness and conductivity of the shell between the tube and the front.

    The PCM-side resistance in `compute_U_i` is that of an annulus growing
    OUTWARD FROM THE TUBE WALL,

        R' = ln(1 + delta/r_e) / (2 pi k),

    so the question this function answers is: which material is that annulus
    made of, and how thick is it?

    It depends on the direction of the phase change, because the front always
    grows away from the tube -- the tube is the driven boundary in both
    half-cycles.

      MELTING (fluid hotter than the PCM). Melting begins at the tube wall and
      the front moves outward, so the shell is LIQUID and its thickness is that
      of the melted area. Heat crosses the melt to reach the remaining solid.

      FREEZING (fluid colder). Solidification also begins at the tube wall, so
      the shell is the FROZEN material and its thickness is that of the solid
      area, A_avail - A_melt. The liquid is displaced outward, beyond the front,
      and is no longer in the conduction path at all.

    All four limits come out right without special-casing:

      fully solid, melting     A_melt = 0            -> delta = 0, no resistance
      fully melted, melting    A_melt = A_avail      -> the full liquid annulus
      fully melted, freezing   A_avail - A_melt = 0  -> delta = 0, no resistance
      fully solid, freezing    A_avail - A_melt = A  -> the full frozen annulus

    The last of those is the one that matters. Until v0.7 the melted thickness
    was used in BOTH directions, which put the conduction path on the wrong
    side of the front during discharge: a segment that had frozen solid was
    given ZERO PCM-side resistance, exactly where the physical resistance is
    largest. The consequence was that the modelled exchanger IMPROVED as it
    froze, when it should degrade. See DN-13 for what that was worth.

    LIMIT OF ANY LUMPED FRONT. One thickness can describe one front. After the
    first half-cycle from a fully solid store the geometry is generally
    three-region -- at cyclic steady state the charge begins with residual
    liquid left in the OUTER part of the cell, so melting produces liquid at the
    tube, solid in the middle and liquid outside. Neither this function nor the
    retired one can represent that; this one is the better approximation, not a
    correct treatment. It is hypothesis H5, stated honestly.
    """
    if case.front_geometry not in ("directional", "melt_side"):
        raise ValueError(f"front_geometry must be 'directional' or "
                         f"'melt_side', not {case.front_geometry!r}")

    A_melt = np.asarray(A_melt, float)
    if case.front_geometry == "melt_side":          # retired, see above
        return (delta_from_area(A_melt, case.r_e, case.num_fins, case.fin_t,
                                case.fin_L),
                np.full(A_melt.shape, case.k_l))

    # ---- two-phase: the shell between the tube and the front --------------
    A_path = A_melt if melting else (A_avail - A_melt)
    delta = delta_from_area(A_path, case.r_e, case.num_fins, case.fin_t,
                            case.fin_L)
    k_path = np.full(A_melt.shape, case.k_l if melting else case.k_s)

    # ---- single phase: there is no front, so there is no side ------------
    # `E` is optional only so that callers which do not have it (none, now)
    # keep the two-phase behaviour; when it is given, the single-phase cells
    # get the bulk resistance of `bulk_shape_factor`, the SAME in both
    # directions, because a one-phase cell has no front and cannot care which
    # way the heat is going.
    if E is not None:
        E = np.asarray(E, float)
        d_bulk = bulk_equivalent_delta(case)
        solid = E < 0.0                              # subcooled solid
        liquid = E > E_lat                           # superheated liquid
        delta = np.where(solid | liquid, d_bulk, delta)
        k_path = np.where(solid, case.k_s,
                          np.where(liquid, case.k_l, k_path))
    return delta, k_path


def mixed_mean_outlet(res):
    """Flow-mixed mean outlet temperature of a march [K].

    The temperature at which the stream returns to the plant, averaged over the
    half-cycle. At constant mass flow and c_p the mixing average reduces to the
    time average, which on a logarithmic grid still has to be weighted by the
    step -- the plain mean of `T_out` is dominated by the early levels and is
    wrong by several kelvin.

    This is the quantity the ORC actually sees, and since v0.6 it is a RESULT
    rather than a prescribed state point (`cycle_state_points`).
    """
    t = np.asarray(res["t"], float)
    dt = np.diff(np.concatenate(([0.0], t)))
    return float(np.sum(np.asarray(res["T_out"], float) * dt) / np.sum(dt))


def unmirror_march(res):
    """Re-index a reversed-flow march result into charge (depth) coordinates.

    Discharge marches from the opposite end of the tube, so its segment index
    runs backwards relative to the charge. `run_cycle` and `simulate_css` hand
    it a mirrored cascade and a mirrored initial state and un-mirror the result
    when they close the loop; anything that *plots* a discharge result has to
    do the same or it draws the well upside down.

    Returns a shallow copy with every per-segment array reversed, including the
    recorded history. `z`, `t` and `t_hist` are coordinates, not fields, and
    are left alone. Scalars are left alone. The cascade in depth coordinates is
    the same for both half-cycles, which is the point: after this call a
    discharge result is plotted against `T_m_seg`, not against `T_m_seg[::-1]`.
    """
    n = len(res["E"])
    out = dict(res)
    for k in ("E", "A_melt", "T_pcm", "delta", "eps_local", "merge_proximity"):
        if k in out and out[k] is not None:
            out[k] = np.asarray(out[k])[::-1].copy()
    h = res.get("history")
    if h:
        hh = {}
        for k, v in h.items():
            v = np.asarray(v)
            # T_fluid has n+1 nodes (inlet plus one per segment); the rest have n
            hh[k] = v[:, ::-1].copy() if v.ndim == 2 and v.shape[1] in (n, n + 1) else v
        out["history"] = hh
    return out


# ==========================================================================
# the formation   (BOREAS v1.0 -- THUMS hypothesis H6 said this was adiabatic)
# ==========================================================================
# H6 made eta_storage identically 1. Not approximately 1: EXACTLY 1, by
# construction, because nothing could leave. That is an indefensible headline
# for a borehole storage paper and it is the first thing a reviewer attacks.
#
# The resistances are wildly unequal and it is worth knowing which one is the
# model. At the design point, per metre of borehole:
#
#     casing steel        0.0003 m K/W     0.1 %
#     cement              0.018            5 %
#     PCM to casing       ~0.12            25 %
#     GROUND              0.32             ~70 %
#
# Everything near the well is a rounding error. The ground resistance is the
# physics, and unlike the others it GROWS WITH TIME as the thermal front
# spreads -- which is why `t_operation_yr` is a Case parameter and why a
# single quoted loss figure is meaningless without the age that goes with it.
#
# THE SIGN IS NOT ALWAYS NEGATIVE. With a 45 K/km gradient the rock at 8000 ft
# is near 128 C, and at the Wilmington gradient of 56.5 K/km it is near 156 C
# -- hotter than the coldest cascade layers. Those segments GAIN heat from the
# formation. A model that assumes loss and takes an absolute value somewhere
# would hide that, so nothing here does.


def ground_resistance(case):
    """Quasi-steady line-source resistance of the formation, per metre.

    Eskilson's form: at time t the disturbance has spread to roughly
    2.5*sqrt(alpha t), and the resistance is the log of that over the borehole
    radius. It is quasi-steady in the sense that the DAILY cycle is far faster
    than the formation's response, so within one charge-discharge the ground
    looks like a fixed resistance to a fixed far-field temperature. The cycle
    time is 10 h and the formation time constant is years; the separation is
    four orders of magnitude, so this is a good approximation rather than a
    convenient one.

    SINGLE WELL, NO FIELD SHIELDING. This is deliberate and it is the
    conservative choice: an isolated well is the worst case, because in a
    field the neighbours hold the ground between them warm and only the
    perimeter leaks. The shielding factor is 29x to 42x, essentially
    independent of the driving temperature and only weakly dependent on rock
    conductivity, so a field result can be obtained by dividing. Nothing here
    should be read as the loss a real cluster would see.
    """
    if case.k_rock <= 0.0:
        return np.inf                      # no formation: adiabatic wall
    alpha = case.k_rock / case.rho_cp_rock
    t = case.t_operation_yr * 365.25 * 24.0 * 3600.0
    arg = 2.5 * np.sqrt(alpha * t) / case.r_bore
    R = np.log(arg) / (2.0 * np.pi * case.k_rock) if arg > 0.0 else -1.0
    if not np.isfinite(R) or R <= 0.0:
        # THE LINE SOURCE HAS BROKEN DOWN, and the branch matters.
        # An earlier version returned 0.0 here, meaning the front "has not
        # cleared the hole yet". Zero RESISTANCE is infinite CONDUCTANCE: it
        # says the rock is a perfect heat sink, which is the opposite of what
        # was meant and the most damaging possible failure of this function.
        # It was caught the first time the adiabatic limit was watched -- a
        # test that drives k_rock to zero and expects eta_storage = 1 got
        # eta_storage = 0.20 and a LARGER loss than the physical case.
        #
        # Infinity is the right fallback: outside its validity the line
        # source is not evidence of a loss, and `validate_case` warns rather
        # than letting a silent number through.
        return np.inf
    return R


def field_radius(case, depth):
    """Radius of the well cluster at a given depth [m].

    Vertical above the kickoff point, then fanning at a constant effective
    half-angle. A real directional programme builds angle over a section and
    then holds it, so the fan is not a straight cone -- but the quantity this
    feeds is a logarithm, and the difference between a cone and a build-and-
    hold is well inside the uncertainty in the spacing data itself.
    """
    if case.n_wells_field is None:
        return None
    R0 = case.B_surface_m * np.sqrt(case.n_wells_field / np.pi)
    spread = np.maximum(0.0, np.asarray(depth, float) - case.z_kickoff_m)
    return R0 + spread * np.tan(np.radians(case.fan_angle_deg))


def well_spacing(case, depth):
    """Centre-to-centre spacing at depth, scaled with the fan radius [m]."""
    Rf = field_radius(case, depth)
    if Rf is None:
        return None
    R0 = case.B_surface_m * np.sqrt(case.n_wells_field / np.pi)
    return case.B_surface_m * Rf / R0


def formation_resistance(case, depth=None):
    """Total series resistance, PCM bulk to undisturbed rock [m K / W].

    Per metre of BOREHOLE, not per metre of tube. The PCM-to-casing term is
    the crude one: n_ft parallel paths from the tube walls out to the casing,

        R_pcm = ln(R_casing_i / r_e) / (2 pi n_ft k_pcm)

    which ignores that the paths crowd together near the wall. It is an
    approximation and it is flagged as one -- but it is a quarter of a total
    that the ground dominates, so refining it would move the answer less than
    the uncertainty in k_rock does.
    """
    R_in = case.D_well / 2.0
    k_pcm = 0.5 * (case.k_s + case.k_l)
    R_pcm = np.log(R_in / case.r_e) / (2.0 * np.pi * case.n_ft * k_pcm)
    R_steel = np.log(case.r_casing_o / R_in) / (2.0 * np.pi * case.k_wall)
    R_cem = np.log(case.r_bore / case.r_casing_o) / (2.0 * np.pi * case.k_cement)
    R_g = ground_resistance(case)
    if not np.isfinite(R_g):
        return np.inf                      # adiabatic; K_loss becomes 0
    R_near = R_pcm + R_steel + R_cem

    # ---- FIELD SHIELDING -------------------------------------------------
    # A well inside a cluster cannot lose what an isolated one would: its
    # neighbours hold the rock between them warm, and once the cells have
    # merged the only heat leaving is across the cluster PERIMETER, shared
    # among N wells. Writing that as a resistance,
    #
    #     R_perimeter(z) = N sqrt(pi alpha t) / (2 pi k R_field(z))
    #
    # -- the semi-infinite-slab flux at the boundary, over the perimeter the
    # field presents, divided by the wells that share it. It carries no
    # driving temperature, so it is a true resistance and drops straight into
    # the series chain.
    #
    # The two regimes are combined by taking the LARGER resistance. That is
    # not a smooth blend and does not pretend to be: a well cannot leak more
    # than if it stood alone, and a field cannot leak more than its perimeter
    # allows, so whichever constraint binds is the one that governs. The
    # crossover is sharp in reality too -- it is the moment the cells merge.
    if depth is None or case.n_wells_field is None:
        return R_near + R_g
    alpha = case.k_rock / case.rho_cp_rock
    t = case.t_operation_yr * 365.25 * 24.0 * 3600.0
    Rf = field_radius(case, depth)
    R_per = (case.n_wells_field * np.sqrt(np.pi * alpha * t)
             / (2.0 * np.pi * case.k_rock * Rf))
    return R_near + np.maximum(R_g, R_per)


def formation_profile(case, n_segments=None, discharge=False):
    """Undisturbed rock temperature [K] at each segment, and its depth [m].

    THE COORDINATE. In BOREAS z is DEPTH, which it could not be in THUMS: the
    hairpin's two legs occupied every depth twice, so one rock temperature
    would have had two PCM segments claiming it. Here the fin tubes only go
    down and the downcomer does not exchange, so the map is one-to-one.

    On CHARGE the water enters the fin tubes at the BOTTOM -- it arrives there
    down the downcomer -- so the solver's z = 0 is the bottom of the well and
    depth decreases with z. On DISCHARGE the flow reverses and the march is
    run mirrored, so the profile is mirrored with it, exactly as T_m_seg is.
    Getting this backwards would put the hot cascade layers against the cold
    shallow rock, which is the THUMS orientation and the thing BOREAS exists
    to avoid -- so it is asserted in `validate_case`, not left to a comment.
    """
    n = n_segments or case.n_segments
    xi = (np.arange(n) + 0.5) / n                 # cell centres, 0 at inlet
    depth = case.L_well * (1.0 - xi)              # charge frame: z=0 is bottom
    T_rock = case.T_surface_C + 273.15 + case.grad_K_per_km / 1000.0 * depth
    if discharge:
        return T_rock[::-1].copy(), depth[::-1].copy()
    return T_rock, depth


def march_h(case, T_inlet, T_m_seg, m_dot, k_wall, times,
            n_segments=None, E0=None, record=False, T_rock=None,
            R_form=None):
    """Segment march on the ENTHALPY state. Replaces `march` from v0.4.

    Identical to `march` in every respect except what is integrated: the state
    is E' [J per m of tube] rather than A_melt [m2], so no heat is ever clipped
    and the PCM temperature is free to leave T_m.

    Recording convention. Every history frame is stamped, in `t_hist`, with the
    instant at which its state arrays (`E`, `A_melt`, `T_pcm`, `delta`) are
    exact. The profile arrays (`T_fluid`, `NTU`, `U_i`, `q_prime`) in the same
    frame are the ones evaluated from that state, i.e. the profile that drives
    the step *beginning* there. There are `len(t)+1` frames: one at t = 0 and
    one at the end of each accepted step, so the final state is in the record.

    Earlier builds recorded the pre-step state but stamped it with the
    post-step time, which on a logarithmic grid put the last frame a quarter of
    a charge behind its label, and appended `E` after the update while the
    other state arrays were from before it. Neither affected any marched or
    reported quantity -- `Q_cum`, `closure`, `eps_local` and the returned state
    all come from the march itself -- but every recorded plot was lagged.
    """
    n_segments = n_segments or case.n_segments
    dz = case.L_tube / n_segments
    r_i, r_e = case.r_i, case.r_e
    A_avail, E_lat, C_s, C_l = pcm_capacities(case)

    # ---- the formation sink ---------------------------------------------
    # K_loss is a conductance per unit TUBE length, so the per-borehole
    # conductance is divided among the n_ft tubes that share each depth.
    # Passing T_rock=None recovers the THUMS adiabatic wall exactly, which is
    # how the v0.16 reproduction test is run.
    if T_rock is None:
        K_loss = np.zeros(n_segments)
        T_rock_arr = np.zeros(n_segments)
    else:
        R_tot = (formation_resistance(case) if R_form is None
                 else np.asarray(R_form, float))
        R_arr = np.full(n_segments, R_tot) if np.ndim(R_tot) == 0 else R_tot
        K_loss = np.where(np.isfinite(R_arr), 1.0 / (R_arr * case.n_ft), 0.0)
        K_loss = np.nan_to_num(K_loss, posinf=0.0, neginf=0.0)
        T_rock_arr = np.asarray(T_rock, float)
    Q_loss_cum = 0.0
    # ---- A CLAUSIUS GATE ON THE LOSS TERM --------------------------------
    # The closure CANNOT catch a flipped loss sign, and this was watched:
    # negating q_loss leaves the residual at 5e-15, because q_prime and
    # Q_loss_cum both carry the error and it cancels algebraically. The
    # closure is self-consistent, which is the DN-26 trap -- a check reading
    # only its own inputs.
    #
    # This one reads the STATE instead. Heat must flow from hot to cold, so
    # q_loss and (T_pcm - T_rock) must share a sign at every segment and
    # every step. The product is accumulated and gated on being non-negative.
    # With the correct convention it is K (dT)^2 and cannot be negative; with
    # the sign flipped it is -K (dT)^2 and cannot be positive. That is the
    # independent input the closure lacks.
    clausius = 0.0

    E = np.zeros(n_segments) if E0 is None else np.array(E0, float)
    E_start = E.copy()
    t_prev, Q_cum = 0.0, 0.0
    ts, Qs, T_outs = [], [], []
    hist = {k: [] for k in ("T_fluid", "T_pcm", "A_melt", "delta", "NTU",
                            "U_i", "q_prime", "E")} if record else None
    t_hist = []

    for t in np.asarray(times, float):
        dt = t - t_prev
        if dt <= 0:
            t_prev = t
            continue

        A, T_pcm = pcm_state(E, T_m_seg, case)
        # `delta` stays the MELT thickness: it is what the H3 merge diagnostic
        # is about. The conduction path is a separate question -- which side of
        # the front the heat has to cross -- and depends on the direction.
        delta = delta_from_area(A, r_e, case.num_fins, case.fin_t, case.fin_L)
        d_melt, k_melt = conduction_shell(case, A, A_avail, True, E, E_lat)
        d_froz, k_froz = conduction_shell(case, A, A_avail, False, E, E_lat)
        T0 = T_inlet
        q_prime = np.zeros(n_segments)
        if record:
            E_before = E.copy()          # the state this frame is stamped at
            T_prof = np.empty(n_segments + 1); T_prof[0] = T_inlet
            NTU_prof = np.empty(n_segments); U_prof = np.empty(n_segments)

        for i in range(n_segments):
            # Which shell the heat crosses depends on which way the front is
            # moving; the sign of the driving difference settles it before K is
            # known, and it has the same sign as q'.
            melting = T0 > T_pcm[i]
            d_path = float(d_melt[i] if melting else d_froz[i])
            k_m = float(k_melt[i] if melting else k_froz[i])
            h_i, cp_d, _, _ = h_internal(case.fluid2, T0, case.P, r_i, m_dot)
            U_i = compute_U_i(h_i, r_i, r_e, k_wall, case.Rf_i, k_m,
                              d_path, case.L_tube, case.fin_t,
                              case.fin_L, case.num_fins)
            NTU = float(np.clip((2 * np.pi * r_i * U_i * dz) / (m_dot * cp_d),
                                -50.0, 50.0))
            # conductance of the segment to the fluid, per unit tube length:
            # q' = K (T0 - T_pcm), from the NTU effectiveness
            K = m_dot * cp_d * (1.0 - np.exp(-NTU)) / dz

            # TWO linear sinks on one node, combined EXACTLY rather than
            # split operator. For any T_pcm,
            #     K (T0 - T) + K_l (T_rock - T) = (K + K_l) (T_eq - T)
            # with T_eq the conductance-weighted mean, so `advance_segment`
            # is called unchanged and remains implicit in T_pcm. Adding the
            # loss explicitly after the step would have been a half-order
            # scheme on a term that can reverse sign down the well.
            K_l = K_loss[i]
            K_eq = K + K_l
            T_eq = (K * T0 + K_l * T_rock_arr[i]) / K_eq if K_eq > 0 else T0
            E[i], q_tot = advance_segment(E[i], T_m_seg[i], T_eq, K_eq, dt,
                                          E_lat, C_s, C_l,
                                          case.sensible_heat)
            # Recover the split. q_tot is what the PCM NET received, so the
            # effective PCM temperature over the step follows from it, and
            # the two streams separate without any further assumption.
            if K_l > 0.0:
                T_eff = T_eq - q_tot / K_eq
                q_loss = K_l * (T_eff - T_rock_arr[i])   # +ve = leaving
                Q_loss_cum += q_loss * dz * dt
                clausius += q_loss * (T_eff - T_rock_arr[i]) * dz * dt
            else:
                q_loss = 0.0
            q_prime[i] = q_tot + q_loss                  # what the TUBE gave
            # the fluid gives up exactly what the PCM took FROM IT -- the
            # formation term is not the fluid's to pay
            T0 = T0 - q_prime[i] * dz / (m_dot * cp_d)
            if record:
                T_prof[i + 1] = T0; NTU_prof[i] = NTU; U_prof[i] = U_i

        if record:
            # stamped at t_prev, where E_before / A / T_pcm / delta are exact
            t_hist.append(t_prev)
            hist["T_fluid"].append(T_prof.copy()); hist["T_pcm"].append(T_pcm.copy())
            hist["A_melt"].append(A.copy());       hist["delta"].append(delta.copy())
            hist["NTU"].append(NTU_prof.copy());   hist["U_i"].append(U_prof.copy())
            hist["q_prime"].append(q_prime.copy()); hist["E"].append(E_before)

        Q_cum += float(np.sum(q_prime) * dz) * dt
        ts.append(t); Qs.append(float(np.sum(q_prime) * dz)); T_outs.append(T0)
        t_prev = t

    if record:
        # the end state, which the loop above never gets to record
        end = segment_profile(case, T_inlet, T_m_seg, m_dot, k_wall, E,
                              n_segments=n_segments)
        t_hist.append(t_prev)
        for k in hist:
            hist[k].append(end[k])

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
        "Q_loss_J": Q_loss_cum,
        # >= 0 always, by the second law. Negative means the loss term is
        # pushing heat uphill, which is the sign error the closure misses.
        "loss_clausius": clausius,
        # THE CLOSURE NOW HAS THREE TERMS. In THUMS it was
        # |Q_cum - dE| / |Q_cum|: everything the tube gave went into storage,
        # because nothing could leave. Here what the tube gives either stays
        # or goes to the rock, so the identity is Q_tube = dE + Q_loss and the
        # residual is what is left over. This is still a bookkeeping identity
        # rather than a gate -- it cannot fail without an arithmetic bug --
        # but it would have caught the sign error of dropping q_loss from
        # q_prime, which is the mistake this split invites.
        "closure": (abs(Q_cum - dE - Q_loss_cum) / abs(Q_cum)
                    if Q_cum else np.nan),
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
        "t_hist": np.array(t_hist) if record else None,
    }


def layer_mean_T_m(case):
    """<T_m>, the volume-weighted mean melting temperature over the cascade.

    The layers are equal in developed length and each owns an equal share of
    the PCM cross-section, so the volume weighting reduces to a plain mean over
    the N_lay melting temperatures:

        T_m,l = T_m_top - l * dT_glide / N_lay,        l = 0 .. N_lay-1

        <T_m> = T_m_top - (N_lay - 1)/(2 N_lay) * dT_glide

    For N_lay = 1 this returns T_m_top, which is what makes the old capacity
    formula the NO-CASCADE limit of the corrected one below.
    """
    return case.T_m - 0.5 * (case.N_lay - 1) / case.N_lay * case.DT_3C_2C


def well_capacity_kJ(case, T_charge_in, T_discharge_in):
    """Energy ONE well can hold over a full cycle -- the capacity criterion.

    CORRECTED 2026-09 (v0.4b). The previous form was

        E_well = rho V_well [ h_m + cp_l (T_charge_in - T_m) ]

    with `case.T_m` the melting temperature of the TOP layer, so the sensible
    span was DT_4C_M = 10 K for the whole store. But DT_4C_M is the approach
    that fixes the charging inlet relative to the FIRST layer only. Applying it
    to all N_lay layers is the no-cascade limit of this function: set
    N_lay = 1 and the two expressions agree identically.

    It was also not a bound in either direction. At the design point the
    model's own end-of-charge state held 4.822 MWh per well against a claimed
    capacity of 4.744 -- mean superheat reached 13.81 K, not 10 K, because
    segments early in a layer sit under fluid far hotter than their own T_m and
    because a nearly saturated store barely cools the fluid at all.

    The corrected form is a true CEILING on the cycle swing, resolved over the
    cascade:

        E_well = rho V_well [ h_m
                            + cp_l <T_3c  - T_m,l>      liquid superheat
                            + cp_s <T_m,l - T_3d> ]     solid subcooling

    Both spans are physical limits rather than design choices: the PCM cannot
    pass the temperature of the fluid entering the field on either half-cycle.
    The solid term belongs here because the store begins each charge from
    wherever the discharge left it, so the capacity that matters is the full
    swing. That term is the one the previous docstring promised through a
    `cycle=True` argument that was never implemented, alongside a
    `T_discharge_in` parameter that was accepted and never used.

    Consequence at the design point: N_capacity 11.573 -> 9.574, so the RATE
    criterion binds instead and N_wells falls 11.573 -> 11.255. That is the
    point of the change. A bound that binds should be suspect, because it
    asserts the design sits at its thermodynamic ceiling. It does not: it sits
    at 84 % of it, whereas the old formula put it at 102 %.
    """
    m_well = case.rho_latent * case.V_well
    E = case.h_m
    if case.sensible_heat:
        T_m_bar = layer_mean_T_m(case)
        # liquid superheat: ceiling is the charging inlet temperature
        E += case.cp_l * max(T_charge_in - T_m_bar, 0.0)
        # solid subcooling: ceiling is the discharging inlet temperature
        if T_discharge_in is not None:
            E += case.cp_s * max(T_m_bar - T_discharge_in, 0.0)
    return m_well * E / 1000.0                      # kJ
