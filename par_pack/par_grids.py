"""The parametric grids, on the Grissom field (R4).

Everything that was swept against the old packing / fan-angle field model
is gone. The field is now the measured one, and the only thing a grid
point can change about it is the storage length.

WHAT THE FIELD DEPENDS ON, and what it does not. R_g,eff(MD) comes out of
the superposition and is a function of the WELL POSITIONS, the
TRAJECTORIES, the storage length L, the bore radius, the rock
conductivity and diffusivity, and the field age. It is NOT a function of
the PCM (T_m, h_m, N_lay, the cascade span), of the tube count, or of the
geothermal gradient -- the gradient sets the rock TEMPERATURE, not the
resistance between the well and it.

So one superposition solve per L serves every point on every map that
shares that length, and the gradient axis rides for free. That is what
makes a 6 x 6 map at R4 affordable at all.
"""
FT = 0.3048

# ---- the L axis, in metres of MEASURED depth ---------------------------
# Capped at the 5000 ft casing shoe, per the paper task's review. Note the
# inconsistency this exposes, flagged in the handoff: the ladder's sweep
# marks a PROVISIONAL shoe at 1250 m MD, while the maps are asked to stop
# at 5000 ft = 1524 m. They are two statements of the same unknown. The
# grid includes 1250 m so the reference case is on it, and ends at 1524 m.
LW_M = (250.0, 500.0, 750.0, 1000.0, 1250.0, 1524.0)
L_REF = 1250.0
L_SHOE_FT = 5000.0

# ---- PCM and cascade ---------------------------------------------------
# T_m capped at 160 C. The hard ceiling is the secondary loop: water
# saturates at 179.9 C at P = 10 bar, which caps the usable T_m near
# 165 C, and the top layer has to sit below the HTHP condenser anyway.
TM = (110.0, 120.0, 130.0, 140.0, 150.0, 160.0)
DG = (25.0, 40.0, 55.0, 70.0, 85.0, 100.0)      # DT_3C_2C, the grading

# ---- latent heat, the new axis ----------------------------------------
# 200 to 500 kJ/kg spans the inorganic salt hydrates and nitrates that are
# candidates in this temperature band; 380 is the reference.
HM = (200e3, 260e3, 320e3, 380e3, 440e3, 500e3)

NL = (3, 5, 7, 9, 12, 15)          # N_lay
NFT = (1, 2, 3, 4, 5, 6)           # fin tubes per well
GR = (25.0, 35.0, 45.0, 56.5, 70.0, 85.0)       # grad_K_per_km


def nseg(nl):
    """n_segments as the nearest multiple of N_lay to 99.

    A layer boundary inside a segment makes validate_case warn, and
    carrying that inconsistency along the whole N_lay axis would be an
    artifact in the one direction the map is measuring.
    """
    return int(nl * max(1, round(99.0 / nl)))


# N_lay is CONFOUNDED with the cascade span: span = DT_cascade*(N-1)/N, so
# sweeping N at fixed DT_3C_2C moves the span too, and "fewer layers is
# better" may be nothing but the span. Pinning DT_cascade = span0*N/(N-1)
# holds the span at the reference 48.889 K and isolates the layer count.
SPAN0 = 55.0 * 8.0 / 9.0


def grids():
    """(name, list of point-dicts). Every point carries its own L_md."""
    G = {}
    G['cascade'] = [dict(L_md=L_REF, T_m_C=t, DT_3C_2C=d, auto_P=True)
                    for t in TM for d in DG]
    G['L_Nlay'] = [dict(L_md=L, N_lay=n, n_segments=nseg(n))
                   for n in NL for L in LW_M]
    G['L_nft'] = [dict(L_md=L, n_ft=f) for f in NFT for L in LW_M]
    G['hm_Tm'] = [dict(L_md=L_REF, h_m=h, T_m_C=t, auto_P=True)
                  for t in TM for h in HM]
    G['L_grad'] = [dict(L_md=L, grad_K_per_km=g) for g in GR for L in LW_M]
    G['Nlay_free'] = [dict(L_md=L_REF, N_lay=n, n_segments=nseg(n))
                      for n in NL]
    G['Nlay_pinned'] = [dict(L_md=L_REF, N_lay=n, n_segments=nseg(n),
                             DT_cascade=SPAN0 * n / (n - 1.0)) for n in NL]
    G['hm_only'] = [dict(L_md=L_REF, h_m=h) for h in HM]
    return G
