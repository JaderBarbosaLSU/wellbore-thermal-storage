ARRAY = dict(n_segments=99, n_wells_field=1000., fan_angle_deg=0.,
             packing='rows', n_islands=4)
# field config WITHOUT the fan angle, for maps that sweep it
ARRAY_NOFAN = dict(n_segments=99, packing='rows', n_islands=4)

TM = (120, 135, 150, 165, 180, 195)
DG = (25, 40, 55, 70, 85, 100)
CASCADE = [{'T_m_C': float(t), 'DT_3C_2C': float(d), 'auto_P': True}
           for t in TM for d in DG]

LW = (1000, 2000, 3000, 4000, 5000, 6000)       # L_ft
NL = (3, 5, 7, 9, 12, 15)                        # N_lay
NFT = (1, 2, 3, 4, 5, 6)                         # n_ft
NW = (10, 30, 100, 300, 1000, 3000)              # n_wells_field
FA = (0, 15, 30, 45, 60, 75)                     # fan_angle_deg
GR = (25, 35, 45, 56.5, 70, 85)                  # grad_K_per_km


def nseg(nl):
    """Nearest multiple of N_lay to 99.

    n_segments must be a multiple of N_lay or layer_map puts a boundary
    inside a segment and validate_case warns. Choosing it per point removes
    the artefact from the whole N_lay axis for free, rather than carrying a
    1-in-99 inconsistency through every cell and hoping it is small.
    """
    return int(nl * max(1, round(99.0 / nl)))


# C1: how capacity splits between well length and layer count
G_LN = [{'L_ft': float(L), 'N_lay': int(n), 'n_segments': nseg(n)}
        for n in NL for L in LW]
# C2: pumping -- more, smaller tubes means more PCM but more pressure drop
G_FT = [{'n_ft': int(f), 'L_ft': float(L)} for L in LW for f in NFT]
# D1: shielding -- well count against how fast the wells diverge
G_WF = [{'n_wells_field': float(w), 'fan_angle_deg': float(a)}
        for a in FA for w in NW]
# D2: hotter rock against more boundary
G_GL = [{'grad_K_per_km': float(g), 'L_ft': float(L)} for L in LW for g in GR]

# n_segments varies per point in G_LN, so it must not also be in `fixed`
ARRAY_NOSEG = {k: v for k, v in ARRAY.items() if k != 'n_segments'}

# N_lay is CONFOUNDED with the cascade span: span = DT_cascade*(N-1)/N, so
# sweeping N at fixed DT_3C_2C changes the span too and the apparent
# "fewer layers is better" may be nothing but the span effect. Pinning
# DT_cascade = span0*N/(N-1) holds the span at the reference 48.889 K and
# isolates the layer count itself.
SPAN0 = 55.0 * 8.0 / 9.0
G_NLPIN = [{'N_lay': int(n), 'n_segments': nseg(n),
            'DT_cascade': SPAN0 * n / (n - 1.0)} for n in NL]
G_NLFREE = [{'N_lay': int(n), 'n_segments': nseg(n)} for n in NL]


# ---- the configuration table of the paper -------------------------------
# Each row names its packing EXPLICITLY. The model default is disc, and a
# row that omits it silently reports the superseded v1.2 geometry.
DISC = dict(packing='disc', n_islands=1)
ROWS4 = dict(packing='rows', n_islands=4)
CONFIGS = [
    dict(label='adiabatic, k_r -> 0', k_rock=1e-9, **DISC),
    dict(label='isolated well', **DISC),
    dict(label='fan, 100 wells, disc', n_wells_field=100., **DISC),
    dict(label='fan, 100 wells, rows', n_wells_field=100., **ROWS4),
    dict(label='fan, 1000 wells, disc', n_wells_field=1000., **DISC),
    dict(label='fan, 1000 wells, rows', n_wells_field=1000., **ROWS4),
    dict(label='fan, 1000, 70 deg, disc', n_wells_field=1000.,
         fan_angle_deg=70., **DISC),
    dict(label='fan, 1000, 70 deg, rows', n_wells_field=1000.,
         fan_angle_deg=70., **ROWS4),
    dict(label='vertical array, disc', n_wells_field=1000.,
         fan_angle_deg=0., **DISC),
    dict(label='vertical array, rows', n_wells_field=1000.,
         fan_angle_deg=0., **ROWS4),
]
RING4 = dict(packing='ring', n_islands=4)
# ring is the double line closed into a loop, losing from BOTH faces. It
# differs from rows only by the two end-caps of each row -- 2*row_gap per
# island -- so at 1000 wells it is 1.6 % less perimeter than rows and at
# 100 wells 14 % less, the gap mattering more against a shorter row.
CONFIGS += [
    dict(label='fan, 100 wells, ring', n_wells_field=100., **RING4),
    dict(label='fan, 1000 wells, ring', n_wells_field=1000., **RING4),
    dict(label='fan, 1000, 70 deg, ring', n_wells_field=1000.,
         fan_angle_deg=70., **RING4),
    dict(label='vertical array, ring', n_wells_field=1000.,
         fan_angle_deg=0., **RING4),
]

CONFIG_FIXED = dict(n_segments=99)
