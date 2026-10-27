"""Publication-quality cycle diagrams and the reference state-point table.

FONTS. Times New Roman is proprietary and is not installed on Colab, on most
Linux images, or in this sandbox. Asking matplotlib for it by name therefore
produces a silent fallback to DejaVu Sans -- the request appears to succeed
and the figure comes out in the wrong face. The chain below asks for Times
first and then for its metric-compatible clones, Nimbus Roman and Liberation
Serif, which ARE present almost everywhere and are visually indistinguishable
from Times at figure scale. `active_serif()` reports which one resolved, so
the fallback is visible rather than silent.

ISOBARS ARE SAMPLED IN ENTHALPY, NOT TEMPERATURE. Inside the dome T and p are
not independent, so CoolProp cannot return a state from (T, p) there and a
temperature-sampled isobar silently breaks at the saturation line -- which is
exactly where every heat-exchanger leg in both cycles lives. Enthalpy
increases monotonically along an isobar through subcooled liquid, the dome and
superheat alike, so sampling in h gives one unbroken path with no branching on
phase. Throttles are drawn the same way, sampled in pressure at constant h,
rather than as the straight lines a T--s sketch usually shows.

CONNECTIVITY IS READ FROM THE MODEL, NOT GUESSED. The ORC is a REHEAT cycle:
3e expands to 5e, is reheated at p_int to 6e, then expands to 4e; 7e--10e are
the regenerative sub-circuit. The HTHP is a two-stage FLASH-ECONOMIZER cycle:
7h throttles to 8h into a flash drum at p_int, the saturated liquid 9h
throttles on to 4h, and the flash vapour 10h joins the second stage at 11h.
That ordering comes from the enthalpies in the state dictionaries, not from
the numbering.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

SERIF = ["Times New Roman", "Nimbus Roman", "Liberation Serif",
         "STIXGeneral", "DejaVu Serif"]

RED, BLUE, NAVY, PURPLE, GREY = "#b03028", "#1b6ca8", "#23356b", "#5b3a8e", "#6a6a6a"

# Annotations in a busy T-s plot land on curves. A white, edgeless patch
# behind the text masks whatever it crosses; without it the approach-
# temperature symbols and the isobar labels are unreadable where they sit
# on a process line.
MASK = dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.3)


def use_paper_style(scale=1.0):
    """Serif, hairline axes, no top/right spine. Call once per session.

    `scale` multiplies every type size, for figures that will be reproduced
    small in a two-column layout.
    """
    s = float(scale)
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": SERIF,
        "mathtext.fontset": "stix",
        "font.size": 10 * s, "axes.labelsize": 11 * s,
        "axes.titlesize": 11 * s,
        "xtick.labelsize": 9.5 * s, "ytick.labelsize": 9.5 * s,
        "legend.fontsize": 8.6 * s,
        "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "lines.linewidth": 1.6, "legend.frameon": False,
        "figure.dpi": 140, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def active_serif():
    """Which font actually resolved. Never assume the first choice won."""
    from matplotlib.font_manager import findfont, FontProperties
    got = findfont(FontProperties(family="serif"))
    norm = lambda x: x.lower().replace(" ", "").replace("-", "")
    for name in SERIF:
        if norm(name) in norm(got):
            return name, got
    return "(fallback)", got


# ----------------------------------------------------------------- helpers
def _dome(ax, fluid, n=320):
    Tc = CP.PropsSI("Tcrit", fluid)
    T = np.linspace(max(CP.PropsSI("Ttriple", fluid) + 1.0, 180.0), Tc - 0.15, n)
    sl = np.array([CP.PropsSI("S", "T", t, "Q", 0, fluid) for t in T]) / 1e3
    sv = np.array([CP.PropsSI("S", "T", t, "Q", 1, fluid) for t in T]) / 1e3
    ax.plot(np.r_[sl, sv[::-1]], np.r_[T, T[::-1]] - 273.15,
            color="0.30", lw=1.3, zorder=1)


def _iso_p(ax, fluid, p, h1, h2, n=160, **kw):
    """Constant-pressure path, sampled in enthalpy so it crosses the dome."""
    H = np.linspace(h1, h2, n) * 1e3
    T = np.array([CP.PropsSI("T", "P", p, "H", h, fluid) for h in H]) - 273.15
    s = np.array([CP.PropsSI("S", "P", p, "H", h, fluid) for h in H]) / 1e3
    ax.plot(s, T, zorder=4, **kw)


def _iso_h(ax, fluid, h, p1, p2, n=120, **kw):
    """Constant-enthalpy path (a throttle), sampled in pressure."""
    P = np.geomspace(p1, p2, n)
    T = np.array([CP.PropsSI("T", "P", p, "H", h * 1e3, fluid) for p in P]) - 273.15
    s = np.array([CP.PropsSI("S", "P", p, "H", h * 1e3, fluid) for p in P]) / 1e3
    ax.plot(s, T, zorder=3, **kw)


def _frame(ax, pts, padx=0.16, pady=0.18):
    """Zoom to the cycle. The full dome dwarfs both cycles otherwise."""
    s = np.array([p[0] for p in pts]); T = np.array([p[1] for p in pts])
    ds = (s.max() - s.min()) or 1.0
    dT = (T.max() - T.min()) or 1.0
    ax.set_xlim(s.min() - padx * ds, s.max() + padx * ds)
    ax.set_ylim(T.min() - pady * dT, T.max() + pady * dT)


def _label(ax, s, T, lab, off, marker="o", fc="white", ec="#1a1a1a"):
    ax.plot(s, T, marker, ms=4.4, mfc=fc, mec=ec, mew=1.2, zorder=6)
    ax.annotate(lab, (s, T), textcoords="offset points", xytext=off,
                fontsize=8.4, color="#1a1a1a", zorder=7, bbox=MASK)


def _axes(ax, title, loc="lower right"):
    ax.set_xlabel(r"specific entropy, $s$  /  kJ kg$^{-1}$ K$^{-1}$")
    ax.set_ylabel(r"temperature, $T$  /  $^{\circ}$C")
    ax.set_title(title, pad=8)
    ax.legend(loc=loc, handlelength=1.7, borderpad=0.2, labelspacing=0.35)


# ------------------------------------------------------------------- ORC
_ORC_OFF = {"1e": (8, -11), "2e": (8, 2), "3e": (-18, 7), "4e": (7, 0),
            "5e": (-20, -11), "6e": (6, 5), "7e": (-24, -13), "8e": (7, 4),
            "9e": (-25, 2), "10e": (7, -11)}


def ts_orc(case, rank, ax=None, zoom=True):
    """T-s diagram of the double-stage reheat ORC at the reference case."""
    f = case.fluid
    if ax is None:
        _, ax = plt.subplots(figsize=(5.7, 4.5))
    _dome(ax, f)
    S = lambda k: rank["s_" + k] / 1e3
    T = lambda k: rank["T_" + k] - 273.15
    h = lambda k: rank["h_" + k]

    _iso_p(ax, f, rank["p_evape"], h("2e"), h("3e"), color=RED, lw=1.8,
           label=r"heat addition, $p_\mathrm{evap}$")
    ax.plot([S("3e"), S("5e")], [T("3e"), T("5e")], color=NAVY, lw=1.8,
            zorder=4, label="turbine / pump")
    _iso_p(ax, f, rank["p_inte"], h("5e"), h("6e"), color=RED, lw=1.8,
           ls=(0, (5, 2)), label=r"reheat, $p_\mathrm{int}$")
    ax.plot([S("6e"), S("4e")], [T("6e"), T("4e")], color=NAVY, lw=1.8, zorder=4)
    _iso_p(ax, f, rank["p_conde"], h("4e"), h("1e"), color=BLUE, lw=1.8,
           label=r"condensation, $p_\mathrm{cond}$")
    ax.plot([S("1e"), S("2e")], [T("1e"), T("2e")], color=NAVY, lw=1.8, zorder=4)

    pts = []
    for k in ("1e", "2e", "3e", "4e", "5e", "6e"):
        _label(ax, S(k), T(k), k, _ORC_OFF[k]); pts.append((S(k), T(k)))
    # 7e--10e differ by hundredths of a degree and would be four overlapping
    # labels on one dot; they are drawn as four markers but annotated once.
    for k in ("7e", "8e", "9e", "10e"):
        ax.plot(S(k), T(k), "s", ms=4.4, mfc="#b9b9c8", mec="#1a1a1a",
                mew=1.2, zorder=6)
        pts.append((S(k), T(k)))
    ax.annotate("7e--10e", (S("8e"), T("8e")), textcoords="offset points",
                xytext=(8, 3), fontsize=8.4, zorder=7)
    ax.plot([], [], "s", ms=4.4, mfc="#b9b9c8", mec="#1a1a1a", mew=1.2,
            label="regenerator, 7e--10e")

    if zoom:
        _frame(ax, pts)
    _axes(ax, rf"(a) ORC, {f}:  $\eta_\mathrm{{ORC}} = {rank['rank_eff']:.4f}$,"
              rf"  $y = {rank['y_frac']:.4f}$", loc="upper left")
    return ax


# ------------------------------------------------------------------ HTHP
_HP_OFF = {"1h": (7, -10), "2h": (7, 1), "3h": (7, 2), "4h": (-19, -11),
           "5h": (-20, 4), "6h": (-20, -11), "7h": (-22, -11), "8h": (-8, -12),
           "9h": (-23, 1), "10h": (-24, -5), "11h": (7, 1), "12h": (7, -10),
           "13h": (7, 2)}


def ts_hthp(case, hp, ax=None, zoom=True):
    """T-s diagram of the two-stage flash-economizer HTHP."""
    f = case.refrig
    if ax is None:
        _, ax = plt.subplots(figsize=(5.7, 4.5))
    _dome(ax, f)
    S = lambda k: hp["s_" + k] / 1e3
    T = lambda k: hp["T_" + k] - 273.15
    h = lambda k: hp["h_" + k]
    pe, pi, pc = hp["p_evaph"], hp["p_inth"], hp["p_condh"]

    _iso_p(ax, f, pe, h("4h"), h("1h"), color=BLUE, lw=1.8,
           label=r"evaporation $+$ superheat, $p_\mathrm{evap}$")
    ax.plot([S("1h"), S("5h")], [T("1h"), T("5h")], color=PURPLE, lw=1.8,
            zorder=4, label="compression")
    ax.plot([S("11h"), S("2h")], [T("11h"), T("2h")], color=PURPLE, lw=1.8, zorder=4)
    _iso_p(ax, f, pc, h("2h"), h("7h"), color=RED, lw=1.8,
           label=r"condensation $+$ subcool, $p_\mathrm{cond}$")
    _iso_h(ax, f, h("7h"), pc, pi, color=GREY, lw=1.3, ls=(0, (2, 2)),
           label=r"throttle, $h =$ const.")
    _iso_h(ax, f, h("9h"), pi, pe, color=GREY, lw=1.3, ls=(0, (2, 2)))
    ax.plot([S("10h"), S("11h")], [T("10h"), T("11h")], color=PURPLE, lw=1.0,
            ls=(0, (1, 2)), zorder=3, label="flash vapour to stage 2")

    # 5h, 6h and 11h are the SAME state in the model -- three labels on one
    # dot. Annotated once rather than stacked.
    pts = []
    for k in [f"{i}h" for i in range(1, 14)]:
        sq = k in ("5h", "6h", "9h", "10h", "11h", "12h")
        if k in ("5h", "6h"):
            ax.plot(S(k), T(k), "s", ms=4.4, mfc="#c3b9d8", mec="#1a1a1a",
                    mew=1.2, zorder=6)
        else:
            _label(ax, S(k), T(k),
                   "5h, 6h, 11h" if k == "11h" else k, _HP_OFF[k],
                   marker="s" if sq else "o", fc="#c3b9d8" if sq else "white")
        pts.append((S(k), T(k)))

    if zoom:
        _frame(ax, pts)
    _axes(ax, rf"(b) HTHP, {f}:  COP $= {hp['hp_cop']:.4f}$,"
              rf"  $x_\mathrm{{8h}} = {hp['x_8h']:.4f}$", loc="upper left")
    return ax


def ts_both(case, rank, hp, figsize=(11.8, 4.7)):
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    ts_orc(case, rank, ax=axes[0])
    ts_hthp(case, hp, ax=axes[1])
    fig.tight_layout()
    return fig, axes


# --------------------------------------------------- both on one axes
_CMB_ORC = {"1e": (-22, -11), "2e": (7, 3), "3e": (-16, 8), "4e": (7, -2),
            "5e": (-19, -13), "6e": (-20, 7)}
_CMB_HP = {"1h": (7, -11), "2h": (-20, 5), "3h": (-17, 7), "4h": (-20, -13),
           "7h": (-21, -12), "8h": (-6, -13), "9h": (-24, 5),
           "10h": (-26, -10), "11h": (7, 0), "12h": (7, -11), "13h": (-21, 3)}


def ts_combined(case, rank, hp, Tdict, ax=None, show_water=True,
                water="level", model=None, isobars=True, show_values=False,
                title=None):
    """Both cycles on one T-s axes with the external water levels.

    This is the coupling diagram: the two cycles never exchange working
    fluid, they exchange water, so the horizontal lines are the only thing
    that makes the pair a single machine rather than two unrelated plots.
    The four levels are the sink, the HTHP evaporator (geothermal source),
    the HTHP condenser inlet (charge loop) and the ORC evaporator inlet
    (the store outlet, T_2d) -- the last being the one quantity the
    borehole model actually sets.
    """
    f = case.fluid
    if ax is None:
        _, ax = plt.subplots(figsize=(7.6, 5.6))
    _dome(ax, f)
    Se = lambda k: rank["s_" + k] / 1e3
    Te = lambda k: rank["T_" + k] - 273.15
    he = lambda k: rank["h_" + k]
    Sh = lambda k: hp["s_" + k] / 1e3
    Th = lambda k: hp["T_" + k] - 273.15
    hh = lambda k: hp["h_" + k]

    pts = []

    # ---- water as REAL GLIDES, duty-mapped onto the working fluid's
    # entropy, with every pinch marked. `water="level"` falls back to the
    # four inlet temperatures as horizontal lines, which is the conventional
    # presentation but hides where each exchanger actually binds.
    if show_water and water == "glide":
        if model is None:
            raise ValueError("water='glide' needs the model module, "
                             "e.g. ts_combined(..., water='glide', model=B)")
        ov = _water_overlay(model, case, rank, hp, Tdict)
        # symbol, not value: the approach temperatures are specified
        # quantities and belong in the nomenclature, not on the drawing
        style = {
            "ORC evaporator": ("#1a1a1a",
                               r"ORC evaporator water, $T_\mathrm{2d}\to T_\mathrm{3d}$",
                               r"$\Delta T_\mathrm{pinch,ORC}$", (9, 7)),
            "ORC condenser": ("#2a6fb0", "sink water, ORC condenser",
                              r"$\Delta T_\mathrm{E,sink}$", (10, -14)),
            "HTHP condenser": ("#0a9396",
                               r"charge water, $T_\mathrm{2c}\to T_\mathrm{3c}$",
                               r"$\Delta T_\mathrm{2H,3C}$", (8, -2)),
            "HTHP evaporator": ("#9c27a0", "source water, HTHP evaporator",
                                r"$\Delta T_\mathrm{pinch,HPE}$", (9, -15)),
        }
        for nm, (col, lab, sym, off) in style.items():
            d = ov[nm]
            txt = f"{lab}  ({d['pinch_K']:.1f} K)" if show_values else lab
            _draw_water(ax, d, col, label=txt,
                        symbol=None if show_values else sym,
                        dx=off[0], dy=off[1])
            pts += list(zip(d["s"], d["T_w"]))

    # ---- water-side levels, drawn first so the cycles sit on top
    if show_water and water == "level":
        levels = [
            (case.T_sink_C, "#2a6fb0", r"sink, ORC condenser"),
            (case.T_source_C, "#9c27a0", r"source, HTHP evaporator"),
            (Tdict["T_2c"] - 273.15, "#0a9396", r"HTHP condenser inlet, $T_\mathrm{2c}$"),
            (Tdict["T_2d"] - 273.15, "#1a1a1a", r"ORC evaporator inlet, $T_\mathrm{2d}$"),
        ]
        for Tl, col, lab in levels:
            ax.axhline(Tl, color=col, lw=1.0, ls=(0, (1, 2.5)), zorder=2,
                       label=f"{lab}  ({Tl:.0f} $^\\circ$C)")
            pts.append((0.0, Tl))

    # ---- ORC, blue
    _iso_p(ax, f, rank["p_evape"], he("2e"), he("3e"), color=NAVY, lw=1.5)
    ax.plot([Se("3e"), Se("5e")], [Te("3e"), Te("5e")], color=NAVY, lw=1.5, zorder=4)
    _iso_p(ax, f, rank["p_inte"], he("5e"), he("6e"), color=NAVY, lw=1.5)
    ax.plot([Se("6e"), Se("4e")], [Te("6e"), Te("4e")], color=NAVY, lw=1.5, zorder=4)
    _iso_p(ax, f, rank["p_conde"], he("4e"), he("1e"), color=NAVY, lw=1.5)
    ax.plot([Se("1e"), Se("2e")], [Te("1e"), Te("2e")], color=NAVY, lw=1.5, zorder=4)

    # ---- HTHP, red
    g = case.refrig
    _iso_p(ax, g, hp["p_evaph"], hh("4h"), hh("1h"), color=RED, lw=1.5)
    ax.plot([Sh("1h"), Sh("5h")], [Th("1h"), Th("5h")], color=RED, lw=1.5, zorder=4)
    ax.plot([Sh("11h"), Sh("2h")], [Th("11h"), Th("2h")], color=RED, lw=1.5, zorder=4)
    _iso_p(ax, g, hp["p_condh"], hh("2h"), hh("7h"), color=RED, lw=1.5)
    _iso_h(ax, g, hh("7h"), hp["p_condh"], hp["p_inth"], color=RED, lw=1.1,
           ls=(0, (3, 2)))
    _iso_h(ax, g, hh("9h"), hp["p_inth"], hp["p_evaph"], color=RED, lw=1.1,
           ls=(0, (3, 2)))
    ax.plot([Sh("10h"), Sh("11h")], [Th("10h"), Th("11h")], color=RED, lw=1.0,
            ls=(0, (1, 2)), zorder=3)

    for k, off in _CMB_ORC.items():
        _label(ax, Se(k), Te(k), k, off, fc=NAVY, ec=NAVY); pts.append((Se(k), Te(k)))
    for k in ("7e", "8e", "9e", "10e"):
        ax.plot(Se(k), Te(k), "o", ms=4.0, mfc=NAVY, mec=NAVY, zorder=6)
        pts.append((Se(k), Te(k)))
    ax.annotate("7e--10e", (Se("8e"), Te("8e")), textcoords="offset points",
                xytext=(-26, 5), fontsize=8.4, color=NAVY, zorder=7)
    for k, off in _CMB_HP.items():
        lab = "5h, 6h, 11h" if k == "11h" else k
        _label(ax, Sh(k), Th(k), lab, off, fc=RED, ec=RED); pts.append((Sh(k), Th(k)))
    for k in ("5h", "6h"):
        ax.plot(Sh(k), Th(k), "o", ms=4.0, mfc=RED, mec=RED, zorder=6)

    ax.plot([], [], "o-", color=NAVY, ms=4.4, lw=1.5, label="ORC (Rankine)")
    ax.plot([], [], "o-", color=RED, ms=4.4, lw=1.5, label="HTHP (heat pump)")

    _frame(ax, pts, padx=0.10, pady=0.12)

    # ---- constant-pressure reference curves, drawn last so the frame is
    # known and they can be clipped to it. The two superheated legs that
    # look like chords -- 13h->1h and 4e to the saturated vapour -- ARE on
    # their isobars; those isobars are straight there to 0.009 K and
    # 0.057 K respectively. Drawing the isobars makes that visible rather
    # than something the reader has to take on trust.
    if isobars:
        y0, y1 = ax.get_ylim()
        levels = sorted({round(p, 3) for p in
                         (rank["p_conde"], rank["p_inte"], rank["p_evape"],
                          hp["p_evaph"], hp["p_inth"], hp["p_condh"])})
        # Group pressures that are within a few per cent: their isobars are
        # one line at this scale, so two labels would sit on top of each
        # other. Here p_int,HP = 3.82 bar and p_evap,ORC = 3.70 bar.
        groups = []
        for p in levels:
            if groups and abs(p - groups[-1][-1]) / p < 0.06:
                groups[-1].append(p)
            else:
                groups.append([p])
        fracs = [0.97, 0.52, 0.86, 0.97, 0.68, 0.40]
        for q, grp in enumerate(groups):
            # only the NUMBER goes in math mode; "bar" set upright
            lab = "$" + r"\,/\,".join(f"{pp / 1e5:.2f}" for pp in grp) + "$ bar"
            for j, pp in enumerate(grp):
                _isobar_bg(ax, f, pp, max(y0, -40.0), min(y1 + 40.0, 230.0),
                           label=(lab if j == 0 else None),
                           label_frac=fracs[q % len(fracs)])
        ax.set_ylim(y0, y1)

    ax.set_xlabel(r"specific entropy, $s$  /  kJ kg$^{-1}$ K$^{-1}$")
    ax.set_ylabel(r"temperature, $T$  /  $^{\circ}$C")
    ax.set_title(title or "Cyclopentane in the Heat Pump and Organic "
                          "Rankine Cycles", pad=9)
    ax.legend(loc="upper left", handlelength=1.9, borderpad=0.25,
              labelspacing=0.35)
    return ax


# ------------------------------------------- profiles along the well
#
# THE COORDINATE. The solver's z is NOT depth. Quoting formation_profile:
# "On CHARGE the water enters the fin tubes at the BOTTOM -- it arrives
# there down the downcomer -- so the solver's z = 0 is the bottom of the
# well and depth decreases with z." So depth = L_well - z, and plotting
# against the raw z and calling it depth puts the well upside down. Use
# `depth_seg`, taken from formation_profile so there is one definition.
#
# WHICH FRAME THE HISTORY IS IN. Both half-cycles come back in the CHARGE
# frame -- index 0 is the BOTTOM for charge and discharge alike -- even
# though the discharge march is run mirrored. Verified, not assumed: the
# test is that the PCM can never sit above its own melting point while it
# is freezing, and at the end of discharge max(T_pcm - T_m) is +0.00 K read
# forward against +26.9 K read reversed. The fluid boundary alone does not
# settle it, because the inlet is at opposite ends in the two half-cycles.
#
# PHYSICALLY. The hottest layer (150 C) sits at the BOTTOM, against the
# hottest rock, which is what keeps its loss down. On charge the water
# descends the insulated downcomer, enters the fin tubes at the bottom and
# cools as it RISES; on discharge it enters at the top, warms as it
# DESCENDS, and returns up the downcomer. Counter-current.

def _time_subset(t_hist, n=9):
    """Roughly log-spaced time levels, always including first and last."""
    t = np.asarray(t_hist, float)
    want = np.unique(np.concatenate(
        ([t[0]], t[0] + (t[-1] - t[0]) * np.logspace(-3, 0, n - 1))))
    return [int(np.argmin(np.abs(t - w))) for w in want]


def _staircase_xy(z, T_m_seg):
    """Explicit staircase polyline, so it can be drawn either way round.

    `ax.step` only ever steps in y as a function of x, which is wrong when
    depth is the ORDINATE. Building the vertices by hand keeps one
    definition of the staircase for both orientations.
    """
    z = np.asarray(z, float)
    T = np.asarray(T_m_seg, float) - 273.15
    dz = z[1] - z[0] if len(z) > 1 else 1.0
    edges = np.concatenate(([z[0] - dz / 2], z + dz / 2))
    return np.repeat(T, 2), np.repeat(edges, 2)[1:-1]


def _staircase(ax, z, T_m_seg, vertical=False, **kw):
    kw = dict(color="#1a1a1a", lw=1.1, ls=(0, (4, 2)), zorder=6) | kw
    T, zz = _staircase_xy(z, T_m_seg)
    ax.plot(*((T, zz) if vertical else (zz, T)), **kw)


def depth_axis(model, case, n=None):
    """(depth of segment centres, depth of fluid nodes) in metres.

    Taken from the model's own formation_profile so the plotting layer
    cannot drift from the solver's idea of which end is the bottom.
    """
    n = n or case.n_segments
    depth_seg = np.asarray(model.formation_profile(case, n,
                                                   discharge=False)[1], float)
    L = case.L_well
    return depth_seg, L - np.linspace(0.0, L, n + 1)


def well_profiles(r, case, field="T_fluid", axes=None, n_times=9,
                  figsize=None, cmap="viridis", ylabel=None,
                  title=None, staircase=False, scale=1.0, offset=0.0,
                  logy=False, vertical=True, depth_seg=None, depth_node=None,
                  overlay=None, overlay_label=None, overlay_ls=(0, (2, 1.6))):
    """One field along the well at several times, charge and discharge.

    `r` must come from simulate_css_corrected(..., record=True).
    `depth_seg` / `depth_node` come from `depth_axis`; without them the raw
    solver coordinate is used, which runs from the BOTTOM of the well and is
    not depth.
    `overlay` names a second field drawn dashed in the same colour per time
    (used for T_pcm beside T_fluid).
    """
    if figsize is None:
        figsize = (8.8, 5.8) if vertical else (11.6, 4.4)
    if axes is None:
        _, axes = plt.subplots(1, 2, figsize=figsize,
                              sharex=vertical, sharey=not vertical)
    cm = plt.get_cmap(cmap)
    qlabel = ylabel or field
    dlabel = r"depth / m"
    for ax, half, lab in ((axes[0], "charge", "(a) charging"),
                          (axes[1], "discharge", "(b) discharging")):
        H, th = r[half]["history"], r[half]["t_hist"]
        zs = depth_seg if depth_seg is not None else r[half]["z"]
        zn = depth_node if depth_node is not None else None
        idx = _time_subset(th, n_times)
        Y = np.asarray(H[field], float)
        O = np.asarray(H[overlay], float) if overlay else None
        for j, i in enumerate(idx):
            col = cm(j / max(1, len(idx) - 1))
            for arr, ls in ((Y, "-"), (O, overlay_ls)):
                if arr is None:
                    continue
                y = arr[i] * scale + offset
                if len(y) == len(zs):
                    zz = zs
                elif zn is not None and len(y) == len(zn):
                    zz = zn
                else:
                    zz = np.linspace(zs[0], zs[-1], len(y))
                ax.plot(*((y, zz) if vertical else (zz, y)), color=col,
                        lw=1.3, ls=ls,
                        label=(f"{th[i] / 60:.1f}" if ls == "-" else None))
        if staircase:
            _staircase(ax, zs, r["T_m_seg"], vertical=vertical)
        ax.set_title(lab, pad=7)
        if vertical:
            ax.set_xlabel(qlabel)
            if logy:
                ax.set_xscale("log")
        else:
            ax.set_xlabel(dlabel)
            if logy:
                ax.set_yscale("log")
    if vertical:
        axes[0].set_ylabel(dlabel)
        zs = depth_seg if depth_seg is not None else r["charge"]["z"]
        for ax in axes:
            ax.set_ylim(max(zs) * 1.02, -0.02 * max(zs))   # 0 at the top
    else:
        axes[0].set_ylabel(qlabel)
    if overlay:
        axes[1].plot([], [], color="0.35", lw=1.3, ls=overlay_ls,
                     label=overlay_label or overlay)
    if staircase:
        axes[1].plot([], [], color="#1a1a1a", lw=1.1, ls=(0, (4, 2)),
                     label=r"$T_{m,k}$")
    leg = axes[1].legend(title="time / min", ncol=2, fontsize=7.4,
                         title_fontsize=7.8, loc="best", handlelength=1.6,
                         labelspacing=0.25, columnspacing=1.0)
    leg._legend_box.align = "left"
    if title:
        axes[0].figure.suptitle(title, y=0.99)
    return axes


def well_integrals(r, case, axes=None, figsize=(11.6, 4.2)):
    """Heat rate integrated along the well, and the melted fraction, vs time.

    The melted fraction is the volume-weighted mean of the local melted
    cross-section, A_melt / A_avail, which is the same quantity the model
    reports as the cycled fraction at the end of each half-cycle.
    """
    if axes is None:
        _, axes = plt.subplots(1, 2, figsize=figsize)
    A_avail = pcm_capacities_A(case)
    for ax, half, lab in ((axes[0], "charge", "(a) charging"),
                          (axes[1], "discharge", "(b) discharging")):
        H, th, z = r[half]["history"], r[half]["t_hist"], r[half]["z"]
        dz = z[1] - z[0]
        q = np.asarray(H["q_prime"], float).sum(axis=1) * dz / 1e3   # kW/well
        eps = np.asarray(H["A_melt"], float).mean(axis=1) / A_avail
        ax.plot(th / 60.0, q, color=RED, lw=1.8)
        ax.set_xlabel("time / min")
        ax.set_ylabel(r"$\dot{Q}$ per well / kW", color=RED)
        ax.tick_params(axis="y", colors=RED)
        axb = ax.twinx()
        axb.plot(th / 60.0, eps, color="#2f6f4f", lw=1.8, ls=(0, (5, 2)))
        axb.set_ylabel("melted fraction of the PCM", color="#2f6f4f")
        axb.tick_params(axis="y", colors="#2f6f4f")
        axb.spines["right"].set_visible(True)
        ax.set_title(lab, pad=7)
    return axes


def _resample(Y, t_src, t_dst):
    """Linear interpolation of a (time, segment) field onto a new time grid.

    The march records on a roughly LOGARITHMIC time grid, which is right for
    static profiles but makes an animation played at constant frame rate
    appear to accelerate. Resampling onto a uniform grid puts the playback
    clock in proportion to physical time; the sub-minute transient is then
    genuinely brief on screen, which is the honest depiction.
    """
    Y = np.asarray(Y, float)
    return np.column_stack([np.interp(t_dst, t_src, Y[:, j])
                            for j in range(Y.shape[1])])


def animate_well(r, case, path, fps=12, frames_per_half=60, dpi=110,
                 writer=None):
    """Animate the well through one full cycle at cyclic steady state.

    Three panels with depth on the ordinate: fluid temperature against the
    melting staircase, local melted fraction, and the overall coefficient.
    Axis limits are fixed over the whole cycle so that nothing jumps between
    the half-cycles and the eye can compare them directly.

    `path` ending .gif uses pillow, .mp4 uses ffmpeg.
    """
    import matplotlib.animation as manim
    import boreas as _b
    A = pcm_capacities_A(case)
    zs, zn = depth_axis(_b, case, len(r["T_m_seg"]))
    seq = []                              # (phase, t_s, Tf, Tp, eps, Ui)
    for half, lab in (("charge", "charging"), ("discharge", "discharging")):
        H, th = r[half]["history"], r[half]["t_hist"]
        tu = np.linspace(th[0], th[-1], frames_per_half)
        Tf = _resample(H["T_fluid"], th, tu) - 273.15
        Tp = _resample(H["T_pcm"], th, tu) - 273.15
        eps = _resample(H["A_melt"], th, tu) / A
        Ui = _resample(H["U_i"], th, tu)
        for k in range(frames_per_half):
            seq.append((lab, tu[k], Tf[k], Tp[k], eps[k], Ui[k]))

    allT = np.concatenate([np.concatenate([s[2], s[3]]) for s in seq])
    allU = np.concatenate([s[5] for s in seq])
    lims = [(allT.min() - 3, allT.max() + 3), (-0.04, 1.04),
            (0.0, allU.max() * 1.05)]
    labs = [r"temperature / $^\circ$C", r"$A_\mathrm{melt}/A_\mathrm{av}$",
            r"$U_i$ / W m$^{-2}$ K$^{-1}$"]
    cols = [NAVY, "#2f6f4f", PURPLE]

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 5.4), sharey=True)
    lines = []
    for ax, lim, lab, col in zip(axes, lims, labs, cols):
        ax.set_xlim(*lim)
        ax.set_xlabel(lab)
        ln, = ax.plot([], [], color=col, lw=2.0)
        lines.append(ln)
    lpcm, = axes[0].plot([], [], color=RED, lw=1.5, ls=(0, (2, 1.6)))
    lines.append(lpcm)
    _staircase(axes[0], zs, r["T_m_seg"], vertical=True)
    axes[0].plot([], [], color=NAVY, lw=2.0, label=r"$T_f$")
    axes[0].plot([], [], color=RED, lw=1.5, ls=(0, (2, 1.6)),
                 label=r"$T_\mathrm{pcm}$")
    axes[0].plot([], [], color="#1a1a1a", lw=1.1, ls=(0, (4, 2)),
                 label=r"$T_{m,k}$")
    axes[0].legend(loc="lower left", fontsize=8.2, handlelength=1.8)
    axes[0].set_ylabel(r"depth / m")
    # Depth increases downward: setting the limits in descending order IS
    # the inversion -- invert_yaxis() on top of this would undo it.
    top = float(np.max(zs))
    axes[0].set_ylim(top * 1.02, -0.02 * top)
    head = fig.suptitle("", y=0.975, fontsize=11)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.11,
                        wspace=0.10)

    def draw(i):
        lab, t, Tf, Tp, eps, Ui = seq[i]
        lines[0].set_data(Tf, zn if len(Tf) == len(zn) else zs)
        lines[1].set_data(eps, zs)
        lines[2].set_data(Ui, zs)
        lines[3].set_data(Tp, zs)
        head.set_text(f"{lab}   $t$ = {t / 60:6.1f} min   "
                      f"mean melted fraction {eps.mean():.3f}")
        return lines + [head]

    anim = manim.FuncAnimation(fig, draw, frames=len(seq), interval=1000 / fps,
                               blit=False)
    if writer is None:
        writer = "pillow" if str(path).endswith(".gif") else "ffmpeg"
    anim.save(str(path), writer=writer, fps=fps, dpi=dpi)
    plt.close(fig)
    return len(seq)


def pcm_capacities_A(case):
    """PCM cross-section available per fin tube [m^2]."""
    import boreas as _b
    return _b.pcm_capacities(case)[0]


# ----------------------------------------------------------- heat maps
def heatmap(rows, xkey, ykey, zkey, xs, ys, ax=None, cmap="viridis",
            title="", cbar_label="", fmt="{:.3f}", mark=None, annotate=True,
            vmin=None, vmax=None, cbar=True, fontsize=7.6,
            best="max", box_ref=True):
    """One metric over a two-parameter grid, with the cells labelled.

    Cells are ANNOTATED with their value. A heat map whose colour scale is
    the only channel invites the reader to misjudge small differences, and
    several of these maps are deliberately flat -- the flatness is the
    result, and it is only legible if the numbers are on the cells.

    Infeasible points (validate_case blocking, or an exception) are hatched
    rather than dropped, so a hole in the map reads as a constraint and not
    as a bug in the plot.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6.2, 4.8))
    m = {(r[xkey], r[ykey]): r for r in rows}
    Z = np.full((len(ys), len(xs)), np.nan)
    bad = np.zeros_like(Z, dtype=bool)
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            r = m.get((float(x), float(y))) or m.get((x, y))
            if r is None or not r.get("feasible", True) or zkey not in r:
                bad[i, j] = True
            else:
                Z[i, j] = r[zkey]
    im = ax.imshow(Z, origin="lower", aspect="auto", cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation="nearest",
                   extent=(-0.5, len(xs) - 0.5, -0.5, len(ys) - 0.5))
    # Which cells get called out. The reference case and the best cell are
    # marked by STYLING THE VALUE, not by drawing a symbol on top of it: an
    # overlaid star or ring hides the number it is pointing at, which is the
    # one thing the reader wants from an annotated map.
    ref_ij = None
    if mark is not None and mark[0] in list(xs) and mark[1] in list(ys):
        ref_ij = (list(ys).index(mark[1]), list(xs).index(mark[0]))
    best_ij = None
    if best and np.isfinite(Z).any():
        pick = np.nanargmax(Z) if best == "max" else np.nanargmin(Z)
        best_ij = np.unravel_index(pick, Z.shape)

    for i in range(len(ys)):
        for j in range(len(xs)):
            if bad[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           facecolor="0.92", edgecolor="0.6",
                                           hatch="///", lw=0.5))
                continue
            if not annotate:
                continue
            lo, hi = np.nanmin(Z), np.nanmax(Z)
            t = 0.5 if hi == lo else (Z[i, j] - lo) / (hi - lo)
            col = "white" if t < 0.55 else "#111111"
            fs, weight = fontsize, "normal"
            if best_ij is not None and (i, j) == tuple(best_ij):
                fs, weight = fontsize * 1.45, "bold"
            elif ref_ij is not None and (i, j) == ref_ij:
                weight = "bold"
            ax.text(j, i, fmt.format(Z[i, j]), ha="center", va="center",
                    fontsize=fs, fontweight=weight, color=col, zorder=5)
    if ref_ij is not None and box_ref:
        ax.add_patch(plt.Rectangle((ref_ij[1] - 0.5, ref_ij[0] - 0.5), 1, 1,
                                   fill=False, edgecolor="#d62828", lw=1.8,
                                   zorder=6))
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels([f"{x:g}" for x in xs])
    ax.set_yticks(range(len(ys)))
    ax.set_yticklabels([f"{y:g}" for y in ys])
    ax.set_title(title, pad=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(True)
    if cbar:
        cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(cbar_label)
        cb.outline.set_linewidth(0.6)
    return ax, Z


def spread(Z):
    """Peak-to-peak as a percentage of the mean -- is this map flat?"""
    f = Z[np.isfinite(Z)]
    return 100.0 * (f.max() - f.min()) / abs(f.mean()) if f.size else np.nan


# --------------------------------------- surface cycles vs T_m
def cycle_vs_Tm(rows_autoP, rows_fixedP, case, axes=None, figsize=(13.2, 4.3)):
    """How the two surface machines trade against storage temperature.

    `rows_autoP` has the water-loop pressure raised to keep the loop liquid
    at each T_m; `rows_fixedP` holds it at the shipped value and goes
    infeasible once the charging inlet reaches saturation. Plotting both is
    the point: the first gives the thermodynamic optimum, the second says
    what it costs to reach it.
    """
    if axes is None:
        _, axes = plt.subplots(1, 3, figsize=figsize)
    ok = [r for r in rows_autoP if r.get("feasible")]
    T = np.array([r["T_m_C"] for r in ok])
    eta = np.array([r["eta_ORC"] for r in ok])
    cop = np.array([r["COP"] for r in ok])
    prod = np.array([r["product"] for r in ok])
    Pb = np.array([r["P_bar"] for r in ok])
    T2d = np.array([r["T_2d_C"] for r in ok])
    T3e = np.array([r["T_3e_C"] for r in ok])
    okf = [r["T_m_C"] for r in rows_fixedP if r.get("feasible")]
    T_lim = max(okf) if okf else None
    i_opt = int(np.argmax(prod))
    T_ref = case.T_m - 273.15

    def shade(ax):
        if T_lim is not None and T_lim < T.max():
            ax.axvspan(T_lim, T.max(), color="0.88", zorder=0)
            # y is in AXES fraction here, so it must be 1.0 -- not
            # ax.get_ylim()[1], which is a data value and would place the
            # text at several hundred percent of the axes height. With
            # savefig.bbox='tight' that silently inflates the saved figure
            # to tens of thousands of pixels instead of raising.
            ax.annotate(f"needs $P>{case.P / 1e5:.0f}$ bar",
                        (0.5 * (T_lim + T.max()), 1.0),
                        ha="center", va="top", fontsize=8.0, color="0.35",
                        xycoords=("data", "axes fraction"), xytext=(0, -6),
                        textcoords="offset points")
        ax.axvline(T_ref, color="0.45", lw=0.9, ls=(0, (1, 2.5)), zorder=1)

    ax = axes[0]
    ax.plot(T, eta, color=NAVY, lw=1.9, label=r"$\eta_\mathrm{ORC}$")
    ax.set_ylabel(r"$\eta_\mathrm{ORC}$", color=NAVY)
    ax.tick_params(axis="y", colors=NAVY)
    axb = ax.twinx()
    axb.plot(T, cop, color=RED, lw=1.9, label="COP")
    axb.set_ylabel("COP", color=RED)
    axb.tick_params(axis="y", colors=RED)
    axb.spines["right"].set_visible(True)
    shade(ax)
    ax.set_title("(a) the two machines pull opposite ways", pad=8)

    ax = axes[1]
    ax.plot(T, prod, color="#2f6f4f", lw=2.0)
    ax.plot([T[i_opt]], [prod[i_opt]], "o", ms=6, mfc="white",
            mec="#2f6f4f", mew=1.6, zorder=5)
    ax.annotate(rf"optimum {T[i_opt]:.0f} $^\circ$C" "\n"
                rf"$\eta_\mathrm{{ORC}}\,$COP$\,={prod[i_opt]:.4f}$",
                (T[i_opt], prod[i_opt]), textcoords="offset points",
                xytext=(10, -26), fontsize=8.4)
    j = int(np.argmin(np.abs(T - T_ref)))
    ax.plot([T[j]], [prod[j]], "s", ms=5, mfc="white", mec="0.35", mew=1.3,
            zorder=5)
    ax.annotate(rf"reference {T_ref:.0f} $^\circ$C, {prod[j]:.4f}"
                "\n" rf"({100 * (prod[i_opt] / prod[j] - 1):+.1f} % to optimum)",
                (T[j], prod[j]), textcoords="offset points", xytext=(-26, -40),
                fontsize=8.4, ha="left")
    shade(ax)
    axb = ax.twinx()
    axb.plot(T, Pb, color="#9c27a0", lw=1.3, ls=(0, (4, 2)))
    axb.set_ylabel("required loop pressure  /  bar", color="#9c27a0")
    axb.tick_params(axis="y", colors="#9c27a0")
    axb.spines["right"].set_visible(True)
    ax.set_ylabel(r"$\eta_\mathrm{ORC}\times$ COP")
    ax.set_title("(b) a shallow optimum, bought with pressure", pad=8)

    ax = axes[2]
    ax.plot(T, T2d - T3e, color="#b06a28", lw=1.9)
    shade(ax)
    ax.set_ylabel(r"$T_\mathrm{2d} - T_\mathrm{3e}$  /  K")
    ax.set_title("(c) the pinch shortfall is not a quirk\nof the reference case",
                 pad=8)
    ax.set_ylim(0, max(60.0, float((T2d - T3e).max()) * 1.15))

    for a in axes:
        a.set_xlabel(r"storage melting temperature, $T_\mathrm{m}$  /  $^{\circ}$C")
    return axes, dict(T_opt_C=float(T[i_opt]), product_opt=float(prod[i_opt]),
                      T_ref_C=float(T_ref), product_ref=float(prod[j]),
                      gain_pct=float(100 * (prod[i_opt] / prod[j] - 1)),
                      P_feasible_limit_C=T_lim)


# ------------------------- water streams drawn ON the T-s diagram
#
# A water stream has its own entropy, on a scale that has nothing to do with
# cyclopentane's, so it cannot simply be plotted on these axes. What CAN be
# plotted, and is what matters for a pinch, is the water temperature at each
# point OF THE EXCHANGER, placed at the entropy the working fluid has at that
# same point. The two streams are tied together by cumulative duty: at a given
# Q the working fluid is at (s, T_wf) and the water is at T_w. Plotting
# (s, T_w) puts the water line directly above (evaporator) or below
# (condenser) the working-fluid path, and the VERTICAL GAP at every abscissa
# is then the real local approach, with the pinch as its narrowest point.
#
# This is a duty mapping, not a thermodynamic path: the water line is not a
# process in cyclopentane's coordinates and its abscissa carries no meaning
# for the water itself. Said plainly in the legend, it is an honest and very
# readable construction; left unsaid it would be a category error.

def _water_overlay(model, case, rank, hp, Tdict, n=400):
    """name -> dict(s, T_wf, T_w, pinch_K, i_pinch, hot) for each exchanger.

    The stream definitions, the water endpoints and the duty bases are taken
    from `exchanger_UA`, so the pinch recovered here is the same number that
    function reports rather than a parallel construction of it.
    """
    F, R, P = case.fluid, case.refrig, case.P
    out = {}

    def pack(name, s, T_wf, T_w, hot):
        dT = (T_wf - T_w) if hot else (T_w - T_wf)
        i = int(np.argmin(dT))
        out[name] = dict(s=s, T_wf=T_wf - 273.15, T_w=T_w - 273.15,
                         pinch_K=float(dT[i]), i_pinch=i, hot=hot)

    def s_of(h, p, fluid):
        return np.array([CP.PropsSI("S", "P", p, "H", x * 1e3, fluid)
                         for x in h]) / 1e3

    def simple(name, h_cold, h_hot, p, fluid, Tw_cold, Tw_hot, hot):
        """Single-pressure working-fluid stream against a water glide."""
        Q, T = model._stream_curve(h_cold, h_hot, p, fluid, n=n)
        Qw, Tw = model._water_curve(Tw_cold, Tw_hot, P, Q[-1])
        h = h_cold + Q
        pack(name, s_of(h, p, fluid), T, np.interp(Q, Qw, Tw), hot)

    # ---- ORC evaporator. The cold side is TWO streams (main + reheat), so
    # the duty at a given point of the main stream includes whatever the
    # reheat stream has already given up at that temperature. Mapping the
    # main stream alone would shift the water line along Q and move the
    # apparent pinch off the bubble point.
    y, p_e, p_i = rank["y_frac"], rank["p_evape"], rank["p_5e"]
    Qc, _ = model._orc_cold_composite(rank, F)
    # The bubble point must be a NODE. The pinch sits exactly on it, so an
    # evenly spaced grid that steps over it reports a slightly loose value
    # (5.06 K instead of 5.00) and the arrow lands next to the corner rather
    # than on it.
    h_bub = CP.PropsSI("H", "P", p_e, "Q", 0, F) / 1e3
    h = np.unique(np.append(np.linspace(rank["h_10e"], rank["h_3e"], n),
                            h_bub))
    T = np.array([CP.PropsSI("T", "P", p_e, "H", x * 1e3, F) for x in h])
    h_rh = np.linspace(rank["h_5e"], rank["h_6e"], n)
    T_rh = np.array([CP.PropsSI("T", "P", p_i, "H", x * 1e3, F) for x in h_rh])
    Q = (h - rank["h_10e"]) + y * (np.interp(T, T_rh, h_rh) - rank["h_5e"])
    Qw, Tw = model._water_curve(Tdict["T_3d"], Tdict["T_2d"], P, Qc[-1])
    pack("ORC evaporator", s_of(h, p_e, F), T, np.interp(Q, Qw, Tw), False)

    simple("ORC condenser", rank["h_1e"], rank["h_4e"], rank["p_conde"], F,
           case.T_sink_C + 273.15,
           case.T_sink_C + case.DT_sink_glide + 273.15, True)
    simple("HTHP condenser", hp["h_3h"], hp["h_2h"], hp["p_condh"], R,
           Tdict["T_2c"], Tdict["T_3c"], True)
    simple("HTHP evaporator", hp["h_4h"], hp["h_13h"], hp["p_evaph"], R,
           Tdict["T_4a"], case.T_source_C + 273.15, False)
    return out


def _draw_water(ax, d, color, label=None, mark_pinch=True, lw=1.6,
                symbol=None, fontsize=8.8, dx=6, dy=-3):
    """Water glide plus its pinch. `symbol` replaces the numeric value."""
    ax.plot(d["s"], d["T_w"], color=color, lw=lw, ls=(0, (6, 2)),
            zorder=5, label=label)
    if mark_pinch:
        i = d["i_pinch"]
        ax.annotate("", xy=(d["s"][i], d["T_w"][i]),
                    xytext=(d["s"][i], d["T_wf"][i]),
                    arrowprops=dict(arrowstyle="<->", color="#1a1a1a", lw=1.0))
        txt = symbol if symbol is not None else rf"${d['pinch_K']:.1f}$ K"
        ax.annotate(txt, (d["s"][i], 0.5 * (d["T_w"][i] + d["T_wf"][i])),
                    textcoords="offset points", xytext=(dx, dy),
                    fontsize=fontsize, zorder=9, bbox=MASK)


def _isobar_bg(ax, fluid, p, T_lo, T_hi, n=240, label=None, lw=0.7,
               color="0.78", fontsize=7.6, zorder=0, label_frac=0.97):
    """A constant-pressure reference curve spanning the frame.

    Sampled in ENTHALPY between two SINGLE-PHASE anchors, so the path runs
    liquid -> dome -> vapour without a (T, p) call ever landing inside the
    two-phase region, where temperature and pressure are not independent.
    """
    h_lo = CP.PropsSI("H", "T", T_lo + 273.15, "P", p, fluid) / 1e3
    h_hi = CP.PropsSI("H", "T", T_hi + 273.15, "P", p, fluid) / 1e3
    H = np.linspace(h_lo, h_hi, n)
    T = np.array([CP.PropsSI("T", "P", p, "H", h * 1e3, fluid)
                  for h in H]) - 273.15
    s = np.array([CP.PropsSI("S", "P", p, "H", h * 1e3, fluid)
                  for h in H]) / 1e3
    ax.plot(s, T, color=color, lw=lw, zorder=zorder)
    if label:
        # Put the label on the last point that is actually INSIDE the frame.
        # Anchoring it at a fixed fraction of the curve puts it off-plot,
        # because the curve is deliberately drawn wider than the axes.
        x0, x1 = ax.get_xlim()
        y0, y1 = sorted(ax.get_ylim())
        inside = np.where((s > x0) & (s < x1) & (T > y0) & (T < y1))[0]
        if len(inside):
            k = inside[int(np.clip(label_frac, 0.0, 1.0) * (len(inside) - 1))]
            ax.annotate(label, (s[k], T[k]), fontsize=fontsize, color="0.45",
                        zorder=8, ha="right", va="bottom",
                        textcoords="offset points", xytext=(-3, 2),
                        bbox=MASK)
    return s, T


def ts_water(model, case, rank, hp, Tdict, axes=None, figsize=(11.8, 8.8)):
    """One panel per exchanger: working fluid and water on the same T-s axes.

    Four panels rather than one plot, because the four exchangers live in
    four different corners of the diagram and overlaying them all produces a
    figure in which no pinch can actually be read.
    """
    if axes is None:
        _, axes = plt.subplots(2, 2, figsize=figsize)
    ov = _water_overlay(model, case, rank, hp, Tdict)
    order = ["ORC evaporator", "ORC condenser",
             "HTHP condenser", "HTHP evaporator"]
    fluid_of = {"ORC evaporator": case.fluid, "ORC condenser": case.fluid,
                "HTHP condenser": case.refrig, "HTHP evaporator": case.refrig}
    tag = dict(zip(order, "abcd"))
    for name, ax in zip(order, np.ravel(axes)):
        d = ov[name]
        _dome(ax, fluid_of[name])
        wf_col = RED if d["hot"] else NAVY
        ax.plot(d["s"], d["T_wf"], color=wf_col, lw=1.9, zorder=4,
                label="working fluid")
        _draw_water(ax, d, "#0a9396", label="water (duty-mapped)")
        pts = list(zip(d["s"], d["T_wf"])) + list(zip(d["s"], d["T_w"]))
        _frame(ax, pts, padx=0.14, pady=0.22)
        ax.set_xlabel(r"specific entropy of the working fluid,"
                      "\n" r"$s$  /  kJ kg$^{-1}$ K$^{-1}$")
        ax.set_ylabel(r"temperature, $T$  /  $^{\circ}$C")
        ax.set_title(rf"({tag[name]}) {name}:  pinch $= "
                     rf"{d['pinch_K']:.2f}$ K", pad=8)
        ax.legend(loc="best", handlelength=1.9, labelspacing=0.3)
    return axes, ov


# ------------------------------------- ORC evaporator composite curves
def _composites(model, case, rank, Tdict):
    """(Qc, Tc, Qh, Th, Q_total) from the MODEL's own constructions.

    `_orc_cold_composite` and `_water_curve` are reused rather than
    reimplemented, so what is plotted is by construction the same pair of
    curves that `orc_pinch` minimises over. A composite curve redrawn by the
    plotting layer is a second implementation of the physics, and the two
    drift apart exactly when it matters.
    """
    Qc, Tc = model._orc_cold_composite(rank, case.fluid)
    Qh, Th = model._water_curve(Tdict["T_3d"], Tdict["T_2d"], case.P, Qc[-1])
    return Qc, Tc, Qh, Th, Qc[-1]


def orc_composite(model, case, rank, Tdict, axes=None, show_hot_end=True,
                  figsize=(11.8, 4.7)):
    """ORC evaporator composite curves and the approach profile.

    The point of the figure is that the pinch is INTERIOR. It sits at the
    onset of boiling, where the working fluid's curve goes flat against a
    water stream that keeps gliding down, and neither end of the exchanger
    shows it. A hot-end approach rule -- the thing everyone writes first --
    is blind to it by construction.

    With `show_hot_end`, the cold composite that the hot-end rule alone would
    permit (T_3e = T_2d - DT_2D_3E) is overlaid dashed. It crosses the water
    curve, so it is not a conservative approximation of the right answer; it
    is a second-law violation, and the efficiency computed from it is
    unreachable at any size of exchanger.
    """
    if axes is None:
        _, axes = plt.subplots(1, 2, figsize=figsize)
    ax, ax2 = axes
    Qc, Tc, Qh, Th, Qtot = _composites(model, case, rank, Tdict)
    pc = 100.0 * np.asarray(Qc) / Qtot
    ph = 100.0 * np.asarray(Qh) / Qtot

    # --- approach profile and the pinch location, on a common Q grid
    q = np.linspace(0.0, 100.0, 4000)
    Tw = np.interp(q, ph, np.asarray(Th) - 273.15)
    Tf = np.interp(q, pc, np.asarray(Tc) - 273.15)
    dT = Tw - Tf
    i = int(np.argmin(dT))
    q_p, dT_p = q[i], dT[i]

    # --- onset of boiling, read off the composite itself.
    #
    # NOT (h_bub - h_10e)/Qtot. That is the main stream's preheat alone, and
    # it lands 5 points of duty too low, because the composite also carries
    # the REHEAT stream (y of the flow, 5e->6e) whose whole duty lies BELOW
    # the bubble temperature and is therefore already banked by the time the
    # main stream starts boiling. Locating the break on the model's own
    # curve avoids re-deriving the assembly and getting it wrong.
    i_b = int(np.searchsorted(np.asarray(Tc), rank["T_3e"] - 1e-9))
    i_b = min(i_b, len(Qc) - 1)
    q_bub = 100.0 * Qc[i_b] / Qtot

    ax.plot(ph, np.asarray(Th) - 273.15, color=RED, lw=1.9,
            label=r"water (hot), $T_\mathrm{2d}\to T_\mathrm{3d}$")
    ax.plot(pc, np.asarray(Tc) - 273.15, color=NAVY, lw=1.9,
            label="ORC working fluid (cold)")

    cross_pct = float("nan")
    if show_hot_end:
        T_he = rank["T_3e_hot_end"]
        rk = model.double_stage_rankine(case.fluid, rank["T_1e"], T_he,
                                        case.eta_turb_s, case.eta_pump_s)
        Qc2, Tc2 = model._orc_cold_composite(rk, case.fluid)
        p2 = 100.0 * np.asarray(Qc2) / Qc2[-1]
        ax.plot(p2, np.asarray(Tc2) - 273.15,
                color="#8a8a8a", lw=1.3, ls=(0, (5, 2)), zorder=3,
                label=r"hot-end rule alone (infeasible)")
        # how much of the exchanger would actually run BACKWARDS
        cross_pct = 100.0 * float(np.mean(
            np.interp(q, p2, np.asarray(Tc2) - 273.15) > Tw))

    ax.axvline(q_bub, color="0.75", lw=0.9, ls=(0, (1, 2.5)), zorder=1)
    ax.annotate("onset of boiling", (q_bub, ax.get_ylim()[0]),
                textcoords="offset points", xytext=(4, 14), rotation=90,
                fontsize=8.2, color="0.35")
    ax.annotate("", xy=(q_p, Tw[i]), xytext=(q_p, Tf[i]),
                arrowprops=dict(arrowstyle="<->", color="#1a1a1a", lw=1.1))
    ax.annotate(rf"pinch $= {dT_p:.2f}$ K", (q_p, 0.5 * (Tw[i] + Tf[i])),
                textcoords="offset points", xytext=(8, -3), fontsize=8.8)
    ax.set_xlabel("cumulative duty  /  % of evaporator load")
    ax.set_ylabel(r"temperature, $T$  /  $^{\circ}$C")
    ax.set_title("(a) ORC evaporator composite curves", pad=8)
    ax.legend(loc="lower right", handlelength=1.9, labelspacing=0.35)

    ax2.plot(q, dT, color="#2f6f4f", lw=1.9)
    ax2.axhline(case.DT_pinch_ORC, color=RED, lw=1.0, ls=(0, (4, 2)),
                label=rf"$\Delta T_\mathrm{{pinch,ORC}} = "
                      rf"{case.DT_pinch_ORC:.0f}$ K")
    ax2.axvline(q_bub, color="0.75", lw=0.9, ls=(0, (1, 2.5)), zorder=1)
    ax2.plot([q_p], [dT_p], "o", ms=5.0, mfc="white", mec="#1a1a1a", mew=1.4,
             zorder=5)
    ax2.annotate(rf"${dT_p:.2f}$ K at ${q_p:.1f}\,\%$", (q_p, dT_p),
                 textcoords="offset points", xytext=(9, 4), fontsize=8.8)
    ax2.set_xlabel("cumulative duty  /  % of evaporator load")
    ax2.set_ylabel(r"approach, $T_\mathrm{water} - T_\mathrm{ORC}$  /  K")
    ax2.set_title("(b) the approach is worst in the interior", pad=8)
    ax2.set_ylim(0.0, max(dT.max() * 1.12, case.DT_pinch_ORC * 1.6))
    ax2.legend(loc="upper right", handlelength=1.9)
    return axes, dict(pinch_K=dT_p, pinch_at_pct=q_p, boiling_at_pct=q_bub,
                      Q_total_kJ_per_kg=Qtot,
                      hot_end_reversed_pct=cross_pct,
                      T_3e_C=rank["T_3e"] - 273.15,
                      T_3e_hot_end_C=rank["T_3e_hot_end"] - 273.15,
                      T_2d_C=Tdict["T_2d"] - 273.15)


# ------------------------------------------------------------- the table
def state_table(rank, hp):
    """State points of both cycles as one DataFrame.

    NOTE h is already kJ/kg in the state dictionaries while s is J/kg/K.
    The mixed units are the model's, not a bug here -- do not rescale h.
    """
    import pandas as pd
    rows = []
    for d, cyc, keys in ((rank, "ORC", [f"{i}e" for i in range(1, 11)]),
                         (hp, "HTHP", [f"{i}h" for i in range(1, 14)])):
        for k in keys:
            rows.append(dict(cycle=cyc, point=k,
                             T_C=d["T_" + k] - 273.15,
                             p_bar=d["p_" + k] / 1e5,
                             h_kJkg=d["h_" + k],
                             s_kJkgK=d["s_" + k] / 1e3))
    return pd.DataFrame(rows).set_index(["cycle", "point"]).round(
        {"T_C": 2, "p_bar": 4, "h_kJkg": 2, "s_kJkgK": 4})


def cycle_summary(rank, hp, res=None):
    """The scalar figures that go in the caption."""
    import pandas as pd
    d = {
        "eta_ORC": rank["rank_eff"],
        "COP": hp["hp_cop"],
        "y, turbine bleed fraction": rank["y_frac"],
        "x_8h, flash quality": hp["x_8h"],
        "eps_R1, IHX-1 effectiveness": hp["epsilon_IHX_1"],
        "eps_R2, IHX-2 effectiveness": hp["epsilon_IHX_2"],
        "ORC evaporator pinch, K": rank["pinch"],
        "ORC T_3e, C": rank["T_3e"] - 273.15,
        "p_evap ORC, bar": rank["p_evape"] / 1e5,
        "p_int  ORC, bar": rank["p_inte"] / 1e5,
        "p_cond ORC, bar": rank["p_conde"] / 1e5,
        "p_evap HTHP, bar": hp["p_evaph"] / 1e5,
        "p_int  HTHP, bar": hp["p_inth"] / 1e5,
        "p_cond HTHP, bar": hp["p_condh"] / 1e5,
    }
    if res is not None:
        for key, lab, scale in (
            ("T_2d_realised", "T_2d realised, C", None),
            ("N_wells", "N wells", 1.0),
            ("m1_ch", "m_dot per tube, charge, kg/s", 1.0),
            ("m1_dc", "m_dot per tube, discharge, kg/s", 1.0),
            ("Q_charge_kJ", "Q charge, GJ", 1e-6),
            ("Q_discharge_kJ", "Q discharge, GJ", 1e-6),
            ("W_in_kJ", "W in, GJ", 1e-6),
            ("W_out_kJ", "W out, GJ", 1e-6),
            ("eta_RTE", "eta_RTE", 1.0),
        ):
            if key in res:
                v = res[key]
                d[lab] = (v - 273.15) if scale is None else v * scale
    return pd.DataFrame({"value": d}).round(5)
