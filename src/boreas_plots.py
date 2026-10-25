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


def use_paper_style():
    """Serif, hairline axes, no top/right spine. Call once per session."""
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": SERIF,
        "mathtext.fontset": "stix",
        "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 11,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5, "legend.fontsize": 8.6,
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
            color="0.45", lw=1.0, zorder=1)


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
                fontsize=8.4, color="#1a1a1a", zorder=7)


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
_CMB_ORC = {"1e": (7, -12), "2e": (7, 2), "3e": (-17, 7), "4e": (7, -2),
            "5e": (-19, -12), "6e": (6, 5)}
_CMB_HP = {"1h": (7, -11), "2h": (7, 2), "3h": (-17, 7), "4h": (-18, -12),
           "7h": (-21, -12), "8h": (-6, -13), "9h": (-22, 2),
           "10h": (-24, -6), "11h": (7, 0), "12h": (7, -11), "13h": (-21, 3)}


def ts_combined(case, rank, hp, Tdict, ax=None, show_water=True,
                water="level", model=None):
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
        style = {"ORC evaporator": ("#1a1a1a", r"ORC evaporator water, "
                                               r"$T_\mathrm{2d}\to T_\mathrm{3d}$"),
                 "ORC condenser": ("#2a6fb0", "sink water, ORC condenser"),
                 "HTHP condenser": ("#0a9396", r"charge water, "
                                               r"$T_\mathrm{2c}\to T_\mathrm{3c}$"),
                 "HTHP evaporator": ("#9c27a0", "source water, HTHP evaporator")}
        for nm, (col, lab) in style.items():
            d = ov[nm]
            _draw_water(ax, d, col, label=f"{lab}  ({d['pinch_K']:.1f} K)")
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
    xl = r"specific entropy, $s$  /  kJ kg$^{-1}$ K$^{-1}$"
    if show_water and water == "glide":
        xl += ("\n" r"(water lines are duty-mapped onto the working fluid's "
               r"$s$; the vertical gap is the local approach)")
    ax.set_xlabel(xl)
    ax.set_ylabel(r"temperature, $T$  /  $^{\circ}$C")
    ax.set_title(rf"{f}:  $\eta_\mathrm{{ORC}} = {rank['rank_eff']:.4f}$,"
                 rf"  COP $= {hp['hp_cop']:.4f}$,"
                 rf"  $y = {rank['y_frac']:.4f}$,"
                 rf"  $x_\mathrm{{8h}} = {hp['x_8h']:.4f}$", pad=8)
    ax.legend(loc="upper left", handlelength=1.9, borderpad=0.25,
              labelspacing=0.35, fontsize=8.2)
    return ax


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


def _draw_water(ax, d, color, label=None, mark_pinch=True, lw=1.6):
    ax.plot(d["s"], d["T_w"], color=color, lw=lw, ls=(0, (6, 2)),
            zorder=5, label=label)
    if mark_pinch:
        i = d["i_pinch"]
        ax.annotate("", xy=(d["s"][i], d["T_w"][i]),
                    xytext=(d["s"][i], d["T_wf"][i]),
                    arrowprops=dict(arrowstyle="<->", color="#1a1a1a", lw=1.0))
        ax.annotate(rf"${d['pinch_K']:.1f}$ K",
                    (d["s"][i], 0.5 * (d["T_w"][i] + d["T_wf"][i])),
                    textcoords="offset points", xytext=(6, -3), fontsize=8.0,
                    zorder=8)


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
