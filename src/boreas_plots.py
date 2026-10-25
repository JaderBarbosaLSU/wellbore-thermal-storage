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


def ts_combined(case, rank, hp, Tdict, ax=None, show_water=True):
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

    # ---- water-side levels, drawn first so the cycles sit on top
    if show_water:
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
    ax.set_xlabel(r"specific entropy, $s$  /  kJ kg$^{-1}$ K$^{-1}$")
    ax.set_ylabel(r"temperature, $T$  /  $^{\circ}$C")
    ax.set_title(rf"{f}:  $\eta_\mathrm{{ORC}} = {rank['rank_eff']:.4f}$,"
                 rf"  COP $= {hp['hp_cop']:.4f}$,"
                 rf"  $y = {rank['y_frac']:.4f}$,"
                 rf"  $x_\mathrm{{8h}} = {hp['x_8h']:.4f}$", pad=8)
    ax.legend(loc="upper left", handlelength=1.9, borderpad=0.25,
              labelspacing=0.35, fontsize=8.2)
    return ax


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
