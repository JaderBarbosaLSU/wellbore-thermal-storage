"""Fig. 1 of the BOREAS paper: components and state points of the P2H2P system,
(a) charging and (b) discharging, drawn in the style of Fig. 4
(BOREAS_thums_well.pdf, well_geometry.py): serif type, the same palette, thin
steel lines, and the borehole drawn as the converted THUMS well -- fin tubes,
an insulated downcomer and a PCM cascade graded from cool (top) to warm
(bottom).

Topology follows the IHTC-18 block diagrams, with the borehole changed from a
hairpin to the downcomer design: on charge the water goes down the downcomer
and up the fin tubes; on discharge it goes down the fin tubes and up the
downcomer. Four diverter valves (VA-VD) route one pump P to either plant.
Active components are filled; idle ones are white with gray outline; idle water
branches are dashed.

Round 25 (2026-11-03, JRB sketch): a loss coil between 3c and 4c (charge) and
between 1d and 2d (discharge) stands for the thermal and pumping losses of the
lines between plant and field; they are not computed yet. The discharge points
are renumbered so that 4c/4d are the well inlets and 1c/1d the well outlets:
the former 1d (pump to VB) is now 4d, and 1d is the well outlet before the coil.
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon, FancyBboxPatch, Arc
from matplotlib.colors import to_rgb

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif'],
                     'mathtext.fontset': 'stix', 'font.size': 7.5,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})

# palette shared with well_geometry.py (Fig. 4)
C = dict(rock='#efe8dc', steel='#2b2b2b', water='#9cc3e6', ins='#c9c3dd',
         fin='#cfcfcf', loop='#1f4e79', on='#d6e6f4', on_dark='#5b7f9e',
         off_edge='#a3a3a3', off_line='#9a9a9a', src='#e07b39', sink='#9cc3e6')
PCM_STOPS = ['#7fb77e', '#e6d36d', '#c9553a']        # cool (top) -> warm (bottom)


def pcm_colors(n=9):
    """Layer colors, index 0 = top (coolest) ... n-1 = bottom (warmest)."""
    st = np.array([to_rgb(c) for c in PCM_STOPS])
    t = np.linspace(0, 1, n)
    out = []
    for x in t:
        s = x * (len(st) - 1); i = min(int(s), len(st) - 2); f = s - i
        out.append(tuple(st[i] * (1 - f) + st[i + 1] * f))
    return out


LW_PIPE, LW_LOOP, LW_EQ = 0.9, 1.25, 0.9
R_DOT = 6


class Diagram:
    def __init__(self, mode):
        self.mode = mode                                  # 'ch' or 'dc'
        self.fig = plt.figure(figsize=(6.6, 3.9))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, 2000); self.ax.set_ylim(1180, 0)
        self.ax.set_aspect('equal'); self.ax.axis('off')
        self.hops = []                                    # (x, y) crossings on horizontals

    # ------------------------------------------------------------ helpers --
    def active(self, plant):
        return (plant == 'hp' and self.mode == 'ch') or (plant == 'orc' and self.mode == 'dc')

    def pipe(self, pts, plant=None, water=False, arrow=True, dashed=False, z=3):
        on = (plant is None) or self.active(plant)
        if dashed:
            col, lw, ls = C['off_line'], 1.0, (0, (4, 3))
        elif water:
            col, lw, ls = (C['loop'] if on else C['off_line']), LW_LOOP, '-'
        else:
            col, lw, ls = (C['steel'] if on else C['off_line']), LW_PIPE, '-'
        pts = np.asarray(pts, float)
        segs = []
        for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
            if y0 == y1:                                    # horizontal: insert hops
                xs = sorted([h[0] for h in self.hops if h[1] == y0 and min(x0, x1) < h[0] < max(x0, x1)])
                if x1 < x0:
                    xs = xs[::-1]
                a = x0
                for xh in xs:
                    sgn = 1 if x1 > x0 else -1
                    segs.append(((a, y0), (xh - sgn * 13, y0)))
                    self.ax.add_patch(Arc((xh, y0), 26, 26, theta1=0, theta2=180, angle=0,
                                          color=col, lw=lw, ls=ls, zorder=z))
                    a = xh + sgn * 13
                segs.append(((a, y0), (x1, y1)))
            else:
                segs.append(((x0, y0), (x1, y1)))
        for k, ((a0, b0), (a1, b1)) in enumerate(segs):
            last = k == len(segs) - 1
            if last and arrow and not dashed:
                self.ax.annotate('', (a1, b1), (a0, b0), zorder=z,
                                 arrowprops=dict(arrowstyle='-|>,head_length=0.45,head_width=0.22',
                                                 lw=lw, color=col, shrinkA=0, shrinkB=0))
            else:
                self.ax.plot([a0, a1], [b0, b1], color=col, lw=lw, ls=ls,
                             solid_capstyle='butt', zorder=z)

    def dot(self, x, y, s, dx=10, dy=-14, plant=None, ha='left'):
        on = (plant is None) or self.active(plant)
        col = C['steel'] if on else C['off_line']
        self.ax.add_patch(Circle((x, y), R_DOT, fc=col, ec='none', zorder=8))
        self.ax.text(x + dx, y + dy, s, ha=ha, va='center', fontsize=7,
                     color=col, zorder=9)

    def box(self, x0, y0, x1, y1, label, plant, fs=7.5, ty=None):
        on = self.active(plant)
        self.ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                                    fc=C['on'] if on else 'white',
                                    ec=C['steel'] if on else C['off_edge'], lw=LW_EQ, zorder=4))
        self.ax.text((x0 + x1) / 2, (y0 + y1) / 2 if ty is None else ty, label, ha='center',
                     va='center', fontsize=fs, color=C['steel'] if on else C['off_edge'], zorder=5)

    def poly(self, pts, label, plant, lx=None, ly=None):
        on = self.active(plant)
        self.ax.add_patch(Polygon(pts, closed=True, fc=C['on'] if on else 'white',
                                  ec=C['steel'] if on else C['off_edge'], lw=LW_EQ, zorder=4))
        p = np.asarray(pts)
        self.ax.text(p[:, 0].mean() if lx is None else lx, p[:, 1].mean() if ly is None else ly,
                     label, ha='center', va='center', fontsize=7.5,
                     color=C['steel'] if on else C['off_edge'], zorder=5)

    def zigzag(self, x0, x1, y, plant, amp=9, n=4):
        xs = np.linspace(x0, x1, 2 * n + 1)
        ys = y + np.array([0] + [(-amp if k % 2 else amp) for k in range(1, 2 * n)] + [0])
        self.ax.plot(xs, ys, color='#8a8a8a', lw=0.8, zorder=6)

    def coil(self, x0, x1, y, act, n=5, r=11):
        """Loss element on a water line (heat and pumping losses between the
        plant and the field): a row of loops drawn over the pipe. Placeholder
        only -- the losses are not computed (round 25)."""
        th = np.linspace(0, 2 * np.pi * n, 400)
        v = (x1 - x0) / (2 * np.pi * n)                  # advance per radian
        xs = x0 + v * th - 0.6 * r * np.sin(th)          # prolate cycloid: loops
        ys = y - 0.8 * r * (1 - np.cos(th)) + 0.8 * r    # centred on the pipe
        col = C['loop'] if act else C['off_line']
        self.ax.plot(xs, ys, color=col, lw=1.0 if act else 0.8, zorder=9,
                     solid_capstyle='round')

    def pump(self, x, y, r, direction, label=None, plant=None, fill=True):
        on = (plant is None) or self.active(plant)
        ec = C['steel'] if on else C['off_edge']
        self.ax.add_patch(Circle((x, y), r, fc=(C['on'] if (on and fill) else 'white'),
                                 ec=ec, lw=LW_EQ, zorder=6))
        d = dict(left=(-1, 0), right=(1, 0), up=(0, -1), down=(0, 1))[direction]
        tip = (x + d[0] * r, y + d[1] * r)
        nrm = (-d[1], d[0])
        b1 = (x - d[0] * 0.15 * r + nrm[0] * r * 0.98, y - d[1] * 0.15 * r + nrm[1] * r * 0.98)
        b2 = (x - d[0] * 0.15 * r - nrm[0] * r * 0.98, y - d[1] * 0.15 * r - nrm[1] * r * 0.98)
        self.ax.plot([b1[0], tip[0], b2[0]], [b1[1], tip[1], b2[1]], color=ec, lw=0.8, zorder=7)
        if label:
            self.ax.text(x - 0.15 * r * d[0], y, label, ha='center', va='center', fontsize=6.8,
                         color=ec, zorder=8)

    def valve_x(self, x, y, r, label, plant, lx=-34):
        on = self.active(plant)
        ec = C['steel'] if on else C['off_edge']
        self.ax.add_patch(Circle((x, y), r, fc=C['on'] if on else 'white', ec=ec, lw=LW_EQ, zorder=6))
        k = r * 0.68
        self.ax.plot([x - k, x + k], [y - k, y + k], color=ec, lw=0.9, zorder=7)
        self.ax.plot([x - k, x + k], [y + k, y - k], color=ec, lw=0.9, zorder=7)
        self.ax.text(x + lx, y, label, ha='right', va='center', fontsize=7.5, color=ec)

    def diverter(self, x, y, ports_ch, ports_dc, r=22):
        """Three-way valve: draws the open path for the current mode."""
        self.ax.add_patch(Circle((x, y), r, fc='white', ec=C['loop'], lw=1.0, zorder=7))
        vec = dict(left=(-1, 0), right=(1, 0), up=(0, -1), down=(0, 1))
        for p in set(ports_ch) | set(ports_dc):
            dx, dy = vec[p]
            if p in (ports_ch if self.mode == 'ch' else ports_dc):
                self.ax.plot([x, x + dx * r], [y, y + dy * r], color=C['loop'], lw=2.0,
                             solid_capstyle='butt', zorder=8)
            else:
                self.ax.plot([x + dx * 0.55 * r, x + dx * r], [y + dy * 0.55 * r, y + dy * r],
                             color=C['off_line'], lw=1.0, zorder=8)
        self.ax.add_patch(Circle((x, y), 3.5, fc=C['loop'], ec='none', zorder=9))

    # ------------------------------------------------------------ borehole --
    def well(self, cx=1100, top=355, bot=1035):
        ax, cols = self.ax, pcm_colors(9)
        ax.add_patch(Rectangle((cx - 108, top - 15), 216, bot - top + 25, fc=C['rock'], ec='none', zorder=1))
        x0, x1 = cx - 80, cx + 80                          # casing bore
        h = (bot - 20 - top) / 9
        for k in range(9):
            ax.add_patch(Rectangle((x0, top + k * h), x1 - x0, h, fc=cols[k], ec='none', zorder=2))
        for xx in (x0, x1):
            ax.plot([xx, xx], [top - 15, bot], color=C['steel'], lw=1.6, zorder=6)
        ax.add_patch(Rectangle((x0, bot - 20), x1 - x0, 20, fc='#555555', ec='none', zorder=3))  # plug
        # downcomer at x = 1055: insulation + water core
        ax.add_patch(Rectangle((1043, top - 15), 24, bot - 32 - top + 15, fc=C['ins'], ec=C['steel'], lw=0.5, zorder=4))
        ax.add_patch(Rectangle((1050, top - 15), 10, bot - 32 - top + 15, fc=C['water'], ec='none', zorder=5))
        # fin tubes (two of n_ft drawn) joined at the top to the port at x = 1145
        for xt in (1108, 1145):
            ax.add_patch(Rectangle((xt - 8, top), 16, bot - 32 - top, fc=C['fin'], ec=C['steel'], lw=0.5, zorder=4))
            ax.add_patch(Rectangle((xt - 3.5, top), 7, bot - 32 - top, fc=C['water'], ec='none', zorder=5))
        ax.add_patch(Rectangle((1100, top - 8), 53, 9, fc=C['steel'], ec='none', zorder=6))       # top manifold
        ax.add_patch(Rectangle((1043, bot - 34), 110, 9, fc=C['steel'], ec='none', zorder=6))     # bottom manifold
        # flow arrows inside the well
        down_dc = self.mode == 'ch'
        for xa, down in ((1055, down_dc), (1108, not down_dc), (1145, not down_dc)):
            ya0, ya1 = (640, 700) if down else (700, 640)
            ax.annotate('', (xa, ya1), (xa, ya0), zorder=7,
                        arrowprops=dict(arrowstyle='-|>,head_length=0.5,head_width=0.25',
                                        lw=0.9, color=C['loop'], shrinkA=0, shrinkB=0))
        # cascade temperatures
        xl = cx - 116
        ax.text(xl, top + 0.5 * h, '$T_m^\\mathrm{bot}$', fontsize=7.5, va='center', ha='right', color='#3d3d3d')
        ax.text(xl, bot - 20 - 0.5 * h, '$T_m^\\mathrm{top}$', fontsize=7.5, va='center', ha='right', color='#3d3d3d')
        ax.annotate('', (xl - 12, bot - 20 - 1.3 * h), (xl - 12, top + 1.3 * h),
                    arrowprops=dict(arrowstyle='->', lw=0.5, color='#3d3d3d'))
        ax.text(xl - 32, (top + bot) / 2, 'melting temperature $T_{m,k}$', fontsize=6.5, va='center', ha='center', rotation=90, color='#3d3d3d')
        ax.text(cx, bot + 50, 'borehole field', ha='center', va='center', fontsize=7.5)
        ax.text(cx, bot + 85, '$N_\\mathrm{w}$ storage wells', ha='center', va='center',
                fontsize=6.5, color='#3d3d3d')

    # ---------------------------------------------------------------- HTHP --
    def hthp(self):
        hp = 'hp'
        # crossings (hops drawn on horizontals)
        self.hops += [(277, 615), (313, 615)]
        self.box(505, 210, 690, 362, 'condenser', hp)
        self.box(253, 280, 333, 430, 'R2', hp, ty=375)
        self.box(377, 280, 457, 430, 'R1', hp, ty=375)
        self.zigzag(258, 328, 318, hp); self.zigzag(382, 452, 318, hp)
        self.box(302, 848, 487, 998, 'evaporator', hp)
        self.poly([(718, 425), (836, 385), (836, 540), (718, 505)], 'C2', hp)
        self.poly([(800, 655), (918, 618), (918, 785), (800, 745)], 'C1', hp)
        # separator
        on = self.active(hp)
        ec = C['steel'] if on else C['off_edge']
        self.ax.add_patch(FancyBboxPatch((100, 567), 93, 128, boxstyle='round,pad=0,rounding_size=14',
                                         fc='white', ec=ec, lw=LW_EQ, zorder=4))
        self.ax.add_patch(FancyBboxPatch((100, 632), 93, 63, boxstyle='round,pad=0,rounding_size=14',
                                         fc=C['on_dark'] if on else '#e6e6e6', ec='none', zorder=5))
        self.ax.text(146.5, 600, 'S', ha='center', va='center', fontsize=7.5, color=ec, zorder=6)
        self.valve_x(148, 473, 32, 'V2', hp)
        self.valve_x(148, 775, 32, 'V1', hp)
        # refrigerant circuit
        P = lambda pts, **k: self.pipe(pts, plant=hp, **k)
        P([(733, 425), (733, 318), (690, 318)])                       # 2h into condenser
        P([(505, 318), (457, 318)], arrow=False); P([(377, 318), (333, 318)], arrow=False)
        P([(253, 318), (148, 318), (148, 441)], arrow=False)          # 7h
        P([(148, 505), (148, 567)])                                    # 8h
        P([(193, 615), (400, 615), (400, 430)])                        # 10h to R1
        P([(435, 430), (435, 615), (818, 615)], arrow=False)           # 11h
        P([(818, 649), (818, 615)], arrow=False)                       # 5h
        P([(818, 615), (818, 540)])                                    # 6h
        P([(148, 695), (148, 743)], arrow=False)                       # 9h
        P([(148, 807), (148, 890), (302, 890)])                        # 4h
        P([(487, 890), (551, 890), (551, 808), (277, 808), (277, 430)])  # 13h to R2
        P([(313, 430), (313, 775), (593, 775), (593, 890), (890, 890), (890, 775)])  # 1h to C1
        for args in [(733, 380, '2h', 10, 0), (481, 318, '3h', 0, 28), (355, 318, '12h', 0, -52),
                     (148, 380, '7h', 12, 0), (148, 523, '8h', -14, 0), (372, 615, '10h', -10, 24),
                     (694, 615, '11h', -14, -18), (818, 637, '5h', -50, 14), (818, 585, '6h', 12, 0),
                     (148, 717, '9h', -14, 0), (228, 890, '4h', -12, 24), (277, 700, '13h', -12, 0),
                     (778, 890, '1h', -12, -20)]:
            x, y, s, dx, dy = args
            self.dot(x, y, s, dx=dx, dy=dy, plant=hp, ha='center' if dx == 0 else ('left' if dx > 0 else 'right'))
        # source stream
        self.pump(797, 1025, 32, 'up', plant=hp)
        self.pipe([(797, 993), (797, 953), (487, 953)], plant=hp, arrow=False)
        self.pipe([(302, 953), (70, 953), (70, 1075)], plant=hp)
        self.pipe([(797, 1075), (797, 1057)], plant=hp, arrow=False)
        self.dot(735, 953, '3a', dx=-12, dy=22, plant=hp, ha='left')
        self.dot(70, 983, '4a', dx=-12, dy=0, plant=hp, ha='right')
        self.ax.plot([30, 900], [1075, 1075], color=C['loop'], lw=0.8)
        self.ax.add_patch(Polygon([(855, 1030), (899, 1030), (877, 1062)], fc=C['src'], ec=C['steel'], lw=0.6))
        for i, w in enumerate((40, 28)):
            self.ax.plot([877 - w / 2, 877 + w / 2], [1075 + 10 * (i + 1)] * 2, color=C['loop'], lw=0.6)
        self.ax.text(877, 1125, '$T_\\mathrm{source}$', ha='center', va='center', fontsize=7.5)
        self.ax.text(480, 700, 'HTHP', ha='center', va='center', fontsize=11,
                     color=C['steel'] if on else C['off_edge'])

    # ----------------------------------------------------------------- ORC --
    def orc(self):
        o = 'orc'
        on = self.active(o)
        self.hops += [(1565, 592)]
        self.box(1487, 195, 1793, 325, '', o)
        self.ax.text(1690, 268, 'evaporator', ha='center', va='center', fontsize=7.5, color=C['steel'] if on else C['off_edge'], zorder=5)
        rc = C['steel'] if on else C['off_edge']
        self.ax.plot([1527, 1527, 1546, 1565, 1565], [325, 268, 296, 268, 325], color='#8a8a8a', lw=1.2, zorder=5)
        self.poly([(1372, 370), (1492, 410), (1492, 490), (1372, 535)], 'T1', o, lx=1425)
        self.poly([(1232, 368), (1350, 408), (1350, 490), (1232, 535)], 'T2', o, lx=1285)
        self.ax.add_patch(Rectangle((1350, 442), 22, 16, fc=rc, ec='none', zorder=4))
        self.box(1273, 756, 1430, 907, 'condenser', o)
        self.box(1598, 668, 1680, 750, 'R', o, ty=735)
        self.zigzag(1603, 1675, 707, o, amp=8, n=3)
        P = lambda pts, **k: self.pipe(pts, plant=o, **k)
        P([(1487, 285), (1453, 285), (1453, 396)])                      # 3e to T1
        P([(1390, 529), (1390, 592), (1638, 592), (1638, 668)])         # 5e to R
        P([(1527, 592), (1527, 325)])                                   # 5e to reheat
        P([(1565, 325), (1565, 672), (1330, 672), (1330, 509)])         # 6e to T2
        P([(1255, 527), (1255, 800), (1273, 800)])                      # 4e
        P([(1430, 800), (1483, 800)], arrow=False)
        P([(1547, 800), (1560, 800), (1560, 707), (1598, 707)])         # 2e to R
        P([(1680, 707), (1950, 707)])                                   # 9e
        P([(1638, 750), (1638, 823), (1728, 823)], arrow=False)          # 7e
        P([(1792, 823), (1950, 823)])                                   # 8e
        P([(1950, 823), (1950, 285), (1793, 285)])                      # 10e
        self.pump(1515, 800, 32, 'right', 'P1', plant=o)
        self.pump(1760, 823, 32, 'right', 'P2', plant=o)
        for args in [(1453, 322, '3e', -12, 0), (1480, 592, '5e', -10, 22), (1330, 610, '6e', 12, 0),
                     (1255, 677, '4e', 12, 0), (1452, 800, '1e', 0, -26), (1560, 763, '2e', 12, 0),
                     (1858, 707, '9e', -12, 22), (1698, 823, '7e', -12, 22), (1860, 823, '8e', -12, 22),
                     (1950, 495, '10e', -12, 0)]:
            x, y, s, dx, dy = args
            self.dot(x, y, s, dx=dx, dy=dy, plant=o, ha='center' if dx == 0 else ('left' if dx > 0 else 'right'))
        # sink stream
        self.pump(1590, 965, 32, 'up', plant=o)
        P([(1590, 933), (1590, 858), (1430, 858)])
        P([(1273, 858), (1250, 858), (1250, 1070)])
        P([(1590, 1070), (1590, 997)], arrow=False)
        self.dot(1550, 858, '3b', dx=-12, dy=22, plant=o)
        self.dot(1250, 890, '4b', dx=-12, dy=0, plant=o, ha='right')
        self.ax.plot([1212, 1925], [1070, 1070], color=C['loop'], lw=0.8)
        self.ax.add_patch(Polygon([(1868, 1025), (1912, 1025), (1890, 1057)], fc=C['sink'], ec=C['steel'], lw=0.6))
        for i, w in enumerate((40, 28)):
            self.ax.plot([1890 - w / 2, 1890 + w / 2], [1070 + 10 * (i + 1)] * 2, color=C['loop'], lw=0.6)
        self.ax.text(1890, 1120, '$T_\\mathrm{sink}$', ha='center', va='center', fontsize=7.5)
        self.ax.text(1770, 563, 'ORC', ha='center', va='center', fontsize=11, color=rc)

    # -------------------------------------------------------- water loop --
    def loop(self):
        ch = self.mode == 'ch'
        VA, VB, VC, VD = (1055, 248), (1145, 180), (455, 180), (1145, 68)
        self.hops += [(1145, 248)]
        W = lambda pts, act, **k: self.pipe(pts, water=True, dashed=not act, **k)
        # shared: pump line, VD -> P -> VC
        W([(VD[0] - 22, 68), (770, 68)], True, arrow=False)
        W([(670, 68), (455, 68), (455, VC[1] - 22)], True)
        self.pump(720, 68, 50, 'left', 'P', fill=True)
        # charging branches
        W([(455, VC[1] + 22), (455, 245), (505, 245)], ch)                       # 2c into condenser
        W([(690, 248), (VA[0] - 22, 248)], ch)                                   # 3c to VA
        W([(VB[0], VB[1] - 22), (VD[0], VD[1] + 22)], ch)                        # fin tubes -> VD (1c)
        # discharging branches
        W([(VC[0] + 22, VC[1]), (VB[0] - 22, VB[1])], not ch)                    # 4d to VB
        W([(VA[0] + 22, 248), (1487, 248)], not ch)                              # 2d to evaporator
        W([(1793, 248), (1985, 248), (1985, 68), (VD[0] + 22, 68)], not ch)      # 3d back to VD
        # wellhead risers (always used)
        W([(VA[0], VA[1] + 22), (VA[0], 340)], True, arrow=False)
        W([(VB[0], VB[1] + 22), (VB[0], 347)], True, arrow=False)
        self.diverter(*VA, ['left', 'down'], ['right', 'down'])
        self.diverter(*VB, ['down', 'up'], ['left', 'down'])
        self.diverter(*VC, ['up', 'down'], ['up', 'right'])
        self.diverter(*VD, ['down', 'left'], ['right', 'left'])
        lab = lambda x, y, s, act, ha='left': self.dot(x, y, s, dx=10 if ha == 'left' else -10, dy=-16,
                                                         plant=None if act else 'none', ha=ha)
        self.dot(590, 68, '2c', dx=0, dy=-20, ha='center', plant=None if ch else 'x')
        self.dot(780, 248, '3c', dx=0, dy=-20, ha='center', plant=None if ch else 'x')
        self.coil(815, 935, 248, ch)                                             # losses 3c -> 4c
        self.dot(970, 248, '4c', dx=0, dy=-20, ha='center', plant=None if ch else 'x')
        self.dot(1145, 125, '1c', dx=12, dy=0, plant=None if ch else 'x')
        self.dot(880, 180, '4d', dx=0, dy=-20, ha='center', plant=None if not ch else 'x')
        self.dot(1205, 248, '1d', dx=0, dy=-20, ha='center', plant=None if not ch else 'x')
        self.coil(1235, 1340, 248, not ch)                                       # losses 1d -> 2d
        self.dot(1380, 248, '2d', dx=0, dy=-20, ha='center', plant=None if not ch else 'x')
        self.dot(1985, 160, '3d', dx=-14, dy=0, plant=None if not ch else 'x', ha='right')
        for (x, y), s in zip((VA, VB, VC, VD), ('VA', 'VB', 'VC', 'VD')):
            self.ax.text(x + 26, y - 26, s, fontsize=6.3, color=C['loop'], ha='left', va='center')

    def draw(self, path):
        self.loop(); self.hthp(); self.orc(); self.well()
        self.fig.savefig(path)
        self.fig.savefig(path.replace('.pdf', '.png'), dpi=200)


if __name__ == '__main__':
    import sys
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    Diagram('ch').draw(f'{out}/BOREAS_layout_ch.pdf')
    Diagram('dc').draw(f'{out}/BOREAS_layout_dc.pdf')
    print('ok')
