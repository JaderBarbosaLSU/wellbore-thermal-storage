import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
"""Master architecture for the DBHE plant — Jader's fork structure.

PROCESS AIR
    mix(T_pc return + OA) -> T_pmc,w_pmc -> fan
      fork 1:  d_pD -> DW -> cooler
               d_pE -> [HX winter duty -> humidifier]
               changeover d_pH selects where the d_pE prong rejoins:
                   route 1 -> at T_pmo, feeding fork 2   (summer: plain bypass)
                   route 2 -> the supply duct downstream of fork 2  (winter)
      merge -> T_pmo,w_pmo
      fork 2:  d_pF -> CHW coil
               d_pG -> IEC
      merge -> T_pi,w_pi -> building

Switchable components carry keys, so topology variants come from this one
layout: draw({'DW','cooler',...}).
"""

C = dict(hot='#e34948', chw='#1baf7a', proc='#2a78d6', scav='#52514e', wct='#eb6834')
INK, INK2, MUTED, OFF = '#0b0b0b', '#52514e', '#9a998f', '#d7d6cf'
SURF = '#fcfcfb'

RET, YE, YP, YS, YH, YC = 122, 158, 250, 395, 505, 660
HI, LO = 216, 296
MX1, JX2, MX2, EJ = 706, 762, 940, 1000


class SVG:
    def __init__(s, w, h):
        s.w, s.h, s.p = w, h, []

    def box(s, x, y, w, h, lab, col, on=True, sub=None, dashed=False, fs=13):
        c, tc = (col, INK) if on else (OFF, MUTED)
        d = ' stroke-dasharray="7 5"' if dashed else ''
        s.p.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="5" fill="{SURF}" '
                   f'stroke="{c}" stroke-width="2.2"{d}/>')
        s.p.append(f'<text x="{x+w/2}" y="{y+h/2+(-5 if sub else 0)}" text-anchor="middle" '
                   f'dominant-baseline="middle" font-family="Helvetica,Arial" font-size="{fs}" '
                   f'font-weight="600" fill="{tc}">{lab}</text>')
        if sub:
            s.p.append(f'<text x="{x+w/2}" y="{y+h/2+11}" text-anchor="middle" '
                       f'font-family="Helvetica,Arial" font-size="9" fill="{MUTED}">{sub}</text>')

    def circ(s, x, y, r, lab, col, on=True):
        c = col if on else OFF
        s.p.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{SURF}" stroke="{c}" stroke-width="2.2"/>')
        s.p.append(f'<text x="{x}" y="{y+1}" text-anchor="middle" dominant-baseline="middle" '
                   f'font-family="Helvetica,Arial" font-size="10" fill="{c}">{lab}</text>')

    def line(s, pts, col, on=True, arrow=True, dashed=False, w=2.2):
        c = col if on else OFF
        d = ' '.join(('M' if i == 0 else 'L') + f'{x},{y}' for i, (x, y) in enumerate(pts))
        da = ' stroke-dasharray="6 4"' if dashed else ''
        hd = f' marker-end="url(#m{c.lstrip("#")})"' if arrow else ''
        s.p.append(f'<path d="{d}" fill="none" stroke="{c}" stroke-width="{w}"{da}{hd}/>')

    def txt(s, x, y, t, size=10, col=None, an='start', wt='400', it=False):
        st = ' font-style="italic"' if it else ''
        s.p.append(f'<text x="{x}" y="{y}" text-anchor="{an}" font-family="Helvetica,Arial" '
                   f'font-size="{size}" font-weight="{wt}" fill="{col or INK2}"{st}>{t}</text>')

    def valve(s, x, y, lab, col, on=True, dy=-11):
        c = col if on else OFF
        s.p.append(f'<circle cx="{x}" cy="{y}" r="6.5" fill="{SURF}" stroke="{c}" stroke-width="1.7"/>')
        s.p.append(f'<path d="M{x-4.6},{y-4.6} L{x+4.6},{y+4.6} M{x+4.6},{y-4.6} L{x-4.6},{y+4.6}" '
                   f'stroke="{c}" stroke-width="1.3"/>')
        s.txt(x, y + dy, lab, 9, c, 'middle', '700')

    def out(s):
        mk = ''.join(f'<marker id="m{c.lstrip("#")}" viewBox="0 0 10 10" refX="9" refY="5" '
                     f'markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">'
                     f'<path d="M0,1 L9,5 L0,9 z" fill="{c}"/></marker>'
                     for c in set(list(C.values()) + [OFF]))
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{s.w}" height="{s.h}" '
                f'viewBox="0 0 {s.w} {s.h}"><defs>{mk}</defs>'
                f'<rect width="{s.w}" height="{s.h}" fill="{SURF}"/>' + ''.join(s.p) + '</svg>')


def draw(active=None, title='Master architecture', subtitle=None, note=None, route=None):
    """route: None = show both changeover routes, 1 = rejoin at T_pmo, 2 = past fork 2."""
    on = (lambda k: True) if active is None else (lambda k: k in active)
    s = SVG(1300, 800)
    P, S, H, W, K = C['proc'], C['scav'], C['hot'], C['chw'], C['wct']
    r1, r2 = route in (None, 1), route in (None, 2)

    s.txt(24, 34, title, 19, INK, wt='700')
    if subtitle:
        s.txt(24, 55, subtitle, 11.5, INK2, it=True)
    lx = 24.0
    for lab, col in (('process air', P), ('scavenging air', S), ('hot water', H),
                     ('chilled water', W), ('cooling-tower water', K)):
        s.p.append(f'<rect x="{lx}" y="76" width="20" height="4" rx="2" fill="{col}"/>')
        s.txt(lx + 26, 83, lab, 9.5, INK2)
        lx += 33 + 6.6 * len(lab)

    # ---------------- tall shared blocks ----------------
    s.box(250, 215, 80, 215, 'DW', P if on('DW') else OFF, on('DW'))
    s.box(620, 215, 80, 215, 'cooler', P if on('cooler') else OFF, on('cooler'), fs=11.5)
    s.box(420, 365, 80, 210, 'HX', H if on('HX') else OFF, on('HX'), fs=12)
    s.box(320, 475, 84, 215, 'ABS', H if on('ABS') else OFF, on('ABS'))

    # ================= PROCESS AIR =================
    s.box(60, YP - 28, 90, 56, 'mix', P, True, 'd_pB / d_pC')
    s.circ(188, YP, 17, 'fan', P)
    s.line([(150, YP), (171, YP)], P)
    s.txt(150, YP - 36, 'T_pmc , w_pmc', 9.5, P)
    # ---- fork 1 ----
    s.line([(205, YP), (228, YP)], P, arrow=False)
    s.line([(228, YP), (228, YE), (400, YE)], P)
    s.line([(228, YP), (250, YP)], P, on('DW'))
    s.valve(238, YP + 30, 'd_pD', P, on('DW'), dy=15)
    s.valve(286, YE, 'd_pE', P, True, dy=-12)
    # prong D
    s.line([(330, YP), (620, YP)], P, on('DW'))
    s.txt(400, YP - 9, 'T_pd , w_pd', 9.5, P)
    # prong E
    s.box(400, YE - 22, 104, 44, 'HX', H if on('HXw') else OFF, on('HXw'), 'winter duty', fs=11.5)
    s.box(556, YE - 22, 104, 44, 'humidifier', W if on('HUM') else OFF, on('HUM'), fs=11)
    s.line([(504, YE), (556, YE)], P, arrow=False)
    s.line([(660, YE), (EJ, YE)], P, arrow=False)
    s.txt(392, YE + 44, 'd_pE prong: bypass in summer, heat + humidify in winter',
          9.5, P, it=True)
    s.valve(MX1, YE + 30, 'd_pH', P, True, dy=15)
    s.line([(MX1, YE), (MX1, YP - 4)], P, r1)
    s.line([(EJ, YE), (EJ, YP - 4)], P, r2)
    if route is None:
        s.txt(MX1 - 10, YE + 52, 'route 1', 9, P, 'end', it=True)
        s.txt(EJ + 9, YE - 9, 'route 2', 9, P, it=True)
    # ---- merge, fork 2 ----
    s.line([(700, YP), (JX2, YP)], P, on('DW') or r1, arrow=False)
    s.txt(706, YP + 24, 'T_pmo , w_pmo', 9.5, P)
    s.box(780, HI - 22, 108, 44, 'CHW coil', W if on('CHW') else OFF, on('CHW'), fs=11.5)
    s.box(780, LO - 22, 108, 44, 'IEC', P if on('IEC') else OFF, on('IEC'), 'M-cycle', fs=11.5)
    s.line([(JX2, YP), (JX2, HI), (780, HI)], P, on('CHW'))
    s.line([(JX2, YP), (JX2, LO), (780, LO)], P, on('IEC'))
    s.valve(JX2, HI - 38, 'd_pF', P, on('CHW'))
    s.valve(JX2, LO + 34, 'd_pG', P, on('IEC'))
    s.line([(888, HI), (MX2, HI), (MX2, LO), (888, LO)], P, arrow=False)
    s.line([(MX2, YP), (1090, YP)], P)
    s.txt(1044, YP - 10, 'T_pi , w_pi', 9.5, P, 'middle')
    s.box(1090, YP - 36, 152, 72, 'BUILDING', P, True, 'comfort setpoint')
    # return / exhaust / OA
    s.line([(1166, YP - 36), (1166, RET), (106, RET), (106, YP - 28)], P)
    s.txt(760, RET - 8, 'return   T_pc , w_pc', 9.5, P)
    s.valve(170, RET, 'd_pA', P)
    s.line([(163, RET), (52, RET)], P)
    s.txt(48, RET - 7, 'exhaust', 9.5, P, 'end')
    s.line([(40, YP + 46), (106, YP + 46), (106, YP + 28)], P)
    s.txt(36, YP + 50, 'OA', 9.5, P, 'end')
    s.txt(1166, YP + 78, 'q_s , q_l   (both signs)', 10, P, 'middle', '700')
    s.line([(1132, YP + 62), (1146, YP + 40)], P, w=1.6)
    s.line([(1200, YP + 62), (1186, YP + 40)], P, w=1.6)

    # ================= SCAVENGING AIR =================
    s.box(840, YS - 28, 100, 56, 'mix', S, on('DW'), 'd_rB / d_rC')
    s.circ(790, YS, 17, 'fan', S, on('DW'))
    s.line([(840, YS), (807, YS)], S, on('DW'))
    s.line([(773, YS), (700, YS)], S, on('DW'))
    s.line([(620, YS), (500, YS)], S, on('DW'))
    s.txt(528, YS - 9, 'T_rc', 9.5, S)
    s.line([(420, YS), (330, YS)], S, on('DW'))
    s.txt(356, YS - 9, 'T_rx', 9.5, S)
    s.line([(250, YS), (168, YS)], S, on('DW'))
    s.txt(162, YS - 9, 'exhaust  T_rd', 9.5, S, 'end')
    s.line([(1000, YS + 44), (940, YS + 44)], S, on('DW'), arrow=False)
    s.line([(940, YS + 44), (940, YS + 28)], S, on('DW'))
    s.txt(1006, YS + 48, 'OA', 9.5, S)
    s.txt(700, YS + 58, 'heat recovery — the process stream preheats the scavenging air',
          9.5, S, it=True)

    # ================= HOT WATER =================
    s.box(60, 475, 74, 96, 'DBHE', H, True)
    s.circ(166, YH, 16, 'pump', H)
    s.box(196, YH - 25, 104, 50, 'accumulator', H, on('ACC'), 'buffer', fs=11)
    s.line([(134, YH), (150, YH)], H)
    s.line([(182, YH), (196, YH)], H)
    s.valve(190, YH - 24, 'v_E', H)
    s.line([(300, YH), (320, YH)], H, on('ABS'))
    s.valve(310, YH - 24, 'v_C', H, on('ABS'))
    s.line([(404, YH), (420, YH)], H, on('HX'))
    s.valve(412, YH - 24, 'v_G', H, on('HX'))
    s.txt(524, YH + 2, 'cascade — generator first (79–100 °C), then the coil (>66 °C)',
          9.5, H, it=True)
    s.line([(478, 575), (478, 752), (95, 752), (95, 571)], H, on('HX'))
    s.txt(232, 746, 'return to well   T_wo', 9.5, H)

    # ================= CHILLED WATER + WCT =================
    s.box(150, YC - 25, 106, 50, 'WCT', K if on('ABS') else OFF, on('ABS'), 'heat rejection', fs=11)
    s.line([(320, YC - 10), (256, YC - 10)], K, on('ABS'))
    s.line([(256, YC + 12), (320, YC + 12)], K, on('ABS'))
    s.txt(178, YC + 46, 'T_ti / T_to', 9.5, K)
    s.circ(432, YC, 16, 'pump', W, on('CHW'))
    s.line([(404, YC), (416, YC)], W, on('CHW'))
    s.line([(448, YC), (516, YC)], W, on('CHW'))
    s.txt(556, YC + 4, 'chilled water  →  CHW coil   (T_cp , T_ch)', 9.5,
          W if on('CHW') else MUTED)
    s.line([(834, HI - 22), (834, HI - 42)], W, on('CHW'), arrow=False, dashed=True)
    s.txt(842, HI - 38, 'from / to ABS evaporator', 9, W if on('CHW') else MUTED, it=True)

    if note:
        s.txt(24, 788, note, 9.5, MUTED, it=True)
    return s.out()


ALL = {'DW', 'cooler', 'HX', 'HXw', 'HUM', 'CHW', 'IEC', 'ABS', 'ACC'}

if __name__ == '__main__':
    import cairosvg
    svg = draw(subtitle='all components present; the d_pE prong is drawn with both changeover routes',
               note='The winter coil on the d_pE prong is the SAME device as the regeneration '
                    'coil, re-ducted — drawn separately to avoid crossing ducts. '
                    'Chilled-water circuit shown as stubs for legibility.')
    open(f'{OUT}/master_architecture.svg', 'w').write(svg)
    cairosvg.svg2png(url=f'{OUT}/master_architecture.svg',
                     write_to=f'{OUT}/master_architecture.png', scale=1.55,
                     background_color='white')
    cairosvg.svg2pdf(url=f'{OUT}/master_architecture.svg',
                     write_to=f'{OUT}/master_architecture.pdf')
    print('written')
