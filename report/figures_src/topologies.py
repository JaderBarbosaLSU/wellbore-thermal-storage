"""Topology variants: one card per operating mode.

Each card shows the same skeleton with only the active path highlighted:
  junction -> fan -> fork1 { HXw -> humidifier | DW -> cooler } -> merge
           -> fork2 { CHW | IEC } -> building
"""
import os, sys
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schematic import SVG, C, INK, INK2, MUTED, OFF, SURF

P, S, H, W, K = C['proc'], C['scav'], C['hot'], C['chw'], C['wct']

MODES = [
 dict(id='M1', quad='cooling + dehumidification', tier='mild',
      on={'DW', 'cooler', 'HX'}, water='≥ 66 °C',
      note='wheel dries; recuperative cooler trims the sensible load'),
 dict(id='M2', quad='cooling + dehumidification', tier='moderate',
      on={'DW', 'cooler', 'HX', 'IEC'}, water='≥ 66 °C',
      note='drying first collapses the wet-bulb, so the IEC bites deep'),
 dict(id='M3', quad='cooling + dehumidification', tier='harsh',
      on={'DW', 'cooler', 'HX', 'ABS', 'CHW'}, water='≥ 79 °C',
      note='chilled water for the sensible peak; wheel still carries the latent'),
 dict(id='M4', quad='cooling + humidification', tier='mild',
      on={'HUM'}, water='none',
      note='adiabatic humidification IS direct evaporative cooling'),
 dict(id='M5', quad='cooling + humidification', tier='harsh',
      on={'ABS', 'CHW', 'HUM'}, water='≥ 79 °C',
      note='dry CHW coil above dew point, then humidify — ORDER PROBLEM, see note'),
 dict(id='M6', quad='heating + dehumidification', tier='—',
      on={'DW', 'cooler', 'HX', 'HXw'}, water='≥ 66 °C',
      note='needs MODULATING d_pD / d_pE: dry one part, heat the other, blend'),
 dict(id='M7', quad='heating + humidification', tier='mild',
      on={'HXw'}, water='≥ supply T',
      note='heating only; latent deficit small enough to ignore'),
 dict(id='M8', quad='heating + humidification', tier='full',
      on={'HXw', 'HUM'}, water='≥ supply T',
      note='the winter workhorse — ~430 kW latent in January'),
]


def card(s, x0, y0, m):
    on = lambda k: k in m['on']
    w, h = 690, 246
    s.p.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" rx="8" fill="{SURF}" '
               f'stroke="#e3e2db" stroke-width="1.6"/>')
    s.txt(x0 + 16, y0 + 26, m['id'], 15, INK, wt='700')
    s.txt(x0 + 46, y0 + 26, m['quad'], 12, INK, wt='600')
    s.txt(x0 + 46, y0 + 42, f"tier: {m['tier']}   ·   geothermal water: {m['water']}", 9.5, INK2)

    yT, yM, yB = y0 + 78, y0 + 118, y0 + 158        # prong E, spine, prong D
    xf, x1, x2, x3, x4 = x0 + 40, x0 + 96, x0 + 300, x0 + 420, x0 + 600

    s.p.append(f'<circle cx="{xf-14}" cy="{yM}" r="3.5" fill="{P}"/>')
    s.circ(xf, yM, 13, 'fan', P)
    s.line([(xf + 13, yM), (x1 - 14, yM)], P, arrow=False)
    # fork 1
    s.line([(x1 - 14, yM), (x1 - 14, yT), (x1, yT)], P, on('HXw') or on('HUM'), arrow=False)
    s.line([(x1 - 14, yM), (x1 - 14, yB), (x1, yB)], P, on('DW'), arrow=False)
    s.box(x1, yT - 15, 88, 30, 'HX winter', H if on('HXw') else OFF, on('HXw'), fs=10)
    s.box(x1 + 100, yT - 15, 92, 30, 'humidifier', W if on('HUM') else OFF, on('HUM'), fs=10)
    s.box(x1, yB - 15, 88, 30, 'DW', P if on('DW') else OFF, on('DW'), fs=10)
    s.box(x1 + 100, yB - 15, 92, 30, 'cooler', P if on('cooler') else OFF, on('cooler'), fs=10)
    s.line([(x1 + 88, yT), (x1 + 100, yT)], P, on('HUM'), arrow=False)
    s.line([(x1 + 88, yB), (x1 + 100, yB)], P, on('cooler'), arrow=False)
    # merge -> fork 2
    s.line([(x1 + 192, yT), (x2, yT), (x2, yM)], P, on('HXw') or on('HUM'), arrow=False)
    s.line([(x1 + 192, yB), (x2, yB), (x2, yM)], P, on('DW'), arrow=False)
    s.p.append(f'<circle cx="{x2}" cy="{yM}" r="3.5" fill="{P}"/>')
    s.line([(x2, yM), (x3 - 14, yM)], P, arrow=False)
    s.txt(x2 + 4, yM - 8, 'T_pmo', 8.5, P)
    s.line([(x3 - 14, yM), (x3 - 14, yT), (x3, yT)], P, on('CHW'), arrow=False)
    s.line([(x3 - 14, yM), (x3 - 14, yB), (x3, yB)], P, on('IEC'), arrow=False)
    s.box(x3, yT - 15, 84, 30, 'CHW', W if on('CHW') else OFF, on('CHW'), fs=10)
    s.box(x3, yB - 15, 84, 30, 'IEC', P if on('IEC') else OFF, on('IEC'), fs=10)
    s.line([(x3 + 84, yT), (x3 + 110, yT), (x3 + 110, yM)], P, on('CHW'), arrow=False)
    s.line([(x3 + 84, yB), (x3 + 110, yB), (x3 + 110, yM)], P, on('IEC'), arrow=False)
    s.line([(x3 + 110, yM), (x4, yM)], P)
    s.box(x4, yM - 19, 74, 38, 'BLDG', P, True, fs=10)
    # water-side subsystems actually needed
    subs = []
    if on('HX') or on('HXw'):
        subs.append(('HX / hot water', H))
    if on('ABS'):
        subs.append(('ABS + WCT', K))
    if on('CHW'):
        subs.append(('chilled water', W))
    if on('DW'):
        subs.append(('scavenging loop', S))
    lx = x0 + 18
    for lab, col in subs:
        s.p.append(f'<rect x="{lx}" y="{y0+196}" width="12" height="4" rx="2" fill="{col}"/>')
        s.txt(lx + 17, y0 + 202, lab, 8.5, INK2)
        lx += 26 + 5.4 * len(lab)
    s.txt(x0 + 18, y0 + 228, m['note'], 9, MUTED, it=True)


def draw():
    s = SVG(1440, 1096)
    s.txt(24, 32, 'Topology variants — one per load quadrant and intensity tier', 18, INK, wt='700')
    s.txt(24, 52, 'quadrants named by the processes they require; the tier says which components '
                  'can deliver them', 11, INK2, it=True)
    for i, m in enumerate(MODES):
        card(s, 24 + (i % 2) * 700, 74 + (i // 2) * 254, m)
    s.txt(24, 1086, 'Greyed components are present in the master architecture but inactive in '
                    'that mode. Fork-1 dampers are binary except in M6.', 9, MUTED, it=True)
    return s.out()


if __name__ == '__main__':
    import cairosvg
    open(f'{OUT}/topology_variants.svg', 'w').write(draw())
    cairosvg.svg2png(url=f'{OUT}/topology_variants.svg',
                     write_to=f'{OUT}/topology_variants.png', scale=1.5,
                     background_color='white')
    cairosvg.svg2pdf(url=f'{OUT}/topology_variants.svg',
                     write_to=f'{OUT}/topology_variants.pdf')
    print('ok')
