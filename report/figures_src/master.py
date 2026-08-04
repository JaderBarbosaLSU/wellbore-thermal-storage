"""Master architecture — DBHE plant, all components present.

Process air : mix -> fan -> fork1 { d_pD: DW -> cooler | d_pE: HX(winter) -> humidifier }
              the d_pE prong ends in two dampers: d_pH rejoins at T_pmo,
              d_pI joins the supply duct downstream of fork 2
              -> T_pmo -> fork2 { d_pF: CHW | d_pG: IEC } -> T_pi -> building
Scavenging  : mix(OA + building exhaust) -> fan -> cooler -> HX -> DW -> exhaust
Hot water   : see water_circuit.py for the same chain drawn larger
"""
import os, sys
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schematic import SVG, C, INK, INK2, MUTED, OFF

P, S, H, W, K = C['proc'], C['scav'], C['hot'], C['chw'], C['wct']
RET, YE, YP, YS, YW, YWB, YWR, YC = 112, 158, 250, 400, 520, 620, 700, 800
HI, LO = 212, 296
TPMO, JX2, MX2, R2 = 985, 1010, 1200, 1230
ALL = {'DW', 'cooler', 'HX', 'HXw', 'HUM', 'CHW', 'IEC', 'ABS', 'ACC'}


def draw(active=None, title='Master architecture', subtitle=None, note=None):
    on = (lambda k: True) if active is None else (lambda k: k in active)
    s = SVG(1490, 880)
    s.txt(24, 34, title, 19, INK, wt='700')
    if subtitle:
        s.txt(24, 55, subtitle, 11.5, INK2, it=True)
    lx = 24.0
    for lab, col in (('process air', P), ('scavenging air', S), ('hot water', H),
                     ('chilled water', W), ('cooling-tower water', K)):
        s.p.append(f'<rect x="{lx}" y="76" width="20" height="4" rx="2" fill="{col}"/>')
        s.txt(lx + 26, 83, lab, 9.5, INK2)
        lx += 33 + 6.6 * len(lab)

    # ---------- tall shared blocks ----------
    s.box(250, 214, 80, 214, 'DW', P if on('DW') else OFF, on('DW'))
    s.box(880, 214, 80, 214, 'cooler', P if on('cooler') else OFF, on('cooler'), fs=11.5)
    s.box(700, 372, 80, 188, 'HX', H if on('HX') else OFF, on('HX'), fs=12)

    # ================= PROCESS AIR =================
    s.circ(186, YP, 17, 'fan', P)
    s.p.append(f'<circle cx="120" cy="{YP}" r="4" fill="{P}"/>')
    s.line([(120, YP), (169, YP)], P)
    s.txt(118, YP - 26, 'T_pmc , w_pmc', 9.5, P, wt='700')
    s.line([(203, YP), (228, YP)], P, arrow=False)
    s.line([(228, YP), (228, YE), (400, YE)], P)
    s.line([(228, YP), (250, YP)], P, on('DW'))
    s.valve(240, YP + 30, 'd_pD', P, on('DW'), dy=15)
    s.valve(300, YE, 'd_pE', P, True)
    s.line([(330, YP), (880, YP)], P, on('DW'))
    s.txt(430, YP - 9, 'T_pd , w_pd', 9.5, P)
    # prong E
    s.box(400, YE - 22, 110, 44, 'HX', H if on('HXw') else OFF, on('HXw'), 'winter duty', fs=11.5)
    s.box(560, YE - 22, 110, 44, 'humidifier', W if on('HUM') else OFF, on('HUM'), fs=11)
    s.line([(510, YE), (560, YE)], P, arrow=False)
    s.line([(670, YE), (R2, YE)], P, arrow=False)
    s.txt(400, YE - 32, 'd_pE prong — plain bypass in summer, heating + humidification in winter',
          9.5, P, it=True)
    s.p.append(f'<circle cx="{TPMO}" cy="{YE}" r="4" fill="{P}"/>')
    s.line([(TPMO, YE), (TPMO, YP - 4)], P)
    s.valve(TPMO, YE + 46, 'd_pH', P, True, dy=15)
    s.valve(1090, YE, 'd_pI', P, True)
    s.line([(R2, YE), (R2, YP - 4)], P)
    # fork 2
    s.line([(960, YP), (JX2, YP)], P, arrow=False)
    s.txt(898, YP + 44, 'T_pmo , w_pmo', 9.5, P)
    s.box(1050, HI - 22, 112, 44, 'CHW coil', W if on('CHW') else OFF, on('CHW'), fs=11.5)
    s.box(1050, LO - 22, 112, 44, 'IEC', P if on('IEC') else OFF, on('IEC'), '1ary / 2ary', fs=11.5)
    s.line([(JX2, YP), (JX2, HI), (1050, HI)], P, on('CHW'))
    s.line([(JX2, YP), (JX2, LO), (1050, LO)], P, on('IEC'))
    s.valve(JX2, HI - 38, 'd_pF', P, on('CHW'))
    s.valve(JX2, LO + 34, 'd_pG', P, on('IEC'))
    s.line([(1162, HI), (MX2, HI), (MX2, LO), (1162, LO)], P, arrow=False)
    s.line([(1108, LO + 22), (1108, LO + 52)], P, on('IEC'), arrow=True)
    s.txt(1116, LO + 52, '2ary discharge  T_iec , w_iec  (mass leaves the supply)', 9,
          P if on('IEC') else MUTED, it=True)
    s.line([(MX2, YP), (1260, YP)], P)
    s.txt(1256, YP + 22, 'T_pi , w_pi', 9.5, P, 'end')
    s.box(1260, YP - 36, 170, 72, 'BUILDING', P, True, 'comfort setpoint')
    # return / exhaust / OA
    s.line([(1345, YP - 36), (1345, RET), (120, RET)], P, arrow=False)
    s.txt(700, RET - 8, 'return   T_pc , w_pc', 9.5, P)
    s.p.append(f'<circle cx="120" cy="{RET}" r="4" fill="{P}"/>')
    # exhaust branch
    s.line([(120, RET), (78, RET)], P)
    s.valve(100, RET, 'd_pA', P)
    s.txt(74, RET - 7, 'exhaust', 9.5, P, 'end')
    s.txt(74, RET + 14, 'T_ea , RH_ea', 9, P, 'end')
    # recirculation branch
    s.line([(120, RET), (120, YP)], P)
    s.valve(120, 180, 'd_pB', P, dy=-13)
    s.txt(132, 196, 'recirculation', 9, P, it=True)
    # outdoor-air branch
    s.line([(74, YP + 58), (120, YP + 58)], P, arrow=False)
    s.line([(120, YP + 58), (120, YP)], P)
    s.valve(98, YP + 58, 'd_pC', P, dy=15)
    s.txt(70, YP + 62, 'OA', 9.5, P, 'end')
    s.txt(1345, YP + 74, 'q_s , q_l   (both signs)', 10, P, 'middle', '700')
    s.line([(1312, YP + 58), (1326, YP + 40)], P, w=1.6)
    s.line([(1378, YP + 58), (1364, YP + 40)], P, w=1.6)

    # ================= SCAVENGING AIR =================
    s.circ(1032, YS, 17, 'fan', S, on('DW'))
    s.p.append(f'<circle cx="1150" cy="{YS}" r="4" fill="{S}"/>')
    s.line([(1150, YS), (1049, YS)], S, on('DW'))
    s.txt(1060, YS - 12, 'T_rm , w_rm', 9.5, S, wt='700')
    s.line([(1015, YS), (960, YS)], S, on('DW'))
    s.line([(880, YS), (780, YS)], S, on('DW'))
    s.txt(800, YS - 9, 'T_rc', 9.5, S)
    s.line([(700, YS), (330, YS)], S, on('DW'))
    s.txt(500, YS - 9, 'T_rx', 9.5, S)
    s.line([(250, YS), (160, YS)], S, on('DW'))
    s.txt(154, YS - 9, 'exhaust  T_rd', 9.5, S, 'end')
    # building-exhaust branch (d_rA) and outdoor-air branch (d_rB)
    s.p.append(f'<circle cx="1345" cy="{RET}" r="4" fill="{P}"/>')
    s.line([(1345, RET), (1462, RET), (1462, YS), (1150, YS)], P, on('DW'))
    s.valve(1300, YS, 'd_rA', S, on('DW'))
    s.txt(1456, YS - 44, 'building exhaust', 9.5, P, 'end')
    s.txt(1456, YS - 32, 'T_ea , RH_ea', 9.5, P, 'end')
    s.line([(1262, YS + 70), (1150, YS + 70)], S, on('DW'), arrow=False)
    s.line([(1150, YS + 70), (1150, YS)], S, on('DW'))
    s.valve(1206, YS + 70, 'd_rB', S, on('DW'), dy=15)
    s.txt(1268, YS + 74, 'OA', 9.5, S)
    s.txt(560, YS + 30, 'heat recovery — the process stream preheats the scavenging air',
          9.5, S, it=True)

    # ================= HOT WATER =================
    s.box(46, 548, 68, 142, 'DBHE', H, True)
    s.circ(150, 600, 17, 'pump', H)
    s.line([(150, YWR), (150, 617)], H)
    s.line([(150, 583), (150, 496)], H, arrow=False)
    s.line([(150, 496), (66, 496)], H, arrow=False)
    s.valve(112, 496, 'v_E', H)
    s.line([(66, 496), (66, 548)], H)
    s.txt(60, 488, 'T_in', 9.5, H, 'end', wt='700')
    s.line([(96, 548), (96, YW)], H)
    s.txt(102, 544, 'T_ret', 9.5, H, wt='700')
    s.line([(150, 496), (150, 474)], H, arrow=False)
    s.line([(150, 474), (150, YW)], H)
    s.valve(150, 474, 'v_F', H, dy=-12)
    s.line([(96, YW), (240, YW)], H, arrow=False)
    s.p.append(f'<circle cx="150" cy="{YW}" r="4" fill="{H}"/>')
    s.txt(176, YW - 10, 'T_wr', 10, H, wt='700')
    # accumulator
    s.p.append(f'<circle cx="240" cy="{YW}" r="4" fill="{H}"/>')
    s.line([(240, YW), (290, YW)], H, on('ACC'))
    s.valve(268, YW, 'v_A', H, on('ACC'))
    s.box(290, YW - 25, 110, 50, 'accumulator', H, on('ACC'), 'buffer', fs=11)
    s.line([(400, YW), (450, YW)], H)
    s.txt(410, YW - 12, 'T_wa', 10, H, wt='700')
    s.line([(240, YW), (240, YWB)], H, arrow=False)
    s.valve(240, 570, 'v_B', H, dy=-12)
    # generator
    s.p.append(f'<circle cx="450" cy="{YW}" r="4" fill="{H}"/>')
    s.line([(450, YW), (486, YW)], H, on('ABS'))
    s.valve(470, YW, 'v_C', H, on('ABS'))
    s.box(486, YW - 26, 118, 52, 'ABS gen', H if on('ABS') else OFF, on('ABS'), fs=11.5)
    s.line([(604, YW), (645, YW)], H, on('ABS'), arrow=False)
    s.txt(610, YW - 12, 'T_wg', 10, H, wt='700')
    s.line([(450, YW), (450, YWB)], H, arrow=False)
    s.line([(240, YWB), (830, YWB)], H, arrow=False)
    s.p.append(f'<circle cx="450" cy="{YWB}" r="4" fill="{H}"/>')
    s.txt(300, YWB + 20, 'T_wma  =  generator-skip branch  +  v_B bypass', 9.5, H, wt='700')
    # v_D temper -> T_wmg -> v_G -> HX
    s.p.append(f'<circle cx="645" cy="{YWB}" r="4" fill="{H}"/>')
    s.line([(645, YWB), (645, YW)], H)
    s.valve(645, 578, 'v_D', H, dy=-12)
    s.p.append(f'<circle cx="645" cy="{YW}" r="4" fill="{H}"/>')
    s.line([(645, YW), (700, YW)], H, on('HX'))
    s.txt(652, YW - 12, 'T_wmg', 10, H, wt='700')
    s.valve(676, YW + 24, 'v_G', H, on('HX'), dy=15)
    s.line([(780, YW), (830, YW)], H, on('HX'), arrow=False)
    s.txt(788, YW - 12, 'T_wx', 10, H, wt='700')
    s.valve(830, YWB, 'v_H', H, dy=16)
    s.line([(830, YWB), (830, YW + 10)], H)
    s.p.append(f'<circle cx="830" cy="{YW}" r="4" fill="{H}"/>')
    s.line([(830, YW), (890, YW), (890, YWR)], H, arrow=False)
    s.txt(846, YW - 12, 'T_wmx', 10, H, wt='700')
    s.line([(890, YWR), (150, YWR)], H)
    s.txt(430, YWR - 10, 'return header   →   pump suction', 9.5, H)

    # ================= CHILLED WATER + WCT (stubs) =================
    s.line([(520, YW + 26), (520, YW + 74)], K, on('ABS'), arrow=False, dashed=True)
    s.txt(470, YW + 88, 'condenser + absorber  →  WCT', 9, K if on('ABS') else MUTED, it=True)
    s.box(300, YC - 25, 116, 50, 'WCT', K if on('ABS') else OFF, on('ABS'), 'heat rejection', fs=11)
    s.txt(424, YC - 4, 'T_ti / T_to  —  from the ABS condenser and absorber', 9,
          K if on('ABS') else MUTED, it=True)
    s.circ(880, YC, 16, 'pump', W, on('CHW'))
    s.line([(910, YC), (980, YC)], W, on('CHW'))
    s.txt(990, YC + 4, 'chilled water  →  CHW coil   (T_cp , T_ch)', 9.5,
          W if on('CHW') else MUTED)
    s.line([(810, YC), (864, YC)], W, on('CHW'), arrow=False)
    s.txt(806, YC + 4, 'ABS evaporator', 9.5, W if on('CHW') else MUTED, 'end', it=True)
    s.line([(1106, HI - 22), (1106, HI - 42)], W, on('CHW'), arrow=False, dashed=True)
    s.txt(1114, HI - 38, 'from / to ABS evaporator', 9, W if on('CHW') else MUTED, it=True)

    if note:
        s.txt(24, 866, note, 9.5, MUTED, it=True)
    return s.out()


if __name__ == '__main__':
    import cairosvg
    svg = draw(subtitle='junctions rather than mixing boxes; d_pA/d_pB/d_pC and d_rA/d_rB on the branches '
                        'themselves',
               note='Water circuit drawn in reduced form — see the dedicated geothermal-water '
                    'figure for the same chain at full size. ABS sub-circuits and the '
                    'chilled-water loop are shown as labelled stubs for legibility.')
    open(f'{OUT}/master_architecture.svg', 'w').write(svg)
    cairosvg.svg2png(url=f'{OUT}/master_architecture.svg',
                     write_to=f'{OUT}/master_architecture.png', scale=1.45,
                     background_color='white')
    cairosvg.svg2pdf(url=f'{OUT}/master_architecture.svg',
                     write_to=f'{OUT}/master_architecture.pdf')
    print('ok')
