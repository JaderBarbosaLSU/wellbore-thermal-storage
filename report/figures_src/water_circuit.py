"""Geothermal water circuit — corrected per Jader's 26 Jul description.

  T_wmx -> pump -> tee:  v_E -> down the annulus  (T_in)
                         v_F -> bypasses the well, rejoins the tubing return
  mix(T_ret, v_F) -> T_wr
  T_wr  -> v_A -> accumulator -> T_wa      |  v_B -> bypass
  T_wa  -> v_C -> ABS generator -> T_wg    |  generator-skip branch
  (generator-skip + v_B) -> T_wma
  T_wg + (T_wma through v_D) -> T_wmg -> v_G -> HX -> T_wx
  T_wx + (T_wma remainder through v_H) -> T_wmx -> back to the pump
"""
import os, sys
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schematic import SVG, C, INK, INK2, MUTED

H = C['hot']


def draw():
    s = SVG(1300, 700)
    s.txt(24, 34, 'Geothermal water circuit', 19, INK, wt='700')
    s.txt(24, 55, 'v_D tempers the generator outlet before the HX; v_H blends the remainder '
                  'after it', 11.5, INK2, it=True)

    YT, YB, YR = 190, 430, 596          # main chain, T_wma bypass, return header

    # ---- well ----
    s.box(70, 300, 84, 250, 'DBHE', H, True)
    for yy in (350, 392, 434, 476, 518):
        s.line([(26, yy), (66, yy)], H, w=1.6)
    s.txt(26, 332, 'q_f', 11, H, wt='700')
    s.line([(96, 250), (96, 300)], H)
    s.txt(90, 242, 'T_in', 10.5, H, 'end', wt='700')
    s.line([(134, 300), (134, YT)], H)
    s.txt(140, 286, 'T_ret', 10.5, H, wt='700')
    s.txt(60, 578, 'annulus down / tubing up', 9.5, MUTED, it=True)

    # ---- pump, v_E, v_F ----
    s.circ(240, 350, 18, 'pump', H)
    s.line([(240, YR), (240, 368)], H)
    s.line([(240, 332), (240, 250)], H, arrow=False)
    s.line([(240, 250), (96, 250)], H, arrow=False)
    s.valve(170, 250, 'v_E', H)
    s.line([(240, 250), (240, 224)], H, arrow=False)
    s.line([(240, 224), (240, YT)], H)
    s.valve(240, 224, 'v_F', H, dy=-12)
    s.txt(256, 228, 'bypasses the well', 9.5, MUTED, it=True)
    s.line([(134, YT), (300, YT)], H, arrow=False)
    s.p.append(f'<circle cx="240" cy="{YT}" r="4" fill="{H}"/>')
    s.txt(268, YT - 10, 'T_wr', 11, H, wt='700')

    # ---- accumulator fork ----
    s.line([(300, YT), (340, YT)], H)
    s.valve(322, YT, 'v_A', H)
    s.box(340, YT - 26, 116, 52, 'accumulator', H, True, 'buffer', fs=12)
    s.line([(456, YT), (516, YT)], H)
    s.txt(470, YT - 12, 'T_wa', 11, H, wt='700')
    s.p.append(f'<circle cx="300" cy="{YT}" r="4" fill="{H}"/>')
    s.line([(300, YT), (300, YB)], H, arrow=False)
    s.valve(300, 300, 'v_B', H, dy=-13)

    # ---- generator fork ----
    s.p.append(f'<circle cx="516" cy="{YT}" r="4" fill="{H}"/>')
    s.line([(516, YT), (572, YT)], H)
    s.valve(548, YT, 'v_C', H)
    s.box(572, YT - 42, 170, 132, '', H, True)
    s.txt(657, YT + 84, 'ABS', 13, INK, 'middle', wt='700')
    s.box(586, YT - 22, 142, 44, 'generator', H, True, fs=11.5)
    s.txt(657, YT + 48, 'condenser / evaporator / absorber', 9, MUTED, 'middle', it=True)
    s.line([(728, YT), (742, YT)], H, arrow=False)
    s.line([(742, YT), (800, YT)], H, arrow=False)
    s.txt(752, YT - 12, 'T_wg', 11, H, wt='700')
    # generator-skip branch joins the v_B stream
    s.line([(516, YT), (516, YB)], H, arrow=False)
    s.line([(300, YB), (860, YB)], H, arrow=False)
    s.p.append(f'<circle cx="516" cy="{YB}" r="4" fill="{H}"/>')
    s.p.append(f'<circle cx="300" cy="{YB}" r="4" fill="{H}"/>')
    s.txt(560, YB + 22, 'T_wma  =  generator-skip branch  +  v_B bypass', 10.5, H, wt='700')

    # ---- v_D temper -> T_wmg -> v_G -> HX ----
    s.p.append(f'<circle cx="800" cy="{YB}" r="4" fill="{H}"/>')
    s.line([(800, YB), (800, YT)], H)
    s.valve(800, 330, 'v_D', H, dy=-12)
    s.txt(812, 334, 'tempers the generator outlet', 9.5, MUTED, it=True)
    s.p.append(f'<circle cx="800" cy="{YT}" r="4" fill="{H}"/>')
    s.line([(800, YT), (884, YT)], H)
    s.txt(806, YT - 26, 'T_wmg', 11, H, wt='700')
    s.valve(866, YT, 'v_G', H)
    s.box(884, YT - 44, 196, 88, 'HX', H, True, fs=13)
    s.txt(982, YT + 12, 'summer: DW regeneration', 9, MUTED, 'middle')
    s.txt(982, YT + 25, 'winter: process-air heating', 9, MUTED, 'middle')

    # ---- T_wx + v_H remainder -> T_wmx ----
    s.line([(1080, YT), (1152, YT)], H, arrow=False)
    s.txt(1090, YT - 12, 'T_wx', 11, H, wt='700')
    s.valve(920, YB, 'v_H', H, dy=16)
    s.line([(860, YB), (1152, YB)], H, arrow=False)
    s.line([(1152, YB), (1152, YT + 12)], H)
    s.p.append(f'<circle cx="1152" cy="{YT}" r="4" fill="{H}"/>')
    s.txt(1164, YT - 12, 'T_wmx', 11, H, wt='700')
    s.line([(1152, YT), (1220, YT), (1220, YR)], H, arrow=False)
    s.line([(1220, YR), (240, YR)], H)
    s.txt(600, YR - 10, 'return header   →   pump suction', 9.5, H)
    s.txt(946, YB - 12, 'remainder of T_wma, blended after the HX', 9.5, MUTED, it=True)
    return s.out()


if __name__ == '__main__':
    import cairosvg
    open(f'{OUT}/water_circuit.svg', 'w').write(draw())
    cairosvg.svg2png(url=f'{OUT}/water_circuit.svg',
                     write_to=f'{OUT}/water_circuit.png', scale=1.7, background_color='white')
    cairosvg.svg2pdf(url=f'{OUT}/water_circuit.svg',
                     write_to=f'{OUT}/water_circuit.pdf')
    print('ok')
