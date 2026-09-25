import sys; sys.path.insert(0,'/tmp/build')
import numpy as np, CoolProp.CoolProp as CP
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import thums_model as T

F = 'cyclopentane'
def Ts(h_kJ, p):  return CP.PropsSI('T','H',h_kJ*1000,'P',p,F), CP.PropsSI('S','H',h_kJ*1000,'P',p,F)/1000.
def isobar(h1, h2, p, n=120):
    hs = np.linspace(h1, h2, n); out = [Ts(h,p) for h in hs]
    return np.array([o[1] for o in out]), np.array([o[0]-273.15 for o in out])
def isentrope(s, p1, p2, n=60):
    ps = np.linspace(p1, p2, n)
    Tt = [CP.PropsSI('T','S',s,'P',p,F)-273.15 for p in ps]
    return np.full(n, s/1000.), np.array(Tt)

def dome(ax):
    Tt = np.linspace(CP.PropsSI('Ttriple',F)+0.5, CP.PropsSI('Tcrit',F)-0.01, 400)
    sl = [CP.PropsSI('S','T',t,'Q',0,F)/1000. for t in Tt]
    sv = [CP.PropsSI('S','T',t,'Q',1,F)/1000. for t in Tt]
    ax.plot(sl, Tt-273.15, color='0.55', lw=1.1, zorder=1)
    ax.plot(sv, Tt-273.15, color='0.55', lw=1.1, zorder=1)
    ax.plot([CP.PropsSI('S','T',CP.PropsSI('Tcrit',F)-0.01,'Q',0,F)/1000.],
            [CP.PropsSI('Tcrit',F)-273.15], 'o', ms=3, color='0.4', zorder=2)

def label(ax, s, t, txt, dx=0.012, dy=3.5, **kw):
    ax.plot([s],[t],'o',ms=4.2,color=kw.get('c','k'),zorder=6)
    ax.annotate(txt,(s,t),textcoords='offset points',xytext=(dx*380,dy),
                fontsize=8.2,zorder=7,color=kw.get('c','k'))

c = T.Case(); rank, hp, Tst = T.cycle_state_points(c)

# =========================== ORC =========================================
fig, ax = plt.subplots(figsize=(7.4, 5.4))
dome(ax)
pe, pi, pc = rank['p_evape'], rank['p_inte'], rank['p_conde']
BLU, ORG, GRN = '#1f5f9e', '#d2691e', '#2e8b57'

s,t = isobar(rank['h_10e'], rank['h_3e'], pe);      ax.plot(s,t,color=BLU,lw=2.4,label='evaporator / reheater (borehole water)')
s,t = isentrope(rank['s_3e'], pe, pi);              ax.plot(s,t,color=ORG,lw=2.0,label='turbines (isentropic)')
s,t = isobar(rank['h_5e'], rank['h_6e'], pi);       ax.plot(s,t,color=BLU,lw=2.4)
s,t = isentrope(rank['s_6e'], pi, pc);              ax.plot(s,t,color=ORG,lw=2.0)
s,t = isobar(rank['h_4e'], rank['h_1e'], pc);       ax.plot(s,t,color=GRN,lw=2.4,label='condenser (sink)')
s,t = isobar(rank['h_7e'], rank['h_9e'], pi);       ax.plot(s,t,color='0.35',lw=1.4,ls='--',label='pumps / regenerator')
ax.plot([rank['s_1e']/1000, rank['s_2e']/1000],
        [rank['T_1e']-273.15, rank['T_2e']-273.15], color='0.35', lw=1.4, ls='--')

for k,(dx,dy) in {'1e':(-.03,-13),'2e':(.004,5),'3e':(.006,5),'4e':(.006,4),
                  '5e':(-.028,6),'6e':(.006,4),'7e':(-.030,-12),'10e':(.004,6)}.items():
    label(ax, rank['s_'+k]/1000, rank['T_'+k]-273.15, k, dx, dy)

# the isothermal plateau that causes the pinch
h_bub = CP.PropsSI('H','P',pe,'Q',0,F)/1000.
s_b = CP.PropsSI('S','P',pe,'Q',0,F)/1000.; s_d = rank['s_3e']/1000.
ax.annotate('', xy=(s_d, rank['T_3e']-273.15+6), xytext=(s_b, rank['T_3e']-273.15+6),
            arrowprops=dict(arrowstyle='<->', color='crimson', lw=1.3))
ax.text(0.5*(s_b+s_d), rank['T_3e']-273.15+9,
        '76%% of the heat input, at ONE temperature\n$T_{3e}=%.1f$ °C  — this is what pinches'%(rank['T_3e']-273.15),
        ha='center', va='bottom', fontsize=8.4, color='crimson')
ax.axhline(Tst['T_2d']-273.15, color='steelblue', ls=':', lw=1.2)
ax.axhline(Tst['T_3d']-273.15, color='steelblue', ls=':', lw=1.2)
ax.text(ax.get_xlim()[1], Tst['T_2d']-273.15, ' water in %.0f °C'%(Tst['T_2d']-273.15),
        fontsize=7.8, color='steelblue', va='bottom', ha='right')
ax.text(ax.get_xlim()[1], Tst['T_3d']-273.15, ' water out %.0f °C'%(Tst['T_3d']-273.15),
        fontsize=7.8, color='steelblue', va='top', ha='right')
ax.set_xlabel('specific entropy  $s$  [kJ kg$^{-1}$ K$^{-1}$]')
ax.set_ylabel('temperature  $T$  [°C]')
ax.set_title('ORC — double-stage with reheat and regeneration, cyclopentane\n'
             r'$\eta_{ORC}=%.4f$,  $y=%.3f$,  v0.9 (pinch-feasible)'%(rank['rank_eff'],rank['y_frac']),
             fontsize=10)
ax.set_xlim(-0.45, 1.75); ax.set_ylim(0, 250)
ax.grid(alpha=.25); ax.legend(fontsize=8, loc='upper left', framealpha=.92)
fig.tight_layout(); fig.savefig('/tmp/ts_orc.png', dpi=175); fig.savefig('/tmp/ts_orc.pdf')
print('ORC figure written')

# =========================== HTHP ========================================
def isenthalp(h_kJ, p1, p2, n=80):
    ps = np.geomspace(p1, p2, n)
    out = [Ts(h_kJ, p) for p in ps]
    return np.array([o[1] for o in out]), np.array([o[0]-273.15 for o in out])

fig, ax = plt.subplots(figsize=(7.4, 5.4))
dome(ax)
pc_, pi_, pe_ = hp['p_condh'], hp['p_inth'], hp['p_evaph']
RED, PUR = '#c02a2a', '#7b3fa0'

s,t = isentrope(hp['s_1h'], pe_, pi_);            ax.plot(s,t,color=ORG,lw=2.2,label='compressors (isentropic)')
s,t = isentrope(hp['s_6h'], pi_, pc_);            ax.plot(s,t,color=ORG,lw=2.2)
s,t = isobar(hp['h_2h'], hp['h_3h'], pc_);        ax.plot(s,t,color=RED,lw=2.4,label='condenser (delivers to the store)')
s,t = isobar(hp['h_3h'], hp['h_7h'], pc_);        ax.plot(s,t,color=PUR,lw=2.2,label='IHX-1 then IHX-2, liquid side')
s,t = isenthalp(hp['h_7h'], pc_, pi_);            ax.plot(s,t,color='0.45',lw=1.5,ls='--',label='expansion valves (isenthalpic)')
s,t = isenthalp(hp['h_9h'], pi_, pe_);            ax.plot(s,t,color='0.45',lw=1.5,ls='--')
s,t = isobar(hp['h_4h'], hp['h_13h'], pe_);       ax.plot(s,t,color=BLU,lw=2.4,label='evaporator (source)')
s,t = isobar(hp['h_13h'], hp['h_1h'], pe_);       ax.plot(s,t,color=GRN,lw=2.2,label='IHX vapour side (suction / flash gas)')
s,t = isobar(hp['h_10h'], hp['h_11h'], pi_);      ax.plot(s,t,color=GRN,lw=2.2)

for k,(dx,dy) in {'1h':(.006,-4),'2h':(.006,4),'3h':(.006,4),'4h':(-.004,-13),
                  '7h':(-.036,-4),'9h':(-.034,-4),'10h':(.006,-12),'11h':(.006,4),
                  '12h':(-.040,6),'13h':(.006,-12),'8h':(.006,-4)}.items():
    label(ax, hp['s_'+k]/1000, hp['T_'+k]-273.15, k, dx, dy)

ax.annotate('flash separator at $p_{int}$:\n$x_{8h}=%.4f$ vapour $\\to$ IHX-1\n$1-x_{8h}=%.4f$ liquid $\\to$ evaporator'
            %(hp['x_8h'],1-hp['x_8h']),
            xy=(hp['s_8h']/1000, hp['T_8h']-273.15), xytext=(0.66, 126),
            arrowprops=dict(arrowstyle='->', color='crimson', lw=1.2),
            fontsize=8.2, color='crimson')
ax.text(0.98, 0.03, 'the $(1-x_{8h})$ split is what v0.8 got wrong\n(it used the evaporator-inlet quality $x_{4h}$)',
        transform=ax.transAxes, fontsize=8.2, color='crimson', va='bottom', ha='right')
ax.annotate('IHX-1 + IHX-2\nsubcool the liquid', xy=(0.5*(hp['s_3h']+hp['s_7h'])/1000, 161),
            xytext=(0.24, 118), arrowprops=dict(arrowstyle='->', color=PUR, lw=1.1),
            fontsize=8.2, color=PUR, ha='center')
ax.set_xlabel('specific entropy  $s$  [kJ kg$^{-1}$ K$^{-1}$]')
ax.set_ylabel('temperature  $T$  [°C]')
ax.set_title('High-temperature heat pump — two stage, flash separator, two regenerators\n'
             r'cyclopentane,  COP $=%.4f$,  $x_{8h}=%.4f$,  v0.9'%(hp['hp_cop'],hp['x_8h']), fontsize=10)
ax.set_xlim(-0.05, 1.75); ax.set_ylim(20, 250)
ax.grid(alpha=.25); ax.legend(fontsize=7.6, loc='upper left', framealpha=.92)
fig.tight_layout(); fig.savefig('/tmp/ts_hthp.png', dpi=175); fig.savefig('/tmp/ts_hthp.pdf')
print('HTHP figure written')
