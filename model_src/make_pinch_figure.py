import sys; sys.path.insert(0,'/tmp/build')
import numpy as np, CoolProp.CoolProp as CP
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import thums_model as T

F='cyclopentane'
c=T.Case(); rank9,hp,Tst = T.cycle_state_points(c)
T_wh, T_wc = Tst['T_2d'], Tst['T_3d']
rank8 = T.double_stage_rankine(c.fluid, Tst['T_1e'], Tst['T_3e_hot_end_rule'])

def water(Qtot, n=300):
    Tw = np.linspace(T_wc, T_wh, n)
    hw = np.array([CP.PropsSI('H','T',t,'P',c.P,'Water')/1000. for t in Tw])
    return (hw-hw[0])*Qtot/(hw[-1]-hw[0]), Tw-273.15

def streams(rank, n=200):
    """the two cold streams, separately, as (Q,T) from their own cold ends"""
    y=rank['y_frac']; pe,pi=rank['p_evape'],rank['p_5e']
    h_bub=CP.PropsSI('H','P',pe,'Q',0,F)/1000.
    hm=np.linspace(rank['h_10e'],rank['h_3e'],n)
    Tm=np.array([CP.PropsSI('T','H',h*1000,'P',pe,F) for h in hm])-273.15
    Qm=hm-rank['h_10e']
    hr=np.linspace(rank['h_5e'],rank['h_6e'],n)
    Tr=np.array([CP.PropsSI('T','H',h*1000,'P',pi,F) for h in hr])-273.15
    Qr=y*(hr-rank['h_5e'])
    return (Qm,Tm),(Qr,Tr)

fig,axes = plt.subplots(1,3,figsize=(15.2,4.9))
BLU,RED,GRN,ORG='#1f5f9e','#c02a2a','#2e8b57','#d2691e'

# ---------- (a) how the cold composite is assembled -----------------------
ax=axes[0]
(Qm,Tm),(Qr,Tr)=streams(rank9)
Qc,Tc=T._orc_cold_composite(rank9,F); Tc=Tc-273.15
ax.plot(Qm,Tm,color=GRN,lw=2.0,ls='--',label='main stream 10e$\\to$3e  (1 kg)')
ax.plot(Qr,Tr,color=ORG,lw=2.0,ls='--',label='reheat 5e$\\to$6e  ($y$ kg)')
ax.plot(Qc,Tc,color=RED,lw=2.8,label='COLD COMPOSITE = sum at each $T$')
ax.annotate('liquid preheat\n(sloping)',xy=(35,85),xytext=(8,118),fontsize=8.2,color=GRN,
            arrowprops=dict(arrowstyle='->',color=GRN,lw=1))
ax.annotate('boiling: FLAT in $(Q,T)$ —\n76 % of the duty at one $T$.\nThis is the whole difficulty.',xy=(250,95.4),xytext=(150,128),
            fontsize=8.4,color=RED,arrowprops=dict(arrowstyle='->',color=RED,lw=1.1))
ax.set_title('(a) building the cold composite',fontsize=10.5)
ax.legend(fontsize=7.8,loc='lower right',framealpha=.95)

# ---------- (b) v0.8 hot-end rule : curves cross -------------------------
ax=axes[1]
Qc8,Tc8=T._orc_cold_composite(rank8,F); Tc8=Tc8-273.15
Qh8,Th8=water(Qc8[-1])
ax.plot(Qh8,Th8,color=BLU,lw=2.6,label='water (hot stream)')
ax.plot(Qc8,Tc8,color=RED,lw=2.6,label='ORC composite (cold)')
Qg=np.linspace(0,Qc8[-1],400)
tw=np.interp(Qg,Qh8,Th8); tc=np.interp(Qg,Qc8,Tc8)
ax.fill_between(Qg,tw,tc,where=tc>tw,color='red',alpha=.16,interpolate=True)
i=np.argmin(tw-tc)
ax.plot([Qg[i],Qg[i]],[tc[i],tw[i]],color='k',lw=1.4)
ax.annotate('pinch $=%+.1f$ K\nat %.0f %% of the duty'%(np.min(tw-tc),100*Qg[i]/Qc8[-1]),
            xy=(Qg[i],0.5*(tw[i]+tc[i])),xytext=(Qg[i]+62,58),fontsize=8.8,fontweight='bold',
            arrowprops=dict(arrowstyle='->',color='k',lw=1.2))
ax.text(0.5,0.95,'heat flows COLD $\\to$ HOT in the shaded region',transform=ax.transAxes,
        ha='center',va='top',fontsize=9,color='darkred',fontweight='bold')
ax.set_title('(b) v0.8 — hot-end rule, $T_{3e}=%.1f$ °C, $\\eta=%.4f$'%(rank8['T_3e']-273.15,rank8['rank_eff']),
             fontsize=10.5,color='darkred')
ax.legend(fontsize=8,loc='lower right')

# ---------- (c) v0.9 pinch-feasible --------------------------------------
ax=axes[2]
Qc9,Tc9=T._orc_cold_composite(rank9,F); Tc9=Tc9-273.15
Qh9,Th9=water(Qc9[-1])
ax.plot(Qh9,Th9,color=BLU,lw=2.6,label='water (hot stream)')
ax.plot(Qc9,Tc9,color=RED,lw=2.6,label='ORC composite (cold)')
Qg=np.linspace(0,Qc9[-1],400)
tw=np.interp(Qg,Qh9,Th9); tc=np.interp(Qg,Qc9,Tc9)
ax.fill_between(Qg,tw,tc,color=BLU,alpha=.10)
i=np.argmin(tw-tc)
ax.plot([Qg[i],Qg[i]],[tc[i],tw[i]],color='k',lw=1.4)
ax.annotate('pinch $=%+.2f$ K\n$=\\Delta T_{pinch,ORC}$'%np.min(tw-tc),
            xy=(Qg[i],0.5*(tw[i]+tc[i])),xytext=(Qg[i]+95,62),fontsize=8.8,fontweight='bold',
            arrowprops=dict(arrowstyle='->',color='k',lw=1.2))
ax.set_title('(c) v0.9 — pinch-feasible, $T_{3e}=%.1f$ °C, $\\eta=%.4f$'%(rank9['T_3e']-273.15,rank9['rank_eff']),
             fontsize=10.5,color='darkgreen')
ax.legend(fontsize=8,loc='lower right')

for ax in axes:
    ax.set_xlabel('cumulative duty  $Q$  [kJ per kg of ORC evaporator flow]')
    ax.set_ylim(40,165); ax.grid(alpha=.25)
axes[0].set_ylabel('temperature  [°C]')
fig.suptitle('The ORC evaporator pinch: why an approach at the hot end is not a constraint',
             fontsize=12, y=1.00)
fig.tight_layout(); fig.savefig('/tmp/pinch.png',dpi=165); fig.savefig('/tmp/pinch.pdf')
print('pinch figure written; pinch v0.8 %.2f K, v0.9 %.2f K'
      %(T.orc_pinch(rank8,T_wh,T_wc,F,c.P), T.orc_pinch(rank9,T_wh,T_wc,F,c.P)))
