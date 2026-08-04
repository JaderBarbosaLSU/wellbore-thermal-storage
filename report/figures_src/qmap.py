import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTED='#8a8983'
S=['#2a78d6','#eb6834','#1baf7a']
p=101325.
def psat(T): return 610.94*np.exp(17.625*T/(T+243.04))
def w_rh(T,RH): pv=RH*psat(T); return 0.622*pv/(p-pv)
T_pc, w_pc = 37.0, w_rh(37,0.50)
CP,hfg=1006.,2.5e6
A=1e6/200.; V=A*3.5; Aw=4*np.sqrt(A)*3.5
Htr=0.5*(2*A+Aw); mi=1.2*V*0.5/3600.; mv=1.2*(1e6*0.030)/3600.
Qint=1e6*0.13; mev=1e6*4.6e-9
def loads(T,RH):
    w=w_rh(T,RH)
    return ((Qint+(Htr+(mi+mv)*CP)*(T-T_pc))/1e3,
            (mev*hfg+(mi+mv)*(w-w_pc)*hfg)/1e3)

# label offsets tuned by hand to avoid collisions
pts=[('Jan mean',12,.75,(10,6),'left'),('Mar mean',18,.74,(10,6),'left'),
     ('Apr mean',21,.74,(10,6),'left'),('Jul mean',28,.78,(-10,-4),'right'),
     ('warm humid night',28,.90,(-10,8),'right'),
     ('summer design humid',31,.75,(12,-20),'left'),
     ('summer design',33,.55,(12,4),'left'),('hot-dry excursion',36,.30,(-10,6),'right'),
     ('cold snap',4,.70,(10,6),'left'),('hot humid extreme',35,.70,(-10,4),'right')]

fig,ax=plt.subplots(figsize=(11.6,8.4),facecolor=SURF)
ax.set_facecolor(SURF)
XL,YL=(-430,240),(-540,240)
ax.add_patch(Rectangle((0,0),XL[1],YL[1],color='#f3f6fa',zorder=0))
ax.add_patch(Rectangle((0,YL[0]),XL[1],-YL[0],color='#fdf4ef',zorder=0))
ax.add_patch(Rectangle((XL[0],0),-XL[0],YL[1],color='#f1faf6',zorder=0))
ax.add_patch(Rectangle((XL[0],YL[0]),-XL[0],-YL[0],color='#f7f7f4',zorder=0))
ax.axhline(0,color='#9a998f',lw=1.6,zorder=2); ax.axvline(0,color='#9a998f',lw=1.6,zorder=2)

for title,modes,x,y,ha,va in [
   ('cooling + dehumidification','M1 mild · M2 moderate · M3 harsh',232,228,'right','top'),
   ('cooling + humidification','M4 mild · M5 harsh',232,-528,'right','bottom'),
   ('heating + dehumidification','M6  (needs modulating d_pD / d_pE)',-420,228,'left','top'),
   ('heating + humidification','M7 mild · M8 full',-420,-528,'left','bottom')]:
    dy=-17 if va=='top' else 17
    ax.text(x,y,title,ha=ha,va=va,fontsize=11.5,fontweight='bold',color=INK,zorder=5)
    ax.text(x,y+dy,modes,ha=ha,va=va,fontsize=9.5,color=INK2,zorder=5)

for lab,T,RH,off,ha in pts:
    qs,ql=loads(T,RH)
    col=S[0] if (qs>0)==(ql>0) else S[1]
    ax.plot(qs,ql,'o',ms=8,mfc=col,mec=SURF,mew=1.8,zorder=4)
    ax.annotate(f'{lab}\n{T:.0f} °C / {RH:.0%}',(qs,ql),xytext=off,
                textcoords='offset points',fontsize=8.2,color=INK2,zorder=5,ha=ha)
ax.plot(0,0,'s',ms=9,mfc=INK,mec=SURF,mew=1.6,zorder=6)
ax.annotate('setpoint  37 °C / 50 % RH',(0,0),xytext=(14,10),textcoords='offset points',
            fontsize=9,fontweight='bold',color=INK,zorder=6)
ax.annotate('M3 and M5 need ≥ 79 °C at the generator —\navailable in year 1, marginal by year 20',
            (232,150),fontsize=9,color=S[1],style='italic',ha='right',va='top',zorder=5)
ax.annotate('M1, M2, M6, M7, M8 need only ≈ 66 °C —\nthe well holds this for the whole horizon',
            (-420,-300),fontsize=9,color=S[2],style='italic',ha='left',va='top',zorder=5)
ax.set_xlim(*XL); ax.set_ylim(*YL)
ax.set_xlabel('sensible load  q_s  [kW]      (positive = cooling required)',fontsize=10.5,color=INK2)
ax.set_ylabel('latent load  q_l  [kW]      (positive = dehumidification required)',fontsize=10.5,color=INK2)
ax.set_title('Load-quadrant map with topology assignment\n'
             'New Orleans conditions against a 37 °C / 50 % RH setpoint, hatchery loads as coded',
             fontsize=13,fontweight='bold',color=INK,loc='left',pad=14)
for sp in ('top','right'): ax.spines[sp].set_visible(False)
for sp in ('left','bottom'): ax.spines[sp].set_color('#d6d5d0')
ax.tick_params(colors=INK2,labelsize=9)
ax.grid(True,color='#e9e8e3',lw=0.7,zorder=1)
fig.text(0.008,0.006,'Loads from a whole-building control volume including forced ventilation '
         '(10 kg/s). The notebook accounts for ventilation through the supply-air balance, so its '
         'q_s differs — reconcile conventions before publication.',fontsize=7.6,color=MUTED)
fig.tight_layout(rect=[0,0.028,1,1])
fig.savefig(f'{OUT}/quadrant_map.png',dpi=150,facecolor=SURF)
fig.savefig(f'{OUT}/quadrant_map.pdf',facecolor=SURF)
print('ok')
