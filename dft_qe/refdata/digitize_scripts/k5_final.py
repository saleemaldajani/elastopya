import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from tmatch import *
from PIL import Image
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, finaldir = sys.argv[1:3]
a=np.asarray(Image.open(img).convert('L')); ink=a>128
Td=make_template(ink,567,59,11,11); Ts=make_template(ink,453,324,11,11)
sc={'diamond':match(ink,Td),'square':match(ink,Ts)}
off={'diamond':(0.5,0.5),'square':(0.0,0.0)}
P={'a':dict(region=(125,785,24,530), x=(109.0,797.5,0,1.0), y=(10.5,540.0,200,-50)),
   'b':dict(region=(112,762,697,1195), x=(101.0,772.5,0,0.2), y=(687.5,1204.0,40,-10))}
conv=lambda pn,x,y:(P[pn]['x'][2]+(x-P[pn]['x'][0])/(P[pn]['x'][1]-P[pn]['x'][0])*(P[pn]['x'][3]-P[pn]['x'][2]),
                    P[pn]['y'][2]+(y-P[pn]['y'][0])/(P[pn]['y'][1]-P[pn]['y'][0])*(P[pn]['y'][3]-P[pn]['y'][2]))
inv=lambda pn,c,E:(P[pn]['x'][0]+(c-P[pn]['x'][2])/(P[pn]['x'][3]-P[pn]['x'][2])*(P[pn]['x'][1]-P[pn]['x'][0]),
                   P[pn]['y'][0]+(E-P[pn]['y'][2])/(P[pn]['y'][3]-P[pn]['y'][2])*(P[pn]['y'][1]-P[pn]['y'][0]))
rows=[]
for pn in ('a','b'):
    allp=[]
    for name in (('diamond','square') if pn=='a' else ('diamond',)):
        for x,y,s in peaks(sc[name],0.55,mind=4,region=P[pn]['region']):
            allp.append([name,x+off[name][0],y+off[name][1],s,'auto'])
    allp.sort(key=lambda p:-p[3]); keep=[]
    for p in allp:
        if all((p[1]-q[1])**2+(p[2]-q[2])**2>5**2 for q in keep): keep.append(p)
    if pn=='a':
        keep=[p for p in keep if conv('a',p[1],p[2])[0]>0.2]   # c<=0.2 taken from the close-up panel (b)
        keep.append(['diamond',453.5,13.0,float(sc['diamond'][13,453]),'manual: on top frame, bbox rows 5-21'])
    else:
        keep.append(['diamond',287.5,1151.5,float(sc['diamond'][1151,287]),'manual: stacked, interior visible rows 1150-1153'])
    for n,x,y,s,src in keep:
        c,E=conv(pn,x,y); rows.append((pn,n,c,E,x,y,s,src))
rows.sort(key=lambda r:(r[0],r[2],r[3]))
lab={'diamond':'Delta_E_f (diamond: Fe-Cr structure in magnetic ground state)','square':'Delta_E_f (square: not well mixed system)'}
fn=f'{finaldir}/klaver2006_fig5_heat_of_formation.csv'
with open(fn,'w') as f:
    f.write("""# 0 K heat of formation Delta_E_f of Fe-Cr structures (magnetic ground state) vs Cr concentration, DIGITIZED from Fig. 5(a) and 5(b) (p. 094435-5)
# Source: T. P. C. Klaver, R. Drautz, M. W. Finnis, Phys. Rev. B 74, 094435 (2006)
# Caption: '(a) The 0 K heat of formation Delta E_f of Fe-Cr structures in their magnetic ground state. The five data points with squares indicate
#   systems that were not very well mixed, i.e., Fe interacted mostly with Fe and Cr mostly with Cr, which explains the low Delta E_f values.
#   (b) A close-up of the Fe-rich part of the concentration range.'  No legend; series named from the caption.
# Definition (Eq. 3): Delta_E_f = [E(Fe_x Cr_1-x) - x n E_ref(Fe) - (1-x) n E_ref(Cr)]/n, references = pure FM bcc Fe and AFM bcc Cr at their
#   equilibrium spacing (for some large Fe-rich cells the Fe reference is a large pure-Fe cell). 65 structures, 2-atom B2 up to 128-atom 4x4x4
#   cells. Method: VASP PAW, GGA-PBE, 400 eV, positions + cell shape + volume fully relaxed (Sec. II). Abstract: 'From 0-12 % Cr a small negative
#   Delta E_f down to -8 meV/atom is found'; mostly positive up to ~100 meV/atom (text) - the plot shows points up to ~199 meV/atom.
#   NOT a random-alloy (SQS/CPA) mixing enthalpy: specific ordered/clustered supercells.
# Panels used: c_Cr <= 0.2 from the close-up panel (b) (all points; more precise); c_Cr > 0.2 from panel (a). Column 'panel' says which.
# Units: c_Cr = Cr atomic fraction (axis 'Cr concentration'); dEf_meV_per_atom in meV/atom.
# Digitization: embedded 1-bit raster (800x1288 px, ~300 dpi, both panels) at native resolution; binary template matching (templates cut from
#   isolated diamond/square markers), score >= 0.55, plus 2 manual points (noted). Calibration (pixel x,y):
#   panel a: x: c=0 -> 109 (frame), ticks every 0.05 = 34.43 px, c=1 -> 797.5; y: 200 -> 10.5 (frame), 0 -> 434.5 (zero line), -50 -> 540.0;
#     ticks every 10 meV = 21.2 px (1 px = 0.47 meV/atom, 0.00145 in c).
#   panel b: x: c=0 -> 101, ticks every 0.01 = 33.57 px, c=0.20 -> 772.5; y: 40 -> 687.5, 0 -> 1100.5 (zero line), -10 -> 1204.0;
#     ticks every 2 meV = 20.66 px (1 px = 0.097 meV/atom, 0.0003 in c).
#   The values cluster on integers in panel (b) (e.g. -2.0, -4.0, 12.0 meV), suggesting the authors plotted rounded values.
# Uncertainty: panel b +/-0.2 meV/atom, +/-0.001 in c; panel a +/-1 meV/atom, +/-0.003 in c (+/-2 meV for stacked/overlapping markers).
# Columns: panel, series, symbol, c_Cr, dEf_meV_per_atom, x_px, y_px, match_score, source
panel,series,symbol,c_Cr,dEf_meV_per_atom,x_px,y_px,match_score,source
""")
    for r in rows: f.write(f'{r[0]},"{lab[r[1]]}",{r[1]},{r[2]:.4f},{r[3]:.2f},{r[4]:.1f},{r[5]:.1f},{r[6]:.2f},"{r[7]}"\n')
fig,axs=plt.subplots(1,2,figsize=(18,12),dpi=100)
for ax,pn,(y0,y1) in zip(axs,('a','b'),((0,560),(680,1215))):
    ax.imshow(255-a,cmap='gray',interpolation='nearest'); ax.set_xlim(90,805); ax.set_ylim(y1,y0)
    for r in rows:
        if r[0]!=pn: continue
        X,Y=inv(pn,r[2],r[3]); ax.plot(X,Y,'+' if r[1]=='diamond' else 'x',color='red' if r[1]=='diamond' else 'blue',ms=11,mew=1.5)
    ax.set_title(f'Klaver 2006 Fig.5({pn}) digitized (data->pixel) over original')
fig.tight_layout(); fig.savefig(f'{finaldir}/klaver2006_fig5_heat_of_formation_overlay.png')
from collections import Counter; print(Counter((r[0],r[1]) for r in rows))
for r in rows: print(r[0],r[1][0],round(r[2],4),round(r[3],1))
