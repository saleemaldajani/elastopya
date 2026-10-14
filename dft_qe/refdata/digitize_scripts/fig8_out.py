import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, pk, finaldir = sys.argv[1:4]
a=load(img).astype(np.uint8); d=pickle.load(open(pk,'rb')); pc,pm=d['pc'],d['pm']
cx=lambda x: np.polyval(pc,x); my=lambda y: np.polyval(pm,y); xc=lambda c:(np.asarray(c)-pc[1])/pc[0]; ym=lambda m:(np.asarray(m)-pm[1])/pm[0]
rows=[]
for x in sorted(d['fe']): rows.append(('M_Fe FM','broken (dashed) black line',cx(x),my(d['fe'][x]),x,d['fe'][x],'line'))
for x in sorted(d['cr']): rows.append(('M_Cr FM','solid black line',cx(x),my(d['cr'][x]),x,d['cr'][x],'line'))
seen=set()
for s in d['sym']:
    key=(s[0],round(s[1]))
    if key in seen: continue
    seen.add(key); rows.append((s[0].split(' (')[0], s[0].split('(')[1].rstrip(')'), cx(s[1]), my(s[2]), s[1], s[2], 'symbol: '+s[6]))
fn=f'{finaldir}/korzhavyi2009_fig8_moments.csv'
with open(fn,'w') as f:
    f.write("""# Spin magnetic moments of Fe and Cr atoms in random bcc Fe-Cr vs Cr content, FM (lines) and paramagnetic DLM (symbols), DIGITIZED from Fig. 8 (p. 054202-7)
# Source: P. A. Korzhavyi, A. V. Ruban, J. Odqvist, J.-O. Nilsson, B. Johansson, Phys. Rev. B 79, 054202 (2009)
# Caption: 'Calculated concentration dependencies of the magnetic moments on the Fe and Cr atoms in the ferromagnetic (lines) and paramagnetic
#   (symbols) states of Fe-Cr alloys.'  Method: EMTO-CPA (LDA self-consistent densities), FM and DLM.
# Sec. III D (p. 054202-7), quoted: 'The results to be presented in this section have been calculated at a fixed lattice parameter of
#   a0 = 2.875 A, which is the average of the room-temperature lattice parameters of bcc Fe and bcc Cr.' (volume dependence stated to be weak)
# IMPORTANT definition: these are NOT sphere-integrated local moments. FM: partial-molar moments from the total moment per atom m_tot(c):
#   M_Fe = m_tot - c dm_tot/dc, M_Cr = m_tot + (1-c) dm_tot/dc (Eqs. 6-8). DLM: M = m_tot/(2 delta) after flipping a small fraction delta of
#   one species' moments (Eq. 9; delta = 0.5 at.% at 5 at.% of that component, 1 at.% at 10 at.%, 2 at.% otherwise).
# Line styles: the LEGEND shows 'M_Fe, FM' as a broken (long-dash, sample ends in a short dash/dot) line and 'M_Cr, FM' as SOLID; the TEXT
#   (p. 054202-7) says 'dashed (dot-dashed) line for Fe (Cr)'. In the plot the Cr line is solid and the Fe line is the broken one (legend wins). DLM: open red circles = M_Fe, DLM; open red squares = M_Cr, DLM (text: M_Cr vanishes in DLM).
# Units: c_Cr in at.% Cr; moment in mu_B (signed; negative = antiparallel to Fe / global magnetization).
# Digitization: embedded raster (700x566 px JPEG, ~300 dpi) at native resolution. Calibration (pixel x,y in that image):
#   x: c=0 -> 85, 20 -> 204, 40 -> 322.5, 60 -> 441, 80 -> 560, 100 -> 678 (resid <0.06 at.%); 1 px = 0.169 at.%.
#   y: 3 -> 10, 2 -> 105, 1 -> 200, 0 -> 295, -1 -> 389.5, -2 -> 484 (resid <0.004 mu_B); 1 px = 0.0105 mu_B.
#   Lines traced column-by-column (centre of ink run). Symbols: centre of the white interior of each open symbol.
# Gaps: the M_Fe FM dashed line is largely HIDDEN beneath the DLM circles (white-filled symbols drawn on top) and the legend box (c~29-63, M<1.85
#   region only; the line itself runs above the box); only dash fragments between circles are digitized (interpolate between them; values agree
#   with the circles within ~0.05 mu_B for c<75). M_Cr FM line is partly hidden under the DLM squares for c > ~60 (where it is ~0).
# Uncertainty: lines +/-0.015 mu_B, +/-0.3 at.%; symbols +/-0.01 mu_B, +/-0.2 at.%.
# Columns: series, style, c_Cr_at_pct, moment_muB, x_px, y_px, kind
series,style,c_Cr_at_pct,moment_muB,x_px,y_px,kind
""")
    for r in rows: f.write(f'{r[0]},{r[1]},{r[2]:.2f},{r[3]:.4f},{r[4]:.1f},{r[5]:.1f},"{r[6]}"\n')
fig,axs=plt.subplots(2,1,figsize=(10,14),dpi=130)
for ax,(x0,x1,y0,y1) in zip(axs,[(60,700,0,500),(60,400,60,440)]):
    ax.imshow(a,interpolation='nearest'); ax.set_xlim(x0,x1); ax.set_ylim(y1,y0)
    for ser,col,m in [('M_Fe FM','cyan','.'),('M_Cr FM','lime','.'),('M_Fe DLM','blue','+'),('M_Cr DLM','blue','x')]:
        pts=[r for r in rows if r[0]==ser]
        ax.plot(xc([r[2] for r in pts]), ym([r[3] for r in pts]), m, color=col, ms=2 if m=='.' else 8, label=ser)
    ax.legend(fontsize=7,loc='lower right')
axs[0].set_title('Korzhavyi 2009 Fig.8: digitized (data->pixel) over original; bottom: zoom')
fig.tight_layout(); fig.savefig(f'{finaldir}/korzhavyi2009_fig8_moments_overlay.png')
from collections import Counter; print(Counter(r[0] for r in rows))
# key numbers
fe=sorted([(r[2],r[3]) for r in rows if r[0]=='M_Fe FM']); cr=sorted([(r[2],r[3]) for r in rows if r[0]=='M_Cr FM'])
for c in (20,35,50,65,95):
    print(c, 'MFe FM', round(float(np.interp(c,*zip(*fe))),3), 'nearest pt dist', round(min(abs(p[0]-c) for p in fe),2), 'MCr FM', round(float(np.interp(c,*zip(*cr))),3))
