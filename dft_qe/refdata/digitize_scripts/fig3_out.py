import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, pk, outdir = sys.argv[1:4]
a = load(img).astype(np.uint8); d = pickle.load(open(pk,'rb'))
pb, pd_ = d['pb'], d['pd']
X0,X100=99.,678.; cx=lambda x:(x-X0)*100/(X100-X0); px=lambda c: X0+np.asarray(c)*(X100-X0)/100
yB=lambda v:(np.asarray(v)-pb[1])/pb[0]; yN=lambda v:(np.asarray(v)-pd_[1])/pd_[0]
hdr_common = """# Source: P. A. Korzhavyi, A. V. Ruban, J. Odqvist, J.-O. Nilsson, B. Johansson, Phys. Rev. B 79, 054202 (2009),
#   "Electronic structure and effective chemical and magnetic exchange interactions in bcc Fe-Cr alloys"
# Calculations (Sec. II A): EMTO method + single-site CPA for random bcc Fe-Cr; paramagnetic state modelled by DLM within CPA;
#   total energies in GGA (PW91, Ref. 39) evaluated on LDA (Perdew-Wang 1992) self-consistent densities (gradient corrections off in the SCF potential).
"""
# ---------- Fig 3(b) ----------
rows = d['b']
with open(f'{outdir}/korzhavyi2009_fig3b_bulk_modulus.csv','w') as f:
    f.write("# Bulk modulus B of bcc Fe-Cr vs Cr content: experimental markers + EMTO-CPA theory points, DIGITIZED from Fig. 3(b) (p. 054202-4)\n")
    f.write(hdr_common)
    f.write("""# Text (p. 054202-4): theory B 'obtained from the calculated equations of state for FM alloys on a static lattice (i.e., without phonon
#   contributions)' i.e. each alloy at its own calculated equilibrium volume (static, 0 K); DLM curve = paramagnetic (DLM) state, same method.
#   Paper: FM B has 'a sharp minimum at 6 at % Cr'; DLM B is monotonous.
# Experimental data as plotted (Refs. 49-52 of the paper; legend shapes): Expt. RT FM = filled black triangles (room temperature, ferromagnetic; Ref. 49 = Speich, Schwoeble, Leslie, Metall. Trans. 3, 2031 (1972));
#   Expt. LT FM = open triangles (77 K, ferromagnetic, Ref. 49); Expt. RT AFM = filled diamond and Expt. LT AFM = open diamond (pure Cr; Ref. 51 Lahteenkorva & Lenkkeri J. Phys. F 11, 767 (1981); Ref. 52 Bolef & de Klerk PR 129, 1063 (1963));
#   Expt. RT PM = filled red squares (RT paramagnetic Cr-rich alloys; Ref. 50 = Lenkkeri, J. Phys. F 10, 611 (1980); the c=100 square is the 'experimental estimate of the RT bulk modulus
#   of paramagnetic Cr', Refs. 50,51). No error bars are plotted.
# Theory: 'Theory FM' = black open circles joined by solid line (legend of panel a); 'Theory DLM' = red open circles joined by dashed line.
#   The plotted lines are straight segments between the circles, so the circle centres define the curves.
#   For c ~ 4-9 at.% the FM circles overlap into an unresolvable chain (and overlap the RT-FM triangles); there the FM curve is given as
#   'envelope estimate' = top edge of the circle chain + ring outer radius (10.5 px), uncertainty ~ +/-3 GPa. Individual FM points between
#   c=3 and c=13 (the paper's minimum near 6 at.%) could not be resolved; the c=50 and c=75 FM circles are half-covered by DLM circles.
# Digitization: embedded raster image (700x1363 px JPEG, ~300 dpi; whole Fig. 3) extracted at native resolution with PyMuPDF.
#   Calibration (pixel x,y in that image): x: c=0 -> x=99, 20->215, 40->331, 60->446, 80->562, 100->678 (major ticks; linear, max resid <0.5 px).
#   y (panel b): B=250 -> y=369, 200 -> 451, 150 -> 532, 100 -> 614 (major ticks; linear fit, resid <0.2 px); 1 px = 0.613 GPa, 0.173 at.%.
#   Markers: centre of symbol bounding box (Grace-style symbols are bbox-centred), read from pixel maps; theory circles: centre of white core.
# Uncertainty: isolated markers +/-1 GPa, +/-0.3 at.%; overlapping markers near c=0-2 and c~10.7 +/-2 GPa, +/-0.5 at.% (assignment of hidden
#   triangles is a best reading); theory circle centres +/-1 GPa (+/-2 GPa for 'partial core' rows); envelope rows +/-3 GPa.
#   Marker counts may be incomplete where symbols are fully hidden (c=0-2 cluster contains at least 2 filled + 1 open triangle).
# Units: c_Cr in at.% Cr; B_GPa in GPa. x_px,y_px = pixel position used (for audit).
# Columns: series, symbol, c_Cr_at_pct, B_GPa, kind, x_px, y_px, note
series,symbol,c_Cr_at_pct,B_GPa,kind,x_px,y_px,note
""")
    for s,m,c,v,x,y,kind,note in rows:
        f.write(f'{s},{m},{c:.2f},{v:.1f},{kind},{x:.1f},{y:.1f},"{note}"\n')
fig,axs=plt.subplots(2,1,figsize=(10,11),dpi=150)
for ax,(y0,y1,x0,x1) in zip(axs,[(329,648,60,700),(440,550,85,250)]):
  ax.imshow(a[y0:y1,x0:x1],extent=(x0-.5,x1-.5,y1-.5,y0-.5),interpolation='nearest')
  sty={'Expt. RT FM':('x','cyan'),'Expt. LT FM':('+','cyan'),'Expt. RT AFM':('x','lime'),'Expt. LT AFM':('+','lime'),'Expt. RT PM':('x','blue'),'Theory FM':('.','orange'),'Theory DLM':('.','lime')}
  for s in sty:
    pts=[r for r in rows if r[0]==s]
    if s.startswith('Theory'):
        pts=sorted(pts,key=lambda r:r[2])
        ax.plot(px([r[2] for r in pts]), yB([r[3] for r in pts]), '-'+sty[s][0], color=sty[s][1], lw=0.8, ms=5, label=s+' (digitized)')
    else:
        ax.plot(px([r[2] for r in pts]), yB([r[3] for r in pts]), sty[s][0], color=sty[s][1], ms=9, mew=1.5, label=s)
  ax.legend(fontsize=6,loc='lower right'); ax.set_xlim(x0,x1); ax.set_ylim(y1,y0)
axs[0].set_title('Korzhavyi 2009 Fig.3(b): digitized data (data->pixel via calibration) over original; bottom: zoom c=0-26')
fig.savefig(f'{outdir}/korzhavyi2009_fig3b_bulk_modulus_overlay.png'); plt.close(fig)
# ---------- Fig 3(d) ----------
rows=d['d']; lx0,lx1=d['legend_x']
with open(f'{outdir}/korzhavyi2009_fig3d_nef.csv','w') as f:
    f.write("# Density of states at the Fermi level N(E_F) of random bcc Fe-Cr vs Cr content, DIGITIZED from Fig. 3(d) (p. 054202-4)\n")
    f.write(hdr_common)
    f.write(f"""# Series (legend): 'Total DLM' = red dashed (paramagnetic DLM, total), 'Total FM' = black solid (ferromagnetic, total, both spins),
#   'Spin down FM' = black dash-dot (FM, spin-down = minority-Fe channel only). Text (p. 054202-5): the slope discontinuity near 7 at.% Cr occurs
#   only in the spin-up channel. Units on axis: N(E_F) (Ry^-1), i.e. states/Ry; per atom presumably (not stated explicitly in caption/text).
# Lattice parameter for these DOS calculations is not stated in the caption or Sec. III C (Sec. III D says only the magnetic-moment results use
#   fixed a0 = 2.875 A); presumably each alloy at its calculated equilibrium volume - NOT verified from the text.
# Digitization: embedded raster image (700x1363 px JPEG, whole Fig. 3) at native resolution; curves traced column-by-column (1 px = 0.173 at.%)
#   by colour (red) / darkness (black) with slope-predicted tracking; y = centre of the line's vertical run in each column.
#   Calibration: x as Fig.3(b): c=0 -> x=99, 20->215, 40->331, 60->446, 80->562, 100->678.
#   y (panel d): N=20 -> y=1048, 10 -> 1166, 0 -> 1284 (bottom frame = 0 tick); minor ticks 25/15/5 at 989/1107/1225 consistent; 1 px = 0.0847 Ry^-1.
# Gaps: dashed/dash-dot lines only sampled where ink is present (dash gaps). The 'Spin down FM' line is HIDDEN behind the legend box for
#   c = {cx(lx0):.1f}-{cx(lx1):.1f} at.% (no data there). Near c > 95 the DLM and FM curves converge (both ~9.8 Ry^-1 at c=100).
# Uncertainty: +/-0.15 Ry^-1 in N (about 1.5-2 px incl. line width ~4 px), +/-0.3 at.% in c.
# Columns: series, c_Cr_at_pct, NEF_states_per_Ry, x_px, y_px
series,c_Cr_at_pct,NEF_states_per_Ry,x_px,y_px
""")
    for s,c,v,x,y in rows: f.write(f'{s},{c:.2f},{v:.3f},{x},{y:.1f}\n')
fig,ax=plt.subplots(figsize=(10,6),dpi=150)
y0,y1,x0,x1=965,1287,60,700
ax.imshow(a[y0:y1,x0:x1],extent=(x0-.5,x1-.5,y1-.5,y0-.5),interpolation='nearest')
for s,col in [('Total DLM','blue'),('Total FM','lime'),('Spin down FM','magenta')]:
    pts=[r for r in rows if r[0]==s]
    ax.plot(px([r[1] for r in pts]), yN([r[2] for r in pts]), '.', color=col, ms=1.0, label=s)
ax.legend(fontsize=7,loc='upper right'); ax.set_title('Korzhavyi 2009 Fig.3(d): digitized data (data->pixel) over original')
fig.savefig(f'{outdir}/korzhavyi2009_fig3d_nef_overlay.png'); plt.close(fig)
