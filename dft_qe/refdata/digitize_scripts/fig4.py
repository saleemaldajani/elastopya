import sys, os, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *; from viterbi import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, outdir, finaldir = sys.argv[1:4]
a = load(img); R_,G_,B_ = a[...,0],a[...,1],a[...,2]
grn=(G_>130)&(R_<160)&(B_<160)&(G_-R_>50)
blk=(a.mean(2)<110)&~grn
excl=np.zeros(blk.shape,bool)
excl[140:292,125:300]=True          # legend box
excl[60:120,240:400]=True           # 'spin up' label
excl[840:915,195:380]=True          # 'spin down' label
excl[:24,:]=True; excl[955:,:]=True # top/bottom frame ticks
excl[:, :97]=True; excl[:, 694:]=True
excl[0:470,97:111]=True; excl[512:,97:111]=True; excl[0:455,683:694]=True; excl[530:,683:694]=True  # left/right frame ticks
for xt in (160,227,295,361,429,496,563,630): excl[478:505, xt-3:xt+4]=True   # ticks on middle axis
excl[489:494,:]=True                # middle axis line (DOS=0)
blk&=~excl; grn&=~excl
xs=list(range(97,694))
# calibration
XT = np.array([227,361.5,496,630.5]); ET=np.array([-0.4,-0.2,0.0,0.2]); pe=np.polyfit(XT,ET,1)
Ex=lambda x: np.polyval(pe,x); xE=lambda e:(np.asarray(e)-pe[1])/pe[0]
YU=np.array([8,129,249.5,370,491.]); VU=np.array([40,30,20,10,0.]); pu=np.polyfit(YU,VU,1)
YD=np.array([491,612,732.5,853,973.5]); VD=np.array([0,10,20,30,40.]); pdn=np.polyfit(YD,VD,1)
print('E fit',pe,'resid',ET-Ex(XT)); print('up',pu,VU-np.polyval(pu,YU)); print('down',pdn,VD-np.polyval(pdn,YD))
def cluster(rr, gap=20):
    rr=sorted(rr); out=[]
    for r in rr:
        if out and r[0]-out[-1][1] <= gap: out[-1]=(out[-1][0], max(out[-1][1],r[1]))
        else: out.append(r)
    return out
names={('black','solid'):'Fe(98%) [Fe98Cr02]',('green','solid'):'Cr(2%) [Fe98Cr02]',('black','dashed'):'Fe(90%) [Fe90Cr10]',('green','dashed'):'Cr(10%) [Fe90Cr10]'}
rows=[]; curves={}
for cname,mask,other in [('black',blk,grn),('green',grn,blk)]:
    for panel,(y0,y1) in [('up',(10,489)),('down',(494,972))]:
        Rr=col_runs(mask, xs, y0, y1)
        if cname=='black':
            for x in range(493,500): Rr[x]=[r for r in Rr[x] if r[1]-r[0]>3]
        occl={x: other[y0:y1,x].any() for x in xs}
        sy = 488 if panel=='up' else 495
        path,c = viterbi_solid(Rr, xs, occl=occl, start_y=sy)
        sol={x:r for x,r in path}
        resid={}
        for x in xs:
            rr=[r for r in Rr[x] if not (x in sol and r==sol[x])]
            if cname=='black' and 493<=x<=499: rr=[]
            resid[x]=cluster(rr)
        dpath,dc = viterbi_solid(resid, xs, G=300, C_gap=2., C_y=0.5, C_dc=0.05)
        dsh={x:r for x,r in dpath}
        conv = (lambda y: np.polyval(pu,y)) if panel=='up' else (lambda y: np.polyval(pdn,y))
        for style,d in [('solid',sol),('dashed',dsh)]:
            pts=[]
            for x in sorted(d):
                s_,e_=d[x]; yc=(s_+e_)/2; ln=e_-s_+1
                pts.append((x,yc,ln))
                rows.append((names[(cname,style)], cname, style, panel, float(Ex(x)), float(conv(yc)), x, yc, ln))
            curves[(cname,style,panel)]=pts
        print(cname,panel,'solid',len(sol),'dashed',len(dsh))
pickle.dump({'rows':rows,'pe':pe,'pu':pu,'pd':pdn}, open(f'{outdir}/fig4.pkl','wb'))
# ---- CSV ----
fn=f'{finaldir}/korzhavyi2009_fig4_pdos_fm.csv'
with open(fn,'w') as f:
    f.write("""# Atom-projected (site-projected) electronic DOS of Fe and Cr in random bcc Fe98Cr02 and Fe90Cr10, FM state, DIGITIZED from Fig. 4 (p. 054202-5)
# Source: P. A. Korzhavyi, A. V. Ruban, J. Odqvist, J.-O. Nilsson, B. Johansson, Phys. Rev. B 79, 054202 (2009)
# Caption: 'Calculated atom-projected densities of electron states in random bcc alloys Fe98Cr02 and Fe90Cr10 in the FM state.'
# Method (Sec. II A): EMTO + CPA, LDA (Perdew-Wang) self-consistent densities. Lattice parameter used for Fig. 4 is NOT stated in the caption or
#   Sec. III C (only Sec. III D, the moments, states a fixed a0 = 2.875 A); do not assume it.
# Legend (line style -> series): black solid = 'Fe (98%)' (Fe site in Fe98Cr02); green solid = 'Cr (2%)' (Cr site in Fe98Cr02);
#   black dashed = 'Fe (90%)' (Fe site in Fe90Cr10); green dashed = 'Cr (10%)' (Cr site in Fe90Cr10).
#   Upper panel 'spin up' (majority channel of Fe); lower panel 'spin down', plotted downward in the paper; here stored as POSITIVE values with spin=down.
#   Per-atom DOS of the given species (the curves are NOT concentration-weighted: Cr(2%) peaks at ~32 Ry^-1).
# Units: E_minus_EF_Ry in Ry; DOS_states_per_Ry in states/Ry (axis label 'Density of states (Ry^-1)'), per atom of that species, per spin.
# Digitization: embedded raster (700x1047 px JPEG, ~300 dpi) at native resolution. Calibration (pixel x,y in that image):
#   x: E=-0.4 -> 227, -0.2 -> 361.5, 0 -> 496, 0.2 -> 630.5 (major ticks; linear, 1 px = 0.00149 Ry); frame spans -0.60..0.30 Ry.
#   y (spin up): 40 -> 8, 30 -> 129, 20 -> 249.5, 10 -> 370, 0 -> 491; y (spin down): 0 -> 491, 10 -> 612, 20 -> 732.5, 30 -> 853, 40 -> 973.5
#   (1 px = 0.083 Ry^-1). Lines separated by colour (black/green) and by style: the solid curve of each colour is the min-cost connected chain of
#   column runs (Viterbi; continuous line), the dashed curve is the remaining ink of that colour. One point per image column (y = centre of the
#   ink run in that column; for steep segments the run is long and its centre is used). 'run_px' = vertical ink extent in that column.
# Gaps / caveats: dashed series have no points in dash gaps and where the dashed line lies on top of / merged with the solid one of the same
#   colour; black curves have no points where hidden under green lines (notably spin-down black for E < ~-0.35 Ry, all curves < ~2 Ry^-1 there)
#   and at the E_F dotted line (x=493-499, black only). Points closer than ~4 px (0.3 Ry^-1) to another same-colour curve may be a blend of both.
#   Sharp peak tips (Fe spin-up ~36 at -0.05 Ry; Fe spin-down ~37 at +0.12 Ry; Cr(2%) spin-down ~31 at -0.14 Ry) are where the run is longest;
#   tip heights from the column-run centre underestimate the true maximum by up to ~1-2 Ry^-1 - see 'run_px' (tip = top/bottom end of the run).
# Uncertainty: +/-0.002 Ry in E; +/-0.3 Ry^-1 in DOS on smooth parts, up to +/-1.5 Ry^-1 on steep flanks/peaks and where curves overlap.
# Columns: series, colour, linestyle, spin, E_minus_EF_Ry, DOS_states_per_Ry, x_px, y_px, run_px
series,colour,linestyle,spin,E_minus_EF_Ry,DOS_states_per_Ry,x_px,y_px,run_px
""")
    for r in rows: f.write(f'{r[0]},{r[1]},{r[2]},{r[3]},{r[4]:.4f},{r[5]:.3f},{r[6]},{r[7]:.1f},{r[8]}\n')
# ---- overlay from CSV values (data -> pixel) ----
import csv
dat=[l for l in open(fn) if not l.startswith('#')]
rd=list(csv.DictReader(dat))
fig,axs=plt.subplots(1,2,figsize=(18,11),dpi=120,gridspec_kw={"width_ratios":[1,1]})
col={('black','solid'):'red',('black','dashed'):'blue',('green','solid'):'magenta',('green','dashed'):'cyan'}
for ax,(x0,x1,yl0,yl1) in zip(axs,[(85,700,980,0),(370,600,960,500)]):
    ax.imshow(a.astype(np.uint8), interpolation='nearest'); ax.set_xlim(x0,x1); ax.set_ylim(yl0,yl1)
    for key,cc in col.items():
        for sp in ('up','down'):
            pts=[r for r in rd if r['colour']==key[0] and r['linestyle']==key[1] and r['spin']==sp]
            xx=xE([float(r['E_minus_EF_Ry']) for r in pts]); v=np.array([float(r['DOS_states_per_Ry']) for r in pts])
            yy=(v-pu[1])/pu[0] if sp=='up' else (v-pdn[1])/pdn[0]
            ax.plot(xx,yy,'.',ms=1.5,color=cc,label=f'{names[key]}' if sp=='up' else None)
axs[0].legend(fontsize=8,loc='lower left'); axs[0].set_title('Korzhavyi 2009 Fig.4: digitized (from CSV, data->pixel) over original'); axs[1].set_title('zoom')
fig.tight_layout(); fig.savefig(f'{finaldir}/korzhavyi2009_fig4_pdos_fm_overlay.png')
from collections import Counter; print(Counter((r[0],r[3]) for r in rows))
