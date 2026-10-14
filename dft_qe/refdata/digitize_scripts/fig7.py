import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *; from curves import runs
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, finaldir = sys.argv[1:3]
a=load(img); g=a.mean(2)
from scipy import ndimage as ndi
D=ndi.gaussian_filter(255-g, 1.2)
def runlen(mask, y, x, axis):
    if not mask[y,x]: return 0
    if axis==0:
        s=y; e=y
        while s>0 and mask[s-1,x]: s-=1
        while e<mask.shape[0]-1 and mask[e+1,x]: e+=1
    else:
        s=x; e=x
        while s>0 and mask[y,s-1]: s-=1
        while e<mask.shape[1]-1 and mask[y,e+1]: e+=1
    return e-s+1
def ridges(x0,x1,y0,y1,thr=120,excl=()):
    m=(D>thr); pts=[]
    def ok(x,y): return not any(ex[0]<=x<=ex[1] and ex[2]<=y<=ex[3] for ex in excl)
    for x in range(x0,x1+1):
        col=D[y0:y1+1,x]
        for i in range(1,len(col)-1):
            if col[i]>thr and col[i]>=col[i-1] and col[i]>col[i+1]:
                y=y0+i
                if runlen(m,y,x,0)<=runlen(m,y,x,1) and ok(x,y):
                    # parabolic subpixel
                    d=(col[i-1]-col[i+1])/(2*(col[i-1]-2*col[i]+col[i+1])+1e-9)
                    pts.append((x, y+float(np.clip(d,-0.5,0.5)), 'col'))
    for y in range(y0,y1+1):
        row=D[y,x0:x1+1]
        for j in range(1,len(row)-1):
            if row[j]>thr and row[j]>=row[j-1] and row[j]>row[j+1]:
                x=x0+j
                if runlen(m,y,x,1)<runlen(m,y,x,0) and ok(x,y):
                    d=(row[j-1]-row[j+1])/(2*(row[j-1]-2*row[j]+row[j+1])+1e-9)
                    pts.append((x+float(np.clip(d,-0.5,0.5)), y, 'row'))
    return pts
guides={'a':{1:[(8,342),(30,338),(60,325),(90,300),(110,270),(125,230),(140,170),(155,100),(171,38)],
             2:[(8,342),(30,355),(60,385),(90,430),(110,475),(131,521)],
             3:[(8,480),(40,472),(70,450),(100,415),(130,375),(160,335),(190,302),(215,280),(240,268)],
             4:[(8,482),(60,482),(90,475),(110,465),(140,435),(170,390),(195,340),(215,300),(230,278),(246,268)]},
        'd':{1:[(8,734),(30,728),(60,712),(90,692),(120,665),(140,640),(152,630)],
             2:[(8,734),(40,760),(70,800),(90,840),(110,890),(130,940),(160,1010),(190,1060),(220,1095),(240,1110)],
             3:[(8,893),(40,880),(70,855),(90,838),(120,795),(150,750),(180,710),(210,680),(235,658),(246,652)],
             4:[(8,895),(40,898),(70,897),(90,893),(110,885),(140,850),(170,800),(200,745),(225,690),(240,660),(246,652)],
             5:[(155,1111),(170,1095),(183,1090),(200,1095),(215,1104),(230,1111)]}}
def segdist(p, poly):
    p=np.array(p); best=1e9
    for (ax,ay),(bx,by) in zip(poly[:-1],poly[1:]):
        A=np.array([ax,ay],float); Bv=np.array([bx,by],float); t=np.clip(np.dot(p-A,Bv-A)/np.dot(Bv-A,Bv-A),0,1)
        best=min(best, np.linalg.norm(p-(A+t*(Bv-A))))
    return best
cal={'a':dict(x=(5.0,248.5), y=np.polyfit([34.5,156.5,278.5,399.5,521.5],[0.2,0.1,0.0,-0.1,-0.2],1)),
     'd':dict(x=(5.0,248.5), y=np.polyfit([626.5,748.5,869.5,991.5,1113.0],[0.2,0.1,0.0,-0.1,-0.2],1))}
for k in cal: print(k, cal[k]['y'])
win={'a':(8,246,37,519,[(15,62,62,118),(8,20,150,163),(8,20,272,285),(236,246,272,285),(8,20,393,406),(236,246,150,163),(236,246,393,406),(48,59,37,48),(97,108,37,48),(145,156,37,48),(194,205,37,48),(48,59,507,519),(97,108,507,519),(145,156,507,519),(194,205,507,519)]),
     'd':(8,246,629,1111,[(8,50,1045,1104),(8,20,742,755),(8,20,863,876),(8,20,985,998),(236,246,742,755),(236,246,863,876),(236,246,985,998),(48,59,629,640),(97,108,629,640),(145,156,629,640),(194,205,629,640),(48,59,1099,1111),(97,108,1099,1111),(145,156,1099,1111),(194,205,1099,1111)])}
out={}
for pn in ('a','d'):
    x0,x1,y0,y1,ex=win[pn]
    pts=ridges(x0,x1,y0,y1,excl=ex)
    rows=[]
    for x,y,kind in pts:
        ds=sorted((segdist((x,y),poly),b) for b,poly in guides[pn].items())
        br=ds[0][1] if ds[0][0]<15 else -1
        amb = 1 if (len(ds)>1 and ds[1][0]<6) else 0
        if pn=='d' and br==5 and x>=198:
            br = 2 if segdist((x,y),guides['d'][2])<5 else -1
        xx=(x-cal[pn]['x'][0])/(cal[pn]['x'][1]-cal[pn]['x'][0]); E=np.polyval(cal[pn]['y'],y)
        rows.append((br,xx,E,x,y,kind,amb))
    out[pn]=rows
    from collections import Counter; print(pn, len(rows), Counter(r[0] for r in rows), 'ambiguous', sum(r[6] for r in rows))
pickle.dump(out, open(finaldir+'/../out/fig7.pkl','wb'))
fn=finaldir+'/korzhavyi2009_fig7_bands_fe.csv'
with open(fn,'w') as f:
    f.write('''# Electronic spectrum (Bloch spectral function maxima, 'band lines') of bcc Fe (0% Cr) along [x00] (Gamma -> H), FM state, DIGITIZED from
#   Fig. 7 panels (a) '0% Cr spin Up' and (d) '0% Cr spin Down' (p. 054202-7)
# Source: P. A. Korzhavyi, A. V. Ruban, J. Odqvist, J.-O. Nilsson, B. Johansson, Phys. Rev. B 79, 054202 (2009)
# Caption: 'Evolution of the electronic spectrum of bcc Fe-Cr alloys with composition.' Grey-scale map of the calculated (EMTO-CPA) Bloch spectral
#   function; dark lines = bands. For 0% Cr this is pure bcc Fe (FM). Spin Up = majority Fe channel. Lattice parameter for Fig. 7 is not stated
#   in the caption or Sec. III C (Sec. III D's a0 = 2.875 A refers to the moment calculations only) - do not assume it.
# Axes: x = reduced wave vector along [x00] in units of 2pi/a (x=0 Gamma, x=1 H); E_minus_EF_Ry = E - E_F in Ry (panel range -0.2..0.2 Ry).
# Branch labels (assigned here, NOT from the paper; by continuity using hand-placed guide polylines, nearest guide within 15 px):
#   panel a (spin up): 1 = band rising steeply from -0.055 Ry at Gamma (leaves frame at E=+0.2 near x~0.68); 2 = band falling from -0.055 Ry at
#     Gamma (leaves frame at -0.2 near x~0.52); 3 = upper of the two bands from -0.165 Ry at Gamma to H at ~+0.01 Ry; 4 = lower one (flat to
#     x~0.3, then rising to H ~+0.01 Ry). Branches 1/2 are degenerate at Gamma, 3/4 at Gamma and at H.
#   panel d (spin down): 1 = band rising from +0.11 Ry at Gamma (leaves frame at +0.2 near x~0.6); 2 = band falling from +0.11 Ry at Gamma to
#     ~-0.2 Ry near H; 3 = upper of the two bands from -0.02 Ry at Gamma to H at ~+0.18 Ry; 4 = lower one (flat ~-0.02 to x~0.4, then rising to H);
#     5 = small band hump at the bottom edge (x~0.63-0.80, top ~-0.19 Ry; partly outside the frame) - low confidence.
#   branch -1 = point not assignable (bottom-edge feature for x>0.79 in panel d). near_crossing=1: another guide within 6 px (degeneracy/crossing),
#     assignment there is uncertain.
# Digitization: embedded raster (1000x1183 px JPEG; whole Fig. 7) at native resolution; darkness map Gaussian-smoothed (sigma 1.2 px);
#   ridge points = local darkness maxima along columns (where the band is flatter than 45 deg) and along rows (steeper bands), threshold
#   darkness > 120/255, sub-pixel parabolic refinement. Tick marks and the panel labels '(a)', '(d)' masked.
#   Calibration (pixel x,y in that image): panel a: x=0 -> 5.0, 0.2 -> 53.5, 0.4 -> 102.5, 0.6 -> 150.5, 0.8 -> 199.5, 1 -> 248.5;
#     E: +0.2 -> 34.5, +0.1 -> 156.5, 0 -> 278.5, -0.1 -> 399.5, -0.2 -> 521.5 (1 px = 0.00082 Ry).
#   panel d: same x; E: +0.2 -> 626.5, +0.1 -> 748.5, 0 -> 869.5, -0.1 -> 991.5, -0.2 -> 1113.0.
# Uncertainty: +/-0.002 Ry in E and +/-0.005 in x on sharp bands; +/-0.005 Ry where bands are broad (flat parts near Gamma at -0.165 Ry (up)
#   and -0.02 Ry (down)) or near crossings.
# Columns: panel, spin, branch, x_2pi_over_a, E_minus_EF_Ry, x_px, y_px, ridge_dir, near_crossing
panel,spin,branch,x_2pi_over_a,E_minus_EF_Ry,x_px,y_px,ridge_dir,near_crossing
''')
    for pn in ('a','d'):
        for r in sorted(out[pn], key=lambda r:(r[0],r[1])):
            f.write(f"{pn},{'up' if pn=='a' else 'down'},{r[0]},{r[1]:.4f},{r[2]:.4f},{r[3]:.1f},{r[4]:.1f},{r[5]},{r[6]}\n")
fig,axs=plt.subplots(1,2,figsize=(14,14),dpi=110)
cols={1:'red',2:'blue',3:'magenta',4:'orange',5:'cyan',-1:'lime'}
for ax,(pn,(yy0,yy1)) in zip(axs,[('a',(30,526)),('d',(622,1118))]):
    ax.imshow(a.astype(np.uint8),interpolation='nearest'); ax.set_xlim(0,255); ax.set_ylim(yy1,yy0)
    for b,c in cols.items():
        r=[q for q in out[pn] if q[0]==b]
        # back-transform from data
        X=[cal[pn]['x'][0]+q[1]*(cal[pn]['x'][1]-cal[pn]['x'][0]) for q in r]; Y=[(q[2]-cal[pn]['y'][1])/cal[pn]['y'][0] for q in r]
        ax.plot(X,Y,'.',ms=2,color=c,label=f'branch {b}')
    ax.legend(fontsize=7,loc='upper left'); ax.set_title(f'Fig.7({pn}) digitized band points (data->pixel) over original')
fig.tight_layout(); fig.savefig(finaldir+'/korzhavyi2009_fig7_bands_fe_overlay.png')
