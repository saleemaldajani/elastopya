import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from tmatch import *
from PIL import Image
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, pk, finaldir = sys.argv[1:4]
a=np.asarray(Image.open(img).convert('L')); res=pickle.load(open(pk,'rb'))
X0,X1=98.0,797.4; Y3,Y0=11.5,334.5
cx=lambda x:(x-X0)/(X1-X0); my=lambda y:(Y0-y)*3/(Y0-Y3)
xc=lambda c:X0+np.asarray(c)*(X1-X0); ym=lambda m:Y0-np.asarray(m)*(Y0-Y3)/3
excl=[(0,108,0,640),(788,800,0,640),(0,800,0,22),(0,800,546,640),(625,660,60,90),(230,265,485,512),(260,350,215,242)]
allp=[]
for name in ('circle','square','diamond'):
    for x,y,s in peaks(res[name], 0.45, mind=3, region=(100,795,14,548)):
        if any(e[0]<=x<=e[1] and e[2]<=y<=e[3] for e in excl): continue
        nm=name
        if nm=='square' and y>535: continue          # bottom-frame tick false positives
        if nm=='square' and y<140 and x>240: nm='circle'   # 2nd Fe circle in the c=0.5 stack matched the square template better
        allp.append([nm,x,y,s,'auto'])
allp.sort(key=lambda p:-p[3]); keep=[]
for p in allp:
    if all((p[1]-q[1])**2+(p[2]-q[2])**2 > 4**2 for q in keep): keep.append(p)
# manual additions in crowded clusters (centre estimated from partially visible symbol interiors; score n/a)
manual=[('circle',176.0,84.6),('circle',162.8,86.8),('circle',188.8,86.7),('circle',214.2,86.7),('circle',229.3,99.1),('circle',214.2,95.8),
        ('circle',155.9,82.2),('circle',142.5,88.0),('circle',129.0,88.0),('circle',110.6,92.0),('circle',117.4,97.9),
        ('circle',273.2,102.0),('circle',273.2,109.2),
        ('square',228.7,156.9),('square',228.6,163.0),('square',214.1,149.8),
        ('diamond',188.1,504.2),('diamond',447.2,342.4),('diamond',184.7,456.0),('diamond',137.5,541.0),('diamond',157.5,541.5)]
for n,x,y in manual:
    if all((x-q[1])**2+(y-q[2])**2 > 4**2 for q in keep): keep.append([n,x,y,float(res[n][int(round(y)),int(round(x))]),'manual'])
lab={'circle':'Fe','diamond':'Cr','square':'average'}
rows=[]
for n,x,y,s,src in sorted(keep,key=lambda p:(p[0],p[1],p[2])):
    q = 'clear' if (src=='auto' and s>=0.6) else 'partial'
    rows.append((lab[n], n, cx(x), my(y), x, y, s, src, q))
fn=f'{finaldir}/klaver2006_fig4_moments.csv'
with open(fn,'w') as f:
    f.write("""# Magnitude of local magnetic moments on Fe atoms, Cr atoms and the average, vs Cr concentration, DIGITIZED from Fig. 4 (p. 094435-5)
# Source: T. P. C. Klaver, R. Drautz, M. W. Finnis, Phys. Rev. B 74, 094435 (2006), 'Magnetism and thermodynamics of defect-free Fe-Cr alloys'
# Caption: 'The magnitude of magnetic moments on atoms in FeCr as a function of concentration. The circles, diamonds, and squares indicate the
#   moments on Fe atoms, Cr atoms, and the average, respectively. The large variation in the magnitude of the Cr moments appears because systems
#   with both clustered and monatomically dispersed Cr atoms are represented. ...' In-plot labels: 'Fe', 'Cr', 'average' (no legend box).
#   Despite 'magnitude', the plotted values are SIGNED (Cr moments mostly negative = antiparallel to Fe).
# Method (Sec. II): VASP, PAW, GGA-PBE, 400 eV cutoff, Monkhorst-Pack k-points; supercells (2-atom B2 ... 3x3x3 and 4x4x4 bcc cells) with atomic
#   positions, cell shape AND volume fully relaxed, collinear magnetic ground state (several initial moment configurations tried).
#   Moments = spin moment integrated in spheres of radius 1.4 A around each site (Fe and Cr); text notes Cr moments are sensitive to this radius
#   and larger than with standard VASP radii. Each marker = one structure (per-structure Fe-site / Cr-site / all-atom average; the paper does not
#   say whether Fe/Cr values are averages over inequivalent sites). Ordered/SQS-like supercells, not random-alloy CPA.
# Units: c_Cr = Cr atomic fraction (0-1, as on the axis); moment_muB in mu_B.
# Digitization: embedded 1-bit raster (800x637 px, ~300 dpi) at native resolution. Calibration (pixel x,y): x: c=0 -> 98 (left frame),
#   c=0.2 -> 237.9, c=1 -> 797.4 (frame; ticks every 0.05 = 34.97 px, checked at 0.2..1.0); y: 3 -> 11.5 (top frame), 0 -> 334.5 (zero line),
#   -2 -> 549.8 (bottom frame); ticks every 0.2 mu_B = 21.5 px. 1 px = 0.00143 in c, 0.0093 mu_B.
#   Marker centres by binary template matching (templates cut from isolated circle/square/diamond markers), score = mean(+-1 agreement);
#   'auto' = score >= 0.45 peak; 'manual' = centre estimated by eye from a partly visible symbol interior in a crowded cluster.
# Completeness: c >= 0.24 markers are essentially complete (stacks at the same c resolved where interiors are visible). For c < 0.2 the Fe
#   circles, average squares and Cr diamonds overlap into dense clusters; many symbols there are fully hidden and NOT recovered. Treat c<0.2 as a
#   representative subset only.
# Uncertainty: 'clear' +/-0.01 mu_B, +/-0.002 in c; 'partial' +/-0.03 mu_B (+/-0.06 for manual cluster points), +/-0.004 in c.
# Columns: series, symbol, c_Cr, moment_muB, x_px, y_px, match_score, source, quality
series,symbol,c_Cr,moment_muB,x_px,y_px,match_score,source,quality
""")
    for r in rows: f.write(f'{r[0]},{r[1]},{r[2]:.4f},{r[3]:.3f},{r[4]:.1f},{r[5]:.1f},{r[6]:.2f},{r[7]},{r[8]}\n')
fig,axs=plt.subplots(2,1,figsize=(14,20),dpi=100)
for ax,(x0,x1,y0,y1) in zip(axs,[(90,800,0,560),(98,300,70,550)]):
    ax.imshow(255-a,cmap='gray',interpolation='nearest'); ax.set_xlim(x0,x1); ax.set_ylim(y1,y0)
    for s,c in [('Fe','red'),('average','blue'),('Cr','green')]:
        for qq,mk in [('clear','+'),('partial','x')]:
            pp=[r for r in rows if r[0]==s and r[8]==qq]
            ax.plot(xc([r[2] for r in pp]), ym([r[3] for r in pp]), mk, color=c, ms=9, mew=1.4, label=f'{s} ({qq})')
    ax.legend(fontsize=7,loc='upper right')
axs[0].set_title('Klaver 2006 Fig.4 digitized (data->pixel) over original; + clear, x partial/manual')
fig.tight_layout(); fig.savefig(f'{finaldir}/klaver2006_fig4_moments_overlay.png')
from collections import Counter; print(Counter((r[0],r[8]) for r in rows))
for s in ('Fe','Cr','average'):
    print(s, sorted([(round(r[2],3),round(r[3],2)) for r in rows if r[0]==s and r[2]>0.24]))
