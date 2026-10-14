import sys, os, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *; from curves import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, outdir = sys.argv[1], sys.argv[2]
a = load(img); R,G,B = a[...,0],a[...,1],a[...,2]
red = (R>150)&(G<120)&(B<120); blk = (a.mean(2)<110)&~red
X0, X100 = 99.0, 678.0
cx = lambda x: (x-X0)*100/(X100-X0); px = lambda c: X0 + np.asarray(c)*(X100-X0)/100
yb = np.array([369,451,532,614.]); vb=np.array([250,200,150,100.]); pb = np.polyfit(yb, vb, 1)
By = lambda y: np.polyval(pb,y); yB = lambda v:(np.asarray(v)-pb[1])/pb[0]
yd = np.array([1048,1166,1284.]); vd=np.array([20,10,0.]); pd_ = np.polyfit(yd, vd, 1)
Ny = lambda y: np.polyval(pd_,y); yN=lambda v:(np.asarray(v)-pd_[1])/pd_[0]
rowsb=[]
mk = [
 ('Expt. RT FM', 'filled triangle', 99.0, 503.0, 'partly overlapped by open triangle; base row 513 x88-110'),
 ('Expt. RT FM', 'filled triangle', 108.5, 500.0,'partly hidden; base row 510 right edge x120'),
 ('Expt. RT FM', 'filled triangle', 131.0, 510.0,'apex merged with theory circles; base row 520 x120-142'),
 ('Expt. RT FM', 'filled triangle', 138.5, 533.5,'isolated; apex row 524 base row 543 x127-150'),
 ('Expt. RT FM', 'filled triangle', 161.0, 517.5,'upper part hidden by open triangle; base row 527 x150-172'),
 ('Expt. LT FM', 'open triangle',   99.0, 499.0, 'interior visible rows 500-505 either side of frame line; base ~509'),
 ('Expt. LT FM', 'open triangle',  161.5, 508.5, 'apex ~499 base 518; interior rows 507-514'),
 ('Expt. RT AFM','filled diamond', 678.0, 507.0, 'at c=100 on right frame; vertices rows 497-517'),
 ('Expt. LT AFM','open diamond',   678.0, 466.0, 'at c=100 on right frame; vertices rows 459-473'),
 ('Expt. RT PM', 'filled red square', 549.5, 472.5, 'red bbox x540-559 y464-481'),
 ('Expt. RT PM', 'filled red square', 615.0, 461.5, 'red bbox x606-624 y452-471'),
 ('Expt. RT PM', 'filled red square', 678.0, 449.5, 'red bbox x669-687 y440-459 (split by frame line); "paramagnetic Cr" estimate'),
]
for s,m,x,y,note in mk: rowsb.append((s,m,cx(x),By(y),x,y,'marker',note))
fm_pts = [(99.0,480.0,'circle c=0; core split by frame line'),(116.0,483.5,'full white core'),
          (174.0,475.5,'full white core'),(185.0,471.5,'full white core'),(196.5,469.0,'full white core'),(208.0,465.5,'partial core'),
          (220.0,462.5,'partial core'),(232.0,461.0,'full white core'),(388.5,428.2,'ring top edge row 417.5 +10.5 (lower half under DLM circle)'),
          (446.0,414.0,'full white core'),(533.0,395.0,'ring top edge row 384.5 +10.5 (lower half under DLM circle)'),
          (562.0,391.5,'full white core'),(619.7,381.0,'core split in two pieces'),(678.0,369.3,'core cut by frame; coincides with DLM')]
dlm_pts = [(99.5,497.0,'mostly hidden by triangles; from visible red arc (+/-3 px)'),(157.0,485.5,'full white core'),(214.5,471.5,'full white core'),
           (272.5,457.5,'full white core'),(388.5,431.5,'full white core'),(532.5,401.5,'full white core'),(678.0,369.3,'coincides with FM circle at c=100')]
for x,y,note in fm_pts: rowsb.append(('Theory FM','open circle (calc. point)',cx(x),By(y),x,y,'circle',note))
for x,y,note in dlm_pts: rowsb.append(('Theory DLM','open red circle (calc. point)',cx(x),By(y),x,y,'circle',note))
# dense FM region c~4-11.5: circles unresolvable; estimate curve as top envelope of black chain + 10.5 px (ring outer radius)
env=[]
for x in range(123, 152, 2):
    ys_ = np.nonzero(blk[455:500, x])[0]
    if len(ys_): env.append((x, 455+ys_[0]+10.5))
for x,y in env: rowsb.append(('Theory FM','envelope estimate (dense overlapping circles)',cx(x),By(y),x,y,'envelope','top edge of overlapping-circle chain + ring radius; +/-3 GPa'))
# ---- panel d ----
xs_all = list(range(101, 677))
d_red = track(red, xs_all, 968, 1282, y_start=1036, win=6, maxrun=12)
d_fm  = track(blk, list(range(102,676)), 1080, 1175, y_start=1112, win=8, maxrun=16)
d_sd  = track(blk, list(range(102,676)), 1195, 1240, y_start=1232, win=5, maxrun=6)
leg = [x for x in range(250,560) if blk[1180:1260,x].mean()>0.8]
lx0, lx1 = min(leg)-1, max(leg)+1
d_sd = {x:y for x,y in d_sd.items() if not (lx0<=x<=lx1)}
rowsd=[]
for name,d in [('Total DLM',d_red),('Total FM',d_fm),('Spin down FM',d_sd)]:
    for x in sorted(d): rowsd.append((name, cx(x), Ny(d[x]), x, d[x]))
pickle.dump({'b':rowsb,'d':rowsd,'pb':pb,'pd':pd_,'legend_x':(lx0,lx1)}, open(f'{outdir}/fig3.pkl','wb'))
print('legend cols', lx0, lx1, 'c', cx(lx0), cx(lx1))
for r in rowsb: print(r[0], r[1][:20], round(r[2],2), round(r[3],1))
from collections import Counter; print(Counter(r[0] for r in rowsd))
