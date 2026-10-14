import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, pk, out, cname, panel = sys.argv[1:6]
xr = tuple(map(int, sys.argv[6:8])) if len(sys.argv)>7 else (90,700)
a = load(img).astype(np.uint8); res=pickle.load(open(pk,'rb'))[(cname,panel)]
y0,y1 = (5,495) if panel=='up' else (487,978)
fig,ax=plt.subplots(figsize=(14,12),dpi=110)
ax.imshow(a, interpolation='nearest'); ax.set_xlim(*xr); ax.set_ylim(y1,y0)
for x,r in res['solid'].items(): ax.plot([x,x],[r[0],r[1]],'-',color='red',lw=0.8)
for x,rr in res['resid'].items():
    for r in rr: ax.plot([x,x],[r[0],r[1]],'-',color='blue',lw=0.8)
ax.set_title(f'{cname} {panel}: red=solid path runs, blue=residual runs'); fig.tight_layout(); fig.savefig(out)
