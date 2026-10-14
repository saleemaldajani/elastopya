import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *; from curves import runs
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, out = sys.argv[1:3]
a=load(img); g=a.mean(2)
def ridge_points(x0,x1,y0,y1,thr=140,maxlen=8):
    pts=[]
    sub=g[y0:y1+1,x0:x1+1]; D=255-sub
    for j in range(sub.shape[1]):
        for s,e in runs(sub[:,j]<thr):
            if e-s+1<=maxlen:
                w=D[s:e+1,j]; yc=(w*np.arange(s,e+1)).sum()/w.sum()
                pts.append((x0+j, y0+yc, 'col', e-s+1))
    for i in range(sub.shape[0]):
        for s,e in runs(sub[i,:]<thr):
            if e-s+1<=maxlen:
                w=D[i,s:e+1]; xc=(w*np.arange(s,e+1)).sum()/w.sum()
                pts.append((x0+xc, y0+i, 'row', e-s+1))
    return pts
res={}
for name,(x0,x1,y0,y1) in {'a':(8,246,37,519),'d':(8,246,629,1111)}.items():
    res[name]=ridge_points(x0,x1,y0,y1)
    print(name, len(res[name]))
pickle.dump(res, open(out,'wb'))
fig,axs=plt.subplots(1,2,figsize=(14,14),dpi=100)
for ax,(name,(y0,y1)) in zip(axs,[('a',(30,526)),('d',(622,1118))]):
    ax.imshow(a.astype(np.uint8),interpolation='nearest'); ax.set_xlim(0,255); ax.set_ylim(y1,y0)
    p=res[name]; ax.plot([q[0] for q in p if q[2]=='col'],[q[1] for q in p if q[2]=='col'],'r.',ms=1)
    ax.plot([q[0] for q in p if q[2]=='row'],[q[1] for q in p if q[2]=='row'],'b.',ms=1)
    ax.set_xticks(range(0,256,10)); ax.set_yticks(range(y0-y0%10,y1,10)); ax.tick_params(labelsize=5); ax.grid(lw=0.2)
fig.tight_layout(); fig.savefig(out.replace('.pkl','.png'))
