import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from tmatch import *
from PIL import Image
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, pk, out, thr = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
a=np.asarray(Image.open(img).convert('L')); res=pickle.load(open(pk,'rb'))
region=(103,792,16,545)
allp=[]
for name,sc in res.items():
    for x,y,s in peaks(sc, thr, mind=4, region=region):
        allp.append([name,x,y,s])
allp.sort(key=lambda p:-p[3]); keep=[]
for p in allp:
    if all((p[1]-q[1])**2+(p[2]-q[2])**2 > 4.5**2 for q in keep): keep.append(p)
print(len(keep)); 
from collections import Counter; print(Counter(p[0] for p in keep))
pickle.dump(keep, open(out,'wb'))
fig,axs=plt.subplots(2,1,figsize=(16,22),dpi=90)
for ax,(x0,x1,y0,y1) in zip(axs,[(90,800,0,560),(95,300,60,250)]):
    ax.imshow(255-a,cmap='gray',interpolation='nearest'); ax.set_xlim(x0,x1); ax.set_ylim(y1,y0)
    for name,c in [('circle','red'),('square','blue'),('diamond','green')]:
        pp=[p for p in keep if p[0]==name]
        ax.plot([p[1] for p in pp],[p[2] for p in pp],'+',color=c,ms=8,mew=1.2)
        for p in pp:
            if x0<p[1]<x1 and y0<p[2]<y1: ax.text(p[1]+2,p[2]-2,f'{p[3]:.2f}',fontsize=5,color=c)
fig.tight_layout(); fig.savefig(out.replace('.pkl','.png'))
