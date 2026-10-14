import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from tmatch import *
from PIL import Image
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, pk, out = sys.argv[1:4]; thr=float(sys.argv[4])
x0,x1,y0,y1 = map(int, sys.argv[5:9]); types=sys.argv[9].split(',')
a=np.asarray(Image.open(img).convert('L')); res=pickle.load(open(pk,'rb'))
allp=[]
for name in types:
    for x,y,s in peaks(res[name], thr, mind=3, region=(x0,x1,y0,y1)):
        allp.append([name,x,y,s])
allp.sort(key=lambda p:-p[3]); keep=[]
for p in allp:
    if all((p[1]-q[1])**2+(p[2]-q[2])**2 > 4**2 for q in keep): keep.append(p)
fig,ax=plt.subplots(figsize=(14,14*(y1-y0+20)/(x1-x0+20)),dpi=100)
ax.imshow(255-a,cmap='gray',interpolation='nearest'); ax.set_xlim(x0-10,x1+10); ax.set_ylim(y1+10,y0-10)
for i,p in enumerate(keep):
    c={'circle':'red','square':'blue','diamond':'green'}[p[0]]
    ax.plot(p[1],p[2],'+',color=c,ms=10,mew=1.5); ax.text(p[1]+1.5,p[2]-1.5,f'{i}:{p[3]:.2f}',fontsize=7,color=c)
ax.set_xticks(range(x0-x0%5,x1+10,5)); ax.set_yticks(range(y0-y0%5,y1+10,5)); ax.tick_params(labelsize=6); ax.grid(lw=0.3,color='orange',alpha=0.5)
fig.tight_layout(); fig.savefig(out)
for i,p in enumerate(keep): print(i,p[0],round(p[1],1),round(p[2],1),round(p[3],2))
