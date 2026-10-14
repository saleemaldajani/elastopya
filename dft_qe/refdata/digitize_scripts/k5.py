import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from tmatch import *
from PIL import Image
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, finaldir, thr = sys.argv[1], sys.argv[2], float(sys.argv[3])
a=np.asarray(Image.open(img).convert('L')); ink=a>128
Td=make_template(ink,567,59,11,11); Ts=make_template(ink,453,324,11,11)
sd=match(ink,Td); ss=match(ink,Ts)
off={'diamond':(0.5,0.5),'square':(0.0,0.0)}
panels={'a':dict(region=(112,795,14,537), x=(109.0,797.5,0,1.0), y=(10.5,540.0,200,-50)),
        'b':dict(region=(104,770,690,1201), x=(101.0,772.5,0,0.2), y=(687.5,1204.0,40,-10))}
rows=[]
for pn,P in panels.items():
    allp=[]
    for name,sc in [('diamond',sd),('square',ss)]:
        if pn=='b' and name=='square': continue
        for x,y,s in peaks(sc,thr,mind=4,region=P['region']):
            allp.append([name,x+off[name][0],y+off[name][1],s])
    allp.sort(key=lambda p:-p[3]); keep=[]
    for p in allp:
        if all((p[1]-q[1])**2+(p[2]-q[2])**2>5**2 for q in keep): keep.append(p)
    for n,x,y,s in keep:
        c=P['x'][2]+(x-P['x'][0])/(P['x'][1]-P['x'][0])*(P['x'][3]-P['x'][2])
        E=P['y'][2]+(y-P['y'][0])/(P['y'][1]-P['y'][0])*(P['y'][3]-P['y'][2])
        rows.append((pn,n,c,E,x,y,s))
pickle.dump(rows,open(finaldir+'/../out/k5.pkl','wb'))
fig,axs=plt.subplots(1,2,figsize=(18,12),dpi=100)
for ax,pn,(y0,y1) in zip(axs,('a','b'),((0,560),(680,1215))):
    ax.imshow(255-a,cmap='gray',interpolation='nearest'); ax.set_xlim(90,805); ax.set_ylim(y1,y0)
    for r in rows:
        if r[0]!=pn: continue
        ax.plot(r[4],r[5],'+' if r[6]>=0.6 else 'x',color='red' if r[1]=='diamond' else 'blue',ms=10,mew=1.3)
        ax.text(r[4]+3,r[5]-3,f'{r[6]:.2f}',fontsize=6,color='m')
fig.tight_layout(); fig.savefig(finaldir+'/../out/k5_probe.png')
from collections import Counter; print(Counter((r[0],r[1]) for r in rows))
for r in sorted(rows): print(r[0],r[1],round(r[2],4),round(r[3],1),round(r[6],2))
