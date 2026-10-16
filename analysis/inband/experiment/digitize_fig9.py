import numpy as np, sys
from PIL import Image
im=np.array(Image.open(sys.argv[1]).convert('RGB')).astype(int)
cols={'blue':(31,119,180),'orange':(255,127,14),'green':(44,160,44),'red':(214,39,40),'purple':(148,103,189)}
x2f=lambda x:(x-168.5)/1.046
y2c=lambda y:(787.5-y)/6.475
for name,c in cols.items():
    m=(np.abs(im-np.array(c)).sum(2)<60)
    m[:330,:]=False  # skip legend
    m[:, :170]=False; m[:,1213:]=False
    prof=[]
    for x in range(170,1213):
        r=np.where(m[:,x])[0]
        prof.append(y2c(r.min()) if len(r) else np.nan)
    prof=np.array(prof)
    # segments of near-constant height
    segs=[];start=0
    for i in range(1,len(prof)+1):
        if i==len(prof) or not (abs(prof[i]-prof[start])<0.6 or (np.isnan(prof[i]) and np.isnan(prof[start]))):
            if i-start>=8 and not np.isnan(prof[start]):
                segs.append((round(x2f(170+start)),round(x2f(170+i-1)),round(float(np.nanmedian(prof[start:i])),1)))
            start=i
    print(name,segs)
