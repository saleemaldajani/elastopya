import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
from scipy import ndimage as ndi
a = load(sys.argv[1]); x0,x1,y0,y1 = map(int, sys.argv[2:6])
R,G,B = a[...,0],a[...,1],a[...,2]
sub = a[y0:y1,x0:x1]
white = sub.mean(2) > 160
lab,n = ndi.label(white)
for i,sl in enumerate(ndi.find_objects(lab)):
    ys,xs = sl
    h = ys.stop-ys.start; w = xs.stop-xs.start
    area = (lab[sl]==i+1).sum()
    if area<4 or h>16 or w>16: continue
    # ring colour: look at pixels in annulus around center
    cy = (ys.start+ys.stop-1)/2+y0; cx=(xs.start+xs.stop-1)/2+x0
    yy,xx = np.mgrid[-11:12,-11:12]; rr=np.hypot(yy,xx); ann=(rr>7)&(rr<10.5)
    patch = a[int(round(cy))-11:int(round(cy))+12, int(round(cx))-11:int(round(cx))+12]
    pr = patch[ann]
    redf = ((pr[:,0]>150)&(pr[:,1]<120)).mean(); blkf = (pr.mean(1)<110).mean()
    print(f'core cx={cx:.1f} cy={cy:.1f} w={w} h={h} area={area} red={redf:.2f} blk={blkf:.2f}')
