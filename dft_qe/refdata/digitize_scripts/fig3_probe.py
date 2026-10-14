import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *; from rings import ring_detect
a = load(sys.argv[1]); R,G,B = a[...,0],a[...,1],a[...,2]
red = (R>150)&(G<120)&(B<120)
blk = (a.mean(2)<110)&~red
for name,m in [('black',blk),('red',red)]:
    pts = ring_detect(m, region=(90,690,330,646), thr=0.55)
    pts.sort()
    print(name, len(pts)); print([(round(x,1),round(y,1),round(s,2)) for x,y,s in pts])
