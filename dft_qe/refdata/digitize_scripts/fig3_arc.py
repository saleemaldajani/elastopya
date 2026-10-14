import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *; from arcfit import *
a = load(sys.argv[1]); R,G,B = a[...,0],a[...,1],a[...,2]
red = (R>150)&(G<120)&(B<120); blk=(a.mean(2)<110)&~red; white = a.mean(2)>170
X0,X100=99.,678.
px=lambda c: X0+c*(X100-X0)/100
# validation on full cores in panel (b)
for c,(x,y) in [(3,(116,483.5)),(13,(174,475.5)),(15,(185,471.5)),(23,(232,461))]:
    print('check', c, (x,y), fit_xy(blk, white, x-3, x+3, y-6, y+6, a0=90, a1=270))
for c in [5,6,7,8,9,11,50,75]:
    x=px(c)
    for (a0,a1) in [(90,200),(60,150)]:
        print(c, round(x,1), (a0,a1), fit_xy(blk, white, x-2.5, x+2.5, 470 if c<20 else 380, 500 if c<20 else 435, a0=a0, a1=a1))
