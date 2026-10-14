import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
a = load(sys.argv[1]); g = a.mean(2)
dark = g < 100
H,W = g.shape
cs = dark.sum(0); rs = dark.sum(1)
print('cols with >30% dark:', [ (i,int(cs[i])) for i in range(W) if cs[i] > 0.3*H])
print('rows with >50% dark:', [ (i,int(rs[i])) for i in range(H) if rs[i] > 0.5*W])
