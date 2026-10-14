import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
a = load(sys.argv[1]); g = a.mean(2); dark = g<110
mode = sys.argv[2]
if mode=='x':  # rows y0..y1, scan columns x0..x1 for dark fraction
    y0,y1,x0,x1 = map(int, sys.argv[3:7])
    f = dark[y0:y1, x0:x1].mean(0)
    print([ (x0+i, round(v,2)) for i,v in enumerate(f) if v>0.6])
else:
    x0,x1,y0,y1 = map(int, sys.argv[3:7])
    f = dark[y0:y1, x0:x1].mean(1)
    print([ (y0+i, round(v,2)) for i,v in enumerate(f) if v>0.6])
