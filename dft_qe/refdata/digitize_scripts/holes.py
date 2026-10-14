import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
def find_holes(ink, region):
    x0,x1,y0,y1 = region
    bg = ~ink
    lab,n = ndi.label(bg)
    sizes = ndi.sum(np.ones_like(lab), lab, range(1,n+1))
    big = np.argmax(sizes)+1
    out=[]
    for i,sl in enumerate(ndi.find_objects(lab)):
        if i+1==big: continue
        ys,xs=sl; h=ys.stop-ys.start; w=xs.stop-xs.start
        cy=(ys.start+ys.stop-1)/2; cx=(xs.start+xs.stop-1)/2
        if not (x0<=cx<=x1 and y0<=cy<=y1): continue
        comp = lab[sl]==i+1
        filled = ndi.binary_fill_holes(comp)
        area=comp.sum(); farea=filled.sum()
        has_dot = farea>area
        out.append(dict(cx=cx,cy=cy,w=w,h=h,area=int(area),fill=farea/(w*h),dot=bool(has_dot),
                        ys=ys.start,ye=ys.stop-1,xs=xs.start,xe=xs.stop-1))
    return out
