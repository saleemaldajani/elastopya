import numpy as np
from scipy import ndimage as ndi
def make_template(ink, cx, cy, hw, hh):
    p = ink[cy-hh:cy+hh+1, cx-hw:cx+hw+1].astype(float)
    return 2*p-1
def match(ink, T):
    I = 2*ink.astype(float)-1
    return ndi.correlate(I, T, mode='constant', cval=-1)/T.size
def peaks(score, thr, mind=5, region=None):
    mx = ndi.maximum_filter(score, size=2*mind+1)
    cand = (score==mx)&(score>thr)
    if region is not None:
        x0,x1,y0,y1=region; m=np.zeros_like(cand); m[y0:y1,x0:x1]=True; cand&=m
    ys,xs=np.nonzero(cand)
    out=[]
    for y,x in zip(ys,xs):
        # subpixel via 3x3 centroid of score
        w=score[y-1:y+2,x-1:x+2]; w=w-w.min()
        yy,xx=np.mgrid[y-1:y+2,x-1:x+2]
        out.append((float((w*xx).sum()/w.sum()), float((w*yy).sum()/w.sum()), float(score[y,x])))
    return out
