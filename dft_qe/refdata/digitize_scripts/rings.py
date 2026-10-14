import numpy as np
def ring_detect(mask, r_out=10.5, r_in=7.0, r_core=4.5, thr=0.6, core_max=0.35, mind=8, region=None):
    """Find centres of open-circle markers in a boolean mask."""
    H,W = mask.shape
    R = int(np.ceil(r_out))+1
    yy,xx = np.mgrid[-R:R+1,-R:R+1]; rr = np.hypot(yy,xx)
    ann = (rr<=r_out)&(rr>=r_in); core = rr<=r_core
    m = mask.astype(float)
    from numpy.lib.stride_tricks import sliding_window_view
    pad = np.pad(m, R)
    win = sliding_window_view(pad, (2*R+1,2*R+1))
    sa = (win*ann).sum((-1,-2))/ann.sum()
    sc = (win*core).sum((-1,-2))/core.sum()
    score = sa - sc
    cand = (sa>thr)&(sc<core_max)
    if region is not None:
        x0,x1,y0,y1 = region; rm = np.zeros_like(cand); rm[y0:y1,x0:x1]=True; cand &= rm
    ys,xs = np.nonzero(cand)
    order = np.argsort(-score[ys,xs])
    picked=[]
    for i in order:
        y,x = ys[i],xs[i]
        if all((y-py)**2+(x-px)**2 > mind**2 for py,px,_ in picked):
            picked.append((y,x,score[y,x]))
    # subpixel: score-weighted centroid in 3x3
    out=[]
    for y,x,s in picked:
        ys2,xs2 = np.mgrid[y-1:y+2, x-1:x+2]
        w = np.clip(score[y-1:y+2,x-1:x+2]-score[y,x]+0.15,0,None)
        out.append(((w*xs2).sum()/w.sum(), (w*ys2).sum()/w.sum(), s))
    return out
