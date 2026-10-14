import numpy as np
def arc_score(mask, white, xc, yc, r_ring=8.75, r_in=4.5, a0=100, a1=230):
    ang = np.deg2rad(np.arange(a0, a1+1, 5))
    xs = xc + r_ring*np.cos(ang); ys = yc - r_ring*np.sin(ang)
    xi = xc + r_in*np.cos(ang); yi = yc - r_in*np.sin(ang)
    s1 = np.mean([mask[int(round(y)), int(round(x))] for x,y in zip(xs,ys)])
    s2 = np.mean([white[int(round(y)), int(round(x))] for x,y in zip(xi,yi)])
    return s1+s2
def fit_y(mask, white, xc, y_lo, y_hi, step=0.25, **kw):
    best=None
    for y in np.arange(y_lo, y_hi, step):
        s = arc_score(mask, white, xc, y, **kw)
        if best is None or s>best[1]: best=(y,s)
    return best
def fit_xy(mask, white, x_lo, x_hi, y_lo, y_hi, step=0.5, **kw):
    best=None
    for x in np.arange(x_lo,x_hi,step):
        for y in np.arange(y_lo,y_hi,step):
            s=arc_score(mask,white,x,y,**kw)
            if best is None or s>best[2]: best=(x,y,s)
    return best
