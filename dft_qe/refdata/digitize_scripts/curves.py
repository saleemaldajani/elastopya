import numpy as np
def runs(col):
    r=[]; inrun=False
    for i,v in enumerate(col):
        if v and not inrun: s=i; inrun=True
        elif not v and inrun: r.append((s,i-1)); inrun=False
    if inrun: r.append((s,len(col)-1))
    return r
def track(mask, x_list, y_lo, y_hi, y_start=None, win=8, maxrun=12, prefer=None, slope_n=6, maxgap=40):
    """Column-wise tracking of a line with slope-extrapolated prediction; returns dict x->y centre."""
    out={}; hist=[]
    for x in x_list:
        rs = [(s+y_lo, e+y_lo) for s,e in runs(mask[y_lo:y_hi, x]) if e-s+1<=maxrun]
        if not rs: continue
        cs = [(s+e)/2 for s,e in rs]
        if not hist:
            if y_start is None: c = cs[0] if prefer!='bottom' else cs[-1]
            else:
                c = min(cs, key=lambda t: abs(t-y_start))
                if abs(c-y_start)>win: continue
        else:
            lx,ly = hist[-1]
            if abs(x-lx)>maxgap: pred=ly
            elif len(hist)>=3:
                hx=np.array([h[0] for h in hist[-slope_n:]]); hy=np.array([h[1] for h in hist[-slope_n:]])
                if np.ptp(hx)>0:
                    p=np.polyfit(hx,hy,1); pred=np.polyval(p,x)
                else: pred=ly
            else: pred=ly
            c = min(cs, key=lambda t: abs(t-pred))
            if abs(c-pred) > win: continue
        out[x]=c; hist.append((x,c))
    return out
