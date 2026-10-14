import numpy as np
from curves import runs
def col_runs(mask, xs, y_lo, y_hi):
    out={}
    for x in xs:
        out[x]=[(s+y_lo,e+y_lo) for s,e in runs(mask[y_lo:y_hi,x])]
    return out
def ivdist(a,b):
    """vertical gap between intervals a=(s,e), b=(s,e); 0 if touching/overlapping"""
    if a[1] < b[0]-1: return b[0]-1-a[1]
    if b[1] < a[0]-1: return a[0]-1-b[1]
    return 0
def viterbi_solid(R, xs, occl=None, G=12, C_gap=30., C_y=5., C_dc=0.02, start_y=None, C_start=1.0):
    """R: dict x->list of runs. Find min-cost chain of runs, favouring connected runs in consecutive columns.
       occl: dict x->bool, True if the column has occluding ink (other colour) -> cheaper gap."""
    xs=[x for x in xs]; idx={x:i for i,x in enumerate(xs)}
    cost=[]; back=[]
    for k,x in enumerate(xs):
        cs=[]; bs=[]
        for r in R[x]:
            c=(r[0]+r[1])/2
            best=(C_gap*k + (C_start*abs(c-start_y) if start_y is not None else 0), None)  # start here
            for dk in range(1, G+1):
                if k-dk<0: break
                xp=xs[k-dk]
                for j,rp in enumerate(R[xp]):
                    if cost[k-dk][j] is None: continue
                    g = dk-1
                    if g>0 and occl is not None:
                        occ = sum(occl.get(xs[k-dd],False) for dd in range(1,dk))
                        gc = C_gap*(g-occ) + 3.0*occ
                    else: gc = C_gap*g
                    cc = cost[k-dk][j] + gc + C_y*ivdist(rp,r) + C_dc*abs(c-(rp[0]+rp[1])/2)
                    if cc<best[0]: best=(cc,(k-dk,j))
            cs.append(best[0]); bs.append(best[1])
        cost.append(cs); back.append(bs)
    # end: allow ending at any column with penalty for unvisited columns
    best=None
    for k in range(len(xs)):
        for j,c in enumerate(cost[k]):
            cc=c+C_gap*(len(xs)-1-k)
            if best is None or cc<best[0]: best=(cc,(k,j))
    path=[]; node=best[1]
    while node is not None:
        k,j=node; path.append((xs[k], R[xs[k]][j])); node=back[k][j]
    return path[::-1], best[0]
