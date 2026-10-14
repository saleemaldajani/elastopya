"""Historical (2026-10-14): rebuilt Fig. 1 with bigger panel-(d) cells while keeping the QE compositions of the
then-current figure, so that newly finished QE runs were not added as a side effect. Superseded by
`make_figures_kdense.py --only-fig1` with the full QE set; kept for traceability.
"""
# Rebuild the Fig. 1 versions with bigger panel-(d) cells only: keep the QE compositions of the current figure
# (20/35/50/65/95% Cr), so the newly finished QE runs are not added to panel (a) as a side effect.
import sys
import numpy as np
sys.path.insert(0, "/Users/saleem/MIT/research/Scripta/mike_final_touches/scripts")
import make_figures as mf
_orig = mf.qe_mixing
def qe_mixing(*a, **k):
    r = _orig(*a, **k)
    if r is None:
        return r
    keep = np.isin(np.round(r["x"], 3), [0.2, 0.35, 0.5, 0.65, 0.95])
    return {kk: v[keep] for kk, v in r.items()}
mf.qe_mixing = qe_mixing
import make_figures_kdense as mk
sys.argv = ["make_figures_kdense.py", "--only-fig1"]
mk.main()
