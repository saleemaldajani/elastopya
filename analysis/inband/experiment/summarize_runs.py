#!/usr/bin/env python3
"""Per-grain highest (P)SAW frequency at a 6.4 um grating from saw_random_orientations.py output.

highest mode per grain = max(saw, psaw) if i_psaw/i_saw >= 2e-3 else saw   (Acta selection rule)
usage: summarize_runs.py run1.csv [run2.csv ...]
"""
import sys
import numpy as np

LAM = 6.4e-6
BAND = 609.0  # MHz, Acta Fig. 6 upper edge (304 SS tensor)

print(f"{'run':28s} {'n':>6s} {'fmin':>6s} {'fmed':>6s} {'fp99':>6s} {'fmax':>6s} {'%>609':>7s} {'%>627':>7s} {'SAWmax':>7s}")
for path in sys.argv[1:]:
    d = np.genfromtxt(path, delimiter=",", names=True)
    saw, psaw, isaw, ipsaw = d["saw_m_s"], d["psaw_m_s"], d["i_saw"], d["i_psaw"]
    ok = np.isfinite(psaw) & np.isfinite(ipsaw) & (ipsaw / isaw >= 2e-3)
    top = np.where(ok, np.fmax(saw, psaw), saw) / LAM / 1e6
    top = top[np.isfinite(top)]
    print(f"{path.split('/')[-1]:28s} {len(top):6d} {top.min():6.0f} {np.median(top):6.0f} "
          f"{np.percentile(top, 99):6.0f} {top.max():6.0f} {100 * np.mean(top > BAND):7.2f} "
          f"{100 * np.mean(top > 627):7.2f} {np.nanmax(saw) / LAM / 1e6:7.0f}")
