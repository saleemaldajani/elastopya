#!/usr/bin/env python3
"""Sensitivity of the FFT-on-MD-cell (continuous composition field) result to the smoothing radius used to define the
local composition (shells = 1: 9 atoms, 2: 15 atoms, 3: 27 atoms as in scripts/phase_compositions.py).
Reports f_iso and C11/C12/C44 errors vs MD phi = 1, five seeds per composition."""
import numpy as np
import validate_md as v
import homog as h

lines = ["x  MD_f  | shells: f (df)  C11 C12 C44 errors"]
err = {1: [], 2: [], 3: []}
for x in v.XS:
    ref = v.dec[(x, 1.0)]["mean"]; fr = h.rayleigh_f(h.cubic(*ref), h.rho_of(x))
    s = f"{x:.2f} {fr:6.1f} |"
    for sh in (1, 2, 3):
        _, cont, _, _ = v.fft_md(x, shells=sh)
        c = np.array(h.cubic_constants(cont)); f = h.rayleigh_f(h.cubic(*c), h.rho_of(x))
        err[sh].append(np.r_[f - fr, c - ref])
        s += f"  s{sh}: {f:6.1f} ({f-fr:+5.1f}) {c[0]-ref[0]:+5.1f} {c[1]-ref[1]:+5.1f} {c[2]-ref[2]:+5.1f}"
    print(s, flush=True); lines.append(s)
for sh, e in err.items():
    e = np.array(e)
    lines.append(f"shells={sh}: RMS df {np.sqrt(np.mean(e[:,0]**2)):.2f} MHz, RMS dC11 {np.sqrt(np.mean(e[:,1]**2)):.2f}, "
                 f"dC12 {np.sqrt(np.mean(e[:,2]**2)):.2f}, dC44 {np.sqrt(np.mean(e[:,3]**2)):.2f} GPa")
print("\n".join(lines[-3:]))
open("md_smoothing_sensitivity.txt", "w").write("\n".join(lines) + "\n")
