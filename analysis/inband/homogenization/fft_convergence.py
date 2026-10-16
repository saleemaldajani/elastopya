#!/usr/bin/env python3
"""Resolution / cell-size / seed convergence of the FFT result for the <100> 'brick' geometry (QE Fe-36Cr phases)."""
import numpy as np
import fecr36_tensors as ft
import geometry as geo
import homog as h

m = ft.models()["qe"]; Ca, Cb = m["Ca"], m["Cb"]
print("n  k0  vox/wavelength  seeds -> cubic C11 C12 C44 (sd)  A  f_iso  | laminate_V3 ref:",
      np.round(h.cubic_constants(h.cubic_project(h.laminate(Ca, Cb, ft.W))), 2))
for n, k0 in ((31, 4.13), (45, 6), (63, 8.4), (45, 4), (45, 9), (63, 6)):
    res = []
    for s in range(3):
        mask = geo.spinodal_mask(n, ft.W, kind="axes", k0_cycles=k0, seed=s)
        fields = [np.where(mask, h.cubic_constants(Cb)[i], h.cubic_constants(Ca)[i]) for i in range(3)]
        res.append(h.cubic_constants(h.fft_homogenize(*fields, tol=1e-7)))
    r = np.array(res); c = r.mean(0)
    print(f"{n:3d} {k0:5.2f} {n/k0:6.1f}  3 -> {c[0]:.2f} {c[1]:.2f} {c[2]:.2f} ({r[:,0].std(ddof=1):.2f},{r[:,1].std(ddof=1):.2f})"
          f"  A {h.zener(*c):.3f}  f {h.rayleigh_f(h.cubic(*c), h.rho_of(0.36)):.1f}", flush=True)
