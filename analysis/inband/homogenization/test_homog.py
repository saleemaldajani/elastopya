#!/usr/bin/env python3
"""Self-checks of homog.py. Run: python3 test_homog.py"""
import numpy as np

import homog as h

Ca = h.cubic(237.7, 121.6, 106.7)     # QE-like alpha (A > 1)
Cb = h.cubic(360.6, 93.5, 91.8)       # QE-like alpha' (A < 1)
w = 0.30
ok = True


def check(name, cond, info=""):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {info}")


# 1. Willis estimate with the layer polarization tensor = exact laminate, for any reference medium
L = h.laminate(Ca, Cb, w, axis=2)
for C0 in (Ca, Cb, h.cubic(300, 100, 80), h.cubic(500, 150, 200)):
    Wl = h.willis([Ca, Cb], [1 - w, w], h.p_direction(C0, [0, 0, 1]), C0)
    check("Willis(layer P) == Backus laminate", np.max(np.abs(Wl - L)) < 1e-8, f"max diff {np.max(np.abs(Wl - L)):.1e}")

# 2. FFT on an explicit 1D layered cell reproduces Backus
n = 45; nb = int(round(w * n)); w_grid = nb / n
mask = np.zeros((3, 3, n), bool); mask[:, :, :nb] = True
f = lambda i: np.where(mask, h.cubic_constants(Cb)[i], h.cubic_constants(Ca)[i])
Cfft = h.fft_homogenize(f(0), f(1), f(2), tol=1e-9)
Lg = h.laminate(Ca, Cb, w_grid, axis=2)
check("FFT layered cell == Backus", np.max(np.abs(Cfft - Lg)) < 0.05, f"max diff {np.max(np.abs(Cfft - Lg)):.3f} GPa")

# 3. ordering Reuss <= scheme <= Voigt (Loewner)
V, R = h.voigt(Ca, Cb, w), h.reuss(Ca, Cb, w)
up, lo = h.hs_bounds_cubic(Ca, Cb, w)
sc = h.self_consistent(Ca, Cb, w)
for name, C in (("laminate", L), ("HS upper", up), ("HS lower", lo), ("SC sphere", sc), ("FFT", Cfft)):
    Ref = h.voigt(Ca, Cb, w_grid) if name == "FFT" else V
    Rr = h.reuss(Ca, Cb, w_grid) if name == "FFT" else R
    e1 = np.linalg.eigvalsh(h.to_mandel(Ref - C)).min(); e2 = np.linalg.eigvalsh(h.to_mandel(C - Rr)).min()
    check(f"Reuss <= {name} <= Voigt", e1 > -1e-6 and e2 > -1e-6, f"(min eig V-C {e1:.2e}, C-R {e2:.2e})")
e = np.linalg.eigvalsh(h.to_mandel(up - lo)).min()
check("HS lower <= HS upper", e > -1e-8, f"{e:.2e}")
e1 = np.linalg.eigvalsh(h.to_mandel(up - sc)).min(); e2 = np.linalg.eigvalsh(h.to_mandel(sc - lo)).min()
check("HS lower <= SC <= HS upper", e1 > -1e-6 and e2 > -1e-6)

# 4. isotropic phases: HS bounds reduce to the classical Hashin-Shtrikman K and G bounds
def iso(K, G):
    return h.cubic(K + 4 * G / 3, K - 2 * G / 3, G)
K1, G1, K2, G2, f2 = 160., 80., 200., 120., 0.4
up, lo = h.hs_bounds_cubic(iso(K1, G1), iso(K2, G2), f2)
Ku = K2 + (1 - f2) / (1 / (K1 - K2) + 3 * f2 / (3 * K2 + 4 * G2))
Gu = G2 + (1 - f2) / (1 / (G1 - G2) + 6 * f2 * (K2 + 2 * G2) / (5 * G2 * (3 * K2 + 4 * G2)))
Kl = K1 + f2 / (1 / (K2 - K1) + 3 * (1 - f2) / (3 * K1 + 4 * G1))
Gl = G1 + f2 / (1 / (G2 - G1) + 6 * (1 - f2) * (K1 + 2 * G1) / (5 * G1 * (3 * K1 + 4 * G1)))
Kup, Gup = (up[0, 0] + 2 * up[0, 1]) / 3, up[3, 3]
Klo, Glo = (lo[0, 0] + 2 * lo[0, 1]) / 3, lo[3, 3]
check("HS upper (isotropic) K,G", abs(Kup - Ku) < 1e-6 and abs(Gup - Gu) < 1e-6, f"{Kup:.4f}/{Ku:.4f} {Gup:.4f}/{Gu:.4f}")
check("HS lower (isotropic) K,G", abs(Klo - Kl) < 1e-6 and abs(Glo - Gl) < 1e-6, f"{Klo:.4f}/{Kl:.4f} {Glo:.4f}/{Gl:.4f}")

# 5. VRH for cubic matches the make_figures formula
c11, c12, c44 = 257.5, 122.05, 113.0
b, g = h.vrh(h.cubic(c11, c12, c44))
gv = (c11 - c12 + 3 * c44) / 5; gr = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
check("VRH cubic", abs(b - (c11 + 2 * c12) / 3) < 1e-9 and abs(g - (gv + gr) / 2) < 1e-9)

# 6. rotation / variants: cubic tensor invariant; laminate variants are axis permutations
check("cubic invariant under variants", all(np.allclose(V_, Ca) for V_ in h.variants(Ca)))
Lx = h.laminate(Ca, Cb, w, axis=0)
check("laminate axis 0 is a variant of axis 2", any(np.allclose(V_, Lx) for V_ in h.variants(L)))
print("ALL PASS" if ok else "SOME CHECKS FAILED")
