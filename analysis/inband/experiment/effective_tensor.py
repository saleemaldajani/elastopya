#!/usr/bin/env python3
"""Two-phase (alpha / alpha') effective cubic tensor for coherent, co-oriented spinodal domains.

Both phases are bcc and share the parent-grain orientation, so the Voigt (uniform strain) and
Reuss (uniform stress) averages of two cubic tensors in the same frame are again cubic.

usage: effective_tensor.py C11a C12a C44a C11b C12b C44b w      (w = volume fraction of phase b)
prints Voigt, Reuss, Hill cubic constants (GPa) and Zener A.
"""
import sys
import numpy as np


def voigt66(c11, c12, c44):
    C = np.zeros((6, 6))
    C[:3, :3] = c12
    np.fill_diagonal(C[:3, :3], c11)
    C[3, 3] = C[4, 4] = C[5, 5] = c44
    return C


def cubic_from66(C):
    return C[0, 0], C[0, 1], C[3, 3]


def mix(a, b, w):
    Ca, Cb = voigt66(*a), voigt66(*b)
    V = (1 - w) * Ca + w * Cb
    R = np.linalg.inv((1 - w) * np.linalg.inv(Ca) + w * np.linalg.inv(Cb))
    H = 0.5 * (V + R)
    return {k: cubic_from66(M) for k, M in (("voigt", V), ("reuss", R), ("hill", H))}


def zener(c11, c12, c44):
    return 2 * c44 / (c11 - c12)


if __name__ == "__main__":
    v = list(map(float, sys.argv[1:8]))
    res = mix(v[0:3], v[3:6], v[6])
    for k, (c11, c12, c44) in res.items():
        print(f"{k:6s} C11={c11:.2f} C12={c12:.2f} C44={c44:.2f} A={zener(c11, c12, c44):.3f}")
