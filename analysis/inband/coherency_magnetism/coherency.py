#!/usr/bin/env python3
"""Part A: acoustoelastic change of the alpha / alpha' tensors under coherency strain, and its effect on the
decomposed Fe-36Cr effective tensor.

Geometry: coherent laminate with modulation (layer normal) along [001]. In-plane strains are set by the misfit
eta and force balance, the out-of-plane strain by sigma_33 = 0:
    eps_par(alpha') - eps_par(alpha) = -eta,   (1-w) Y_a eps_par(a) + w Y_a' eps_par(a') = 0,
    eps_perp = -2 C12/C11 eps_par,  Y = C11 + C12 - 2 C12^2/C11 (biaxial modulus of (001)).
Effective wave coefficients of a homogeneously pre-strained crystal in the current (spatial) frame, to first order
(Thurston & Brugger 1964; Thurston 1965), with rho = rho0 (1 - tr eps):
    a_ijkl = (1 - tr eps) C_ijkl + C_ijklmn eps_mn + eps_ip C_pjkl + eps_jp C_ipkl + eps_kp C_ijpl + eps_lp C_ijkp
             + delta_ik sigma_jl
C_ijklmn: Brugger TOECs (C111, C112, C123, C144, C155, C456).
Decomposed medium: Voigt average of a (and rho) over the two phases; the delta*sigma term then vanishes because the
mean stress is zero. Three ways of placing the laminate in a grain:
    'cubic'  : average over the three <100> variants (cubic again) -> 10k-grain SAW with the standard script
    'tetra'  : one variant, layer normal along crystal [001]       -> top grains with the general-tensor solver
Usage: coherency.py [table|run|refine]
"""
import sys
from itertools import product

import numpy as np

import common as cm
import toec_values as tv

X0 = 0.36
W = (X0 - cm.XA) / (cm.XB - cm.XA)
VO = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (1, 2): 3, (2, 1): 3, (0, 2): 4, (2, 0): 4, (0, 1): 5, (1, 0): 5}


def toec_voigt(c111, c112, c123, c144, c155, c456):
    """6x6x6 Brugger TOEC array for cubic symmetry (fully symmetric in its three Voigt indices)."""
    T = np.zeros((6, 6, 6))
    for I, J, K in product(range(6), repeat=3):
        s = sorted((I, J, K))
        n = sum(i < 3 for i in s)
        if n == 3:
            a = s
            if a[0] == a[1] == a[2]:
                T[I, J, K] = c111
            elif a[0] == a[1] or a[1] == a[2]:
                T[I, J, K] = c112
            else:
                T[I, J, K] = c123
        elif n == 1:  # one normal index i, two shear indices
            i, p, q = s[0], s[1], s[2]
            if p == q:
                T[I, J, K] = c144 if p - 3 == i else c155
        elif n == 0:
            if s == [3, 4, 5]:
                T[I, J, K] = c456
    return T


def toec_tensor(T):
    C6 = np.zeros((3,) * 6)
    for idx in product(range(3), repeat=6):
        C6[idx] = T[VO[idx[0:2]], VO[idx[2:4]], VO[idx[4:6]]]
    return C6


def cubic(c11, c12, c44):
    return cm.sd.cubic_cijkl(c11, c12, c44)


def strained(C2, C6, eps):
    """Spatial wave coefficients a_ijkl and stress for a pre-strain eps (3x3), all in GPa."""
    sig = np.einsum("ijkl,kl->ij", C2, eps)
    a = (1 - np.trace(eps)) * C2 + np.einsum("ijklmn,mn->ijkl", C6, eps)
    a += np.einsum("ip,pjkl->ijkl", eps, C2) + np.einsum("jp,ipkl->ijkl", eps, C2)
    a += np.einsum("kp,ijpl->ijkl", eps, C2) + np.einsum("lp,ijkp->ijkl", eps, C2)
    a_sym = a.copy()
    a += np.einsum("ik,jl->ijkl", np.eye(3), sig)
    return a, a_sym, sig


def biaxial(c11, c12):
    return c11 + c12 - 2 * c12 ** 2 / c11


def laminate_strains(Ca, Cb, eta, w=W):
    Ya, Yb = biaxial(Ca[0], Ca[1]), biaxial(Cb[0], Cb[1])
    ea = w * Yb * eta / ((1 - w) * Ya + w * Yb)
    eb = ea - eta
    E = []
    for C, e in ((Ca, ea), (Cb, eb)):
        E.append(np.diag([e, e, -2 * C[1] / C[0] * e]))
    return E


def rot_variants(a):
    """Average of a_ijkl over the three <100> laminate variants (cyclic axis permutations)."""
    out = np.zeros_like(a)
    for perm in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        P = np.eye(3)[list(perm)]
        out += np.einsum("pi,qj,rk,sl,ijkl->pqrs", P, P, P, P, a)
    return out / 3


def to_voigt(a):
    M = np.zeros((6, 6))
    for (i, j), I in VO.items():
        for (k, l), J in VO.items():
            M[I, J] = a[i, j, k, l]
    return M


def cubic_of(a):
    M = to_voigt(a)
    return np.array([M[[0, 1, 2], [0, 1, 2]].mean(), M[[0, 0, 1], [1, 2, 2]].mean(), M[[3, 4, 5], [3, 4, 5]].mean()])


def case(Ca, Cb, Ta, Tb, eta, w=W):
    """Returns dict with per-phase strains, per-phase changes and the decomposed effective tensors."""
    C2a, C2b = cubic(*Ca), cubic(*Cb)
    C6a, C6b = toec_tensor(toec_voigt(*Ta)), toec_tensor(toec_voigt(*Tb))
    ea, eb = laminate_strains(Ca, Cb, eta, w)
    aa, aa_s, sa = strained(C2a, C6a, ea)
    ab, ab_s, sb = strained(C2b, C6b, eb)
    ra, rb = cm.rho(cm.XA) * (1 - np.trace(ea)), cm.rho(cm.XB) * (1 - np.trace(eb))
    a_dec = (1 - w) * aa + w * ab               # sigma term cancels: (1-w) sa + w sb = 0
    a0 = (1 - w) * C2a + w * C2b
    r0 = (1 - w) * cm.rho(cm.XA) + w * cm.rho(cm.XB)
    r_dec = (1 - w) * ra + w * rb
    return dict(ea=ea, eb=eb, sa=sa, sb=sb, aa=aa_s, ab=ab_s, ra=ra, rb=rb, a_dec=a_dec, a0=a0, r0=r0, r_dec=r_dec,
                stress_balance=(1 - w) * sa + w * sb)
