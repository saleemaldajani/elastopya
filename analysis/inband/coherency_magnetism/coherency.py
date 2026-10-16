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


# ---------------------------------------------------------------------------------------------------------------
def fe_sets():
    """alpha-phase TOEC sets: Powell-Skove C111/C112/C123 + C144+2C155 from measured dC44/dP, three splits."""
    out = {}
    for src, dcdp in tv.FE_DCDP.items():
        K = tv.combos_from_dcdp(tv.FE_C_GB, dcdp)
        for d3 in (-1000.0, 0.0, 1000.0):
            c155 = (K[2] - d3) / 3
            out[(src, d3)] = (tv.FE_PS["c111"], tv.FE_PS["c112"], tv.FE_PS["c123"], c155 + d3, c155, 0.0)
    return out


def cr_combos():
    cr_af = (cm.qe_elem("Cr_afm_k16x16x16_e60"), 2.8588542556, cm.qe_elem("Cr_afm_k16x16x16_e60_aexp"), 2.8846)
    cr_nm = (cm.qe_elem("Cr_nm_k16x16x16_e60"), 2.843459079, cm.qe_elem("Cr_nm_k16x16x16_e60_aexp"), 2.8846)
    return {"IXS (Rodiere 2026) + RT Cr": tv.combos_from_dcdp(np.array([350.0, 67.8, 100.8]), tv.cr_ixs_dcdp()),
            "QE AF Cr volume deriv.": tv.combos_from_volume(*cr_af),
            "QE NM Cr volume deriv.": tv.combos_from_volume(*cr_nm)}


def fe_combos():
    return {src: tv.combos_from_dcdp(tv.FE_C_GB, d) for src, d in tv.FE_DCDP.items()}


def sweep(Ca, Cb, eta):
    """All combinations: Fe source x Cr source x deviatoric splits. Returns list of dicts."""
    rows = []
    Kfe, Kcr = fe_combos(), cr_combos()
    d1_fe = tv.FE_PS["c111"] - tv.FE_PS["c112"]
    for (fs, K_a), (cs, K_b) in product(Kfe.items(), Kcr.items()):
        for d3a, d3b, d1b in product((-1000.0, 0.0, 1000.0), (-1000.0, 0.0, 1000.0), (d1_fe, 0.0)):
            Ta = (tv.FE_PS["c111"], tv.FE_PS["c112"], tv.FE_PS["c123"], *tv.toec_from_combos(K_a, d1_fe, d3a)[3:5], 0.0)
            Tb = tv.toec_from_combos(K_b, d1b, d3b)
            r = case(Ca, Cb, Ta, Tb, eta)
            rows.append(dict(fe=fs, cr=cs, d3a=d3a, d3b=d3b, d1b=d1b, Ta=Ta, Tb=Tb, **r))
    return rows


def describe(r, Ca, Cb):
    """Per-phase changes (cubic part, single variant shown as C11/C33, C12/C13, C44/C66) and decomposed changes."""
    out = {}
    for ph, a, C in (("alpha", r["aa"], Ca), ("alpha'", r["ab"], Cb)):
        M = to_voigt(a)
        out[ph] = dict(dC11=M[0, 0] - C[0], dC33=M[2, 2] - C[0], dC12=M[0, 1] - C[1], dC13=M[0, 2] - C[1],
                       dC44=M[3, 3] - C[2], dC66=M[5, 5] - C[2], cub=cubic_of(rot_variants(a)) - C)
    out["dec_cubic"] = cubic_of(rot_variants(r["a_dec"])) - cubic_of(r["a0"])
    Md = to_voigt(r["a_dec"] - r["a0"])
    out["dec_tetra"] = dict(dC11=Md[0, 0], dC33=Md[2, 2], dC12=Md[0, 1], dC13=Md[0, 2], dC44=Md[3, 3], dC66=Md[5, 5])
    return out
