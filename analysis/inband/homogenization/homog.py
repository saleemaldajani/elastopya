#!/usr/bin/env python3
"""Effective elastic tensors of a coherent two-phase mixture of aligned cubic crystals.

All stiffnesses are 6x6 Voigt matrices (order 11,22,33,23,13,12; engineering shear strains), GPa,
unless a function says Mandel. Both phases share the same crystal axes (coherent spinodal decomposition
inside one bcc grain), so every scheme below is a homogenization of an *aligned* two-phase crystal.

Schemes
  voigt, reuss, hill                      : uniform strain / uniform stress / their mean
  laminate(Ca, Cb, w, axis)               : exact static (Backus/Postma) layered medium, welded interfaces
  cubic_project(C), cubic_project_reuss(C): average of the three <100> variants of a tensor (Voigt / Reuss)
  willis(Cs, ws, P, C0)                   : Hashin-Shtrikman/Willis estimate for a reference medium C0 and
                                            the polarization tensor P of the two-point statistics
  p_direction(C0, n), p_sphere(C0), p_100(C0): P for a laminate normal n, isotropic statistics, and statistics
                                            concentrated on the three <100> wave vectors (Cahn modulations)
  hs_bounds_cubic, self_consistent        : HS upper/lower bounds and self-consistent estimate
  fft_homogenize                          : Moulinec-Suquet FFT full-field homogenization on a voxel grid
"""
import numpy as np

VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
W = np.array([1, 1, 1, np.sqrt(2), np.sqrt(2), np.sqrt(2)])
I6 = np.eye(6)


# ------------------------------------------------------------------ basic tensor tools
def cubic(c11, c12, c44):
    C = np.zeros((6, 6))
    C[:3, :3] = c12
    np.fill_diagonal(C[:3, :3], c11)
    C[3, 3] = C[4, 4] = C[5, 5] = c44
    return C


def to_mandel(Cv):
    return Cv * np.outer(W, W)


def from_mandel(Cm):
    return Cm / np.outer(W, W)


def voigt_to_full(Cv):
    C = np.zeros((3, 3, 3, 3))
    for I, (i, j) in enumerate(VOIGT):
        for J, (k, l) in enumerate(VOIGT):
            for a, b in ((i, j), (j, i)):
                for c, d in ((k, l), (l, k)):
                    C[a, b, c, d] = Cv[I, J]
    return C


def full_to_voigt(C):
    return np.array([[C[i, j, k, l] for (k, l) in VOIGT] for (i, j) in VOIGT])


def rotate(Cv, R):
    """C'_{pqrs} = R_pi R_qj R_rk R_sl C_ijkl."""
    return full_to_voigt(np.einsum("pi,qj,rk,sl,ijkl->pqrs", R, R, R, R, voigt_to_full(Cv)))


def cubic_constants(Cv):
    """(C11, C12, C44) of the cubic projection (= average over the three axis permutations)."""
    return (np.trace(Cv[:3, :3]) / 3, (Cv[0, 1] + Cv[0, 2] + Cv[1, 2]) / 3, np.trace(Cv[3:, 3:]) / 3)


def zener(c11, c12, c44):
    return 2 * c44 / (c11 - c12)


# ------------------------------------------------------------------ simple bounds
def voigt(Ca, Cb, w):
    return (1 - w) * Ca + w * Cb


def reuss(Ca, Cb, w):
    return np.linalg.inv((1 - w) * np.linalg.inv(Ca) + w * np.linalg.inv(Cb))


def hill(Ca, Cb, w):
    return 0.5 * (voigt(Ca, Cb, w) + reuss(Ca, Cb, w))


# ------------------------------------------------------------------ exact laminate (Backus 1962 / Postma 1955)
def laminate(Ca, Cb, w, axis=2):
    """Layers normal to x_axis, volume fraction w of phase b, perfectly bonded (welded) interfaces.
    Continuous across interfaces: in-plane strains (p) and normal tractions (n). With
      n = {axis-normal, two transverse shears}, p = the rest,
      C*_nn = <C_nn^-1>^-1,  C*_np = C*_nn <C_nn^-1 C_np>,
      C*_pp = <C_pp - C_pn C_nn^-1 C_np> + <C_pn C_nn^-1> C*_nn <C_nn^-1 C_np>.
    Exact in the quasi-static (layer thickness << wavelength) limit."""
    nidx = {0: [0, 4, 5], 1: [1, 3, 5], 2: [2, 3, 4]}[axis]
    pidx = [i for i in range(6) if i not in nidx]
    avg = lambda f: (1 - w) * f(Ca) + w * f(Cb)
    nn = lambda C: C[np.ix_(nidx, nidx)]
    npm = lambda C: C[np.ix_(nidx, pidx)]
    pn = lambda C: C[np.ix_(pidx, nidx)]
    pp = lambda C: C[np.ix_(pidx, pidx)]
    inv = np.linalg.inv
    Cnn = inv(avg(lambda C: inv(nn(C))))
    A = avg(lambda C: inv(nn(C)) @ npm(C))          # <Cnn^-1 Cnp>
    B = avg(lambda C: pn(C) @ inv(nn(C)))           # <Cpn Cnn^-1>
    Cpp = avg(lambda C: pp(C) - pn(C) @ inv(nn(C)) @ npm(C)) + B @ Cnn @ A
    out = np.zeros((6, 6))
    out[np.ix_(nidx, nidx)] = Cnn
    out[np.ix_(nidx, pidx)] = Cnn @ A
    out[np.ix_(pidx, nidx)] = B @ Cnn
    out[np.ix_(pidx, pidx)] = Cpp
    return 0.5 * (out + out.T)


PERMS = [np.eye(3), np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0.]]), np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0.]])]


def variants(Cv):
    """The three <100> variants of a tensor (cyclic permutations of the axes)."""
    return [rotate(Cv, R) for R in PERMS]


def cubic_project(Cv):
    """Equal-weight Voigt (iso-strain) average of the three <100> variants = cubic projection."""
    return cubic(*cubic_constants(Cv))


def cubic_project_reuss(Cv):
    """Equal-weight Reuss (iso-stress) average of the three variants."""
    S = np.mean([np.linalg.inv(V) for V in variants(Cv)], axis=0)
    return np.linalg.inv(S)


# ------------------------------------------------------------------ Hashin-Shtrikman / Willis
def _p_of_n(C0full, n):
    K = np.einsum("ijkl,j,l->ik", C0full, n, n)
    N = np.linalg.inv(K)
    P = 0.25 * (np.einsum("j,ik,l->ijkl", n, N, n) + np.einsum("i,jk,l->ijkl", n, N, n)
                + np.einsum("j,il,k->ijkl", n, N, n) + np.einsum("i,jl,k->ijkl", n, N, n))
    return P


def _full_to_mandel(T):
    return full_to_voigt(T) * np.outer(W, W)


def p_direction(C0v, n):
    """Polarization tensor (Mandel) of an infinite layer with unit normal n (laminate statistics)."""
    n = np.asarray(n, float) / np.linalg.norm(n)
    return _full_to_mandel(_p_of_n(voigt_to_full(C0v), n))


def p_100(C0v):
    """Statistics whose two-point spectrum lies on the three <100> wave vectors (e.g. c = f(x)+g(y)+h(z))."""
    return np.mean([p_direction(C0v, n) for n in np.eye(3)], axis=0)


def p_sphere(C0v, nth=48, nph=96):
    """Isotropic (spherical) two-point statistics: P = (1/4pi) \\int P(n) dOmega, Gauss-Legendre x uniform."""
    C0 = voigt_to_full(C0v)
    mu, wt = np.polynomial.legendre.leggauss(nth)
    phis = (np.arange(nph) + 0.5) * 2 * np.pi / nph
    P = np.zeros((3, 3, 3, 3))
    for m, wm in zip(mu, wt):
        s = np.sqrt(1 - m * m)
        for ph in phis:
            P += wm * _p_of_n(C0, np.array([s * np.cos(ph), s * np.sin(ph), m]))
    P *= (2 * np.pi / nph) / (4 * np.pi)
    return _full_to_mandel(P)


def willis(Cs, ws, Pm, C0v):
    """C* = <C_r A_r><A_r>^-1, A_r = [I + P (C_r - C0)]^-1 (all in Mandel). Returns Voigt."""
    C0m = to_mandel(C0v)
    num = np.zeros((6, 6)); den = np.zeros((6, 6))
    for C, w in zip(Cs, ws):
        Cm = to_mandel(C)
        A = np.linalg.inv(I6 + Pm @ (Cm - C0m))
        num += w * Cm @ A
        den += w * A
    Cm = num @ np.linalg.inv(den)
    return from_mandel(0.5 * (Cm + Cm.T))


def _cubic_eig(C):
    c11, c12, c44 = cubic_constants(C)
    return np.array([c11 + 2 * c12, c11 - c12, c44])   # 3K, 2C', C44 (cubic eigen-moduli)


def _cubic_from_eig(e):
    k3, c2, c44 = e
    return cubic((k3 + 2 * c2) / 3, (k3 - c2) / 3, c44)


def hs_bounds_cubic(Ca, Cb, w, pfun=p_sphere):
    """HS bounds for aligned cubic phases: aligned cubic tensors share eigenspaces, so the reference
    C0 = cubic(max eigen-moduli) satisfies C0 >= C_r (upper bound) and min eigen-moduli gives C0 <= C_r (lower)."""
    ea, eb = _cubic_eig(Ca), _cubic_eig(Cb)
    out = []
    for e0 in (np.maximum(ea, eb), np.minimum(ea, eb)):
        C0 = _cubic_from_eig(e0)
        out.append(willis([Ca, Cb], [1 - w, w], pfun(C0), C0))
    return out  # [upper, lower]


def self_consistent(Ca, Cb, w, pfun=p_sphere, tol=1e-9, itmax=200):
    C0 = voigt(Ca, Cb, w)
    for _ in range(itmax):
        Cn = willis([Ca, Cb], [1 - w, w], pfun(C0), C0)
        Cn = cubic_project(Cn) if pfun is not None else Cn
        if np.max(np.abs(Cn - C0)) < tol:
            return Cn
        C0 = Cn
    return C0


# ------------------------------------------------------------------ FFT full-field homogenization
def fft_homogenize(c11, c12, c44, tol=1e-7, itmax=2000, verbose=False):
    """Moulinec-Suquet basic scheme. c11, c12, c44: 3D arrays (GPa), locally cubic with common axes.
    Periodic cell, unit voxels. Returns the 6x6 effective Voigt stiffness (columns from 6 macroscopic strains).
    Odd grid sizes avoid the Nyquist-frequency ambiguity."""
    shape = c11.shape
    lam0 = 0.5 * (c12.min() + c12.max())
    # isotropic reference: mu0 from the extreme shear-type eigen-moduli
    smin = min(c44.min(), (0.5 * (c11 - c12)).min()); smax = max(c44.max(), (0.5 * (c11 - c12)).max())
    mu0 = 0.5 * (smin + smax)
    kmin = ((c11 + 2 * c12) / 3).min(); kmax = ((c11 + 2 * c12) / 3).max()
    lam0 = 0.5 * (kmin + kmax) - 2 * mu0 / 3
    freqs = [np.fft.fftfreq(n) for n in shape]
    k = np.stack(np.meshgrid(*freqs, indexing="ij"))            # (3, ...)
    k2 = (k ** 2).sum(0); k2[0, 0, 0] = 1.0
    nvec = k / np.sqrt(k2)
    a1, a2 = 1.0 / mu0, 1.0 / mu0 - 1.0 / (lam0 + 2 * mu0)

    def stress(e):  # e: (3,3,...) symmetric strain
        tr = e[0, 0] + e[1, 1] + e[2, 2]
        s = 2 * c44 * e
        for i in range(3):
            s[i, i] = c12 * tr + (c11 - c12) * e[i, i]
        return s

    def gamma(t_hat):  # Gamma0 : tau in Fourier space
        v = np.einsum("ij...,j...->i...", t_hat, nvec)
        nv = np.einsum("i...,i...->...", nvec, v)
        wv = a1 * v - a2 * nvec * nv
        g = 0.5 * (nvec[:, None] * wv[None, :] + nvec[None, :] * wv[:, None])
        g[:, :, 0, 0, 0] = 0
        return g

    Ceff = np.zeros((6, 6))
    for J, (p, q) in enumerate(VOIGT):
        E = np.zeros((3, 3))
        E[p, q] = E[q, p] = 1.0 if p == q else 0.5           # engineering shear strain gamma = 1
        e = np.broadcast_to(E[:, :, None, None, None], (3, 3) + shape).copy()
        for it in range(itmax):
            s = stress(e)
            s_hat = np.fft.fftn(s, axes=(2, 3, 4))
            # equilibrium error
            div = np.einsum("ij...,j...->i...", s_hat, k)
            err = np.sqrt((np.abs(div) ** 2).sum() / k2.size) / np.abs(s_hat[:, :, 0, 0, 0]).max()
            if err < tol:
                break
            e_hat = np.fft.fftn(e, axes=(2, 3, 4)) - gamma(s_hat)
            e_hat[:, :, 0, 0, 0] = E * k2.size
            e = np.fft.ifftn(e_hat, axes=(2, 3, 4)).real
        if verbose:
            print(f"  load {J}: {it} iterations, err {err:.1e}")
        smean = stress(e).mean(axis=(2, 3, 4))
        Ceff[:, J] = [smean[i, j] for (i, j) in VOIGT]
    return 0.5 * (Ceff + Ceff.T)


# ------------------------------------------------------------------ isotropic averages and Rayleigh estimate
def vrh(Cv):
    S = np.linalg.inv(Cv)
    kv = (Cv[0, 0] + Cv[1, 1] + Cv[2, 2] + 2 * (Cv[0, 1] + Cv[0, 2] + Cv[1, 2])) / 9
    gv = (Cv[0, 0] + Cv[1, 1] + Cv[2, 2] - (Cv[0, 1] + Cv[0, 2] + Cv[1, 2]) + 3 * (Cv[3, 3] + Cv[4, 4] + Cv[5, 5])) / 15
    kr = 1 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * (S[0, 1] + S[0, 2] + S[1, 2]))
    gr = 15 / (4 * (S[0, 0] + S[1, 1] + S[2, 2]) - 4 * (S[0, 1] + S[0, 2] + S[1, 2]) + 3 * (S[3, 3] + S[4, 4] + S[5, 5]))
    return 0.5 * (kv + kr), 0.5 * (gv + gr)


LAMBDA = 6.4e-6
RHO_FE, RHO_CR = 7870.0, 7200.0


def rho_of(x):
    return x * RHO_CR + (1 - x) * RHO_FE


def rayleigh_f(Cv, rho):
    """Isotropic Rayleigh frequency (MHz) at 6.4 um, same formula as scripts/make_figures.py rayleigh_f,
    with Hill B and G (identical to make_figures for cubic tensors)."""
    b, g = vrh(Cv)
    e = 9 * b * g / (3 * b + g)
    nu = (3 * b - e) / (6 * b)
    vr = (0.874 + 0.196 * nu - 0.042 * nu ** 2 - 0.055 * nu ** 3) * np.sqrt(e * 1e9 / (2 * (1 + nu) * rho))
    return vr / LAMBDA / 1e6
