#!/usr/bin/env python3
"""Python port of the Du & Zhao (2017) surface Green's-function SAW/PSAW solver (getSAW.m).

For a cubic crystal with surface given by Euler angles and in-plane direction `deg`, the
surface displacement response G33 is evaluated over a slowness grid; resonances (sign flips of
Re G33) are the SAW and pseudo-SAW modes. This reproduces getSAW.m (MATLAB, provided with the
Du et al. code) line by line, vectorised over the slowness grid:

  * k = (k0, n k0), n = 1..sampling, lambda = 7 um, T from `sampling`, w = w0 (1 + 1e-6 i)
  * partial waves: roots p of det(F p^2 + i B p + N) = 0 with Re p > 0
  * traction-free surface, unit normal load -> amplitudes a, G33 = sum_r a_r A_r3
  * H_L_Peak: indices where diff(Re G33) < 0; SAW = last, PSAW = second-to-last
Velocities returned in m/s (SAW, PSAW).
"""

import numpy as np


def cubic_cijkl(c11, c12, c44):
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    if i == j and k == l:
                        C[i, j, k, l] = c11 if i == k else c12
                    elif (i == k and j == l) or (i == l and j == k):
                        C[i, j, k, l] = c44
    return C


def euler2matrix(a, b, r):
    Rza = np.array([[np.cos(a), np.sin(a), 0], [-np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rxb = np.array([[1, 0, 0], [0, np.cos(b), np.sin(b)], [0, -np.sin(b), np.cos(b)]])
    Rzr = np.array([[np.cos(r), np.sin(r), 0], [-np.sin(r), np.cos(r), 0], [0, 0, 1]])
    return Rzr @ Rxb @ Rza


def matrix2euler(g):
    """Matrix2Euler.m (returns radians)."""
    psid = np.degrees(np.arccos(g[2, 2]))
    s = np.sin(np.radians(psid))
    phy1 = np.degrees(np.arctan2(g[2, 0] / s, -g[2, 1] / s)) % 360
    phy2 = np.degrees(np.arctan2(g[0, 2] / s, g[1, 2] / s)) % 360
    return np.radians([phy1, psid, phy2])


def euler_100_surface(normal):
    """Euler angles for a surface whose normal lies in the (100) zone, as in the MATLAB database:
    columns [x' y' z'] = [cross([1 0 0], n), [1 0 0], n]."""
    n = np.asarray(normal, float); n /= np.linalg.norm(n)
    x = np.array([1.0, 0, 0])
    g = np.column_stack([np.cross(x, n), x, n])
    if abs(n[2] - 1) < 1e-12:
        return np.zeros(3)
    return matrix2euler(g)


def get_saw(C, rho, euler, deg, sampling=4000, psaw=True, return_intensity=False):
    T = {40000: 20e-14, 4000: 20e-13, 400: 20e-12}.get(sampling, 20e-12)
    lam = 7e-6
    k0 = 2 * np.pi / lam
    w0 = 2 * np.pi / T
    w = w0 + 1j * 1e-6 * w0
    d = np.radians(deg)
    MM = np.array([[np.cos(d), np.cos(np.pi / 2 - d), 0], [np.cos(np.pi / 2 + d), np.cos(d), 0], [0, 0, 1]])
    a = (euler2matrix(*euler) @ MM).T
    Cr = np.einsum("pi,qj,rk,sl,ijkl->pqrs", a, a, a, a, C)
    ny = np.arange(1, sampling + 1)
    kx = k0 * np.ones_like(ny, float)
    ky = k0 * ny
    K = np.stack([kx, ky], 1)                                   # (n,2)
    F = Cr[:, 2, 2, :]                                          # F(i,j) = C(i,3,3,j)
    Cu3 = Cr[:, :2, 2, :]                                       # C(i,u,3,l)
    C3v = Cr[:, 2, :2, :]                                       # C(i,3,v,l)
    B = -(np.einsum("nu,iul->nil", K, Cu3) + np.einsum("nv,ivl->nil", K, C3v))
    M = 1j * B
    Cuv = Cr[:, :2, :2, :]                                      # C(i,u,v,l)
    N = rho * w ** 2 * np.eye(3)[None] - np.einsum("nu,nv,iuvl->nil", K, K, Cuv)
    # roots of det(F p^2 + M p + N): quadratic eigenproblem, linearised (6x6 companion)
    Finv = np.linalg.inv(F)
    top = np.concatenate([np.zeros((len(ny), 3, 3)), np.broadcast_to(np.eye(3), (len(ny), 3, 3))], 2)
    bot = np.concatenate([-Finv[None] @ N, -Finv[None] @ M], 2)
    p_all = np.linalg.eigvals(np.concatenate([top, bot], 1))   # (n,6)
    # MATLAB: keep roots with real part > 0, in roots() order; take the first three
    order = np.argsort(-p_all.real, axis=1)
    p = np.take_along_axis(p_all, order[:, :3], 1)              # (n,3) Re p > 0
    S = F[None, None] * p[..., None, None] ** 2 + M[:, None] * p[..., None, None] + N[:, None]
    _, _, Vh = np.linalg.svd(S)
    A = Vh[..., -1, :].conj()                                    # null vector, (n,3 roots,3 comps)
    R = np.einsum("il,nr,nrl->nir", F, p, A)                   # R(i,r) = C(i,3,3,l) p_r A(r,l)
    I = np.einsum("iul,nu,nrl->nir", Cr[:, 2, :2, :], K, A)    # I(i,r) = C(i,3,u,l) k_u A(r,l)
    Comb = -R + 1j * I
    amp = np.linalg.solve(Comb, np.broadcast_to(np.array([0, 0, 1.0 + 0j]), (len(ny), 3))[..., None])[..., 0]
    G33 = np.einsum("nr,nr->n", amp, A[:, :, 2])
    y = G33.real
    drops = np.nonzero(np.diff(y) < 0)[0] + 1                   # MATLAB 1-based indices from find()
    if len(drops) == 0:
        return np.array([np.nan, np.nan])
    idx = [drops[-1], drops[-2]] if psaw and len(drops) > 1 else [drops[-1]]
    num = 1 + np.array(idx)
    slowness = num * k0 / w.real
    v = 1.0 / slowness
    v = v if len(v) == 2 else np.array([v[0], np.nan])
    if not return_intensity:
        return v
    # getSAW.m intensity block: drop height from the preceding local maximum (MATLAB 1-based)
    inten = []
    for jj in idx:
        loop = jj + 1                       # pos_ind (1-based) = YYnew + 1
        while loop - 2 >= 1 and not (y[loop - 3] < y[loop - 2]):
            loop -= 1
        neg = loop - 1
        inten.append(y[neg - 1] - y[jj])    # ynew(neg_ind) - ynew(pos_ind), converted to 0-based
    inten = inten + [np.nan] * (2 - len(inten))
    return v, np.array(inten)


SURFACES_100 = {  # name: normal (as in the MATLAB study, 19 surfaces from (001) to (010))
    "001": (0, 0, 1), "0111": (0, 1, 11), "016": (0, 1, 6), "014": (0, 1, 4), "038": (0, 3, 8),
    "012": (0, 1, 2), "035": (0, 3, 5), "034": (0, 3, 4), "056": (0, 5, 6), "011": (0, 1, 1),
    "054": (0, 5, 4), "032": (0, 3, 2), "053": (0, 5, 3), "021": (0, 2, 1), "031": (0, 3, 1),
    "041": (0, 4, 1), "061": (0, 6, 1), "0121": (0, 12, 1), "010": (0, 1, 0),
}


def spectrum(c11, c12, c44, rho, surfaces=SURFACES_100, degs=range(181), sampling=4000):
    """Return {surface: array (2, len(degs))} of (SAW, PSAW) speeds in m/s; C in Pa."""
    C = cubic_cijkl(c11, c12, c44)
    out = {}
    for name, n in surfaces.items():
        e = euler_100_surface(n)
        out[name] = np.array([get_saw(C, rho, e, dg, sampling) for dg in degs]).T
    return out
