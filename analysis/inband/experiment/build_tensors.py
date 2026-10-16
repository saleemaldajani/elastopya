#!/usr/bin/env python3
"""Measured-input cubic tensors for random and decomposed Fe-Cr (writes tensors.csv here).

Measured inputs (room temperature, GPa):
  Fe  : C11 231.4, C12 134.6, C44 116.4  Rotter & Smith (1966) as tabulated in Ledbetter & Reed,
        J. Phys. Chem. Ref. Data 2, 531 (1973), Table 5 (2.314/1.346/1.164 x10^12 dyn/cm^2); rho 7870.
  Cr  : C11 350.0, C12 67.8, C44 100.8   Bolef & de Klerk, Phys. Rev. 129, 1063 (1963) (antiferromagnetic,
        T_N ~ 311 K; value relayed by the single-crystal literature search, not re-read here); rho 7200.
  Fe-36.3Cr single grain, homogenized 1000 C (Rae et al., arXiv:2601.06659, Table 3):
        C' = 56.7 +/- 0.3 GPa, C44 = 113 +/- 0.9 GPa  (second homogenization: 55.0, 112.8)
        aged 500 C 10 h: C' = 59.8 +/- 0.3, C44 = 111.3 +/- 0.7;  1 h: 58.2, 114.7.
  B is not measured by TGS; taken linear between Fe and Cr (B varies only 167 -> 162 GPa).

Schemes
  vegard : every constant linear in x between Fe and Cr (alloy-blind; Voigt mix == random by construction)
  anchored : B, C', C44 piecewise-linear through Fe, measured Fe-36Cr (homogenized), Cr
"""
import csv
from pathlib import Path

import numpy as np

from effective_tensor import mix, zener

HERE = Path(__file__).resolve().parent
FE = np.array([231.4, 134.6, 116.4]); CR = np.array([350.0, 67.8, 100.8])
RHO = lambda x: 7870 + (7200 - 7870) * x  # noqa: E731
XA, XB, X0 = 0.14, 0.88, 0.36
W = (X0 - XA) / (XB - XA)


def bcc(c):  # (C11, C12, C44) -> (B, C', C44)
    return np.array([(c[0] + 2 * c[1]) / 3, (c[0] - c[1]) / 2, c[2]])


def cub(p):  # (B, C', C44) -> (C11, C12, C44)
    return np.array([p[0] + 4 * p[1] / 3, p[0] - 2 * p[1] / 3, p[2]])


def vegard(x):
    return (1 - x) * FE + x * CR


B36 = (1 - X0) * bcc(FE)[0] + X0 * bcc(CR)[0]
M36 = {"hom": cub([B36, 56.7, 113.0]), "hom2": cub([B36, 55.0, 112.8]),
       "aged1h": cub([B36, 58.2, 114.7]), "aged10h": cub([B36, 59.8, 111.3])}
NODES = np.array([0.0, X0, 1.0])
PN = np.array([bcc(FE), bcc(M36["hom"]), bcc(CR)])


def anchored(x):
    return cub([np.interp(x, NODES, PN[:, i]) for i in range(3)])


rows = []
def add(tag, c, rho, note):  # noqa: E302
    rows.append(dict(tag=tag, c11=round(c[0], 3), c12=round(c[1], 3), c44=round(c[2], 3), rho=round(rho, 1),
                     cprime=round((c[0] - c[1]) / 2, 2), A=round(zener(*c), 3), note=note))


add("fe_pure", FE, 7870, "measured")
add("cr_pure", CR, 7200, "measured (AF)")
for k, c in M36.items():
    add(f"fe36_meas_{k}", c, RHO(X0), "Rae et al. C', C44 + linear B")
for name, f in (("vegard", vegard), ("anch", anchored)):
    a, b = f(XA), f(XB)
    add(f"{name}_alpha_x014", a, RHO(XA), name)
    add(f"{name}_alphap_x088", b, RHO(XB), name)
    for k, c in mix(a, b, W).items():
        add(f"{name}_decomp_{k}", np.array(c), RHO(X0), f"{name}, w={W:.4f}")
# partial decomposition matched to APT amplitude of Fe-36Cr 500C/10h (Zhou et al. 2013, via Rae et al.)
xa, xb = 0.258, 0.459; w2 = (X0 - xa) / (xb - xa)
for k, c in mix(anchored(xa), anchored(xb), w2).items():
    add(f"anch_partialAPT_{k}", np.array(c), RHO(X0), f"anchored, 25.8/45.9 at.%, w={w2:.3f}")

with open(HERE / "tensors.csv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
    wr.writeheader(); wr.writerows(rows)
for r in rows:
    print(f"{r['tag']:26s} {r['c11']:8.2f} {r['c12']:8.2f} {r['c44']:8.2f} rho={r['rho']:7.1f} C'={r['cprime']:6.2f} A={r['A']:.3f}")
