#!/usr/bin/env python3
"""Part B: alpha' tensors for alternative magnetic states and the decomposed Fe-36Cr Voigt tensors.

alpha  (x=0.14): QE aexp SQS set, linear interpolation (0.10-0.20).
alpha' (x=0.88) alternatives
  AF   : QE aexp set as used in the paper (interpolates the x=0.80 Cr-AF cell and x=0.95)
  NM1  : QE aexp set without the magnetically ordered x=0.80 cell (0.50 -> 0.95; x=0.95 cell has |m_Cr|<=0.27 muB)
  NM2  : QE aexp x=0.50 -> non-magnetic Cr at the experimental lattice (x=1)
  NM3  : AF + full Cr (NM - AF) shift at a_exp (upper bound of removing the Cr moment)
  MEAS*: measured tensors (literature), added in measured_tensors.py
Decomposed Fe-36Cr: C = (1-w) C_alpha + w C_alpha', w = (0.36-0.14)/(0.88-0.14); rho linear in x.
Usage: magnetic_tensors.py [run]   (writes tensors_magnetic.csv; 'run' also runs the 10k-grain SAW)
"""
import sys

import numpy as np

import common as cm

# Measured RT Cr (antiferromagnetic at 295 K, T_N = 311 K): Bolef & de Klerk, Phys. Rev. 129, 1063 (1963),
# as tabulated by Ledbetter & Kim (2001) (values supplied by the coordinator; see REPORT.md for verification status).
CR_MEAS = np.array([350.0, 67.8, 100.8])


def emto(x):
    """EMTO-CPA (Razumovskiy et al. 2011, digitised): C11 = B + 4C'/3, C12 = B - 2C'/3, C44."""
    d = np.loadtxt(cm.ROOT / "data" / "sheet_raw_razumovskiy.csv", delimiter=",", skiprows=1)
    xs, cp, c44, b = d.T
    xs = np.clip(xs, 0, 1)
    C = np.c_[b + 4 * cp / 3, b - 2 * cp / 3, c44]
    return np.array([np.interp(x, xs, C[:, i]) for i in range(3)])


X0 = 0.36
W = (X0 - cm.XA) / (cm.XB - cm.XA)


def build():
    q = cm.qe_aexp()
    cr_af, cr_nm = cm.qe_elem("Cr_afm_k16x16x16_e60_aexp"), cm.qe_elem("Cr_nm_k16x16x16_e60_aexp")
    qn = dict(q); qn[1.0] = tuple(cr_nm)
    alpha = cm.interp(q, cm.XA)
    ap = {
        "AF": cm.interp(q, cm.XB),
        "NM1": cm.interp(q, cm.XB, [x for x in sorted(q) if abs(x - 0.80) > 1e-6]),
        "NM2": cm.interp(qn, cm.XB, [0.5, 1.0]),
    }
    ap["NM3"] = ap["AF"] + (cr_nm - cr_af)
    ap["MEAS"] = CR_MEAS.copy()
    ap["MEASFE"] = CR_MEAS + emto(cm.XB) - emto(1.0)
    rows = [("alpha_qe", cm.XA, *alpha, cm.rho(cm.XA))]
    for k, C in ap.items():
        rows.append((f"aprime_{k}", cm.XB, *C, cm.rho(cm.XB)))
        Cd = (1 - W) * alpha + W * C
        rows.append((f"dec36_{k}", X0, *Cd, (1 - W) * cm.rho(cm.XA) + W * cm.rho(cm.XB)))
    return rows


def write(rows, name):
    with open(cm.HERE / name, "w") as f:
        f.write("tag,x_cr,c11,c12,c44,rho\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]:.3f},{r[2]:.4f},{r[3]:.4f},{r[4]:.4f},{r[5]:.2f}\n")


if __name__ == "__main__":
    rows = build()
    write(rows, "tensors_magnetic.csv")
    for r in rows:
        print(r[0], np.round(r[2:5], 2), round(r[5], 1))
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        for r in rows:
            cm.run10k(r[0], *r[2:6])
            print("done", r[0], cm.summary(cm.OUT / f"{r[0]}.csv"), flush=True)
