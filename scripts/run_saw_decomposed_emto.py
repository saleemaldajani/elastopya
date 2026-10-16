#!/usr/bin/env python3
"""Decomposed Fe-Cr from the EMTO-CPA random-alloy moduli (Razumovskiy et al. 2011, digitized C', C44, B), the same
effective lever-rule (Voigt) construction as run_saw_decomposed_qe.py, alpha 0.14 / alpha' 0.88 (CALPHAD, 748 K).
Writes data/saw/litdec_x*.csv (10^4 grains) and landscape_litdec_tensors.csv. 2026-10-15."""
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402

SAW = mf.DATA / "saw"; XA, XB = 0.14, 0.88


def emto():
    d = mf.read_csv("sheet_raw_razumovskiy.csv")
    x, cp, c44, b = d["x_cr"], d["cprime_gpa"], d["c44_gpa"], d["bulk_gpa"]
    o = np.argsort(x)
    return x[o], np.column_stack([b + 4 * cp / 3, b - 2 * cp / 3, c44])[o]


def main():
    w_ = sys.argv[1] if len(sys.argv) > 1 else "10"
    x, C = emto()
    Ca, Cb = (np.array([np.interp(v, x, C[:, i]) for i in range(3)]) for v in (XA, XB))
    lines = ["tag,x_cr,c11,c12,c44,rho"]
    for x0 in (0.20, 0.30, 0.36, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80):
        w = (x0 - XA) / (XB - XA); Cd = (1 - w) * Ca + w * Cb
        rho = x0 * mf.RHO_CR + (1 - x0) * mf.RHO_FE
        lines.append(f"litdec_x{x0:.2f},{x0},{Cd[0]:.4f},{Cd[1]:.4f},{Cd[2]:.4f},{rho:.1f}")
        out = SAW / f"litdec_x{x0:.2f}.csv"
        if not out.exists():
            subprocess.run([sys.executable, str(HERE / "saw_random_orientations.py"), "run", *(f"{c}" for c in Cd), f"{rho}",
                            str(SAW / "angles_10k.csv"), str(out), "10000", w_], check=True)
            print("done", x0, flush=True)
    (SAW / "landscape_litdec_tensors.csv").write_text("\n".join(lines) + "\n")
    print("ALL_DONE", flush=True)


if __name__ == "__main__":
    main()
