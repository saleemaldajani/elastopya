#!/usr/bin/env python3
"""Effective tensors of decomposed Fe-Cr from the QE random-alloy tensors, and their 10^4-grain (P)SAW (2026-10-15).

Spinodal domains (nm) are far smaller than the 6.4 um SAW wavelength, so the wave sees an effective medium. Inside
the miscibility gap a decomposed alloy of composition x is taken as the lever-rule (Voigt) average of the alpha (0.14)
and alpha' (0.88) tensors (CALPHAD, 748 K), each interpolated linearly in composition from the QE set:
    C_dec(x) = (1 - w) C(0.14) + w C(0.88),   w = (x - 0.14) / (0.88 - 0.14).
Check against MD, where the decomposed cell exists: the same construction from MD random tensors reproduces the MD
phi = 1 tensors (C11 within 0.1 GPa) and 80-90 % of their frequency shift (Fe-36Cr: 451 -> 469 vs 473 MHz).
Writes data/saw/qeaexpdec_x*.csv and qedec_x*.csv (new files only) and landscape_dec_tensors.csv.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402

SAW = mf.DATA / "saw"
RUNS = HERE.parent / "dft_qe" / "engaging" / "runs"
XA, XB = 0.14, 0.88


def series(suf):
    rows = {}
    pats = [f"sqs40_x*_k442_e60{suf}"] + (["sqs40_x0.80_k442_e60_aexp_crafm"] if suf else ["sqs40_x0.80_k442_e60_crafm2"])
    for pat in pats:
        for d in sorted(RUNS.glob(pat)):
            if not suf and d.name.split("e60")[1] and "crafm2" not in d.name:
                continue
            j = d / "result.json"
            if j.exists():
                r = json.loads(j.read_text()); rows[round(r["x_cr"], 3)] = (r["C11"], r["C12"], r["C44"])
    x = np.array(sorted(rows)); C = np.array([rows[v] for v in x])
    return x, C


def main():
    w_ = sys.argv[1] if len(sys.argv) > 1 else "10"
    lines = ["tag,x_cr,c11,c12,c44,rho,x_alpha,x_alpha_prime"]
    todo = []
    for suf, pre in (("_aexp", "qeaexpdec"), ("", "qedec")):
        x, C = series(suf)
        Ca, Cb = (np.array([np.interp(v, x, C[:, i]) for i in range(3)]) for v in (XA, XB))
        for x0 in x[(x > XA) & (x < XB)]:
            w = (x0 - XA) / (XB - XA); Cd = (1 - w) * Ca + w * Cb
            rho = x0 * mf.RHO_CR + (1 - x0) * mf.RHO_FE
            lines.append(f"{pre}_x{x0:.2f},{x0},{Cd[0]:.4f},{Cd[1]:.4f},{Cd[2]:.4f},{rho:.1f},{XA},{XB}")
            if not (SAW / f"{pre}_x{x0:.2f}.csv").exists():
                todo.append((pre, x0, Cd, rho))
    (SAW / "landscape_dec_tensors.csv").write_text("\n".join(lines) + "\n")
    print(len(todo), "runs to do", flush=True)
    for pre, x0, Cd, rho in todo:
        subprocess.run([sys.executable, str(HERE / "saw_random_orientations.py"), "run", *(f"{c}" for c in Cd), f"{rho}",
                        str(SAW / "angles_10k.csv"), str(SAW / f"{pre}_x{x0:.2f}.csv"), "10000", w_], check=True)
        print("done", pre, x0, flush=True)
    print("ALL_DONE", flush=True)


if __name__ == "__main__":
    main()
