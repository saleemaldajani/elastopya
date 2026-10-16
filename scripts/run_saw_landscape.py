#!/usr/bin/env python3
"""10^4-grain (P)SAW speeds for the composition-resolved SAW landscape (SM figure and mechanism panel g, v10).

Writes, only where missing (existing files are never overwritten):
  data/saw/qeaexp_x{x}.csv   QE 40-atom SQS (and 2-atom pure elements) at the experimental lattice parameter
  data/saw/md5r_x{x}.csv     MD 25^3 cell, random alloy (phi = 0.5), tensor averaged over the five MD alloys
  data/saw/md5d_x{x}.csv     MD 25^3 cell, fully decomposed (phi = 1), same averaging
The older md_/mddec_ files (single MD alloy, 2026-10-12) are kept for the figures that use them.
Densities: linear between bcc Fe and Cr, as everywhere else.

  run_saw_landscape.py [workers]
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


def jobs():
    out = []
    for d in sorted(RUNS.glob("sqs40_x*_k442_e60_aexp")) + sorted(RUNS.glob("sqs40_x0.80_k442_e60_aexp_crafm")):
        j = d / "result.json"
        if j.exists():
            r = json.loads(j.read_text()); out.append(("qeaexp", r["x_cr"], r["C11"], r["C12"], r["C44"]))
    for x, n in ((0.0, "Fe_fm_k16x16x16_e60_aexp"), (1.0, "Cr_afm_k16x16x16_e60_aexp")):
        j = RUNS / n / "result.json"
        if j.exists():
            r = json.loads(j.read_text()); out.append(("qeaexp", x, r["C11"], r["C12"], r["C44"]))
    dd = mf.read_csv("md_decomposition.csv")
    for x in sorted(set(np.round(dd["x_cr"], 4))):
        for ph, pre in ((0.5, "md5r"), (1.0, "md5d")):
            k = np.isclose(dd["x_cr"], x) & np.isclose(dd["phi"], ph)
            out.append((pre, float(x), float(dd["c11"][k].mean()), float(dd["c12"][k].mean()), float(dd["c44"][k].mean())))
    return out


def main():
    w = sys.argv[1] if len(sys.argv) > 1 else "8"
    todo = [(p, x, c11, c12, c44) for p, x, c11, c12, c44 in jobs() if not (SAW / f"{p}_x{x:.2f}.csv").exists()]
    with open(SAW / "landscape_tensors.csv", "w") as fh:
        fh.write("tag,x_cr,c11,c12,c44,rho\n")
        for p, x, c11, c12, c44 in jobs():
            fh.write(f"{p}_x{x:.2f},{x},{c11:.4f},{c12:.4f},{c44:.4f},{x * mf.RHO_CR + (1 - x) * mf.RHO_FE:.1f}\n")
    print(len(todo), "runs to do", flush=True)
    for p, x, c11, c12, c44 in todo:
        rho = x * mf.RHO_CR + (1 - x) * mf.RHO_FE
        subprocess.run([sys.executable, str(HERE / "saw_random_orientations.py"), "run", f"{c11}", f"{c12}", f"{c44}", f"{rho}",
                        str(SAW / "angles_10k.csv"), str(SAW / f"{p}_x{x:.2f}.csv"), "10000", w], check=True)
        print("done", p, x, flush=True)
    print("ALL_DONE", flush=True)


if __name__ == "__main__":
    main()
