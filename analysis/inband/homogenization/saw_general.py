#!/usr/bin/env python3
"""(P)SAW speeds over the 10^4 random orientations for a GENERAL (any symmetry) stiffness tensor.

Same physics/settings as scripts/saw_random_orientations.py (Du & Zhao Green's-function solver, local copy saw_du.py,
in-plane angle 0 for each Euler triple), but the tensor is read from a 6x6 Voigt matrix (GPa) in a .npy file, so a
tetragonal (single-variant laminate) tensor is handled exactly: every grain is the tensor rotated by its Euler angles.

Usage: saw_general.py <C6x6.npy> <rho> <angles.csv> <out.csv> [n] [workers]
"""
import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import saw_du as sd  # noqa: E402  (local, unmodified copy of scripts/saw_du.py)
import homog as h  # noqa: E402


def one(euler, C, rho):
    v, inten = sd.get_saw(C, rho, euler, 0, return_intensity=True)
    return np.r_[v, inten]


def main():
    Cv = np.load(sys.argv[1]); rho = float(sys.argv[2])
    angles = np.loadtxt(sys.argv[3], delimiter=",")
    n = int(sys.argv[5]) if len(sys.argv) > 5 else len(angles)
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 8
    C = h.voigt_to_full(Cv) * 1e9
    with Pool(workers) as p:
        out = np.array(p.map(partial(one, C=C, rho=rho), list(angles[:n]), chunksize=8))
    np.savetxt(sys.argv[4], out, delimiter=",", header="saw_m_s,psaw_m_s,i_saw,i_psaw", comments="")
    f = out[:, :2] / 6.4e-6 / 1e6
    print(f"{sys.argv[4]}: n={n} SAW {np.nanmin(f[:,0]):.0f}-{np.nanmax(f[:,0]):.0f} PSAW {np.nanmin(f[:,1]):.0f}-"
          f"{np.nanmax(f[:,1]):.0f} MHz", flush=True)


if __name__ == "__main__":
    main()
