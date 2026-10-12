#!/usr/bin/env python3
"""(P)SAW speeds over random crystal orientations with the Du & Zhao method (saw_du.get_saw).

Usage:
  saw_random_orientations.py validate <angles.csv> <speeds.csv> [n]   # compare with stored MATLAB run
  saw_random_orientations.py run <C11> <C12> <C44> <rho> <angles.csv> <out.csv> [n] [workers]
C in GPa, rho in kg/m^3. Each orientation is evaluated at in-plane angle deg = 0 (as in the
10k-orientation MATLAB run behind Acta Fig. 6). Output: SAW and PSAW speed (m/s) per orientation.
"""
import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import saw_du as sd  # noqa: E402


def one(euler, C, rho):
    v, inten = sd.get_saw(C, rho, euler, 0, return_intensity=True)
    return np.r_[v, inten]


def speeds(c11, c12, c44, rho, angles, workers=6):
    C = sd.cubic_cijkl(c11 * 1e9, c12 * 1e9, c44 * 1e9)
    with Pool(workers) as p:
        return np.array(p.map(partial(one, C=C, rho=rho), list(angles), chunksize=8))


def main():
    mode = sys.argv[1]
    if mode == "validate":
        angles = np.loadtxt(sys.argv[2], delimiter=",")
        ref = np.loadtxt(sys.argv[3])
        n = int(sys.argv[4]) if len(sys.argv) > 4 else 120
        mine = speeds(209, 133, 121, 7950, angles[:n])[:, :2]
        for name, r in (("interleaved", ref.reshape(-1, 2)[:n]),
                        ("blocked", np.stack([ref[:len(ref) // 2], ref[len(ref) // 2:]], 1)[:n])):
            rel = np.abs(mine - r) / r
            print(f"{name:12s} SAW median {np.nanmedian(rel[:, 0]):.4f} max {np.nanmax(rel[:, 0]):.4f} | "
                  f"PSAW median {np.nanmedian(rel[:, 1]):.4f} max {np.nanmax(rel[:, 1]):.4f}")
        print(np.round(np.c_[mine[:5], ref.reshape(-1, 2)[:5]]))
    elif mode == "run":
        c11, c12, c44, rho = map(float, sys.argv[2:6])
        angles = np.loadtxt(sys.argv[6], delimiter=",")
        n = int(sys.argv[8]) if len(sys.argv) > 8 else len(angles)
        w = int(sys.argv[9]) if len(sys.argv) > 9 else 6
        out = speeds(c11, c12, c44, rho, angles[:n], w)
        np.savetxt(sys.argv[7], out, delimiter=",", header="saw_m_s,psaw_m_s,i_saw,i_psaw", comments="")
        f = out[:, :2] / 6.4e-6 / 1e6
        print(f"n={n} f(6.4um) MHz: SAW {np.nanmin(f[:,0]):.0f}-{np.nanmax(f[:,0]):.0f}  "
              f"PSAW {np.nanmin(f[:,1]):.0f}-{np.nanmax(f[:,1]):.0f}  all max {np.nanmax(f):.0f}")


if __name__ == "__main__":
    main()
