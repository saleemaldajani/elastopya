#!/usr/bin/env python3
"""Finite-displacement (phonopy) phonons of pure bcc Fe and Cr from the 40-atom (2x2x5 conventional) reference
supercells computed on Engaging (phonopy/<case>/phonopy_disp.yaml + d000.out, pw.x forces, same PAW, cutoffs,
smearing and k-mesh as the SQS elastic runs). Written because the Fe DFPT (ph.x) frequencies are ~60% too high while a
frozen-phonon check with identical settings agrees with neutron data.

  phonopy_bcc_bands.py <case> [<case> ...]       e.g. Fe_ref Cr_ref

For each case: force constants in the 40-atom cell, mapped onto the 1-atom bcc primitive cell; dispersion along
Gamma-H-N-Gamma-P-H (Cartesian, units 2 pi/a) and a 40^3-mesh DOS. Output in phonopy_derived/<case>/:
freq.gp (path distance in 2 pi/a, frequencies in cm^-1, as matdyn) and dos.dat (cm^-1, DOS), so the DFPT plotting
code reads them unchanged. Needs phonopy (run with a Python that has it).
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC, OUT = HERE / "engaging" / "phonopy", HERE / "phonopy_derived"
THZ_CM1 = 33.35641
S2, S3 = np.sqrt(2), np.sqrt(3)
PATH = [(0, 0, 0), (1, 0, 0), (0.5, 0.5, 0), (0, 0, 0), (0.5, 0.5, 0.5), (1, 0, 0)]   # G H N G P H, 2 pi/a
NSEG = 60


def forces(out_file):
    """Final forces (Ry/bohr -> eV/A) from a pw.x output."""
    txt = Path(out_file).read_text().splitlines()
    i = max(n for n, l in enumerate(txt) if "Forces acting on atoms" in l)
    f = []
    for l in txt[i + 2:]:
        if "force =" not in l:
            if f:
                break
            continue
        f.append([float(v) for v in l.split("=")[1].split()])
    return np.array(f) * 13.605693 / 0.529177


def run(case):
    import phonopy
    from phonopy import Phonopy
    d = SRC / case
    ph = phonopy.load(d / "phonopy_disp.yaml", produce_fc=False, log_level=0)
    n = len(ph.supercells_with_displacements)
    ph.forces = np.array([forces(d / f"d{i:03d}.out") for i in range(n)])
    ph.produce_force_constants(); ph.symmetrize_force_constants()
    uc = ph.unitcell
    L = np.array(uc.cell)                                   # rows: 2a, 2a, 5a (orthogonal)
    a = L[0, 0] / 2
    Ap = 0.5 * a * np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], float)   # rows: bcc primitive vectors
    P = np.linalg.inv(L.T) @ Ap.T                            # columns: primitive vectors in unit-cell coordinates
    pp = Phonopy(uc, supercell_matrix=np.eye(3, dtype=int), primitive_matrix=P, log_level=0)
    pp.force_constants = ph.force_constants
    qc, dist, x = [], [], 0.0
    for p0, p1 in zip(PATH[:-1], PATH[1:]):
        p0, p1 = np.array(p0, float), np.array(p1, float)
        for t in np.linspace(0, 1, NSEG, endpoint=False):
            q = p0 + t * (p1 - p0)
            if qc:
                x += np.linalg.norm(q - qc[-1])
            qc.append(q); dist.append(x)
    qc.append(np.array(PATH[-1], float)); x += np.linalg.norm(qc[-1] - qc[-2]); dist.append(x)
    qred = np.array(qc) @ (Ap.T / a)                         # reduced coordinates: q . a_i (units 2 pi)
    pp.run_qpoints(qred)
    fr = pp.get_qpoints_dict()["frequencies"]                # THz, [nq, 3]
    o = OUT / case; o.mkdir(parents=True, exist_ok=True)
    np.savetxt(o / "freq.gp", np.column_stack([dist, fr * THZ_CM1]), fmt="%.6f")
    pp.run_mesh([40, 40, 40]); pp.run_total_dos(sigma=0.05, freq_min=0, freq_max=12, freq_pitch=0.01)
    dd = pp.get_total_dos_dict()
    np.savetxt(o / "dos.dat", np.column_stack([dd["frequency_points"] * THZ_CM1, dd["total_dos"]]), fmt="%.6e")
    i13 = int(np.argmin(abs(np.array(dist) - 1 / 3)))
    print(f"{case}: a = {a:.4f} A, min freq {fr.min():.3f} THz, at zeta=1/3 on G-H: {np.round(fr[i13], 3)} THz")


if __name__ == "__main__":
    for c in sys.argv[1:]:
        run(c)
