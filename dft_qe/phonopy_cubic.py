#!/usr/bin/env python3
"""Finite-displacement phonons of FM bcc Fe in a cubic 3x3x3 conventional supercell (54 atoms).

The Fe DFPT (ph.x) frequencies are ~60% too high, while a frozen-phonon check with identical settings agrees with
neutron data; the earlier phonopy reference used the 2x2x5 SQS-shaped cell, whose tetragonal image pattern makes the
symmetrised force constants inconsistent (spurious imaginary branch). A cubic supercell removes that problem.
Same physics as the SQS elastic runs: PBE, pslibrary PAW, 60/480 Ry, MV 0.02 Ry, ferromagnetic start.

  phonopy_cubic.py prepare <a_A> <dir>      phonopy_disp.yaml + scf input(s) d000.in ... (k = 5x5x5)
  phonopy_cubic.py collect <dir>            freq.gp / dos.dat in matdyn format (via phonopy_bcc_bands)
"""
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
os.environ.setdefault("QE_PSEUDO", "/orcd/pool/004/sdajani/fecr_qe/pseudo")
sys.path.insert(0, str(HERE))
N, K = 3, (5, 5, 5)


def prepare(a, d):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    import elastic_sqs_qe as e
    a = float(a); d = Path(d); d.mkdir(parents=True, exist_ok=True)
    uc = PhonopyAtoms(symbols=["Fe", "Fe"], cell=np.eye(3) * a, scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]])
    ph = Phonopy(uc, supercell_matrix=np.eye(3, dtype=int) * N, primitive_matrix="auto", log_level=0)
    ph.generate_displacements(distance=0.01)
    ph.save(d / "phonopy_disp.yaml")
    for i, sc in enumerate(ph.supercells_with_displacements):
        e.write_input(d / f"d{i:03d}.in", np.array(sc.cell), list(sc.symbols), np.array(sc.scaled_positions), K, 60, "scf")
    print(f"{d}: a = {a} A, {len(ph.supercell)} atoms, {len(ph.supercells_with_displacements)} displacement(s)")


def collect(d):
    import phonopy
    import phonopy_bcc_bands as b
    d = Path(d)
    ph = phonopy.load(d / "phonopy_disp.yaml", produce_fc=False, log_level=0)
    ph.forces = np.array([b.forces(d / f"d{i:03d}.out") for i in range(len(ph.supercells_with_displacements))])
    ph.produce_force_constants(); ph.symmetrize_force_constants()
    a = float(np.array(ph.unitcell.cell)[0, 0])
    Ap = 0.5 * a * np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], float)
    qc, dist, x = [], [], 0.0
    for p0, p1 in zip(b.PATH[:-1], b.PATH[1:]):
        p0, p1 = np.array(p0, float), np.array(p1, float)
        for t in np.linspace(0, 1, b.NSEG, endpoint=False):
            q = p0 + t * (p1 - p0)
            if qc:
                x += np.linalg.norm(q - qc[-1])
            qc.append(q); dist.append(x)
    qc.append(np.array(b.PATH[-1], float)); x += np.linalg.norm(qc[-1] - qc[-2]); dist.append(x)
    prim = ph.primitive
    Ap_found = np.array(prim.cell)
    qred = np.array(qc) @ (Ap_found.T / a)                     # q . a_i in units of 2 pi (a_i of phonopy's primitive)
    ph.run_qpoints(qred)
    fr = ph.get_qpoints_dict()["frequencies"]
    out = HERE / "phonopy_derived" / d.name; out.mkdir(parents=True, exist_ok=True)
    np.savetxt(out / "freq.gp", np.column_stack([dist, fr * b.THZ_CM1]), fmt="%.6f")
    ph.run_mesh([40, 40, 40]); ph.run_total_dos(sigma=0.05, freq_min=0, freq_max=12, freq_pitch=0.01)
    dos = ph.total_dos
    np.savetxt(out / "dos.dat", np.column_stack([dos.frequency_points * b.THZ_CM1, dos.dos]), fmt="%.6e")
    i13 = int(np.argmin(abs(np.array(dist) - 1 / 3)))
    print(f"{d.name}: min {fr.min():.3f} THz; zeta=1/3 on G-H {np.round(fr[i13], 3)}; H {np.round(fr[b.NSEG], 3)} THz")


if __name__ == "__main__":
    {"prepare": lambda: prepare(*sys.argv[2:4]), "collect": lambda: collect(sys.argv[2])}[sys.argv[1]]()
