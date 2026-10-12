#!/usr/bin/env python3
"""Cubic elastic constants of bcc Fe or Cr with Quantum ESPRESSO (stress-strain, central differences).

2-atom conventional bcc cell. Steps: (1) vc-relax the lattice parameter at the chosen k-mesh,
(2) apply strains e_xx = +-d (gives C11, C12) and engineering shear g_yz = +-d (gives C44) and
read the stress tensor. No internal relaxation is needed for elemental bcc (inversion symmetry).

Usage: elastic_qe.py <element Fe|Cr> <kx> <ky> <kz> <ecutwfc_Ry> <outdir> [--mag fm|afm|nm] [--np N]
Writes <outdir>/result.json with a0, C11, C12, C44, magnetic moment, and wall time.
"""
import argparse
import json
import re
import subprocess
import time
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-saleem-MIT-research-Scripta/37f8fe43-8e5f-4d11-883b-ebb44fa03560/scratchpad")
PW = SCR / "qesrc/q-e-qe-7.5/build/bin/pw.x"
MPIRUN = SCR / "qebuild/bin/mpirun"
PSEUDO = SCR / "pseudo"
UPF = {"Fe": ("55.845", "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF"), "Cr": ("51.996", "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF")}
RY_BOHR3_TO_GPA = 14710.507848  # 1 Ry/bohr^3 in GPa (QE prints kbar; we read kbar)


def write_input(path, el, cell, k, ecut, mag, calc):
    mass, upf = UPF[el]
    two_species = mag == "afm"
    species = f"{el}1 {mass} {upf}\n" + (f"{el}2 {mass} {upf}\n" if two_species else "")
    mag_lines = {"fm": f"starting_magnetization(1) = 0.4\n",
                 "afm": f"starting_magnetization(1) = 0.5\nstarting_magnetization(2) = -0.5\n",
                 "nm": ""}[mag]
    nspin = 1 if mag == "nm" else 2
    second = f"{el}2" if two_species else f"{el}1"
    txt = f"""&control
  calculation = '{calc}', prefix = 'x', outdir = './tmp', pseudo_dir = '{PSEUDO}',
  tstress = .true., tprnfor = .true., etot_conv_thr = 1e-6, forc_conv_thr = 1e-4
/
&system
  ibrav = 0, nat = 2, ntyp = {2 if two_species else 1}, ecutwfc = {ecut}, ecutrho = {8 * ecut},
  occupations = 'smearing', smearing = 'mv', degauss = 0.02, nspin = {nspin}
  {mag_lines}/
&electrons
  conv_thr = 1e-9, mixing_beta = 0.2, mixing_mode = 'local-TF', electron_maxstep = 300
/
&ions
/
&cell
  press_conv_thr = 0.1
/
ATOMIC_SPECIES
{species}CELL_PARAMETERS angstrom
{cell[0][0]:.10f} {cell[0][1]:.10f} {cell[0][2]:.10f}
{cell[1][0]:.10f} {cell[1][1]:.10f} {cell[1][2]:.10f}
{cell[2][0]:.10f} {cell[2][1]:.10f} {cell[2][2]:.10f}
ATOMIC_POSITIONS crystal
{el}1 0.0 0.0 0.0
{second} 0.5 0.5 0.5
K_POINTS automatic
{k[0]} {k[1]} {k[2]} 0 0 0
"""
    path.write_text(txt)


def run(d, name, np_):
    out = d / f"{name}.out"
    subprocess.run([str(MPIRUN), "-np", str(np_), str(PW), "-nk", str(min(np_, 4)), "-in", f"{name}.in"],
                   cwd=d, stdout=open(out, "w"), stderr=subprocess.STDOUT, check=True,
                   env={"TMPDIR": "/tmp", "HOME": str(Path.home()), "PATH": "/usr/bin:/bin", "OMP_NUM_THREADS": "1"})
    return out.read_text()


def stress_kbar(text):
    if 'convergence NOT achieved' in text or 'total   stress' not in text:
        raise RuntimeError('SCF not converged or no stress printed')
    blk = text.split("total   stress")[-1].splitlines()[1:4]
    return np.array([[float(v) for v in l.split()[3:6]] for l in blk])  # kbar columns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("el"); ap.add_argument("k", nargs=3, type=int); ap.add_argument("ecut", type=float)
    ap.add_argument("outdir"); ap.add_argument("--mag", default="fm"); ap.add_argument("--np", type=int, default=12)
    ap.add_argument("--strain", type=float, default=0.005)
    a = ap.parse_args()
    d = Path(a.outdir); d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    a_guess = 2.83 if a.el == "Fe" else 2.87
    write_input(d / "relax.in", a.el, np.eye(3) * a_guess, a.k, a.ecut, a.mag, "vc-relax")
    txt = (d / "relax.out").read_text() if (d / "relax.out").exists() and "End final" in (d / "relax.out").read_text() \
        else run(d, "relax", a.np)
    cells = re.findall(r"CELL_PARAMETERS \(angstrom\)\n(.*)\n(.*)\n(.*)", txt)
    vecs = np.array([[float(v) for v in r.split()] for r in cells[-1]])
    a0 = float(np.cbrt(abs(np.linalg.det(vecs))))  # cubic cell edge
    mag = re.findall(r"total magnetization\s+=\s+([-\d.]+)", txt)
    res = {"element": a.el, "k": a.k, "ecutwfc_Ry": a.ecut, "mag": a.mag, "a0_A": a0,
           "total_magnetization_muB_per_cell": float(mag[-1]) if mag else 0.0}
    d0 = a.strain
    sig = {}
    for tag, F in {"xx": np.diag([1 + d0, 1, 1]), "yz": np.array([[1, 0, 0], [0, 1, d0 / 2], [0, d0 / 2, 1]])}.items():
        for sgn in (1, -1):
            Fs = np.eye(3) + sgn * (F - np.eye(3))
            cell = (np.eye(3) * a0) @ Fs.T
            name = f"s_{tag}_{'p' if sgn > 0 else 'm'}"
            write_input(d / f"{name}.in", a.el, cell, a.k, a.ecut, a.mag, "scf")
            sig[(tag, sgn)] = -stress_kbar(run(d, name, a.np)) * 0.1  # kbar -> GPa, QE sign: pressure>0
    dxx = (sig[("xx", 1)] - sig[("xx", -1)]) / (2 * d0)
    dyz = (sig[("yz", 1)] - sig[("yz", -1)]) / (2 * d0)
    res.update(C11=dxx[0, 0], C12=0.5 * (dxx[1, 1] + dxx[2, 2]), C44=dyz[1, 2], wall_s=time.time() - t0)
    (d / "result.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: (round(v, 2) if isinstance(v, float) else v) for k, v in res.items()}))


if __name__ == "__main__":
    main()
