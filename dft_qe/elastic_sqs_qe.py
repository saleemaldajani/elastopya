#!/usr/bin/env python3
"""Elastic constants of a 40-atom Fe-Cr SQS with Quantum ESPRESSO (stress-strain).

Steps: (1) vc-relax of the three box lengths (cell_dofree='xyz') with relaxed ions,
(2) strains e_xx = +-d and engineering shear g_yz = +-d applied to the relaxed cell, ions relaxed
at each strain ('relax'), (3) central differences of the stress:
C11 = dsig_xx/de_xx, C12 = mean(dsig_yy, dsig_zz)/de_xx, C44 = dsig_yz/dg_yz.
Ferromagnetic start (Fe 2.5, Cr 1.0 muB); PBE-PAW (pslibrary); MV smearing 0.02 Ry.

Usage: elastic_sqs_qe.py <POSCAR> <x_cr> <kx> <ky> <kz> <ecut_Ry> <outdir> [--np N] [--only-scf]
"""
import argparse
import json
import os
import sys
import re
import subprocess
import time
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-saleem-MIT-research-Scripta/37f8fe43-8e5f-4d11-883b-ebb44fa03560/scratchpad")
# Local defaults; on a cluster set QE_PW, QE_MPIRUN, QE_PSEUDO (and optionally QE_NK, the number of
# k-point pools). With QE_MPIRUN set, the job's own environment (modules, Slurm) is passed to MPI.
PW = Path(os.environ.get("QE_PW", SCR / "qesrc/q-e-qe-7.5/build/bin/pw.x"))
MPIRUN = Path(os.environ.get("QE_MPIRUN", SCR / "qebuild/bin/mpirun"))
PSEUDO = Path(os.environ.get("QE_PSEUDO", SCR / "pseudo"))


def read_poscar(path):
    L = Path(path).read_text().split("\n")
    scale = float(L[1]); cell = np.array([[float(v) for v in L[i].split()] for i in (2, 3, 4)]) * scale
    species, counts = L[5].split(), [int(c) for c in L[6].split()]
    frac = np.array([[float(v) for v in L[8 + i].split()[:3]] for i in range(sum(counts))])
    symbols = [s for s, c in zip(species, counts) for _ in range(c)]
    return cell, symbols, frac


def write_input(path, cell, symbols, frac, k, ecut, calc):
    sp = sorted(set(symbols))
    mass = {"Fe": 55.845, "Cr": 51.996}
    upf = {"Fe": "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF", "Cr": "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF"}
    mag = {"Fe": 0.6, "Cr": float(os.environ.get("QE_MAG_CR", "0.25"))}   # QE_MAG_CR < 0: Cr antiparallel to Fe
    lines = [f"""&control
  calculation = '{calc}', prefix = 'sqs', outdir = './tmp', pseudo_dir = '{PSEUDO}',
  tstress = .true., tprnfor = .true., etot_conv_thr = 1e-5, forc_conv_thr = 1e-3, nstep = 60
/
&system
  ibrav = 0, nat = {len(symbols)}, ntyp = {len(sp)}, ecutwfc = {ecut}, ecutrho = {8 * ecut},
  occupations = 'smearing', smearing = 'mv', degauss = 0.02, nspin = 2,
""" + "".join(f"  starting_magnetization({i + 1}) = {mag[s]}\n" for i, s in enumerate(sp)) + """/
&electrons
  conv_thr = 1e-8, mixing_beta = """ + os.environ.get("QE_MIXING_BETA", "0.2") + """, mixing_mode = 'local-TF', electron_maxstep = 300"""
  + (", mixing_ndim = " + os.environ["QE_MIXING_NDIM"] if "QE_MIXING_NDIM" in os.environ else "")
  + (",\n  startingpot = 'file', startingwfc = 'file'" if os.environ.get("QE_START_FROM_FILE") else "") + """
/
&ions
/
&cell
  cell_dofree = 'xyz', press_conv_thr = 0.2
/
ATOMIC_SPECIES
""" + "".join(f"{s} {mass[s]} {upf[s]}\n" for s in sp) + "CELL_PARAMETERS angstrom\n"
             + "".join(" ".join(f"{v:.10f}" for v in row) + "\n" for row in cell)
             + "ATOMIC_POSITIONS crystal\n"
             + "".join(f"{s} {f[0]:.10f} {f[1]:.10f} {f[2]:.10f}\n" for s, f in zip(symbols, frac))
             + f"K_POINTS automatic\n{k[0]} {k[1]} {k[2]} 0 0 0\n"]
    path.write_text("".join(lines))


def run(d, name, np_):
    out = d / f"{name}.out"
    if out.exists() and "JOB DONE" in out.read_text():
        return out.read_text()
    nk = int(os.environ.get("QE_NK", min(np_, 4)))
    env = ({**os.environ, "OMP_NUM_THREADS": "1"} if "QE_MPIRUN" in os.environ else
           {"TMPDIR": "/tmp", "HOME": str(Path.home()), "PATH": "/usr/bin:/bin", "OMP_NUM_THREADS": "1"})
    subprocess.run([str(MPIRUN), "-np", str(np_), str(PW), "-nk", str(nk), "-in", f"{name}.in"],
                   cwd=d, stdout=open(out, "w"), stderr=subprocess.STDOUT, check=True, env=env)
    return out.read_text()


def final_structure(text, cell0, frac0):
    cells = re.findall(r"CELL_PARAMETERS \(angstrom\)\n(.*)\n(.*)\n(.*)", text)
    cell = np.array([[float(v) for v in r.split()] for r in cells[-1]]) if cells else cell0
    blocks = text.split("ATOMIC_POSITIONS (crystal)")
    if len(blocks) > 1:
        rows = blocks[-1].strip().split("\n")[: len(frac0)]
        frac = np.array([[float(v) for v in r.split()[1:4]] for r in rows])
    else:
        frac = frac0
    return cell, frac


def stress_gpa(text):
    if "convergence NOT achieved" in text or "total   stress" not in text:
        raise RuntimeError("SCF not converged or no stress")
    blk = text.split("total   stress")[-1].splitlines()[1:4]
    return -np.array([[float(v) for v in l.split()[3:6]] for l in blk]) * 0.1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar"); ap.add_argument("x", type=float); ap.add_argument("k", nargs=3, type=int)
    ap.add_argument("ecut", type=float); ap.add_argument("outdir")
    ap.add_argument("--np", type=int, default=12); ap.add_argument("--strain", type=float, default=0.005)
    ap.add_argument("--only-scf", action="store_true")
    a = ap.parse_args()
    d = Path(a.outdir); d.mkdir(parents=True, exist_ok=True)
    cell, sym, frac = read_poscar(a.poscar)
    a0 = 2.83 * (1 - a.x) + 2.87 * a.x                      # Vegard start from PBE Fe/Cr
    cell = cell / cell[0, 0] * 2 * a0
    t0 = time.time()
    if a.only_scf:
        write_input(d / "bench.in", cell, sym, frac, a.k, a.ecut, "scf")
        txt = run(d, "bench", a.np)
        print(json.dumps({"wall_s": round(time.time() - t0), "scf_iter": len(re.findall("iteration #", txt))}))
        return
    write_input(d / "relax.in", cell, sym, frac, a.k, a.ecut, "vc-relax")
    cell, frac = final_structure(run(d, "relax", a.np), cell, frac)
    # Guard: an under-converged cutoff makes vc-relax collapse the cell (16-atom test at 40 Ry: -12 %).
    # PBE volumes are ~11.3 (Fe) to ~11.6 (Cr) A^3/atom; stop before the expensive strain steps.
    v_at = abs(np.linalg.det(cell)) / len(sym)
    if not 10.8 < v_at < 12.2:
        sys.exit(f"relaxed volume {v_at:.2f} A^3/atom is unphysical; raise ecut (Fe PAW needs >= 60 Ry)")
    sig = {}
    for tag, E in {"xx": np.diag([a.strain, 0, 0]),
                   "yz": np.array([[0, 0, 0], [0, 0, a.strain / 2], [0, a.strain / 2, 0]])}.items():
        for sgn in (1, -1):
            c = cell @ (np.eye(3) + sgn * E).T
            name = f"s_{tag}_{'p' if sgn > 0 else 'm'}"
            write_input(d / f"{name}.in", c, sym, frac, a.k, a.ecut, "relax")
            sig[(tag, sgn)] = stress_gpa(run(d, name, a.np))
    dxx = (sig[("xx", 1)] - sig[("xx", -1)]) / (2 * a.strain)
    dyz = (sig[("yz", 1)] - sig[("yz", -1)]) / (2 * a.strain)
    mags = re.findall(r"total magnetization\s+=\s+([-\d.]+)", (d / "relax.out").read_text())
    res = {"x_cr": a.x, "k": a.k, "ecut_Ry": a.ecut, "volume_A3": float(abs(np.linalg.det(cell))),
           "moment_muB_per_atom": float(mags[-1]) / len(sym) if mags else None,
           "C11": dxx[0, 0], "C12": 0.5 * (dxx[1, 1] + dxx[2, 2]), "C44": dyz[1, 2],
           "wall_s": time.time() - t0}
    (d / "result.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: (round(v, 2) if isinstance(v, float) else v) for k, v in res.items()}))


if __name__ == "__main__":
    main()
