#!/usr/bin/env python3
"""Phonon DOS and vibrational entropy of mixing of the 40-atom Fe-Cr SQS by finite displacements (phonopy +
Quantum ESPRESSO), for comparison with the neutron-scattering results of Fultz et al., PRB 52, 3280 (1995).

The force constants are those of the 40-atom cell itself (supercell = the SQS, 2x2x5 conventional bcc cells);
pure Fe (FM) and pure Cr (NM) are computed in the same cell, so that the vibrational entropy of mixing
  dS_vib(x, T) = S_vib(alloy) - (1 - x) S_vib(Fe) - x S_vib(Cr)
is evaluated with identical Brillouin-zone sampling for all three. Same pseudopotentials, cutoff, smearing and
k-mesh as the elastic constants (elastic_sqs_qe.py); SCF converged to 1e-10 Ry.

  phonopy_sqs.py prepare <poscar> <relax_dir> <outdir>       relaxed alloy (relax.out of elastic_sqs_qe.py)
  phonopy_sqs.py prepare-pure <poscar> <Fe|Cr> <a_A> <outdir> the same cell filled with one element, lattice a
  phonopy_sqs.py collect <outdir>                             forces -> force constants -> DOS, S_vib(T), F_vib(T)
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import elastic_sqs_qe as e  # noqa: E402

RYBOHR_EVA = 13.605693123 / 0.529177210903
DIST = 0.01            # displacement (A)
MESH = [12, 12, 5]     # q-mesh in the reciprocal cell of the 40-atom cell
T_GRID = dict(t_min=0, t_max=1600, t_step=10)


def _phonopy(cell, sym, frac):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    return Phonopy(PhonopyAtoms(symbols=list(sym), cell=cell, scaled_positions=frac),
                   supercell_matrix=np.eye(3, dtype=int), primitive_matrix="P")   # the 40-atom cell itself


def _write_displacements(ph, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    ph.generate_displacements(distance=DIST, is_plusminus=False)
    ph.save(out / "phonopy_disp.yaml")
    for i, sc in enumerate(ph.supercells_with_displacements):
        f = out / f"d{i:03d}.in"
        e.write_input(f, sc.cell, sc.symbols, sc.scaled_positions, (4, 4, 2), 60.0, "scf")
        s = f.read_text().replace("conv_thr = 1e-8", "conv_thr = 1e-10")
        s = re.sub(r"&ions\s*/\s*", "", s); s = re.sub(r"&cell[^/]*/\s*", "", s)
        f.write_text(s.replace("prefix = 'sqs'", f"prefix = 'd{i:03d}'"))
    print(out, "displacements:", len(ph.supercells_with_displacements))


def prepare(poscar, relax_dir, out):
    cell0, sym, frac0 = e.read_poscar(poscar)
    rin = Path(relax_dir, "relax.in").read_text().split("CELL_PARAMETERS angstrom\n")[1].split("\n")[:3]
    cell0 = np.array([[float(v) for v in r.split()] for r in rin])
    text = Path(relax_dir, "relax.out").read_text()
    if "End final coordinates" not in text:
        sys.exit(f"{relax_dir}: relaxation not finished")
    cell, frac = e.final_structure(text, cell0, frac0)
    _write_displacements(_phonopy(cell, sym, frac), out)


def prepare_pure(poscar, element, a, out):
    _, sym, frac = e.read_poscar(poscar)
    n = np.rint(np.diag(e.read_poscar(poscar)[0]) / e.read_poscar(poscar)[0][0, 0] * 2).astype(int)  # 2,2,5
    cell = np.diag(n * float(a))
    ideal = np.round(frac * 2 * n) / (2 * n)            # ideal bcc sites of the 2x2x5 cell
    _write_displacements(_phonopy(cell, [element] * len(sym), ideal), out)


def forces(out_file):
    t = Path(out_file).read_text()
    if "JOB DONE" not in t or "convergence NOT achieved" in t:
        raise RuntimeError(f"{out_file}: not converged")
    blk = t.split("Forces acting on atoms (cartesian axes, Ry/au):")[-1]
    f = [[float(v) for v in m] for m in re.findall(r"atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", blk)]
    return np.array(f) * RYBOHR_EVA


def collect(out):
    import phonopy
    out = Path(out)
    ph = phonopy.load(out / "phonopy_disp.yaml", produce_fc=False, log_level=0)
    n = len(ph.supercells_with_displacements)
    ph.forces = np.array([forces(out / f"d{i:03d}.out") for i in range(n)])
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    ph.save(out / "phonopy_params.yaml", settings={"force_constants": True})
    ph.run_mesh(MESH)
    ph.run_total_dos(sigma=0.15, freq_min=-1, freq_max=12, freq_pitch=0.02)
    d = ph.get_total_dos_dict()
    ph.run_thermal_properties(**T_GRID)
    tp = ph.get_thermal_properties_dict()
    nat = len(ph.supercell)
    R = 8.314462618
    fmin = float(np.min(ph.mesh.frequencies))
    np.savez(out / "phonons.npz", freq_thz=d["frequency_points"], dos_per_atom=d["total_dos"] / nat,
             T=tp["temperatures"], S_kB_per_atom=tp["entropy"] / (nat * R),
             F_eV_per_atom=tp["free_energy"] / (nat * 96.485332), min_freq_thz=fmin)
    print(out.name, "min freq (THz):", round(fmin, 3), " S_vib(300 K) kB/atom:",
          round(float(np.interp(300, tp["temperatures"], tp["entropy"])) / (nat * R), 3))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "prepare":
        prepare(*sys.argv[2:5])
    elif cmd == "prepare-pure":
        prepare_pure(*sys.argv[2:6])
    elif cmd == "collect":
        collect(sys.argv[2])
