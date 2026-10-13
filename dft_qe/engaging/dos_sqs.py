#!/usr/bin/env python3
"""Spin- and atom-projected densities of states of a relaxed 40-atom SQS (cf. Korzhavyi et al., PRB 79,
054202 (2009), Fig. 4).

  dos_sqs.py prepare <poscar> <x> <relax_dir> <out_dir>   scf.in (4x4x2), nscf.in (6x6x3, tetrahedra), projwfc.in
Uses the final (relaxed) structure of <relax_dir>/relax.out and the same settings as elastic_sqs_qe.py.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import elastic_sqs_qe as e  # noqa: E402

VAL = {"Fe": 16, "Cr": 14}


def prepare(poscar, x, relax_dir, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    cell0, sym, frac0 = e.read_poscar(poscar)
    cell, frac = e.final_structure(Path(relax_dir, "relax.out").read_text(), cell0, frac0)
    e.write_input(out / "scf.in", cell, sym, frac, (4, 4, 2), 60.0, "scf")
    import re
    s = (out / "scf.in").read_text().replace("prefix = 'sqs'", "prefix = 'dos'")
    s = re.sub(r"&ions\s*/\s*", "", s)                     # relaxation namelists are not used here
    s = re.sub(r"&cell[^/]*/\s*", "", s)
    (out / "scf.in").write_text(s)
    import os
    # empty-band margin: 1.2x occupied reaches ~+1.8 eV above E_F; QE_NBND_FACTOR=1.6 extends the DOS further
    nbnd = int(sum(VAL[t] for t in sym) / 2 * float(os.environ.get("QE_NBND_FACTOR", "1.2"))) + 8
    n = s.replace("calculation = 'scf'", "calculation = 'nscf'")
    n = n.replace("occupations = 'smearing', smearing = 'mv', degauss = 0.02",
                  f"occupations = 'tetrahedra_opt', nbnd = {nbnd}")
    n = n.replace("K_POINTS automatic\n4 4 2 0 0 0", "K_POINTS automatic\n6 6 3 0 0 0")
    n = n.replace("electron_maxstep = 300", "electron_maxstep = 300, diagonalization = 'cg'")   # Davidson failed (Fe-20Cr)
    (out / "nscf.in").write_text(n)
    # tetrahedron DOS from the nscf run (an explicit degauss here gave all-zero PDOS with tetrahedra_opt)
    (out / "projwfc.in").write_text("&projwfc\n  prefix = 'dos', outdir = './tmp', filpdos = 'pdosT',\n"
                                    "  DeltaE = 0.01\n/\n")
    print(out, "nbnd", nbnd)


if __name__ == "__main__":
    if sys.argv[1] == "prepare":
        prepare(*sys.argv[2:6])
