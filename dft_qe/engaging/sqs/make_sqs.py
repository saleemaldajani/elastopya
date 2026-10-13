#!/usr/bin/env python3
"""40-atom bcc Fe-Cr SQS (2x2x5 conventional cells) with sqsgenerator.

Pair correlations are optimised over the first three neighbour shells (weights 1, 1/2, 1/3).
Writes x{X}.vasp (POSCAR format) and x{X}_sro.txt with the Warren-Cowley parameters.
Usage: make_sqs.py <out_dir> [iterations] [n_Cr,n_Cr,...]  (default 8,20,38)
"""
import json
import sys
from pathlib import Path

import sqsgenerator as sqs
from ase.io import write

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ITER = int(sys.argv[2]) if len(sys.argv) > 2 else 20_000_000
A0 = 2.87          # placeholder lattice parameter; relaxed in DFT
N = 40
N_CR = [int(v) for v in sys.argv[3].split(",")] if len(sys.argv) > 3 else [8, 20, 38]
for n_cr in N_CR:
    cfg = {
        "structure": {"lattice": [[A0, 0, 0], [0, A0, 0], [0, 0, A0]],
                      "coords": [[0, 0, 0], [0.5, 0.5, 0.5]], "species": ["Fe", "Fe"],
                      "supercell": [2, 2, 5]},
        "iterations": ITER,
        "shell_weights": {1: 1.0, 2: 0.5, 3: 1 / 3},
        "composition": {"Fe": N - n_cr, "Cr": n_cr},
        "max_results_per_objective": 1,
    }
    pack = sqs.optimize(cfg)
    best = pack.best()
    st = best.structure() if callable(best.structure) else best.structure
    atoms = sqs.to_ase(st)
    x = n_cr / N
    write(OUT / f"x{x:.2f}.vasp", atoms, format="vasp", direct=True, sort=True)
    with open(OUT / f"x{x:.2f}_sro.txt", "w") as fh:
        fh.write(f"objective {best.objective}\n")
        sro = best.sro() if callable(best.sro) else best.sro
        fh.write(f"sro {sro!r}\n")
    print(f"x={x:.2f}  n_Cr={n_cr}  objective={best.objective:.4g}")
