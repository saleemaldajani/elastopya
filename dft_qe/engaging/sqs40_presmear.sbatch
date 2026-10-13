#!/bin/bash
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Two-stage start for magnetically frustrated SQS (Fe-80Cr): (1) one SCF with degauss 0.05 Ry from the same
# starting structure, (2) the standard elastic workflow (degauss 0.02 Ry, unchanged physics) starting from that
# charge density and wavefunctions. Usage: sbatch --export=ALL,X=0.80 -J qe_sqs40_x0.80_presmear sqs40_presmear.sbatch
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
export QE_NK=20
D=runs/sqs40_x${X}_k442_e60_presmear
mkdir -p $D
if ! grep -q "JOB DONE" $D/scf_presmear.out 2>/dev/null; then
  $PY -I - "$X" "$D" <<'PYEOF'
import re, sys
from pathlib import Path
sys.path.insert(0, ".")
import elastic_sqs_qe as e
x, d = sys.argv[1], Path(sys.argv[2])
cell, sym, frac = e.read_poscar(f"sqs/x{x}.vasp")
e.write_input(d / "scf_presmear.in", cell, sym, frac, (4, 4, 2), 60.0, "scf")
s = (d / "scf_presmear.in").read_text().replace("degauss = 0.02", "degauss = 0.05").replace("conv_thr = 1e-8", "conv_thr = 1e-7")
s = re.sub(r"&ions\s*/\s*", "", s); s = re.sub(r"&cell[^/]*/\s*", "", s)
(d / "scf_presmear.in").write_text(s)
PYEOF
  (cd $D && $QE_MPIRUN -np 60 $QE_PW -nk 20 -in scf_presmear.in > scf_presmear.out)
fi
export QE_START_FROM_FILE=1 QE_MIXING_BETA=0.1 QE_MIXING_NDIM=16
$PY -I elastic_sqs_qe.py sqs/x${X}.vasp ${X} 4 4 2 60 $D --np 60
