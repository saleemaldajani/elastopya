#!/bin/bash
# Elastic constants of random (homogeneous) bcc Fe(1-x)Cr(x) alloys with the Waldron/Acta
# LAMMPS stiffness package, at the compositions of the DFT calculations.
#
# Identical to the Acta stiffness step (same potential.mod, displace.mod, stiffness.in,
# 50x50x50 bcc cell, a = 2.866 A, strain 1e-6, CG box/ionic relaxation) except that the
# structure is a random alloy created with `set type/fraction` instead of a dumped MC
# structure. stiffness.in has one added line that prints the relaxed box.
#
# Usage: run_random_alloys.sh <lmp_mpi> <nproc> <x_cr> [seed ...]
set -euo pipefail
LMP=$1; NP=$2; X=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd)
# Seeds: the ten `set type/fraction` seeds of the Acta alloys (ExtraInformation.docx)
SEEDS=${*:-"23609378 67101604 8085385 38507825 63689744 20009178 67502748 70189207 83842524 90137588"}
for SEED in $SEEDS; do
  D="$HERE/runs/x${X}_s${SEED}"
  mkdir -p "$D"
  cp "$HERE"/template/{displace.mod,potential.mod,stiffness.in,FeCr_d.eam.alloy.txt,FeCr_s.eam.fs.txt} "$D/"
  { cat "$HERE/template/init.mod.head"
    cat <<EOF
region SYSTEM prism 0 \${L} 0 \${L} 0 \${L} 0.0 0.0 0.0 units lattice
create_box 2 SYSTEM
create_atoms 1 region SYSTEM
set group all type/fraction 2 ${X} ${SEED}
mass 1 55.85 #Fe
mass 2 51.996 #Cr
EOF
  } > "$D/init.mod"
  (cd "$D" && TMPDIR=/tmp "$(dirname "$LMP")/mpirun" -np "$NP" "$LMP" -in stiffness.in -log log.lammps > stdout.txt 2>&1)
  # pure elements: the seed has no effect, one run is enough
  if [ "$X" = "0" ] || [ "$X" = "1" ] || [ "$X" = "0.0" ] || [ "$X" = "1.0" ]; then break; fi
done
