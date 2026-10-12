#!/bin/bash
# Lawrence Waldron's spinodal-decomposition workflow (Acta Mater. 246 (2023) 118552, Sec. 3.7)
# at an arbitrary composition and cell size.
#
#  1. generate: generatefinal.in -- random bcc alloy, Monte Carlo atom swaps at 10 K driven by
#     the artificial LJ unlike-pair repulsion (epsilon = 10, sigma = 10, cutoff 2.866 A).
#     fix atom/swap 50 3000 as in the original; run length scaled so that the number of
#     attempted swaps PER ATOM equals the original (40,000 steps for a 50^3 cell).
#  2. convert: as convert10..50.in, with the original seeds. A fraction c of the Cr atoms is
#     turned into Fe and a fraction c*N_Cr/N_Fe of the Fe atoms into Cr, which keeps the
#     composition fixed (identical to the original at 50/50, where N_Cr = N_Fe).
#     c = 2(1-x)(1-phi), so phi = 1 is the decomposed structure and phi = 0.5 is random
#     (the original c = 0.1..0.5 at x = 0.5).
#  3. stiffness: the Acta stiffness package (random_alloys/template), reading each structure.
#
# Usage: run_decomposition.sh <lmp_mpi> <nproc> <x_cr> <L> [alloy 1-10]
# Format-only differences: atom_style atomic and write_data/read_data instead of xyz dumps.
set -euo pipefail
LMP=$1; NP=$2; X=$3; L=$4; A=${5:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
TPL="$HERE/../random_alloys/template"
MPIRUN="$(dirname "$LMP")/mpirun"
SET_SEEDS=(23609378 67101604 8085385 38507825 63689744 20009178 67502748 70189207 83842524 90137588)
SWAP_SEEDS=(50941252 63770018 25875305 97705767 83263562 75547858 63061468 48086192 62931845 58363333)
CONV_A=(68860890 75832863 28717628 15965553 11347432)   # G22 (Cr -> Fe), conversions 10..50
CONV_B=(62081045 89205392 75026443 76192003 3721942)    # G11 (Fe -> Cr)
STEPS=$(python3 -c "print(int(round(40000*($L/50)**3/50))*50)")
D="$HERE/runs/x${X}_L${L}_a${A}"; mkdir -p "$D"; cd "$D"
run() { TMPDIR=/tmp "$MPIRUN" -np "$NP" "$LMP" -in "$1" -log "${1%.in}.log" > "${1%.in}.out" 2>&1; }

cat > gen.in <<EOF
units metal
dimension 3
boundary p p p
atom_style atomic
lattice bcc 2.866
variable L equal $L
region SYSTEM prism 0 \${L} 0 \${L} 0 \${L} 0.0 0.0 0.0 units lattice
create_box 2 SYSTEM
create_atoms 1 region SYSTEM
mass 1 55.85
mass 2 51.996
set group all type/fraction 2 $X ${SET_SEEDS[$((A-1))]}
pair_style lj/cut 2.866
pair_coeff 1 1 0.0 0.0 2.866
pair_coeff 2 2 0.0 0.0 2.866
pair_coeff 1 2 10 10 2.866
fix 3 all atom/swap 50 3000 ${SWAP_SEEDS[$((A-1))]} 10 types 1 2
thermo_style custom step pe
thermo 500
timestep 0.01
run $STEPS
write_data phi1.000.data nocoeff
EOF
[ -f phi1.000.data ] || run gen.in

PHIS=(0.900 0.800 0.700 0.600 0.500)
for k in 0 1 2 3 4; do
  PHI=${PHIS[$k]}
  C=$(python3 -c "print(2*(1-$X)*(1-$PHI))")
  cat > conv$PHI.in <<EOF
units metal
atom_style atomic
read_data phi1.000.data
group G22 type 2
group G11 type 1
variable fB equal $C*count(G22)/count(G11)
set group G22 type/fraction 1 $C ${CONV_A[$k]}
set group G11 type/fraction 2 \${fB} ${CONV_B[$k]}
write_data phi$PHI.data
EOF
  [ -f phi$PHI.data ] || run conv$PHI.in
done

for PHI in 1.000 0.900 0.800 0.700 0.600 0.500; do
  S="stiff_phi$PHI"; mkdir -p "$S"
  cp "$TPL"/{displace.mod,potential.mod,stiffness.in,FeCr_d.eam.alloy.txt,FeCr_s.eam.fs.txt} "$S/"
  { grep -v -E "^(lattice|variable L|region|create_box|read_dump)" "$TPL/init.mod.head"
    echo "read_data ../phi$PHI.data"; } > "$S/init.mod"
  (cd "$S" && [ -f log.lammps ] && grep -q "Total wall" log.lammps) || \
    (cd "$S" && TMPDIR=/tmp "$MPIRUN" -np "$NP" "$LMP" -in stiffness.in -log log.lammps > stdout.txt 2>&1)
done
echo "done x=$X L=$L alloy=$A"
