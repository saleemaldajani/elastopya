#!/bin/bash
#SBATCH -J qe_dos_wide
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mem=200G
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Wider-window DOS: NSCF with 1.6x occupied bands (reusing the converged SCF in dos/x${X}), projwfc -> pdosW.*
# Usage: sbatch [--dependency=...] --export=ALL,X=0.95 dos_wide.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/dos/x${X}
source ../../qe75_env.sh
grep -q "JOB DONE" scf.out || { echo "no converged SCF"; exit 1; }
if [ ! -f nscf_wide.in ]; then
  occ=$(grep -o "nbnd = [0-9]*" nscf.in | awk '{print $3}')
  new=$(python3 -c "print(int(($occ - 8) / 1.2 * 1.6) + 8)")
  sed "s/nbnd = $occ/nbnd = $new/" nscf.in | sed "/diagonalization/!s/electron_maxstep = 300/electron_maxstep = 300, diagonalization = 'cg'/" > nscf_wide.in
  printf "&projwfc\n  prefix = 'dos', outdir = './tmp', filpdos = 'pdosW',\n  DeltaE = 0.01\n/\n" > projwfc_wide.in
fi
grep "nbnd" nscf_wide.in
grep -q "JOB DONE" nscf_wide.out 2>/dev/null || $QE_MPIRUN -np 60 $QE_PW -nk 20 -in nscf_wide.in > nscf_wide.out
grep -q "JOB DONE" nscf_wide.out || { echo NSCF_FAILED; exit 1; }
$QE_MPIRUN -np 16 $(dirname $QE_PW)/projwfc.x -nk 1 -in projwfc_wide.in > projwfc_wide.out
ls pdosW.pdos_tot && echo DOS_WIDE_OK
