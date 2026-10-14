#!/bin/bash
#SBATCH -J qe_unfold
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 2
#SBATCH --ntasks-per-node=64
#SBATCH --mem=250G
#SBATCH -t 2-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# As unfold.sbatch, for the band step that did not finish in 24 h with CG on one node (unfold/x0.00, 510 SC k-points
# x 2 spins, 392 bands): two nodes (k-point pools), the 2-day limit, and ParO diagonalisation (PPCG was removed from QE in Dec. 2024), which converges the same
# eigenproblem; if ParO fails the step is rerun with the original CG input. Same physics as unfold.sbatch.
# Usage: sbatch --export=ALL,X=0.05 unfold2.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/unfold/x${X}
source ../../qe75_env.sh
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 128 $QE_PW -nk 16 -in scf.in > scf.out
[ -f bands_cg.in ] || cp bands.in bands_cg.in
sed "s/diagonalization = 'cg'/diagonalization = 'paro'/" bands_cg.in > bands_paro.in
if ! grep -q "JOB DONE" bands.out 2>/dev/null; then
  $QE_MPIRUN -np 128 $QE_PW -nk 32 -in bands_paro.in > bands.out
  if ! grep -q "JOB DONE" bands.out; then
    mv bands.out bands_paro_failed.out
    $QE_MPIRUN -np 128 $QE_PW -nk 32 -in bands_cg.in > bands.out
  fi
fi
grep -q "JOB DONE" bands.out && [ ! -f unfolded_dw.npz ] && ../../venv/bin/python ../../unfold_fs.py unfold . > unfold.log 2>&1
tail -3 unfold.log
