#!/bin/bash
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=250G
#SBATCH -t 12:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Finite-displacement Fe phonons in a cubic 3x3x3 conventional supercell (phonopy_cubic.py prepare/collect).
# Usage: sbatch --export=ALL,D=Fe_cub phcubic.sbatch      (run from phonopy_cubic/)
cd /orcd/pool/004/sdajani/fecr_qe/phonopy_cubic/$D
source ../../qe75_env.sh
for f in d*.in; do o=${f%.in}.out; grep -q "JOB DONE" $o 2>/dev/null || $QE_MPIRUN -np 64 $QE_PW -nk 8 -in $f > $o; done
grep -l "convergence NOT" d*.out && echo "SCF NOT CONVERGED"
echo TASK_DONE
