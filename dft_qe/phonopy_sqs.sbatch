#!/bin/bash
#SBATCH -J qe_phonopy
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=120G
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%A_%a.out
# Finite-displacement SCFs of phonopy_sqs.py. Array task i runs displacements i, i+STRIDE, i+2*STRIDE, ...
# Usage: sbatch --array=0-<STRIDE-1> --export=ALL,D=phonopy/x0.30,STRIDE=20 phonopy_sqs.sbatch
# Finished displacements (JOB DONE, converged) are skipped when requeued. After all tasks: phonopy_sqs.py collect $D
cd /orcd/pool/004/sdajani/fecr_qe/$D
source ../../qe75_env.sh
N=$(ls d[0-9][0-9][0-9].in | wc -l)
for ((i=SLURM_ARRAY_TASK_ID; i<N; i+=STRIDE)); do
  f=$(printf "d%03d" $i)
  if grep -q "JOB DONE" $f.out 2>/dev/null && ! grep -q "convergence NOT achieved" $f.out; then continue; fi
  $QE_MPIRUN -np 32 $QE_PW -nk 8 -in $f.in > $f.out
  rm -rf tmp/$f.save tmp/$f.wfc* 2>/dev/null
done
echo TASK_DONE
