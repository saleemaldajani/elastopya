#!/bin/bash
#SBATCH -J qe_ebs
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=250G
#SBATCH -t 2-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Unfolded band structure of a 40-atom SQS along Gamma-H-N-Gamma-P-H, reusing the SCF of unfold/x${X}.
# Usage: sbatch --dependency=afterok:<unfold job> --export=ALL,X=0.05 bandpath_sqs.sbatch
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
D=bandpath/sqs_x${X}
[ -f $D/bands.in ] || venv/bin/python bands_path.py sqs-prepare unfold/x${X} $D
cd $D
# ParO first (with ParO the unfolding band step took ~1.5 h instead of >24 h with CG); CG if ParO fails
[ -f bands_cg.in ] || cp bands.in bands_cg.in
sed "s/diagonalization = 'cg'/diagonalization = 'paro'/" bands_cg.in > bands_paro.in
if ! grep -q "JOB DONE" bands.out 2>/dev/null; then
  $QE_MPIRUN -np 64 $QE_PW -nk 16 -in bands_paro.in > bands.out
  grep -q "JOB DONE" bands.out || { mv bands.out bands_paro_failed.out; $QE_MPIRUN -np 64 $QE_PW -nk 16 -in bands_cg.in > bands.out; }
fi
[ -f ebs_dw.npz ] || ../../venv/bin/python ../../bands_path.py sqs-unfold . > unfold.log 2>&1
ls ebs_*.npz && echo EBS_OK
