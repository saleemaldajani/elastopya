#!/bin/bash
#SBATCH -J qe_ebs
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=250G
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Unfolded band structure of a 40-atom SQS along Gamma-H-N-Gamma-P-H, reusing the SCF of unfold/x${X}.
# Usage: sbatch --dependency=afterok:<unfold job> --export=ALL,X=0.05 bandpath_sqs.sbatch
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
D=bandpath/sqs_x${X}
[ -f $D/bands.in ] || venv/bin/python bands_path.py sqs-prepare unfold/x${X} $D
cd $D
grep -q "JOB DONE" bands.out 2>/dev/null || $QE_MPIRUN -np 64 $QE_PW -nk 16 -in bands.in > bands.out
[ -f ebs_dw.npz ] || ../../venv/bin/python ../../bands_path.py sqs-unfold . > unfold.log 2>&1
ls ebs_*.npz && echo EBS_OK
