#!/bin/bash
#SBATCH -J qe_sqs40
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 80
#SBATCH --mem=240G
#SBATCH -t 12:00:00
#SBATCH -o slurm-%x-%j.out
# 40-atom SQS elastic constants, same settings as the local runs (60 Ry, k 4x4x2).
# Usage: sbatch --export=ALL,X=0.20 sqs40.sbatch      (X = 0.20, 0.50 or 0.95)
# 20 k-points x 2 spins: 20 pools of 4 ranks. Finished steps (JOB DONE) are skipped on resubmission.
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
export QE_NK=20
$PY -I elastic_sqs_qe.py sqs/x${X}.vasp ${X} 4 4 2 60 runs/sqs40_x${X}_k442_e60 --np 80
