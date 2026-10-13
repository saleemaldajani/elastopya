#!/bin/bash
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Same as sqs40.sbatch, written to runs/sqs40_x${X}_k442_e60_${TAG}; used for alternative magnetic starting
# configurations (QE_MAG_CR) and SCF mixing (QE_MIXING_BETA, QE_MIXING_NDIM). Same physics settings otherwise.
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
export QE_NK=20
$PY -I elastic_sqs_qe.py sqs/x${X}.vasp ${X} 4 4 2 60 runs/sqs40_x${X}_k442_e60_${TAG} --np 60
