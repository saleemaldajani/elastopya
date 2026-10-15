#!/bin/bash
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mem=200G
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# One strain step of a finished relaxation (elastic_sqs_qe.py with QE_STRAIN_ONLY=$STEP), so the four steps run in
# parallel. Submit with the same QE_* settings as the relaxation job, e.g.
#   sbatch --export=ALL,X=0.80,TAG=crafm2,STEP=s_xx_p,QE_MAG_CR=-0.3,QE_MIXING_BETA=0.1,QE_MIXING_NDIM=16,QE_MAXSTEP=600,QE_UPSCALE=1 strain_par.sbatch
# then run sqs40_variant.sbatch once more (afterok) to assemble result.json without recomputing.
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
export QE_NK=20 QE_STRAIN_ONLY=$STEP
$PY -I elastic_sqs_qe.py sqs/x${X}.vasp ${X} 4 4 2 60 runs/sqs40_x${X}_k442_e60_${TAG} --np 60
