#!/bin/bash
#SBATCH -J qe_bcc2
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o slurm-%x-%j.out
# Pure bcc Fe (FM) or Cr (commensurate AFM), 2-atom cell, 60 Ry.
# Usage: sbatch --export=ALL,EL=Fe,MAG=fm,K="8 8 10" bcc2.sbatch
#   K = "8 8 10": same k-point density as the 40-atom SQS (4x4x2 on 2x2x5 cells); K = "16 16 16": converged.
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
export QE_NK=16
TAG=$(echo $K | tr ' ' 'x')
$PY -I elastic_qe.py $EL $K 60 runs/${EL}_${MAG}_k${TAG}_e60 --mag $MAG --np 32
