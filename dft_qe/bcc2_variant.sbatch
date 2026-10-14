#!/bin/bash
#SBATCH -J qe_bcc2
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o slurm-%x-%j.out
# As bcc2.sbatch, written to runs/${EL}_${MAG}_k${TAG}_e60_${SUFFIX}, for variants set through the environment
# (e.g. QE_FIXED_A=exp: experimental lattice parameter). Same physics settings otherwise.
# Usage: sbatch --export=ALL,EL=Fe,MAG=fm,K="16 16 16",SUFFIX=aexp,QE_FIXED_A=exp bcc2_variant.sbatch
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
export QE_NK=16
TAG=$(echo $K | tr ' ' 'x')
$PY -I elastic_qe.py $EL $K 60 runs/${EL}_${MAG}_k${TAG}_e60_${SUFFIX} --mag $MAG --np 32
