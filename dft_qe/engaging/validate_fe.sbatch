#!/bin/bash
#SBATCH -J qe_fe_check
#SBATCH -p mit_quicktest
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=32G
#SBATCH -t 00:15:00
#SBATCH -o slurm-%x-%j.out
# Validation: bcc Fe, k 12^3, 60 Ry. Local QE 7.5 gave a = 2.827 A, C11/C12/C44 = 311.0/159.2/114.1 GPa.
cd /orcd/pool/004/sdajani/fecr_qe
source qe_env.sh
export QE_NK=16
$PY -I elastic_qe.py Fe 12 12 12 60 runs/Fe_k12_e60 --np 32
