#!/bin/bash
#SBATCH -J saw_qe
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Orientation-resolved (P)SAW over 10,000 grains (Du & Zhao method, saw_random_orientations.py) for a QE tensor.
# Usage: sbatch --export=ALL,TAG=qe_x0.20,C11=..,C12=..,C44=..,RHO=.. saw_qe.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/saw
../venv/bin/python saw_random_orientations.py run $C11 $C12 $C44 $RHO angles_10k.csv ${TAG}.csv 10000 32
