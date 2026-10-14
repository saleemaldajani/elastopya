#!/bin/bash
#SBATCH -J qe_frozen_Fe
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o slurm-%x-%j.out
# Frozen-phonon check of the Fe transverse mode at q = (0,0,1/3) 2pi/a (6-atom 1x1x3 cubic cell), +-0.01 A.
cd /orcd/pool/004/sdajani/fecr_qe/frozen_Fe
source ../qe75_env.sh
for t in p m; do $QE_MPIRUN -np 32 $QE_PW -nk 16 -in fz_$t.in > fz_$t.out; done
grep -A8 "Forces acting" fz_p.out | head -9
