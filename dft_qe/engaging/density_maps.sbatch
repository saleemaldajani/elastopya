#!/bin/bash
#SBATCH -J qe_densmap
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Spin and deformation-charge maps of the SQS (density_maps.py). Usage: sbatch --export=ALL,DIRS="x0.20 x0.50" density_maps.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/dos
source ../qe75_env.sh
B=$(dirname $QE_PW)
for d in $DIRS; do (cd $d && for k in spin drho; do $QE_MPIRUN -np 16 $B/pp.x -in pp_$k.in > pp_$k.out; done; ls spin_110.dat drho_110.dat && echo $d MAPS_OK); done
