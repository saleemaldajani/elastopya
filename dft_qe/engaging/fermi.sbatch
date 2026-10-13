#!/bin/bash
#SBATCH -J qe_fermi
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o slurm-%x-%j.out
# Spin-resolved Fermi-surface cross-sections (fermi_slices.py): SCF, then bands on a k_z = 0 grid.
# Usage: sbatch --export=ALL,CASE=Fe_fm fermi.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/fermi/$CASE
source ../../qe75_env.sh
$QE_MPIRUN -np 32 $QE_PW -nk 8 -in scf.in > scf.out
$QE_MPIRUN -np 32 $QE_PW -nk 16 -in bands.in > bands.out
grep -E "the Fermi energy|the spin up/dw Fermi|JOB DONE" scf.out bands.out | tail -3
