#!/bin/bash
#SBATCH -J qe_phdos
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o slurm-%x-%j.out
# Phonon DOS from existing DFPT force constants (phonons.py dos-input). Usage: sbatch --export=ALL,DIRS="Fe_fm Cr_afm" phonon_dos.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/phonon
source ../qe75_env.sh
for c in $DIRS; do (cd $c && $QE_MPIRUN -np 1 $(dirname $QE_PW)/matdyn.x -in matdyn_dos.in > matdyn_dos.out && ls dos.dat && echo $c DOS_OK); done
