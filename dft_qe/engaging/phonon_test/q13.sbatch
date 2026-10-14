#!/bin/bash
#SBATCH -J qe_phtest
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH -o slurm-%x-%j.out
cd /orcd/pool/004/sdajani/fecr_qe/phonon_test/Cr_nm_q13
source ../../qe75_env.sh
B=$(dirname $QE_PW)
$QE_MPIRUN -np 32 $QE_PW -nk 8 -in scf.in > scf.out
$QE_MPIRUN -np 32 $B/ph.x -nk 8 -in ph.in > ph.out
grep "freq (" ph.out
