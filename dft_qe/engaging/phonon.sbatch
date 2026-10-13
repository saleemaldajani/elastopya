#!/bin/bash
#SBATCH -J qe_phonon
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=96G
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# DFPT phonon dispersion (phonons.py inputs). Usage: sbatch --dependency=afterok:<build> --export=ALL,CASE=Fe_fm phonon.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/phonon/$CASE
source ../../qe75_env.sh
B=$(dirname $QE_PW)
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 48 $QE_PW -nk 12 -in scf.in > scf.out
grep -q "JOB DONE" ph.out 2>/dev/null || $QE_MPIRUN -np 48 $B/ph.x -nk 12 -in ph.in > ph.out
$B/q2r.x -in q2r.in > q2r.out && $B/matdyn.x -in matdyn.in > matdyn.out && ls freq.gp && echo PHONON_OK
