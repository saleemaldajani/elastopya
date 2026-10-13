#!/bin/bash
#SBATCH -J qe_phpost
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# q2r.x + matdyn.x after a finished ph.x run (MPI binaries must go through mpirun on Engaging).
# Usage: sbatch [--dependency=afterany:<phonon job>] --export=ALL,CASE=Fe_fm phonon_post.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/phonon/$CASE
source ../../qe75_env.sh
B=$(dirname $QE_PW)
grep -q "JOB DONE" ph.out || { echo "ph.x not finished"; exit 1; }
$QE_MPIRUN -np 1 $B/q2r.x -in q2r.in > q2r.out && $QE_MPIRUN -np 1 $B/matdyn.x -in matdyn.in > matdyn.out
ls freq.gp && echo PHONON_OK
