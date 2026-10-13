#!/bin/bash
#SBATCH -J qe_dos
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mem=200G
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Spin/atom-projected DOS of the relaxed 40-atom SQS (dos_sqs.py). Runs after the elastic job of the same x.
# Usage: sbatch --dependency=afterany:<sqs job id> --export=ALL,X=0.20 dos.sbatch
# Optional: TAG=aexp reads runs/sqs40_x${X}_k442_e60_${TAG} and writes dos/x${X}_${TAG}; NBND=1.6 widens the window (pdosW)
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
D=dos/x${X}${TAG:+_$TAG}
RUN=runs/sqs40_x${X}_k442_e60${TAG:+_$TAG}
[ -f $D/scf.in ] || QE_NBND_FACTOR=${NBND:-1.2} $PY -I dos_sqs.py prepare sqs/x${X}.vasp ${X} $RUN $D
PROJ=$(dirname $QE_PW)/projwfc.x
cd $D
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 60 $QE_PW -nk 20 -in scf.in > scf.out
grep -q "JOB DONE" nscf.out 2>/dev/null || $QE_MPIRUN -np 60 $QE_PW -nk 20 -in nscf.in > nscf.out
grep -q "JOB DONE" nscf.out || { echo NSCF_FAILED; exit 1; }
grep -q "JOB DONE" projwfc.out 2>/dev/null || $QE_MPIRUN -np 16 $PROJ -nk 1 -in projwfc.in > projwfc.out
[ -n "$NBND" ] && for f in pdosT.*; do mv "$f" "pdosW${f#pdosT}"; done
ls pdos?.pdos_tot && echo DOS_OK
