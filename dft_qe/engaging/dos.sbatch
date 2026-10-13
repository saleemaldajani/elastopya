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
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
D=dos/x${X}
[ -f $D/scf.in ] || $PY -I dos_sqs.py prepare sqs/x${X}.vasp ${X} runs/sqs40_x${X}_k442_e60 $D
PROJ=$(dirname $QE_PW)/projwfc.x
cd $D
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 60 $QE_PW -nk 20 -in scf.in > scf.out
grep -q "JOB DONE" nscf.out 2>/dev/null || $QE_MPIRUN -np 60 $QE_PW -nk 20 -in nscf.in > nscf.out
grep -q "JOB DONE" projwfc.out 2>/dev/null || $QE_MPIRUN -np 60 $PROJ -nk 20 -in projwfc.in > projwfc.out
ls pdos.pdos_tot && echo DOS_OK
