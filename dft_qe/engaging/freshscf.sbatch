#!/bin/bash
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mem=100G
#SBATCH -t 12:00:00
#SBATCH --requeue
#SBATCH --exclude=node2701
#SBATCH -o slurm-%x-%j.out
# Single fresh-basis SCF at a relaxed structure (2026-10-29): energy of the Fe-80Cr AF-branch structure in the same
# basis convention as the low-moment run, whose closing check SCF converged. Usage: sbatch --export=ALL,D=<run dir> freshscf.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/runs/$D
source ../../qe75_env.sh
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 60 $QE_PW -nk 20 -in scf.in > scf.out
grep "^!" scf.out | tail -1
