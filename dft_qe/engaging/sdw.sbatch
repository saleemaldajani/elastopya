#!/bin/bash
#SBATCH -J qe_sdw
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=200G
#SBATCH -t 2-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Cr SDW/CDW supercell (sdw_cr.py). PP_NP sets the pp.x ranks: 16 segfaulted on the large SDW cell (2026-10-14), use 4.
# Cr SDW/CDW supercell (sdw_cr.py). Usage: sbatch -J qe_sdw_<kind> --export=ALL,K=sdw sdw.sbatch   (K = nm | afm | sdw)
cd /orcd/pool/004/sdajani/fecr_qe/sdw/$K
source ../../qe75_env.sh
B=$(dirname $QE_PW)
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 64 $QE_PW -nk 16 -in scf.in > scf.out
grep -q "convergence NOT achieved" scf.out && { echo SCF_NOT_CONVERGED; exit 1; }
$QE_MPIRUN -np ${PP_NP:-16} $B/pp.x -in pp.in > pp.out && $QE_MPIRUN -np 1 $B/average.x < avg.in > avg.out && ls avg.dat && echo SDW_OK
