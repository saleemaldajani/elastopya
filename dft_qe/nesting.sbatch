#!/bin/bash
#SBATCH -J qe_nesting
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Cr nesting function (nesting.py): scf (24^3) -> bands on the full grid (nosym) -> FFT autocorrelation.
cd /orcd/pool/004/sdajani/fecr_qe/nesting/Cr_nm
source ../../qe75_env.sh
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 32 $QE_PW -nk 8 -in scf.in > scf.out
grep -q "JOB DONE" bands.out 2>/dev/null || $QE_MPIRUN -np 32 $QE_PW -nk 32 -in bands.in > bands.out
cp tmp/fs.save/data-file-schema.xml bands_data.xml
../../venv/bin/python -I ../../nesting.py analyse . > analyse.log 2>&1; cat analyse.log
