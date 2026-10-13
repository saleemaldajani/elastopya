#!/bin/bash
#SBATCH -J qe_bandpath
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Band structure along Gamma-H-N-Gamma-P-H (bands_path.py pure-inputs). Usage: sbatch --export=ALL,CASE=Fe_fm bandpath.sbatch
cd /orcd/pool/004/sdajani/fecr_qe/bandpath/$CASE
source ../../qe75_env.sh
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 32 $QE_PW -nk 8 -in scf.in > scf.out
grep -q "JOB DONE" bands.out 2>/dev/null || $QE_MPIRUN -np 32 $QE_PW -nk 16 -in bands.in > bands.out
cp tmp/fs.save/data-file-schema.xml bands_data.xml && echo BANDS_OK
