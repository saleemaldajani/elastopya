#!/bin/bash
#SBATCH -J qe75_ph
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# QE 7.5 phonon codes (ph.x, q2r.x, matdyn.x) in the existing build
source /etc/profile >/dev/null 2>&1
module load gcc/12.2.0 openmpi/4.1.4 openblas/0.3.26 fftw/3.3.10 cmake/3.27.9
cd /orcd/pool/004/sdajani/fecr_qe/qe75/q-e-qe-7.5/build
make -j16 ph
ls bin/ph.x bin/q2r.x bin/matdyn.x && echo PH_OK
