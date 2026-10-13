#!/bin/bash
#SBATCH -J qe75_build
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o slurm-%x-%j.out
# Quantum ESPRESSO 7.5 with the same options as the local build (MPI, no OpenMP/ScaLAPACK, OpenBLAS, FFTW3)
source /etc/profile >/dev/null 2>&1
module load gcc/12.2.0 openmpi/4.1.4 openblas/0.3.26 fftw/3.3.10 cmake/3.27.9
module list
cd /orcd/pool/004/sdajani/fecr_qe/qe75/q-e-qe-7.5
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc \
  -DQE_ENABLE_MPI=ON -DQE_ENABLE_OPENMP=OFF -DQE_ENABLE_SCALAPACK=OFF -DQE_ENABLE_LIBXC=OFF \
  -DQE_FFTW_VENDOR=FFTW3 -DBLA_VENDOR=OpenBLAS -DQE_ENABLE_TEST=OFF
make -j16 pw
ls -la bin/pw.x && echo BUILD_OK
