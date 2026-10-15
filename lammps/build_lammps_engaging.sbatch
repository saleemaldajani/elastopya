#!/bin/bash
#SBATCH -J build_lammps
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o slurm-%x-%j.out
# LAMMPS 22 Jul 2025 update 6 -- the version of the local MD decomposition runs -- built on Engaging with the MANYBODY
# (EAM) and MC (fix atom/swap) packages. The cluster module lammps/20230802.3 lacks MC, which the spinodal-decomposition
# workflow (run_decomposition.sh) needs for its Monte Carlo atom swaps.
set -e
cd /orcd/pool/004/sdajani/fecr_md
source /etc/profile >/dev/null 2>&1
module load community-modules; module load cmake/3.27.9 gcc/12.2.0 openmpi/4.1.4
[ -d lammps-src ] || git clone --depth 1 --branch stable_22Jul2025_update6 https://github.com/lammps/lammps.git lammps-src
mkdir -p lammps-src/build && cd lammps-src/build
cmake ../cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MPI=on -DPKG_MANYBODY=on -DPKG_MC=on -DLAMMPS_MACHINE=mpi
make -j 16
ls -la lmp_mpi && ./lmp_mpi -h | grep -E "^Large-scale|MANYBODY|^MC" | head -5
echo BUILD_OK
