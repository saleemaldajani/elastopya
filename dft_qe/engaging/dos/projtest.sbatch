#!/bin/bash
#SBATCH -J qe_projtest
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH -o slurm-%x-%j.out
cd /orcd/pool/004/sdajani/fecr_qe/dos/x0.95
source ../../qe75_env.sh
$QE_MPIRUN -np 16 $(dirname $QE_PW)/projwfc.x -nk 1 -in projwfc_tetra.in > projwfc_tetra.out
python3 -c "import numpy as np; a=np.loadtxt(\"pdosT.pdos_tot\"); print(\"max dosup\", a[:,1].max(), \"max dosdw\", a[:,2].max())"
