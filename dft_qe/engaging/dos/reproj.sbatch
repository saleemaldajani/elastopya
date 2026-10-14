#!/bin/bash
#SBATCH -J qe_reproj
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Re-run projwfc (tetrahedron DOS, no pools) for DOS runs whose first projection used a Gaussian width (all-zero PDOS)
cd /orcd/pool/004/sdajani/fecr_qe/dos
source ../qe75_env.sh
for x in 0.35 0.50 0.65; do
  (cd x$x && $QE_MPIRUN -np 16 $(dirname $QE_PW)/projwfc.x -nk 1 -in projwfc.in > projwfc.out)
  python3 -c "import numpy as np; a=np.loadtxt(\"x$x/pdosT.pdos_tot\"); print(\"x$x max dosup\", a[:,1].max())"
done
