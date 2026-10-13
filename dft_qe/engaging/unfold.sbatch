#!/bin/bash
#SBATCH -J qe_unfold
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=250G
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH -o slurm-%x-%j.out
# Unfolded Fermi surface of a 40-atom Fe-rich SQS (unfold_fs.py): scf -> bands at SC K-points -> BandUPpy.
# Usage: sbatch --export=ALL,X=0.05 unfold.sbatch     Finished steps are skipped when requeued.
cd /orcd/pool/004/sdajani/fecr_qe/unfold/x${X}
source ../../qe75_env.sh
grep -q "JOB DONE" scf.out 2>/dev/null || $QE_MPIRUN -np 64 $QE_PW -nk 8 -in scf.in > scf.out
grep -q "JOB DONE" bands.out 2>/dev/null || $QE_MPIRUN -np 64 $QE_PW -nk 16 -in bands.in > bands.out
[ -f unfolded_dw.npz ] || ../../venv/bin/python ../../unfold_fs.py unfold . > unfold.log 2>&1
tail -3 unfold.log
