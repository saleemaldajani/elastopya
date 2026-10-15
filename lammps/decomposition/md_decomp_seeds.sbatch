#!/bin/bash
#SBATCH -J md_decomp
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH --requeue
#SBATCH -o logs/slurm-%x-%A_%a.out
# Extra seeds (alloys 2-5) for the MD spinodal-decomposition series (run_decomposition.sh, 25^3 cells), so that
# Fig. 3 carries the seed-to-seed spread. One array task per line of tasks_seeds.txt ("x alloy").
# Task 0 repeats alloy 1 at x = 0.36 (already computed locally with LAMMPS 22 Jul 2025) into check/ as a
# cross-version consistency check of this cluster's LAMMPS 2 Aug 2023. Tasks 1-20: the five compositions of
# Fig. 3b; tasks 21-48: the other seven compositions of Fig. 3a.
cd /orcd/pool/004/sdajani/fecr_md/decomposition
source /etc/profile >/dev/null 2>&1
module load community-modules; module load lammps/20230802.3
export OMPI_MCA_opal_warn_on_missing_libcuda=0
export MPIRUN=/orcd/software/core/001/spack/pkg/openmpi/4.1.4/zahpnmk/bin/mpirun
read X A < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tasks_seeds.txt)
if [ "$SLURM_ARRAY_TASK_ID" = "0" ]; then          # same driver, output under check/ (template path unchanged)
  sed 's|^D="$HERE/runs/|D="$HERE/check/|' run_decomposition.sh > run_check.sh && chmod +x run_check.sh
  ./run_check.sh $(which lmp) 16 $X 25 $A
else
  ./run_decomposition.sh $(which lmp) 16 $X 25 $A
fi
echo TASK_DONE x=$X alloy=$A
