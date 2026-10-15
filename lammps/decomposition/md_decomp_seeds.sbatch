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
# consistency check of the cluster build (same version as the local runs). Tasks 1-20: the five compositions of
# Fig. 3b; tasks 21-48: the other seven compositions of Fig. 3a.
cd /orcd/pool/004/sdajani/fecr_md/decomposition
source /etc/profile >/dev/null 2>&1
module load community-modules; module load gcc/12.2.0 openmpi/4.1.4
LMP=/orcd/pool/004/sdajani/fecr_md/lammps-src/build/lmp_mpi   # LAMMPS 22Jul2025 u6 + MC (build_lammps_engaging.sbatch)
export OMPI_MCA_opal_warn_on_missing_libcuda=0
export MPIRUN=$(which mpirun)
set -e
read X A < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tasks_seeds.txt)
if [ "$SLURM_ARRAY_TASK_ID" = "0" ]; then          # same driver, output under check/ (template path unchanged)
  sed 's|^D="$HERE/runs/|D="$HERE/check/|' run_decomposition.sh > run_check.sh && chmod +x run_check.sh
  ./run_check.sh $LMP 16 $X 25 $A
else
  ./run_decomposition.sh $LMP 16 $X 25 $A
fi
D=$( [ "$SLURM_ARRAY_TASK_ID" = "0" ] && echo check || echo runs )/x${X}_L25_a${A}
grep -q "C44all =" $D/stiff_phi0.500/log.lammps && echo TASK_DONE x=$X alloy=$A || { echo "TASK_FAILED x=$X alloy=$A"; exit 1; }
