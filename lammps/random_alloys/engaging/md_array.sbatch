#!/bin/bash
#SBATCH -J md_rand
#SBATCH -p mit_preemptable
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --requeue
#SBATCH -o logs/slurm-%x-%A_%a.out
# MD (LAMMPS, Bonny 2011 EAM) elastic constants of random bcc Fe-Cr alloys, Fe-rich grid.
# One array task per line of tasks.txt ("x seed"); task 0 repeats a local run (x = 0.1) as a check.
cd /orcd/pool/004/sdajani/fecr_md/random_alloys
source /etc/profile >/dev/null 2>&1
module load community-modules; module load lammps/20230802.3
export OMPI_MCA_opal_warn_on_missing_libcuda=0
export MPIRUN=/orcd/software/core/001/spack/pkg/openmpi/4.1.4/zahpnmk/bin/mpirun
read X SEED < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tasks.txt)
if [ "$SLURM_ARRAY_TASK_ID" = "0" ]; then mkdir -p check && cp -r template check/ && cp run_random_alloys.sh check/ && cd check; fi
./run_random_alloys.sh $(which lmp) 32 $X $SEED
grep -h "C11all =" runs/x${X}_s${SEED}/log.lammps | tail -1
