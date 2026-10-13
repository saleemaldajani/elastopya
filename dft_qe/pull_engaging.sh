#!/bin/bash
# Copy the Engaging QE work (inputs, outputs, job scripts, Slurm logs) into dft_qe/engaging/.
# Wavefunction scratch (tmp/) and the QE 7.5 source/build (qe75/) are not copied.
# Needs a live `ssh engaging` connection (see ~/.ssh/config).
set -euo pipefail
cd "$(dirname "$0")"
REMOTE=/orcd/pool/004/sdajani/fecr_qe
mkdir -p engaging
ssh engaging "cd $REMOTE && tar czf - --exclude='*/tmp' --exclude='*/tmp/*' --exclude='qe75' --exclude='*.wfc*' \
  --exclude='*.save' *.py *.sh *.sbatch sqs pseudo/*.UPF runs slurm-*.out 2>/dev/null" | tar xzf - -C engaging
rm -f engaging/pseudo/*.UPF; rmdir engaging/pseudo 2>/dev/null || true   # pseudopotentials are not redistributed
find engaging -name "*.out" -size +20M -print -delete
du -sh engaging
