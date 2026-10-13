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
# Fermi-surface band data (the only files needed from tmp/): fermi/<case>/tmp/fs.save/data-file-schema.xml -> fermi/<case>/fs_data.xml
for c in $(ssh engaging "cd $REMOTE/fermi 2>/dev/null && ls -d */ 2>/dev/null" | tr -d /); do
  if ssh engaging "test -f $REMOTE/fermi/$c/tmp/fs.save/data-file-schema.xml && grep -q 'JOB DONE' $REMOTE/fermi/$c/bands.out 2>/dev/null"; then
    mkdir -p engaging/fermi/$c && scp -q engaging:$REMOTE/fermi/$c/tmp/fs.save/data-file-schema.xml engaging/fermi/$c/fs_data.xml
  fi
done
ssh engaging "cd $REMOTE && tar czf - --exclude='*/tmp' fermi fermi.sbatch 2>/dev/null" | tar xzf - -C engaging
rm -f engaging/pseudo/*.UPF; rmdir engaging/pseudo 2>/dev/null || true   # pseudopotentials are not redistributed
find engaging -name "*.out" -size +20M -print -delete
du -sh engaging
