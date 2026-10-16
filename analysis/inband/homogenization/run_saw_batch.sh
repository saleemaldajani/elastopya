#!/bin/bash
# Usage: run_saw_batch.sh <workers> <tag1> <tag2> ...   (tag = <model>_<scheme>, tensor saw_runs/<tag>.npy)
# rho is taken from fecr36_tensors.csv. Output saw_runs/<tag>.csv (10^4 grains). Skips existing outputs.
cd "$(dirname "$0")"
W=$1; shift
for tag in "$@"; do
  out=saw_runs/$tag.csv
  [ -f "$out" ] && { echo "skip $tag"; continue; }
  model=${tag%%_*}; scheme=${tag#*_}
  rho=$(python3 -c "import csv;print([r['rho'] for r in csv.DictReader(open('fecr36_tensors.csv')) if r['model']=='$model' and r['scheme']=='$scheme'][0])")
  python3 saw_general.py saw_runs/$tag.npy $rho ../../../data/saw/angles_10k.csv $out 10000 $W
done
echo BATCH_DONE
