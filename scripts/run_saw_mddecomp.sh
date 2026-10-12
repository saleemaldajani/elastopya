#!/bin/bash
# Orientation-resolved (P)SAW for the fully decomposed (phi = 1) MD structures: their effective cubic
# tensors (data/md_decomposition.csv) through the same Du-method 10k-grain calculation.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}
$PY scripts/collect_md_decomposition.py > /dev/null
$PY -I -c "
import csv
rows=[r for r in csv.DictReader(open('data/md_decomposition.csv')) if float(r['phi'])==1.0]
with open('data/saw/mddec_tensors.csv','w') as f:
    for r in rows:
        x=float(r['x_cr']); f.write(f\"mddec_x{x:.2f},{x},{float(r['c11']):.4f},{float(r['c12']):.4f},{float(r['c44']):.4f}\n\")
"
while IFS=, read -r tag x c11 c12 c44; do
  rho=$(python3 -c "print(7870*(1-$x)+7200*$x)")
  [ -f data/saw/$tag.csv ] || $PY scripts/saw_random_orientations.py run "$c11" "$c12" "$c44" "$rho" data/saw/angles_10k.csv data/saw/$tag.csv 10000 12 | sed "s/^/$tag: /"
done < data/saw/mddec_tensors.csv
echo MDDEC_SAW_DONE
