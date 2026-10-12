#!/bin/bash
# Orientation-resolved (P)SAW spectra (Du & Zhao method, 10k random orientations as in Acta Fig. 6)
# for the Fe-Cr tensors used in the paper. Density: linear Fe (7870) - Cr (7200) interpolation.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}
ANG=data/saw/angles_10k.csv
OUT=data/saw
run() { # tag C11 C12 C44 x
  local rho; rho=$(python3 -c "print(7870*(1-$5)+7200*$5)")
  [ -f "$OUT/$1.csv" ] || $PY scripts/saw_random_orientations.py run "$2" "$3" "$4" "$rho" "$ANG" "$OUT/$1.csv" 10000 8 | sed "s/^/$1: /"
}
# DFT (sheet 'Analyzed data' D,G,J)
run dft_x0.00 324.434353 169.838953 127.7242 0
run dft_x0.20 400.025537 164.57738 120.574293 0.2
run dft_x0.50 347.46247 186.14903 123.354977 0.5
run dft_x0.95 499.051533 145.051597 99.27015 0.95
run dft_x1.00 568.728963 158.14884 100.913353 1
# Literature (Razumovskiy et al., aligned rows of the sheet 'Raw data')
while IFS=, read -r tag x c11 c12 c44; do run "$tag" "$c11" "$c12" "$c44" "$x"; done < data/saw/lit_tensors.csv
# MD random alloys (this work)
while IFS=, read -r tag x c11 c12 c44; do run "$tag" "$c11" "$c12" "$c44" "$x"; done < data/saw/md_tensors.csv
echo SAW_FECR_DONE
