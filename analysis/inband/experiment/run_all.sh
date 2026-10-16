#!/bin/zsh
# Runs the 10k-orientation SAW sweeps for tensors.csv rows not yet done (6 workers). Skips finished outputs.
cd /Users/saleem/MIT/research/Scripta/mike_final_touches
E=analysis/inband/experiment; A=data/saw/angles_10k.csv; S=scripts/saw_random_orientations.py
for tag in fe36_meas_hom anch_decomp_voigt fe36_meas_aged10h anch_decomp_hill anch_decomp_reuss fe_pure cr_pure \
           vegard_decomp_reuss anch_partialAPT_voigt anch_alphap_x088; do
  [[ -s $E/runs/$tag.csv ]] && continue
  args=(${=$(awk -F, -v t=$tag '$1==t{print $2,$3,$4,$5}' $E/tensors.csv)})
  python3 $S run $args $A $E/runs/$tag.tmp.csv 10000 6 > $E/runs/$tag.log 2>&1 && mv $E/runs/$tag.tmp.csv $E/runs/$tag.csv
  echo "$tag $args :: $(tail -1 $E/runs/$tag.log)" >> $E/runs/progress.txt
done
echo ALLDONE >> $E/runs/progress.txt
