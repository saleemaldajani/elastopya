#!/bin/bash
# Overnight watcher: every 20 min pull Engaging results into dft_qe/engaging/ and push them to the repo;
# resubmit 40-atom SQS jobs that hit the 12 h limit (the driver resumes from finished steps), max 2 times each.
cd "$(dirname "$0")"
REMOTE=/orcd/pool/004/sdajani/fecr_qe
CNT="${TMPDIR:-/tmp}/engaging_resubmits"; mkdir -p "$CNT"
while true; do
  ./pull_engaging.sh > /dev/null 2>&1 && (cd .. && python3 sync_repo.py "Engaging QE results (automatic pull)" 2>&1 | tail -1)
  for line in $(ssh engaging "sacct -u \$USER -S today -n -X -o JobName%30,State%12 | awk '\$2==\"TIMEOUT\" && \$1 ~ /qe_sqs40_x/ {print \$1}'" 2>/dev/null | sort -u); do
    X=${line#qe_sqs40_x}
    if ssh engaging "test -f $REMOTE/runs/sqs40_x${X}_k442_e60/result.json"; then continue; fi
    if ssh engaging "squeue -h -u \$USER -n $line | grep -q ."; then continue; fi
    n=$(cat "$CNT/$X" 2>/dev/null || echo 0)
    if [ "$n" -lt 2 ]; then
      ssh engaging "cd $REMOTE && sbatch --export=ALL,X=$X -J $line sqs40.sbatch" && echo $((n+1)) > "$CNT/$X"
      echo "$(date '+%H:%M') resubmitted x=$X (attempt $((n+1)))"
    fi
  done
  left=$(ssh engaging "squeue -h -u \$USER | wc -l" 2>/dev/null || echo 1)
  echo "$(date '+%H:%M') jobs left on Engaging: $left"
  if [ "$left" = "0" ]; then ./pull_engaging.sh > /dev/null 2>&1; (cd .. && python3 sync_repo.py "Engaging QE results (final pull)" | tail -1); echo ENGAGING_ALL_DONE; break; fi
  sleep 1200
done
