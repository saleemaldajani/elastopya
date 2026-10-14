#!/bin/bash
# Exit (waking Claude) when a key Engaging result appears or a qe job fails; checks every 10 min, up to 12 h.
KEYS="runs/sqs40_x0.05_k442_e60/result.json runs/sqs40_x0.10_k442_e60/result.json runs/sqs40_x0.70_k442_e60/result.json runs/sqs40_x0.05_k442_e60_aexp/result.json runs/sqs40_x0.10_k442_e60_aexp/result.json runs/sqs40_x0.65_k442_e60_aexp/result.json runs/sqs40_x0.80_k442_e60_aexp_crafm/result.json runs/sqs40_x0.80_k442_e60_crafm2/result.json phonopy_cubic/Fe_cub/d000.out phonopy_cubic/Fe_cub_aexp/d000.out"
have0=$(ssh engaging "cd /orcd/pool/004/sdajani/fecr_qe; for f in $KEYS; do [ -f \$f ] && grep -q . \$f 2>/dev/null && { [[ \$f == *.out ]] && ! grep -q 'JOB DONE' \$f && continue; echo \$f; }; done; true" 2>/dev/null)
start=$(date +%s)
while [ $(( $(date +%s) - start )) -lt 43200 ]; do
  sleep 600
  now=$(ssh engaging "cd /orcd/pool/004/sdajani/fecr_qe; for f in $KEYS; do [ -f \$f ] && { [[ \$f == *.out ]] && ! grep -q 'JOB DONE' \$f && continue; echo \$f; }; done; true" 2>/dev/null) || { echo "CONNECTION PROBLEM $(date)"; exit 0; }
  new=$(comm -13 <(echo "$have0" | sort) <(echo "$now" | sort))
  fail=$(ssh engaging "sacct -u \$USER -S \$(date -d '-10 min' +%Y-%m-%dT%H:%M) -X -n -o JobName%26,State%14 | grep -E 'FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL' | grep qe_" 2>/dev/null)
  use=$(ssh engaging "df --output=pcent /orcd/pool/004/sdajani | tail -1 | tr -dc 0-9" 2>/dev/null)
  if [ -n "$use" ] && [ "$use" -ge 90 ]; then echo "DISK ${use}% full"; date; exit 0; fi
  if [ -n "$new" ] || [ -n "$fail" ]; then echo "NEW: $new"; echo "FAILED: $fail"; date; exit 0; fi
done
echo "no key results in 12 h"
