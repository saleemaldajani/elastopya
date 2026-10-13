#!/bin/bash
# Detached guardian for the laptop QE runs (started in its own session so Claude restarts cannot stop it).
# 1) every 2 min: resume (SIGCONT) any stopped QE / MPI / driver / keep-alive / watcher process;
# 2) when the running Fe-50Cr driver (PID given as $1) exits with a result, run Fe-20Cr then Fe-95Cr.
DRIVER_PID=$1
cd "$(dirname "$0")"
PY=/private/tmp/claude-501/-Users-saleem-MIT-research-Scripta/37f8fe43-8e5f-4d11-883b-ebb44fa03560/scratchpad/pdfenv/bin/python
resume_stopped() {
  ps -ax -o pid=,stat=,command= | awk '$2 ~ /^T/' | grep -E "pw\.x|prterun|mpirun|elastic_sqs_qe|caffeinate|keepalive_engaging|watch_engaging" \
    | awk '{print $1}' | while read p; do kill -CONT "$p" 2>/dev/null && echo "$(date '+%H:%M') resumed $p"; done
}
while kill -0 "$DRIVER_PID" 2>/dev/null; do resume_stopped; sleep 120; done
if [ ! -f runs/sqs40_x0.50_k442_e60/result.json ]; then echo "$(date '+%H:%M') x0.50 produced no result; queue not started"; exit 1; fi
for X in 0.20 0.95; do
  START=$(date +%s)
  caffeinate -i $PY -I elastic_sqs_qe.py sqs/x$X.vasp $X 4 4 2 60 runs/sqs40_x${X}_k442_e60 --np 12 &
  QP=$!
  while kill -0 $QP 2>/dev/null; do resume_stopped; sleep 120; done
  wait $QP || echo "$(date '+%H:%M') FAILED x=$X"
  echo "$(date '+%H:%M') QE40 x$X elapsed $(( $(date +%s)-START )) s"
done
echo "$(date '+%H:%M') QE40 QUEUE DONE"
