#!/bin/bash
# Keep the shared `ssh engaging` connection busy so it never idles out; log its state.
# Never prompts: if the shared connection is gone, BatchMode fails immediately and the state is LOST.
LOG=${1:-/tmp/engaging_keepalive.log}
while true; do
  if ssh -o BatchMode=yes -o ConnectTimeout=20 engaging true 2>/dev/null; then
    echo "$(date '+%Y-%m-%d %H:%M') OK" > "$LOG"
  else
    echo "$(date '+%Y-%m-%d %H:%M') LOST" > "$LOG"
  fi
  sleep 240
done
