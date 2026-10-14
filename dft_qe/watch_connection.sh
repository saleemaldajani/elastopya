#!/bin/bash
# Keep the shared `ssh engaging` connection alive and the Mac awake while Engaging work is in progress.
#  - every 120 s: a no-op command over the shared connection (never prompts: BatchMode), state logged to $LOG
#  - caffeinate -i for as long as this script runs: no idle system sleep (the display may still sleep)
#  - on the first OK -> LOST transition: one macOS notification asking for `ssh engaging` + Duo (it cannot log in itself)
# Stop it with:  pkill -f watch_connection.sh   (caffeinate exits with it)
LOG=${1:-/tmp/engaging_keepalive.log}
caffeinate -i -w $$ &
state=OK
while true; do
  if ssh -o BatchMode=yes -o ConnectTimeout=20 engaging true 2>/dev/null; then
    new=OK
  else
    new=LOST
  fi
  echo "$(date '+%Y-%m-%d %H:%M') $new" > "$LOG"
  if [ "$new" = LOST ] && [ "$state" = OK ]; then
    osascript -e 'display notification "Run ssh engaging in Terminal and approve Duo (or ask Claude to open the login window)." with title "Engaging connection lost"' 2>/dev/null
  fi
  state=$new
  sleep 120
done
