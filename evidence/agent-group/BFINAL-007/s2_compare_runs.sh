#!/bin/bash
# BFINAL-007 S2: byte-identical comparison of the three scratch runs.
# Usage: bash s2_compare_runs.sh <dir_pre> <dir_off> <dir_on>
set -u
PRE=$1; OFF=$2; ON=$3
echo "== pre vs off (production outputs must be byte-identical) =="
diff -r --exclude='Log.stageB7.*' "$PRE" "$OFF" > /tmp/s2_pre_off.diff 2>&1
echo "diff -r pre off (excluding Logs) exit=$?  ($(wc -l < /tmp/s2_pre_off.diff) diff lines)"
head -40 /tmp/s2_pre_off.diff
echo
echo "== off vs on (production outputs must be byte-identical; probe writes only to evidence dir) =="
diff -r --exclude='Log.stageB7.*' "$OFF" "$ON" > /tmp/s2_off_on.diff 2>&1
echo "diff -r off on (excluding Logs) exit=$?  ($(wc -l < /tmp/s2_off_on.diff) diff lines)"
head -40 /tmp/s2_off_on.diff
echo
echo "== Log comparison: pre vs off (must be byte-identical) =="
if cmp -s "$PRE/Log.stageB7.pre" "$OFF/Log.stageB7.pre"; then
  echo "Log.pre == Log.off : IDENTICAL"
else
  echo "Log.pre != Log.off : DIFFERENT"
  diff "$PRE/Log.stageB7.pre" "$OFF/Log.stageB7.pre" | head -20
fi
echo
echo "== Log comparison: off vs on (expected: only stageB7 probe Info lines differ) =="
diff "$OFF/Log.stageB7.pre" "$ON/Log.stageB7.pre" > /tmp/s2_off_on_log.diff 2>&1 || true
echo "log diff lines: $(wc -l < /tmp/s2_off_on_log.diff)"
head -30 /tmp/s2_off_on_log.diff
echo
echo "== RC check =="
grep -h "MTO_RC" "$PRE/Log.stageB7.pre" "$OFF/Log.stageB7.pre" "$ON/Log.stageB7.pre" 2>/dev/null
echo "S2_COMPARE_DONE"
