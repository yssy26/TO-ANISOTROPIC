#!/bin/bash
echo "=== tail log.b4export ==="
tail -20 /home/ys/b2_case_smoke/log.b4export
echo "=== key milestone count ==="
grep -cE 'ExplicitJToracle' /home/ys/b2_case_smoke/log.b4export
echo "=== process ==="
pgrep -a MTO_HF | head -2
