#!/bin/bash
cd /home/ys/b2_case_smoke
echo "=== explicit solution load / residual ==="
grep -nE 'Loaded explicit solution|relRes=|true relative residual|ExplicitJToracle' log.b4import | head -12
