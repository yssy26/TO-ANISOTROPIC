#!/bin/bash
cd /home/ys/b2_case_smoke
echo "=== ExplicitJToracle / dot-test / export lines ==="
grep -nE 'ExplicitJToracle|Reduced cold-flow operator transpose|BlockDot|dot test|FatalError|Loaded explicit' log.b4export | head -40
echo "=== exported files ==="
ls -la explicitJT.mtx explicitRhs*.mtx 2>&1
echo "=== explicitJT.mtx header (first 3 lines) ==="
head -3 explicitJT.mtx 2>&1
echo "=== explicitRhs_pressureDrop.mtx header ==="
head -2 explicitRhs_pressureDrop.mtx 2>&1
