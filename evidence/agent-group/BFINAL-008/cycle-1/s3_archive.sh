#!/bin/bash
# BFINAL-008 S3 — archive run-1 (probe) + run-2 (production) artifacts + logs
set -u
ART=/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-008/cycle-1/artifacts
S3="$ART/s3"
mkdir -p "$S3"
cd /home/ys/dsH

# run logs
cp -f b8_scratch_r1/Log.stageB8.r1.txt "$S3/" 2>/dev/null
cp -f b8_scratch_prod/Log.stageB8.prod.txt "$S3/" 2>/dev/null

# corrected explicit operator exports (run-1 probe path)
cp -f b8_scratch_r1/explicitJT.mtx "$S3/explicitJT.mtx" 2>/dev/null
cp -f b8_scratch_r1/explicitRhs_pressureDrop.mtx "$S3/" 2>/dev/null
cp -f b8_scratch_r1/explicitRhs_thermalCoupling.mtx "$S3/" 2>/dev/null

# P6 (stageB5 G1/G2) + P7 (stageB6 R_x) raw exports from run-1
cp -f b8_scratch_r1/stageB5_*.mtx "$S3/" 2>/dev/null
cp -f b8_scratch_r1/stageB6_*.mtx "$S3/" 2>/dev/null

# stageB8 probe exports are already in $ART (written by the probe itself);
# refresh the sha256 manifest for everything under $ART.
cd "$ART"
sha256sum $(find . -type f | sort) > sha256_artifacts.txt 2>/dev/null
echo "=== artifact inventory (S3) ==="
ls -la "$S3" | head -60
echo "=== total sha256 manifest ==="
cat sha256_artifacts.txt | head -80
