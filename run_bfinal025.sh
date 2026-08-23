#!/bin/bash
# BFINAL-025 case runner: launches MTO_HF in background from the case dir.
# Usage: bash run_bfinal025.sh <caseDir> <logName>
#
# Protocol mirrors the B24 gauge run (Log.verify_b24_gauge.txt Exec line =
# /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF):
#   - clean PATH first, then source /opt/openfoam7/etc/bashrc
#   - FOAM_USER_APPBIN points at the freshly built binary
#   - unset FOAM_SIGFPE (adjoint NaN/Inf must not abort)
set -o pipefail
CASE="${1:?usage: run_bfinal025.sh <caseDir> <logName>}"
LOG="${2:?usage: run_bfinal025.sh <caseDir> <logName>}"

export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
unset FOAM_SIGFPE

BIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF
test -x "$BIN" || { echo "ERROR: binary not found: $BIN" >&2; exit 1; }

cd "$CASE"
nohup "$BIN" > "$LOG" 2>&1 &
echo "RUN_PID=$!"
