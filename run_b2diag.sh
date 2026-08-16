#!/bin/bash
# B2/B3 diagnostic run: non-fatal adjoint + gradient-convergence checkpoint.

# Clean PATH first: drop inherited Windows entries (some contain spaces and
# break OpenFOAM's unquoted `export PATH=...:$PATH` in etc/bashrc).
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
unset FOAM_SIGFPE

cd /home/ys/b2_case_smoke || { echo "ERROR: cd failed"; exit 1; }

rm -rf 1 stageB2 adjointCheckpoint_*.tsv optimization_log.dat
rm -f log.b2diag

echo "=== B2/B3 diagnostic: $(date) ==="
echo "MTO_HF = $(command -v MTO_HF)"

MTO_HF > log.b2diag 2>&1
rc=$?
echo "RUN_EXIT=$rc"
echo "--- tail log.b2diag ---"
tail -60 log.b2diag
