#!/bin/bash
# Build the fixed solver, then run a mini-FD (no flow adjoints) to validate
# the frozen-primal oracle and obtain reference FD_J / FD_gDP.

# 1) Clean PATH: drop inherited Windows entries (spaces break OpenFOAM bashrc).
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
unset FOAM_SIGFPE

# 2) Build
cd /home/ys/TO-ANISOTROPIC/src || exit 1
wmake > /home/ys/TO-ANISOTROPIC/src/validation_stage_b2diag_build2.log 2>&1
build_rc=$?
if [ "$build_rc" -ne 0 ]; then
    echo "BUILD FAILED (rc=$build_rc)"
    tail -50 /home/ys/TO-ANISOTROPIC/src/validation_stage_b2diag_build2.log
    exit 1
fi
echo "BUILD PASSED"

# 3) Configure case: solveFlowAdjoints false -> skip flow adjoints (FD only)
cd /home/ys/b2_case_smoke || exit 1
sed -i 's/^solveFlowAdjoints.*/solveFlowAdjoints                 false;/' constant/optProperties
echo "--- solveFlowAdjoints ---"
grep -n 'solveFlowAdjoints' constant/optProperties

# 4) Clean stale products and run
rm -rf 1 stageB2 adjointCheckpoint_*.tsv optimization_log.dat
rm -f log.minifd
echo "=== mini-FD run: $(date) ==="
MTO_HF > log.minifd 2>&1
echo "RUN_EXIT=$?"
echo "--- tail log.minifd ---"
tail -50 log.minifd
