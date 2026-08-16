#!/bin/bash
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
unset FOAM_SIGFPE

cd /home/ys/b2_case_smoke || exit 1

# Backup current config
cp constant/optProperties constant/optProperties.b2diag_run_backup

# Configure for B4 export (full explicit J^T + RHS, oracle verification)
sed -i 's/^solveFlowAdjoints.*/solveFlowAdjoints                 true;/' constant/optProperties
sed -i 's/^stageB2Enabled.*/stageB2Enabled false;/' constant/optProperties
sed -i 's/^stageB4JacobianProbe.*/stageB4JacobianProbe true;/' constant/optProperties
sed -i 's/^discreteExportOnly.*/discreteExportOnly true;/' constant/optProperties
sed -i 's/^discreteUseExplicitSolution.*/discreteUseExplicitSolution false;/' constant/optProperties
sed -i 's#^discreteExplicitSolutionFile.*#discreteExplicitSolutionFile "/home/ys/b2_case_smoke/explicitSol";#' constant/optProperties
sed -i 's#^discreteExplicitMatrixFile.*#discreteExplicitMatrixFile "/home/ys/b2_case_smoke/explicitJT.mtx";#' constant/optProperties

echo "=== reconfigured key lines ==="
grep -nE 'solveFlowAdjoints|stageB2Enabled|stageB4JacobianProbe|discreteExportOnly|discreteUseExplicitSolution|discreteExplicitSolutionFile|discreteExplicitMatrixFile|mmaUpdateEnabled' constant/optProperties

# Clean stale exports
rm -f explicitJT.mtx explicitRhs*.mtx explicitSol*.mtx
rm -rf 1 stageB2 adjointCheckpoint_*.tsv
rm -f log.b4export

echo "=== B4 export run: $(date) ==="
MTO_HF > log.b4export 2>&1
echo "RUN_EXIT=$?"
echo "=== ExplicitJToracle / export / dot-test lines ==="
grep -nE 'ExplicitJToracle|dot-test|Reduced cold-flow|BlockDot|FatalError|Loaded explicit|Export|explicitJT' log.b4export | head -50
echo "=== tail ==="
tail -25 log.b4export
