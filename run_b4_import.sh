#!/bin/bash
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
unset FOAM_SIGFPE

cd /home/ys/b2_case_smoke || exit 1

# Configure for B4 import (load exact explicit solution, verify true residual,
# compute gsensPressureDrop via the main-loop sensitivity chain)
sed -i 's/^solveFlowAdjoints.*/solveFlowAdjoints                 true;/' constant/optProperties
sed -i 's/^stageB2Enabled.*/stageB2Enabled false;/' constant/optProperties
sed -i 's/^stageB4JacobianProbe.*/stageB4JacobianProbe true;/' constant/optProperties
sed -i 's/^discreteExportOnly.*/discreteExportOnly true;/' constant/optProperties
sed -i 's/^discreteUseExplicitSolution.*/discreteUseExplicitSolution true;/' constant/optProperties
# add the includes-deviatoric gate switch (required by validateMmaUnlockGate)
grep -q 'discreteExplicitSolutionIncludesDeviatoric' constant/optProperties \
  || sed -i '/^discreteUseExplicitSolution/a discreteExplicitSolutionIncludesDeviatoric true;' constant/optProperties
sed -i 's/^discreteExplicitSolutionIncludesDeviatoric.*/discreteExplicitSolutionIncludesDeviatoric true;/' constant/optProperties
sed -i 's#^discreteExplicitSolutionFile.*#discreteExplicitSolutionFile "/home/ys/b2_case_smoke/explicitSol";#' constant/optProperties

echo "=== reconfigured ==="
grep -nE 'solveFlowAdjoints|stageB2Enabled|stageB4JacobianProbe|discreteExportOnly|discreteUseExplicitSolution|discreteExplicitSolutionIncludesDeviatoric|discreteExplicitSolutionFile' constant/optProperties

rm -rf 1 stageB2 adjointCheckpoint_*.tsv
rm -f log.b4import

echo "=== B4 import run: $(date) ==="
MTO_HF > log.b4import 2>&1
echo "RUN_EXIT=$?"
echo "=== explicit solution load / residual lines ==="
grep -nE 'Loaded explicit solution|relRes=|ExplicitJToracle|FatalError|PRODRESID|true relative residual' log.b4import | head -40
echo "=== tail ==="
tail -25 log.b4import
