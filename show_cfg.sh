#!/bin/bash
cd /home/ys/b2_case_smoke
echo "=== CURRENT optProperties (key lines) ==="
grep -nE 'mmaUpdateEnabled|stageBEnabled|stageB2Enabled|stageB4JacobianProbe|discreteExportOnly|discreteUseExplicitSolution|discreteExplicitSolutionIncludesDeviatoric|solveFlowAdjoints|freezeTurbulenceForValidation|gradientValidated|ransDirectionValidated|discreteCheckpointGradient|discreteProdConvergeFatal|discreteExplicitSolutionFile|discreteExplicitMatrixFile|discreteFlowAdjointTolerance|discreteFlowAdjointRestart|discreteFlowAdjointMaxIter' constant/optProperties
echo "=== backups available ==="
ls -la constant/optProperties.*
