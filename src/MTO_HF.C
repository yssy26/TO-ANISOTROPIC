#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "MMA/MMA.h"
#include "diff.c"
#include "cellSet.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvOptions.H"
    #include "createFields.H"
    #include "readTransportProperties.H"
    #include "readThermalProperties.H"
    #include "createFrozenHotRegionFields.H"
    #include "createFrozenTurbulenceFields.H"
    #include "validateDiscreteObjectiveDerivatives.H"
    #include "validateRegionDefinitions.H"
    #include "validateFrozenHotCase.H"
    #include "initContinuityErrs.H"
    #include "opt_initialization.H"

    while (simple.loop(runTime))
    {
        #include "update.H"
        #include "updateFrozenTurbulenceFields.H"
        Info << "MTO_HF: Finished update.H" << endl << flush;

        if (!freezeColdFlowForValidation)
        {
            Info << "MTO_HF: Starting NS.H" << endl << flush;
            #include "NS.H"
            Info << "MTO_HF: Finished NS.H" << endl << flush;
        }
        else
        {
            Info<< "MTO_HF: cold U/p/phi frozen for thermal validation"
                << endl << flush;
        }

        #include "updateThermalFlux.H"

        Info << "MTO_HF: Starting HeatTransfer.H" << endl << flush;
        #include "HeatTransfer.H"
        Info << "MTO_HF: Finished HeatTransfer.H" << endl << flush;

        // Refresh objective derivative based on current mass flux before
        // solving the adjoint (mass-flow-weighted temperature depends on phi).
        #include "computeObjective.H"

        Info << "MTO_HF: Starting AdjHeatTransfer.H" << endl << flush;
        #include "AdjHeatTransfer.H"
        Info << "MTO_HF: Finished AdjHeatTransfer.H" << endl << flush;

        if (!freezeColdFlowForValidation && solveFlowAdjoints)
        {
            Info << "MTO_HF: Starting AdjNS_HT.H" << endl << flush;
            #include "AdjNS_HT.H"
            Info << "MTO_HF: Finished AdjNS_HT.H" << endl << flush;
        }

        if
        (
            !freezeColdFlowForValidation
         && solveFlowAdjoints
         && objectiveMode == "legacyComposite"
        )
        {
            Info << "MTO_HF: Starting AdjNS_FF.H" << endl << flush;
            #include "AdjNS_FF.H"
            Info << "MTO_HF: Finished AdjNS_FF.H" << endl << flush;
        }

        if (!freezeColdFlowForValidation && solveFlowAdjoints)
        {
            Info << "MTO_HF: Starting AdjNS_PD.H" << endl << flush;
            #include "AdjNS_PD.H"
            Info << "MTO_HF: Finished AdjNS_PD.H" << endl << flush;
        }

        Info << "MTO_HF: Starting costfunction.H" << endl << flush;
        #include "costfunction.H"
        #include "validateOptimizationState.H"
        Info << "MTO_HF: Finished costfunction.H" << endl << flush;

        #include "writeOptimizationState.H"

        Info << "MTO_HF: Starting sensitivity.H" << endl << flush;
        #include "sensitivity.H"
        Info << "MTO_HF: Finished sensitivity.H" << endl << flush;

        // Phase J: Save optimizer state for restart.
        // Always save, even on finalSolvedIteration, so that the state is
        // available for restart. The saved state reflects the current
        // solved design (pre-MMA-update on finalSolvedIteration, post-MMA
        // on normal iterations).
        #include "saveOptimizerState.H"

        #include "validateCommon.H"
        #include "validateGradientChain.H"
        #include "validateFrozenGradient.H"
        #include "validateSSTDirection.H"
        // Strict replacement for the legacy sign-only Gate 6. This version
        // removes the duplicate pressure normalization, checks gradient
        // magnitude, enforces a minimum number of valid directions and
        // verifies full-SST repeatability before unlocking gradientValidated.
        #include "validateGate6SSTStrict.H"
    }

    #include "finalize.H"
    return 0;
}
