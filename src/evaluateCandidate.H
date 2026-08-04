//-----------------------------------------------------------------------------
// evaluateCandidate.H
// Phase I: Full-SST candidate evaluation and acceptance check.
//
// This file is included from sensitivity.H when enableSSTAcceptanceCheck
// is true.  At this point the MMA candidate has already been written to the
// design fields (x, xp, xh) and the filter/projection chain has been
// recomputed via filter_x.H.
//
// Steps:
//   1. Recompute material properties (alpha, DT, Q) from candidate xh
//   2. Update frozen-turbulence derived fields for warm start
//   3. Solve full SST primal (NS.H + HeatTransfer.H)
//   4. Compute candidate objective and constraints
//   5. Check convergence and acceptance criteria
//   6. Set candidateAccepted = true/false
//
// The variable `candidateAccepted` (bool) is expected to be in scope and
// is set by this file.
//-----------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// 1. Recompute material properties from candidate xh
// ---------------------------------------------------------------------------
#include "updateMaterialProperties.H"

// ---------------------------------------------------------------------------
// 2. Update frozen-turbulence derived fields (warm start with candidate xh)
// ---------------------------------------------------------------------------
#include "updateFrozenTurbulenceFields.H"

// ---------------------------------------------------------------------------
// 3. Full SST primal solve
// ---------------------------------------------------------------------------
// Temporarily disable turbulence freezing so NS.H calls turbulence->correct().
const Switch freezeTurbBackup = freezeTurbulenceForValidation;
freezeTurbulenceForValidation = false;

Info<< "\n=== Transaction: evaluating candidate with full SST ===\n"
    << endl;

#include "NS.H"
#include "updateThermalFlux.H"
#include "HeatTransfer.H"

freezeTurbulenceForValidation = freezeTurbBackup;

// ---------------------------------------------------------------------------
// 4–6. Compute candidate metrics, check health, and decide acceptance.
//      Uses do-while(false) + break to allow early exit without goto.
// ---------------------------------------------------------------------------
scalar candidateMeanT = 0.0;
scalar candidateVolumeConstraint = 0.0;
scalar candidatePressureDropPa = 0.0;
scalar candidatePressureDropConstraint = 0.0;
scalar candidateObjective = 0.0;
scalar candidatePowerDiss = 0.0;
scalar candidateMassImbalance = 0.0;
scalar candidateThermalBalanceError = 0.0;
bool candidateSSTConverged = !ransFlowModel || turbulenceConverged;
word candidateRejectionReason("");

do
{
    // --- Thermal metric (MeanT) ---
    if (thermalObjectiveType == "legacyMeanTemperature")
    {
        candidateMeanT =
            fvc::domainIntegrate(T*thermalObjectiveMask).value()
           /Foam::max(thermalObjectiveVolume, SMALL);
    }
    else  // maximizeColdOutletTemperature
    {
        scalar mFlow = sum(phi.boundaryField()[coldOutletPatchID]);
        scalar mT = sum
        (
            phi.boundaryField()[coldOutletPatchID]
           *T.boundaryField()[coldOutletPatchID]
        );
        reduce(mFlow, sumOp<scalar>());
        reduce(mT, sumOp<scalar>());

        if (mag(mFlow) < SMALL)
        {
            Info<< "Transaction: REJECT - cold outlet mass flow near zero."
                << endl;
            candidateAccepted = false;
            candidateRejectionReason = "ColdOutletMassFlowNearZero";
            break;
        }
        candidateMeanT = mT/mFlow;
    }

    // --- Volume constraint ---
    const scalar candidateDesignFluidFraction =
        fvc::domainIntegrate(xh*designMask).value()
       /Foam::max(designVolume, SMALL);
    const scalar candidateDesignSolidFraction =
        1.0 - candidateDesignFluidFraction;
    candidateVolumeConstraint =
        solidFractionMin - candidateDesignSolidFraction;

    // --- Pressure drop ---
    if (nObjPatch >= 2)
    {
        const label inletPatchI = conPatchList[0];
        const label outletPatchI = conPatchList[1];

        scalar inletArea =
            sum(mesh.magSf().boundaryField()[inletPatchI]);
        scalar outletArea =
            sum(mesh.magSf().boundaryField()[outletPatchI]);
        scalar inletPIntegral = sum
        (
            p.boundaryField()[inletPatchI]
           *mesh.magSf().boundaryField()[inletPatchI]
        );
        scalar outletPIntegral = sum
        (
            p.boundaryField()[outletPatchI]
           *mesh.magSf().boundaryField()[outletPatchI]
        );

        reduce(inletArea, sumOp<scalar>());
        reduce(outletArea, sumOp<scalar>());
        reduce(inletPIntegral, sumOp<scalar>());
        reduce(outletPIntegral, sumOp<scalar>());

        candidatePressureDropPa =
            rhoFluid
           *(
                inletPIntegral/Foam::max(inletArea, SMALL)
              - outletPIntegral/Foam::max(outletArea, SMALL)
            );
    }
    candidatePressureDropConstraint =
        candidatePressureDropPa/Foam::max(pressureDropMaxPa, SMALL) - 1.0;

    // --- Normalized objective ---
    const scalar candidateNormalizedThermal =
        thermalObjectiveType == "maximizeColdOutletTemperature"
      ? -candidateMeanT/Foam::max(thermalReferenceTemperature, SMALL)
      :  candidateMeanT/Foam::max(thermalReferenceTemperature, SMALL);

    {
        scalar rawPowerDiss = 0.0;
        for (label patchI = 0; patchI < nObjPatch; ++patchI)
        {
            rawPowerDiss -= sum
            (
                phi.boundaryField()[conPatchList[patchI]]
               *
               (
                    p.boundaryField()[conPatchList[patchI]]
                  + 0.5*magSqr(U.boundaryField()[conPatchList[patchI]])
                )
            );
        }
        reduce(rawPowerDiss, sumOp<scalar>());
        candidatePowerDiss = rawPowerDiss/Foam::max(PowerDiss0, SMALL);
    }

    candidateObjective =
        objectiveMode == "outletTemperature"
      ? candidateNormalizedThermal
      : (1.0 - w1)*candidateNormalizedThermal + w1*candidatePowerDiss;

    // --- Transaction test: artificially worsen candidate objective ---
    // When transactionTestMode is active and this is the first attempt,
    // inflate the candidate objective to guarantee rejection.  This
    // simulates a candidate that would be rejected by the acceptance
    // criterion, allowing full rollback verification.
    if (transactionTestMode && candidateAttempt == 1)
    {
        const scalar inflator = mag(objectiveTxBackup)*0.1 + 0.001;
        candidateObjective += inflator;
        Info<< "  TX_TEST: Artificially inflated candidate objective by "
            << inflator << " to force rejection" << endl;

        // [v2 Gate 5] The feasibility-based acceptance logic may ACCEPT a
        // more-feasible candidate unconditionally (no objective check),
        // which at pressure-infeasible design points defeats the inflation.
        // The transaction test must deterministically exercise the rollback
        // path, so force the rejection of the first attempt.
        candidateAccepted = false;
        candidateRejectionReason = "TransactionTestForcedRejection";
        break;
    }

    // --- Mass balance ---
    {
        scalar coldIn = -sum(phi.boundaryField()[conPatchList[0]]);
        scalar coldOut = sum(phi.boundaryField()[conPatchList[1]]);
        reduce(coldIn, sumOp<scalar>());
        reduce(coldOut, sumOp<scalar>());
        candidateMassImbalance =
            mag(coldOut - coldIn)
           /Foam::max
            (
                Foam::max(mag(coldIn), mag(coldOut)),
                SMALL
            );
    }

    // --- Thermal balance ---
    if (frozenHotRegionEnabled)
    {
        tmp<surfaceScalarField> tCandFlux = TEqn.flux();
        scalar coldEnthalpyIn = -sum
        (
            tCandFlux().boundaryField()[coldInletPatchID]
        );
        scalar coldEnthalpyOut = sum
        (
            tCandFlux().boundaryField()[coldOutletPatchID]
        );
        scalar hotEnthalpyIn = -sum
        (
            tCandFlux().boundaryField()[hotInletPatchID]
        );
        scalar hotEnthalpyOut = sum
        (
            tCandFlux().boundaryField()[hotOutletPatchID]
        );

        reduce(coldEnthalpyIn, sumOp<scalar>());
        reduce(coldEnthalpyOut, sumOp<scalar>());
        reduce(hotEnthalpyIn, sumOp<scalar>());
        reduce(hotEnthalpyOut, sumOp<scalar>());

        scalar QcoldCand =
            rhoc.value()*(coldEnthalpyOut - coldEnthalpyIn);
        scalar QhotCand =
            rhoc.value()*(hotEnthalpyIn - hotEnthalpyOut);
        candidateThermalBalanceError =
            mag(QhotCand - QcoldCand)
           /Foam::max(max(mag(QhotCand), mag(QcoldCand)), SMALL);
    }

    // --- Hard failure checks ---
    const scalar massBalanceTolerance =
        optProperties.lookupOrDefault<scalar>("massBalanceTolerance", 5e-3);
    const scalar thermalBalanceToleranceTx =
        optProperties.lookupOrDefault<scalar>("thermalBalanceTolerance", 1e-2);

    if (!candidateSSTConverged)
    {
        Info<< "Transaction: REJECT - SST did not converge." << endl;
        candidateAccepted = false;
        candidateRejectionReason = "SSTNotConverged";
        break;
    }

    if (candidateMassImbalance > massBalanceTolerance)
    {
        Info<< "Transaction: REJECT - mass imbalance "
            << candidateMassImbalance << " > " << massBalanceTolerance
            << endl;
        candidateAccepted = false;
        candidateRejectionReason = "MassImbalance";
        break;
    }

    if
    (
        frozenHotRegionEnabled
     && candidateThermalBalanceError > thermalBalanceToleranceTx
    )
    {
        Info<< "Transaction: REJECT - thermal balance error "
            << candidateThermalBalanceError << " > "
            << thermalBalanceToleranceTx << endl;
        candidateAccepted = false;
        candidateRejectionReason = "ThermalBalanceError";
        break;
    }

    // --- Acceptance criteria ---
    const scalar acceptedTheta =
        sqr(Foam::max(0.0, volumeConstraintTxBackup))
      + sqr(Foam::max(0.0, pressureDropTxBackup));

    const scalar candidateTheta =
        sqr(Foam::max(0.0, candidateVolumeConstraint))
      + sqr(Foam::max(0.0, candidatePressureDropConstraint));

    const scalar objectiveTolerance =
        mag(objectiveTxBackup)*relativeObjectiveTolerance;

    Info<< "Transaction comparison:" << endl;
    Info<< "  Accepted: J=" << objectiveTxBackup
        << ", V=" << volumeConstraintTxBackup
        << ", dP=" << pressureDropTxBackup
        << ", theta=" << acceptedTheta
        << ", MeanT=" << meanTTxBackup << endl;
    Info<< "  Candidate: J=" << candidateObjective
        << ", V=" << candidateVolumeConstraint
        << ", dP=" << candidatePressureDropConstraint
        << ", theta=" << candidateTheta
        << ", MeanT=" << candidateMeanT << endl;
    Info<< "  Health: SST=" << (candidateSSTConverged ? "OK" : "FAIL")
        << ", massImb=" << candidateMassImbalance
        << ", thermalBalance=" << candidateThermalBalanceError << endl;

    if (acceptedTheta > absoluteConstraintTolerance)
    {
        // Accepted design is infeasible
        if (candidateTheta < acceptedTheta - absoluteConstraintTolerance)
        {
            candidateAccepted = true;
            Info<< "  Decision: ACCEPT (candidate more feasible)" << endl;
        }
        else if
        (
            candidateTheta > acceptedTheta + absoluteConstraintTolerance
         && candidateObjective > objectiveTxBackup
        )
        {
            candidateAccepted = false;
            candidateRejectionReason = "ConstraintsAndObjectiveWorse";
            Info<< "  Decision: REJECT (constraints and objective worse)"
                << endl;
        }
        else
        {
            candidateAccepted =
                (candidateObjective <= objectiveTxBackup + objectiveTolerance);
            if (!candidateAccepted)
            {
                candidateRejectionReason = "InfeasibleSimilarObjectiveWorse";
            }
            Info<< "  Decision: "
                << (candidateAccepted ? "ACCEPT" : "REJECT")
                << " (similar feasibility, objective comparison)" << endl;
        }
    }
    else
    {
        // Accepted design is feasible
        if (candidateTheta > absoluteConstraintTolerance)
        {
            candidateAccepted =
                (candidateObjective
                 < objectiveTxBackup - 2.0*objectiveTolerance);
            if (!candidateAccepted)
            {
                candidateRejectionReason = "CandidateInfeasibleObjectiveWorse";
            }
            Info<< "  Decision: "
                << (candidateAccepted ? "ACCEPT" : "REJECT")
                << " (candidate infeasible, accepted feasible)" << endl;
        }
        else
        {
            candidateAccepted =
                (candidateObjective <= objectiveTxBackup + objectiveTolerance);
            if (!candidateAccepted)
            {
                candidateRejectionReason = "BothFeasibleObjectiveWorse";
            }
            Info<< "  Decision: "
                << (candidateAccepted ? "ACCEPT" : "REJECT")
                << " (both feasible, objective comparison)" << endl;
        }
    }
} while (false);

// Update global objective/constraint values for accepted candidate
if (candidateAccepted)
{
    Objective = candidateObjective;
    V = candidateVolumeConstraint;
    PressureDropPa = candidatePressureDropPa;
    MeanT = candidateMeanT;
    PowerDiss = candidatePowerDiss;
}

Info<< "=== Transaction: "
    << (candidateAccepted ? "CANDIDATE ACCEPTED" : "CANDIDATE REJECTED")
    << " ===\n" << endl;

// Phase K: Write candidate-trial row to comprehensive CSV log
csvCandidateTrial = candidateAttempt;
csvAcceptedFlag = candidateAccepted ? 1 : 0;
rejectionReason = candidateRejectionReason;
#include "writeCSVLog.H"
