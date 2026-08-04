//-----------------------------------------------------------------------------
// writeCSVLog.H
// Phase K: Comprehensive CSV logging for topology optimization.
//
// This file is self-contained: it computes all quantities from the current
// field state (U, p, T, phi, xh, x) and global diagnostic variables, then
// writes a single row to optimization_history.csv.
//
// It is included from two call sites:
//   1. costfunction.H  — accepted iteration row (csvCandidateTrial = 0)
//   2. evaluateCandidate.H — candidate trial row (csvCandidateTrial > 0)
//
// Before including, the caller must set:
//   csvCandidateTrial  — 0 for accepted, 1+ for candidate attempt
//   rejectionReason    — word describing rejection (empty if accepted)
//   csvAcceptedFlag    — 1 if this row is accepted, 0 if rejected candidate
//
// CSV write requirements (Plan §13.4):
//   - First run writes header;
//   - Restart appends without duplicate header;
//   - Candidate trials and accepted iterations have different identifiers;
//   - Existing human-readable logs (optimization_log.dat, mma_step_log.dat)
//     are preserved.
//
// IMPORTANT: All reduce() calls MUST be outside Pstream::master() blocks
// because reduce() is a collective MPI operation requiring ALL processes.
//-----------------------------------------------------------------------------

if (enableCSVLog)
{
    // --- Compute quantities from current field state -----------------------
    // All processes participate in these computations (including reduce calls).

    // Cold outlet mixed (mass-flow-weighted) temperature
    scalar csvMeanT = 0.0;
    if (thermalObjectiveType == "maximizeColdOutletTemperature")
    {
        scalar mFlow = sum(phi.boundaryField()[coldOutletPatchID]);
        scalar mT = sum
        (
            phi.boundaryField()[coldOutletPatchID]
           *T.boundaryField()[coldOutletPatchID]
        );
        reduce(mFlow, sumOp<scalar>());
        reduce(mT, sumOp<scalar>());
        csvMeanT = mag(mFlow) > SMALL ? mT/mFlow : 0.0;
    }
    else
    {
        csvMeanT =
            fvc::domainIntegrate(T*thermalObjectiveMask).value()
           /Foam::max(thermalObjectiveVolume, SMALL);
    }

    // Volume constraint
    const scalar csvDesignFluidFraction =
        fvc::domainIntegrate(xh*designMask).value()
       /Foam::max(designVolume, SMALL);
    const scalar csvDesignSolidFraction = 1.0 - csvDesignFluidFraction;
    const scalar csvVolumeConstraint =
        solidFractionMin - csvDesignSolidFraction;

    // Pressure drop
    scalar csvPressureDropPa = 0.0;
    if (nObjPatch >= 2)
    {
        const label inletPatchI = conPatchList[0];
        const label outletPatchI = conPatchList[1];
        scalar inletArea = sum(mesh.magSf().boundaryField()[inletPatchI]);
        scalar outletArea = sum(mesh.magSf().boundaryField()[outletPatchI]);
        scalar inletPInt = sum
        (
            p.boundaryField()[inletPatchI]
           *mesh.magSf().boundaryField()[inletPatchI]
        );
        scalar outletPInt = sum
        (
            p.boundaryField()[outletPatchI]
           *mesh.magSf().boundaryField()[outletPatchI]
        );
        reduce(inletArea, sumOp<scalar>());
        reduce(outletArea, sumOp<scalar>());
        reduce(inletPInt, sumOp<scalar>());
        reduce(outletPInt, sumOp<scalar>());
        csvPressureDropPa = rhoFluid
           *(
                inletPInt/Foam::max(inletArea, SMALL)
              - outletPInt/Foam::max(outletArea, SMALL)
            );
    }
    const scalar csvPressureDropConstraint =
        csvPressureDropPa/Foam::max(pressureDropMaxPa, SMALL) - 1.0;

    // Constraint violation (sum of squared positive violations)
    const scalar csvConstraintViolation =
        sqr(Foam::max(0.0, csvVolumeConstraint))
      + sqr(Foam::max(0.0, csvPressureDropConstraint));

    // Mass flows
    scalar csvColdInletFlow = -sum(phi.boundaryField()[conPatchList[0]]);
    scalar csvColdOutletFlow = sum(phi.boundaryField()[conPatchList[1]]);
    reduce(csvColdInletFlow, sumOp<scalar>());
    reduce(csvColdOutletFlow, sumOp<scalar>());
    const scalar csvColdMassImbalance =
        mag(csvColdOutletFlow - csvColdInletFlow)
       /Foam::max
        (
            Foam::max(mag(csvColdInletFlow), mag(csvColdOutletFlow)),
            SMALL
        );

    // Thermal balance (NaN when frozen-hot region is not enabled, per Plan §13.2)
    scalar csvQhotW = std::numeric_limits<double>::quiet_NaN();
    scalar csvQcoldW = std::numeric_limits<double>::quiet_NaN();
    scalar csvThermalBalanceError = std::numeric_limits<double>::quiet_NaN();
    if (frozenHotRegionEnabled)
    {
        tmp<surfaceScalarField> tCsvThermalFlux = TEqn.flux();
        const surfaceScalarField& csvThermalFlux = tCsvThermalFlux();

        scalar coldEnthIn = -sum
        (
            csvThermalFlux.boundaryField()[coldInletPatchID]
        );
        scalar coldEnthOut = sum
        (
            csvThermalFlux.boundaryField()[coldOutletPatchID]
        );
        scalar hotEnthIn = -sum
        (
            csvThermalFlux.boundaryField()[hotInletPatchID]
        );
        scalar hotEnthOut = sum
        (
            csvThermalFlux.boundaryField()[hotOutletPatchID]
        );
        reduce(coldEnthIn, sumOp<scalar>());
        reduce(coldEnthOut, sumOp<scalar>());
        reduce(hotEnthIn, sumOp<scalar>());
        reduce(hotEnthOut, sumOp<scalar>());

        csvQcoldW = rhoc.value()*(coldEnthOut - coldEnthIn);
        csvQhotW = rhoc.value()*(hotEnthIn - hotEnthOut);
        csvThermalBalanceError =
            mag(csvQhotW - csvQcoldW)
           /Foam::max(max(mag(csvQhotW), mag(csvQcoldW)), SMALL);
    }

    // Temperature statistics
    scalar csvMaxT = gMax(T.primitiveField());
    scalar csvAvgT =
        fvc::domainIntegrate(T).value()
       /Foam::max(gSum(mesh.V()), SMALL);

    // Gray volume fraction: design-region volume where 0.05 < xh < 0.95
    scalar grayVolume = 0.0;
    scalar totalDesignVolume = 0.0;
    forAll(xh.primitiveField(), celli)
    {
        if (designMask[celli] > 0.5)
        {
            const scalar v = mesh.V()[celli];
            totalDesignVolume += v;
            if (xh[celli] > 0.05 && xh[celli] < 0.95)
            {
                grayVolume += v;
            }
        }
    }
    reduce(grayVolume, sumOp<scalar>());
    reduce(totalDesignVolume, sumOp<scalar>());
    const scalar csvGrayVolumeFraction =
        grayVolume/Foam::max(totalDesignVolume, SMALL);

    // Design variable statistics (active design cells only)
    scalar csvMinDesign = GREAT;
    scalar csvMaxDesign = -GREAT;
    forAll(activeDesignCells, activeI)
    {
        const scalar xv = x[activeDesignCells[activeI]];
        csvMinDesign = Foam::min(csvMinDesign, xv);
        csvMaxDesign = Foam::max(csvMaxDesign, xv);
    }
    reduce(csvMinDesign, minOp<scalar>());
    reduce(csvMaxDesign, maxOp<scalar>());

    // Normalized objective
    const scalar csvNormalizedThermal =
        thermalObjectiveType == "maximizeColdOutletTemperature"
      ? -csvMeanT/Foam::max(thermalReferenceTemperature, SMALL)
      :  csvMeanT/Foam::max(thermalReferenceTemperature, SMALL);

    scalar csvRawPowerDiss = 0.0;
    for (label patchI = 0; patchI < nObjPatch; ++patchI)
    {
        csvRawPowerDiss -= sum
        (
            phi.boundaryField()[conPatchList[patchI]]
           *
           (
                p.boundaryField()[conPatchList[patchI]]
              + 0.5*magSqr(U.boundaryField()[conPatchList[patchI]])
            )
        );
    }
    reduce(csvRawPowerDiss, sumOp<scalar>());
    const scalar csvPowerDiss =
        csvRawPowerDiss/Foam::max(PowerDiss0, SMALL);

    const scalar csvObjective =
        objectiveMode == "outletTemperature"
      ? csvNormalizedThermal
      : (1.0 - w1)*csvNormalizedThermal + w1*csvPowerDiss;

    // SST relative changes: NaN for non-RANS flow models (Plan §13.2)
    const scalar csvRelChangeK =
        ransFlowModel
      ? lastRelativeChangeK
      : std::numeric_limits<double>::quiet_NaN();
    const scalar csvRelChangeOmega =
        ransFlowModel
      ? lastRelativeChangeOmega
      : std::numeric_limits<double>::quiet_NaN();

    // Adjoint residuals: NaN when adjoints are not solved (Plan §13.2)
    const scalar csvAdjResU =
        (!freezeColdFlowForValidation && solveFlowAdjoints)
      ? adjointResidualU
      : std::numeric_limits<double>::quiet_NaN();
    const scalar csvAdjResP =
        (!freezeColdFlowForValidation && solveFlowAdjoints)
      ? adjointResidualP
      : std::numeric_limits<double>::quiet_NaN();
    const scalar csvAdjResT = adjointResidualT;

    // --- Write CSV row (master only) --------------------------------------
    // File I/O is done by master only. All reduce() calls are already done.

    if (Pstream::master())
    {
        // Check if file exists (for header writing / restart detection)
        const bool csvFileExists = Foam::isFile(csvLogPath);

        // Use std::ofstream for portable append mode
        std::ofstream csvFile(csvLogPath.c_str(), std::ios::app);
        if (!csvFile.is_open())
        {
            FatalErrorInFunction
                << "Cannot open " << csvLogPath << " for CSV writing."
                << abort(FatalError);
        }
        csvFile.precision(16);

        if (!csvFileExists)
        {
            // First run: write header
            csvFile
                << "Iteration"
                << ",AcceptedIteration"
                << ",CandidateTrial"
                << ",Objective"
                << ",ColdOutletMixedTemperature"
                << ",VolumeConstraint"
                << ",PressureDropConstraint"
                << ",ConstraintViolation"
                << ",DesignSolidFraction"
                << ",ColdPressureDropPa"
                << ",ColdInletMassFlow"
                << ",ColdOutletMassFlow"
                << ",ColdMassImbalance"
                << ",QhotW"
                << ",QcoldW"
                << ",ThermalBalanceError"
                << ",MaxTemperature"
                << ",AverageTemperature"
                << ",ForwardInitialResidualU"
                << ",ForwardFinalResidualU"
                << ",ForwardInitialResidualP"
                << ",ForwardFinalResidualP"
                << ",ForwardInitialResidualT"
                << ",ForwardFinalResidualT"
                << ",RelativeChangeK"
                << ",RelativeChangeOmega"
                << ",AdjointResidualU"
                << ",AdjointResidualP"
                << ",AdjointResidualT"
                << ",MoveLimit"
                << ",CandidateAccepted"
                << ",RejectionReason"
                << ",GrayVolumeFraction"
                << ",MinDesignVariable"
                << ",MaxDesignVariable"
                << ",PrimalOuterIterations"
                << ",ElapsedTime"
                << "\n" << flush;
        }

        csvFile
            << opt << ","
            << csvAcceptedIteration << ","
            << csvCandidateTrial << ","
            << csvObjective << ","
            << csvMeanT << ","
            << csvVolumeConstraint << ","
            << csvPressureDropConstraint << ","
            << csvConstraintViolation << ","
            << csvDesignSolidFraction << ","
            << csvPressureDropPa << ","
            << csvColdInletFlow << ","
            << csvColdOutletFlow << ","
            << csvColdMassImbalance << ","
            << csvQhotW << ","
            << csvQcoldW << ","
            << csvThermalBalanceError << ","
            << csvMaxT << ","
            << csvAvgT << ","
            << forwardInitialResidualU << ","
            << forwardFinalResidualU << ","
            << forwardInitialResidualP << ","
            << forwardFinalResidualP << ","
            << forwardInitialResidualT << ","
            << forwardFinalResidualT << ","
            << csvRelChangeK << ","
            << csvRelChangeOmega << ","
            << csvAdjResU << ","
            << csvAdjResP << ","
            << csvAdjResT << ","
            << currentMoveLimit << ","
            << csvAcceptedFlag << ","
            << rejectionReason << ","
            << csvGrayVolumeFraction << ","
            << csvMinDesign << ","
            << csvMaxDesign << ","
            << forwardPrimalOuterIterations << ","
            << runTime.elapsedClockTime()
            << "\n" << flush;

        csvFile.close();

        Info<< "CSV log written: " << csvLogPath
            << " (trial=" << csvCandidateTrial
            << ", accepted=" << csvAcceptedFlag << ")" << endl << flush;
    }
}
