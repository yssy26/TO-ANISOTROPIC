Info<< "cost function" << endl << flush;

if (thermalObjectiveType == "legacyMeanTemperature")
{
    MeanT =
        fvc::domainIntegrate(T*thermalObjectiveMask).value()
       /Foam::max(thermalObjectiveVolume, SMALL);
}
else if (thermalObjectiveType == "maximizeColdOutletTemperature")
{
    // Mass-flow-weighted mixed outlet temperature:
    //   Tmix = sum(phi_out * T_out) / sum(phi_out)
    scalar coldOutletMassFlow = sum
    (
        phi.boundaryField()[coldOutletPatchID]
    );
    scalar coldOutletMassFlowTemperatureIntegral = sum
    (
        phi.boundaryField()[coldOutletPatchID]
       *T.boundaryField()[coldOutletPatchID]
    );
    reduce(coldOutletMassFlow, sumOp<scalar>());
    reduce
    (
        coldOutletMassFlowTemperatureIntegral,
        sumOp<scalar>()
    );

    if (mag(coldOutletMassFlow) < SMALL)
    {
        FatalErrorInFunction
            << "Cold outlet mass flow is near zero ("
            << coldOutletMassFlow << "). Cannot compute "
            << "mass-flow-weighted outlet temperature."
            << abort(FatalError);
    }

    MeanT =
        coldOutletMassFlowTemperatureIntegral
       /coldOutletMassFlow;
}
else
{
    FatalErrorInFunction
        << "Unsupported thermalObjectiveType '" << thermalObjectiveType
        << "'." << abort(FatalError);
}

const scalar designFluidFraction =
    fvc::domainIntegrate(xh*designMask).value()
   /Foam::max(designVolume, SMALL);
const scalar designSolidFraction = 1.0 - designFluidFraction;
V = solidFractionMin - designSolidFraction;

scalar rawPowerDiss = 0.0;
for (i = 0; i < nObjPatch; ++i)
{
    rawPowerDiss -= sum
    (
        phi.boundaryField()[conPatchList[i]]
       *
        (
            p.boundaryField()[conPatchList[i]]
          + 0.5*magSqr(U.boundaryField()[conPatchList[i]])
        )
    );
}
reduce(rawPowerDiss, sumOp<scalar>());
PowerDiss = rawPowerDiss/Foam::max(PowerDiss0, SMALL);

const scalar normalizedThermalObjective =
    thermalObjectiveType == "maximizeColdOutletTemperature"
  ? -MeanT/Foam::max(thermalReferenceTemperature, SMALL)
  :  MeanT/Foam::max(thermalReferenceTemperature, SMALL);

Objective =
    objectiveMode == "outletTemperature"
  ? normalizedThermalObjective
  : (1.0 - w1)*normalizedThermalObjective + w1*PowerDiss;

scalar pInletAvg = 0.0;
scalar pOutletAvg = 0.0;
if (nObjPatch >= 2)
{
    const label inletPatchI = conPatchList[0];
    const label outletPatchI = conPatchList[1];

    scalar inletArea = sum(mesh.magSf().boundaryField()[inletPatchI]);
    scalar outletArea = sum(mesh.magSf().boundaryField()[outletPatchI]);
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

    pInletAvg = inletPIntegral/Foam::max(inletArea, SMALL);
    pOutletAvg = outletPIntegral/Foam::max(outletArea, SMALL);
}
else
{
    FatalErrorInFunction
        << "Pressure-drop constraint requires inlet and outlet patches in "
        << "constraintFunctionDict_flow." << abort(FatalError);
}

PressureDropPa = rhoFluid*(pInletAvg - pOutletAvg);
const scalar pressureDropConstraint =
    PressureDropPa/Foam::max(pressureDropMaxPa, SMALL) - 1.0;

scalar coldInletFlow = -sum(phi.boundaryField()[conPatchList[0]]);
scalar coldOutletFlow = sum(phi.boundaryField()[conPatchList[1]]);
reduce(coldInletFlow, sumOp<scalar>());
reduce(coldOutletFlow, sumOp<scalar>());
const scalar coldMassImbalance =
    mag(coldOutletFlow - coldInletFlow)
   /Foam::max
    (
        Foam::max(mag(coldInletFlow), mag(coldOutletFlow)),
        SMALL
    );

scalar QcoldW = 0.0;
scalar QhotW = 0.0;
scalar thermalBalanceError = 0.0;
scalar maxColdVelocityInBlockedRegion = 0.0;
scalar maxSolvedColdVelocity = 0.0;
scalar relativeColdLeakage = 0.0;

if (frozenHotRegionEnabled)
{
    // Use the flux of the assembled temperature matrix so this diagnostic
    // matches the selected convection, laplacian and non-orthogonal schemes.
    tmp<surfaceScalarField> tThermalEquationFlux = TEqn.flux();
    const surfaceScalarField& thermalEquationFlux =
        tThermalEquationFlux();

    scalar coldEnthalpyIn = -sum
    (
        thermalEquationFlux.boundaryField()[coldInletPatchID]
    );
    scalar coldEnthalpyOut = sum
    (
        thermalEquationFlux.boundaryField()[coldOutletPatchID]
    );
    scalar hotEnthalpyIn = -sum
    (
        thermalEquationFlux.boundaryField()[hotInletPatchID]
    );
    scalar hotEnthalpyOut = sum
    (
        thermalEquationFlux.boundaryField()[hotOutletPatchID]
    );

    reduce(coldEnthalpyIn, sumOp<scalar>());
    reduce(coldEnthalpyOut, sumOp<scalar>());
    reduce(hotEnthalpyIn, sumOp<scalar>());
    reduce(hotEnthalpyOut, sumOp<scalar>());

    QcoldW = rhoc.value()*(coldEnthalpyOut - coldEnthalpyIn);
    QhotW = rhoc.value()*(hotEnthalpyIn - hotEnthalpyOut);
    thermalBalanceError =
        mag(QhotW - QcoldW)
       /Foam::max(max(mag(QhotW), mag(QcoldW)), SMALL);

    forAll(hotRegionCells, hotI)
    {
        maxColdVelocityInBlockedRegion = Foam::max
        (
            maxColdVelocityInBlockedRegion,
            mag(U[hotRegionCells[hotI]])
        );
    }
    forAll(wallRegionCells, wallI)
    {
        maxColdVelocityInBlockedRegion = Foam::max
        (
            maxColdVelocityInBlockedRegion,
            mag(U[wallRegionCells[wallI]])
        );
    }
    reduce(maxColdVelocityInBlockedRegion, maxOp<scalar>());
    forAll(U.primitiveField(), celli)
    {
        maxSolvedColdVelocity = Foam::max
        (
            maxSolvedColdVelocity,
            mag(U[celli])
        );
    }
    reduce(maxSolvedColdVelocity, maxOp<scalar>());
    relativeColdLeakage =
        maxColdVelocityInBlockedRegion
       /Foam::max(maxSolvedColdVelocity, SMALL);

    Info<< "Frozen-hot diagnostics: Qhot=" << QhotW << " W"
        << ", Qcold=" << QcoldW << " W"
        << ", balanceError=" << thermalBalanceError
        << ", max|Ucold| in hot/wall="
        << maxColdVelocityInBlockedRegion
        << ", relativeColdLeakage=" << relativeColdLeakage << endl << flush;

    const scalar thermalBalanceTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "thermalBalanceTolerance",
            1e-2
        );
    const scalar coldLeakageRelativeTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "coldLeakageRelativeTolerance",
            1e-3
        );

    if
    (
        mmaUpdateEnabled
     &&
        (
            thermalBalanceError > thermalBalanceTolerance
         || relativeColdLeakage > coldLeakageRelativeTolerance
        )
    )
    {
        FatalErrorInFunction
            << "Forward-health gate failed before MMA: thermalBalanceError="
            << thermalBalanceError << " (limit "
            << thermalBalanceTolerance << "), relativeColdLeakage="
            << relativeColdLeakage << " (limit "
            << coldLeakageRelativeTolerance << ")." << abort(FatalError);
    }
}

Info<< "Thermal objective type=" << thermalObjectiveType
    << ", objective mode=" << objectiveMode
    << ", thermal metric=" << MeanT
    << ", normalized thermal term=" << normalizedThermalObjective
    << ", total objective=" << Objective << endl << flush;

Info<< "Forward residuals: U(init=" << forwardInitialResidualU
    << ",fin=" << forwardFinalResidualU << ")"
    << ", P(init=" << forwardInitialResidualP
    << ",fin=" << forwardFinalResidualP << ")"
    << ", T(init=" << forwardInitialResidualT
    << ",fin=" << forwardFinalResidualT << ")"
    << ", primalIters=" << forwardPrimalOuterIterations << endl << flush;

if (Pstream::master())
{
    static Foam::OFstream optimizationLogFile
    (
        runTime.path()/"optimization_log.dat"
    );

    if (!optimizationLogFile.good())
    {
        FatalErrorInFunction
            << "Unable to open optimization_log.dat."
            << abort(FatalError);
    }

    static bool optimizationLogHeaderWritten = false;
    optimizationLogFile.precision(16);
    if (!optimizationLogHeaderWritten)
    {
        optimizationLogFile
            << "# Iteration\tVolumeConstraint\tPressureDropPa"
            << "\tPressureDropConstraint\tDesignSolidFraction"
            << "\tPowerDiss\tThermalMetric\tNormalizedThermalObjective"
            << "\tObjective\tQhotW\tQcoldW\tThermalBalanceError"
            << "\tMaxColdVelocityBlocked\tRelativeColdLeakage"
            << "\tColdInletFlow"
            << "\tColdOutletFlow\tColdMassImbalance"
            << "\tFwdInitResU\tFwdFinResU"
            << "\tFwdInitResP\tFwdFinResP"
            << "\tFwdInitResT\tFwdFinResT"
            << "\tPrimalIters\tElapsedTime" << endl << flush;
        optimizationLogHeaderWritten = true;
    }

    optimizationLogFile
        << opt << "\t"
        << V << "\t"
        << PressureDropPa << "\t"
        << pressureDropConstraint << "\t"
        << designSolidFraction << "\t"
        << PowerDiss << "\t"
        << MeanT << "\t"
        << normalizedThermalObjective << "\t"
        << Objective << "\t"
        << QhotW << "\t"
        << QcoldW << "\t"
        << thermalBalanceError << "\t"
        << maxColdVelocityInBlockedRegion << "\t"
        << relativeColdLeakage << "\t"
        << coldInletFlow << "\t"
        << coldOutletFlow << "\t"
        << coldMassImbalance << "\t"
        << forwardInitialResidualU << "\t"
        << forwardFinalResidualU << "\t"
        << forwardInitialResidualP << "\t"
        << forwardFinalResidualP << "\t"
        << forwardInitialResidualT << "\t"
        << forwardFinalResidualT << "\t"
        << forwardPrimalOuterIterations << "\t"
        << runTime.elapsedClockTime() << endl << flush;
}

// Phase K: Write accepted-iteration row to comprehensive CSV log
csvAcceptedIteration = opt;
csvCandidateTrial = 0;
csvAcceptedFlag = 1;
rejectionReason = word("");
Info<< "Phase K: Writing CSV log (enableCSVLog=" << enableCSVLog << ")..." << endl << flush;
#include "writeCSVLog.H"
