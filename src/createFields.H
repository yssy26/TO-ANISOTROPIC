Info<< "Reading optimization parameters\n" << endl;
IOdictionary optProperties
(
    IOobject
    (
        "optProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar raa0(readScalar(optProperties.lookup("raa0")));
scalar mma_init(readScalar(optProperties.lookup("mma_init")));
scalar mma_dec(readScalar(optProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(optProperties.lookup("mma_inc")));
scalar movlim(readScalar(optProperties.lookup("movlim")));
scalar voluse(readScalar(optProperties.lookup("voluse")));
scalar solidFractionMin
(
    optProperties.lookupOrDefault<scalar>("solidFractionMin", 1.0 - voluse)
);
scalar PowerDiss0(readScalar(optProperties.lookup("PowerDiss0")));
scalar MeanT0(readScalar(optProperties.lookup("MeanT0")));
scalar rhoFluid(optProperties.lookupOrDefault<scalar>("rhoFluid", 1.0));
scalar pressureDropMaxPa
(
    optProperties.lookupOrDefault<scalar>("pressureDropMaxPa", 40.0)
);
scalar filterR(readScalar(optProperties.lookup("filterR")));
scalar w1(readScalar(optProperties.lookup("w1")));
word objectiveMode
(
    optProperties.lookupOrDefault<word>
    (
        "objectiveMode",
        "outletTemperature"
    )
);
word adjointMode
(
    optProperties.lookupOrDefault<word>
    (
        "adjointMode",
        "discrete"
    )
);
word flowModel
(
    optProperties.lookupOrDefault<word>("flowModel", "laminar")
);
Switch frozenTurbulenceAdjoint
(
    optProperties.lookupOrDefault<Switch>
    (
        "frozenTurbulenceAdjoint",
        false
    )
);
Switch ransDirectionValidated
(
    optProperties.lookupOrDefault<Switch>("ransDirectionValidated", false)
);Switch freezeTurbulenceForValidation
(
    optProperties.lookupOrDefault<Switch>
    (
        "freezeTurbulenceForValidation",
        false
    )
);
scalar turbulentPrandtl
(
    optProperties.lookupOrDefault<scalar>("turbulentPrandtl", 0.85)
);
scalar turbulenceDampingExponent
(
    optProperties.lookupOrDefault<scalar>
    (
        "turbulenceDampingExponent",
        3.0
    )
);
label turbulenceMinCorrectors
(
    optProperties.lookupOrDefault<label>("turbulenceMinCorrectors", 40)
);
label turbulenceMaxCorrectors
(
    optProperties.lookupOrDefault<label>("turbulenceMaxCorrectors", 200)
);
scalar turbulenceConvergenceTolerance
(
    optProperties.lookupOrDefault<scalar>
    (
        "turbulenceConvergenceTolerance",
        1e-6
    )
);

if
(
    objectiveMode != "outletTemperature"
 && objectiveMode != "legacyComposite"
)
{
    FatalErrorInFunction
        << "objectiveMode must be 'outletTemperature' or "
        << "'legacyComposite', received '" << objectiveMode << "'."
        << abort(FatalError);
}

if (adjointMode != "discrete" && adjointMode != "continuous")
{
    FatalErrorInFunction
        << "adjointMode must be 'discrete' or 'continuous', received '"
        << adjointMode << "'." << abort(FatalError);
}
if (flowModel != "laminar" && flowModel != "incompressibleRANSFrozen")
{
    FatalErrorInFunction
        << "flowModel must be 'laminar' or 'incompressibleRANSFrozen', "
        << "received '" << flowModel << "'." << abort(FatalError);
}
if (flowModel == "incompressibleRANSFrozen" && !frozenTurbulenceAdjoint)
{
    FatalErrorInFunction
        << "flowModel=incompressibleRANSFrozen requires "
        << "frozenTurbulenceAdjoint=true." << abort(FatalError);
}if
(
    freezeTurbulenceForValidation
 && flowModel != "incompressibleRANSFrozen"
)
{
    FatalErrorInFunction
        << "freezeTurbulenceForValidation is available only for "
        << "flowModel=incompressibleRANSFrozen." << abort(FatalError);
}
if
(
    turbulentPrandtl <= SMALL
 || turbulenceDampingExponent <= SMALL
 || turbulenceMinCorrectors < 1
 || turbulenceMaxCorrectors < turbulenceMinCorrectors
 || turbulenceConvergenceTolerance <= SMALL
)
{
    FatalErrorInFunction
        << "Invalid frozen-RANS controls: require turbulentPrandtl>0, "
        << "turbulenceDampingExponent>0, 1<=turbulenceMinCorrectors<="
        << "turbulenceMaxCorrectors and "
        << "turbulenceConvergenceTolerance>0." << abort(FatalError);
}
scalar validationProjectionBeta
(
    optProperties.lookupOrDefault<scalar>("validationProjectionBeta", 0.0)
);
scalar projectionEtaTolerance
(
    optProperties.lookupOrDefault<scalar>("projectionEtaTolerance", 1e-10)
);
Switch recomputeDesignMappingWhenMmaDisabled
(
    optProperties.lookupOrDefault<Switch>
    (
        "recomputeDesignMappingWhenMmaDisabled",
        true
    )
);
Switch freezeColdFlowForValidation
(
    optProperties.lookupOrDefault<Switch>
    (
        "freezeColdFlowForValidation",
        false
    )
);
Switch solveFlowAdjoints
(
    optProperties.lookupOrDefault<Switch>("solveFlowAdjoints", true)
);
Switch directionOnlyOptimizationApproved
(
    optProperties.lookupOrDefault<Switch>
    (
        "directionOnlyOptimizationApproved",
        false
    )
);
Switch gradientValidated
(
    optProperties.lookupOrDefault<Switch>("gradientValidated", false)
);
Switch directionOnlyFeasibilityFirst
(
    optProperties.lookupOrDefault<Switch>
    (
        "directionOnlyFeasibilityFirst",
        false
    )
);
Switch continuationEnabled
(
    optProperties.lookupOrDefault<Switch>("continuationEnabled", false)
);
scalar objectiveGradientScale
(
    optProperties.lookupOrDefault<scalar>("objectiveGradientScale", 1.0)
);
scalar pressureGradientScale
(
    optProperties.lookupOrDefault<scalar>("pressureGradientScale", 1.0)
);

if (objectiveGradientScale <= 0.0 || pressureGradientScale <= 0.0)
{
    FatalErrorInFunction
        << "objectiveGradientScale and pressureGradientScale must be positive."
        << abort(FatalError);
}

if
(
    validationProjectionBeta < 0.0
 || projectionEtaTolerance <= SMALL
 || projectionEtaTolerance >= 0.1
)
{
    FatalErrorInFunction
        << "Require validationProjectionBeta >= 0 and "
        << "0 < projectionEtaTolerance < 0.1."
        << abort(FatalError);
}

// ---------------------------------------------------------------------------
// Phase H: Transaction mechanism configuration parameters
// ---------------------------------------------------------------------------
Switch enableSSTAcceptanceCheck
(
    optProperties.lookupOrDefault<Switch>("enableSSTAcceptanceCheck", false)
);
scalar initialMoveLimit
(
    optProperties.lookupOrDefault<scalar>("initialMoveLimit", movlim)
);
scalar minimumMoveLimit
(
    optProperties.lookupOrDefault<scalar>("minimumMoveLimit", 0.01)
);
scalar maximumMoveLimit
(
    optProperties.lookupOrDefault<scalar>("maximumMoveLimit", 0.5)
);
scalar moveLimitReduction
(
    optProperties.lookupOrDefault<scalar>("moveLimitReduction", 0.5)
);
scalar moveLimitGrowth
(
    optProperties.lookupOrDefault<scalar>("moveLimitGrowth", 1.2)
);
label maxCandidateRetries
(
    optProperties.lookupOrDefault<label>("maxCandidateRetries", 3)
);
scalar relativeObjectiveTolerance
(
    optProperties.lookupOrDefault<scalar>("relativeObjectiveTolerance", 0.01)
);
scalar absoluteConstraintTolerance
(
    optProperties.lookupOrDefault<scalar>("absoluteConstraintTolerance", 0.01)
);
Switch enableCSVLog
(
    optProperties.lookupOrDefault<Switch>("enableCSVLog", false)
);
Switch enableOptimizerRestart
(
    optProperties.lookupOrDefault<Switch>("enableOptimizerRestart", false)
);

// ---------------------------------------------------------------------------
// Transaction test mode: forces first candidate rejection by amplifying
// the step, then verifies full rollback and retry acceptance.
// ---------------------------------------------------------------------------
Switch transactionTestMode
(
    optProperties.lookupOrDefault<Switch>("transactionTestMode", false)
);
scalar transactionTestAmplification
(
    optProperties.lookupOrDefault<scalar>("transactionTestAmplification", 5.0)
);

if
(
    initialMoveLimit <= 0.0
 || minimumMoveLimit <= 0.0
 || maximumMoveLimit < minimumMoveLimit
 || moveLimitReduction <= 0.0 || moveLimitReduction >= 1.0
 || moveLimitGrowth < 1.0
 || maxCandidateRetries < 0
 || relativeObjectiveTolerance < 0.0
 || absoluteConstraintTolerance < 0.0
)
{
    FatalErrorInFunction
        << "Invalid transaction mechanism parameters. Require: "
        << "initialMoveLimit>0, minimumMoveLimit>0, "
        << "maximumMoveLimit>=minimumMoveLimit, "
        << "0<moveLimitReduction<1, moveLimitGrowth>=1, "
        << "maxCandidateRetries>=0, tolerances>=0."
        << abort(FatalError);
}

// enableSSTAcceptanceCheck requires the frozen-RANS flow model so that
// full SST re-convergence can be performed during candidate evaluation.
if (enableSSTAcceptanceCheck && flowModel != "incompressibleRANSFrozen")
{
    FatalErrorInFunction
        << "enableSSTAcceptanceCheck=true requires "
        << "flowModel=incompressibleRANSFrozen."
        << abort(FatalError);
}

// The transaction mechanism evaluates candidates with full SST.  If the
// accepted primal was solved with frozen turbulence, the objective
// comparison would be inconsistent.  Warn (do not abort) so that users
// can still run mixed configurations for debugging.
if (enableSSTAcceptanceCheck && freezeTurbulenceForValidation)
{
    WarningInFunction
        << "enableSSTAcceptanceCheck=true with "
        << "freezeTurbulenceForValidation=true: the accepted design's "
        << "objective is based on frozen turbulence while candidates use "
        << "full SST.  Set freezeTurbulenceForValidation=false for "
        << "consistent comparisons." << endl;
}

// Runtime move limit (adapts during optimization)
scalar currentMoveLimit = initialMoveLimit;
label consecutiveAccepts = 0;

Info<< "Transaction mechanism configuration:" << endl;
Info<< "  enableSSTAcceptanceCheck=" << enableSSTAcceptanceCheck << endl;
Info<< "  initialMoveLimit=" << initialMoveLimit << endl;
Info<< "  minimumMoveLimit=" << minimumMoveLimit << endl;
Info<< "  maximumMoveLimit=" << maximumMoveLimit << endl;
Info<< "  moveLimitReduction=" << moveLimitReduction << endl;
Info<< "  moveLimitGrowth=" << moveLimitGrowth << endl;
Info<< "  maxCandidateRetries=" << maxCandidateRetries << endl;
Info<< "  relativeObjectiveTolerance=" << relativeObjectiveTolerance << endl;
Info<< "  absoluteConstraintTolerance=" << absoluteConstraintTolerance << endl;
Info<< "  enableCSVLog=" << enableCSVLog << endl;
Info<< "  enableOptimizerRestart=" << enableOptimizerRestart << endl;

dimensionedScalar b
(
    "b",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    1.0
);

Switch solid_area(optProperties.lookup("solid_area"));
Switch fluid_area(optProperties.lookup("fluid_area"));
Switch test_area(optProperties.lookup("test_area"));
label flowInnerIter
(
    optProperties.lookupOrDefault<label>("flowInnerIter", 40)
);
label adjointInnerIter
(
    optProperties.lookupOrDefault<label>("adjointInnerIter", 40)
);
scalar vini(readScalar(optProperties.lookup("vini")));

Info<< "Reading pseudo-density fields\n" << endl;
volScalarField xp
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("xp", dimless, vini),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField xh
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("xh", dimless, vini),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField x
(
    IOobject
    (
        "x",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    xp
);

volScalarField drho
(
    IOobject
    (
        "drho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("drho", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField designMask
(
    IOobject
    (
        "designMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("designMask", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Creating heat source field Q\n" << endl;
volScalarField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "Q",
        dimensionSet(0, 0, -1, 1, 0, 0, 0),
        0.0
    )
);

labelList cells_solid;
labelList cells_fluid;

if (solid_area)
{
    const word zoneName("solid_area");
    const label zoneID = mesh.cellZones().findZoneID(zoneName);
    if (zoneID < 0)
    {
        FatalErrorInFunction
            << "cellZone '" << zoneName << "' was requested but not found."
            << abort(FatalError);
    }

    cells_solid = mesh.cellZones()[zoneID];
    forAll(cells_solid, zoneI)
    {
        const label celli = cells_solid[zoneI];
        x[celli] = 0.0;
        xp[celli] = 0.0;
        xh[celli] = 0.0;
        designMask[celli] = 0.0;
    }
}

if (fluid_area)
{
    const word zoneName("fluid_area");
    const label zoneID = mesh.cellZones().findZoneID(zoneName);
    if (zoneID < 0)
    {
        FatalErrorInFunction
            << "cellZone '" << zoneName << "' was requested but not found."
            << abort(FatalError);
    }

    cells_fluid = mesh.cellZones()[zoneID];
    forAll(cells_fluid, zoneI)
    {
        const label celli = cells_fluid[zoneI];
        x[celli] = 0.999;
        xp[celli] = 0.999;
        xh[celli] = 0.999;
        designMask[celli] = 0.0;
    }
}

if (test_area)
{
    const word setName("test2_area");
    IOobject testSetIO
    (
        setName,
        mesh.time().constant(),
        fileName("polyMesh")/"sets",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    cellSet testCells(testSetIO);
    const labelList testCellLabels(testCells.toc());

    if (testCellLabels.empty())
    {
        FatalErrorInFunction
            << "cellSet '" << setName << "' is empty."
            << abort(FatalError);
    }

    forAll(testCellLabels, setI)
    {
        const label celli = testCellLabels[setI];
        x[celli] = 0.4;
        xp[celli] = 0.4;
        xh[celli] = 0.4;
    }
}

x.correctBoundaryConditions();
xp.correctBoundaryConditions();
xh.correctBoundaryConditions();
designMask.correctBoundaryConditions();


volScalarField gsensVol
(
    IOobject
    (
        "gsensVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gsensVol", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gsenshVol
(
    IOobject
    (
        "gsenshVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gsenshVol", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

scalar designVolume = fvc::domainIntegrate(designMask).value();
if (designVolume <= SMALL)
{
    FatalErrorInFunction
        << "Design volume is zero or negative."
        << abort(FatalError);
}

const dictionary& conFunc =
    mesh.solutionDict().subDict("constraintFunctionDict_flow");
const label nObjPatch =
    conFunc.lookupOrDefault<label>("numberConstraintPatches", 2);
const wordList conPatchNames(conFunc.lookup("constraintPatchesNames"));

if (nObjPatch <= 0 || conPatchNames.size() != nObjPatch)
{
    FatalErrorInFunction
        << "numberConstraintPatches=" << nObjPatch
        << " but constraintPatchesNames has size "
        << conPatchNames.size() << '.' << abort(FatalError);
}

labelList conPatchList(nObjPatch, -1);
forAll(conPatchNames, patchI)
{
    conPatchList[patchI] =
        mesh.boundaryMesh().findPatchID(conPatchNames[patchI]);

    if (conPatchList[patchI] < 0)
    {
        FatalErrorInFunction
            << "Constraint patch '" << conPatchNames[patchI]
            << "' was not found." << abort(FatalError);
    }
}
