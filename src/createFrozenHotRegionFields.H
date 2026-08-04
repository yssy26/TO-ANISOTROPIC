Info<< "Reading optional frozen-hot-region configuration\n" << endl;

Switch frozenHotRegionEnabled(false);
Switch strictFrozenHotValidation(true);
Switch mmaUpdateEnabled(true);
word hotZoneName("hotFluidZone");
word wallZoneName("separatorWallZone");
word hotVelocityFieldName("UHotFrozen");
word hotInletPatchName("hotInlet");
word hotOutletPatchName("hotOutlet");
word coldInletPatchName("inlet");
word coldOutletPatchName("outlet");
word thermalObjectiveType("legacyMeanTemperature");

scalar hotMassBalanceTolerance(1e-6);
scalar hotDivergenceTolerance(1e-8);
scalar thermalReferenceTemperature(Foam::max(MeanT0, SMALL));

dimensionedScalar hotRhoCp
(
    "hotRhoCp",
    rhoc.dimensions(),
    rhoc.value()
);

dimensionedScalar hotK
(
    "hotK",
    kf.dimensions(),
    kf.value()
);

// wallK removed: the separator wall now uses the anisotropic solidK tensor
// via the xh=0 interpolation path (K_s/(rho*c)), so no separate scalar is
// needed.

dimensionedScalar coldBlockAlpha
(
    "coldBlockAlpha",
    alpha.dimensions(),
    alphaMAX.value()
);

if (thermalProperties.found("frozenHotRegion"))
{
    const dictionary& frozenHotDict =
        thermalProperties.subDict("frozenHotRegion");

    frozenHotRegionEnabled =
        frozenHotDict.lookupOrDefault<Switch>("enabled", false);
    strictFrozenHotValidation =
        frozenHotDict.lookupOrDefault<Switch>("strictValidation", true);
    hotZoneName =
        frozenHotDict.lookupOrDefault<word>("hotZone", hotZoneName);
    wallZoneName =
        frozenHotDict.lookupOrDefault<word>("wallZone", wallZoneName);
    hotVelocityFieldName =
        frozenHotDict.lookupOrDefault<word>
        (
            "hotVelocityField",
            hotVelocityFieldName
        );
    hotInletPatchName =
        frozenHotDict.lookupOrDefault<word>
        (
            "hotInletPatch",
            hotInletPatchName
        );
    hotOutletPatchName =
        frozenHotDict.lookupOrDefault<word>
        (
            "hotOutletPatch",
            hotOutletPatchName
        );
    coldInletPatchName =
        frozenHotDict.lookupOrDefault<word>
        (
            "coldInletPatch",
            coldInletPatchName
        );
    coldOutletPatchName =
        frozenHotDict.lookupOrDefault<word>
        (
            "coldOutletPatch",
            coldOutletPatchName
        );
    hotMassBalanceTolerance =
        frozenHotDict.lookupOrDefault<scalar>
        (
            "hotMassBalanceTolerance",
            hotMassBalanceTolerance
        );
    hotDivergenceTolerance =
        frozenHotDict.lookupOrDefault<scalar>
        (
            "hotDivergenceTolerance",
            hotDivergenceTolerance
        );

    if (frozenHotRegionEnabled)
    {
        hotRhoCp = dimensionedScalar(frozenHotDict.lookup("hotRhoCp"));
        hotK = dimensionedScalar(frozenHotDict.lookup("hotK"));
        // wallK no longer read: separator wall uses solidK tensor via xh=0.
        coldBlockAlpha =
            dimensionedScalar(frozenHotDict.lookup("coldBlockAlpha"));
    }
}

mmaUpdateEnabled =
    optProperties.lookupOrDefault<Switch>("mmaUpdateEnabled", true);
thermalObjectiveType =
    optProperties.lookupOrDefault<word>
    (
        "thermalObjectiveType",
        frozenHotRegionEnabled
      ? word("maximizeColdOutletTemperature")
      : word("legacyMeanTemperature")
    );
thermalReferenceTemperature =
    optProperties.lookupOrDefault<scalar>
    (
        "thermalReferenceTemperature",
        Foam::max(MeanT0, SMALL)
    );

if (thermalReferenceTemperature <= SMALL)
{
    FatalErrorInFunction
        << "thermalReferenceTemperature must be positive."
        << abort(FatalError);
}

if
(
    frozenHotRegionEnabled
 && thermalObjectiveType != "maximizeColdOutletTemperature"
)
{
    FatalErrorInFunction
        << "Frozen-hot mode supports thermalObjectiveType "
        << "'maximizeColdOutletTemperature' only. Received '"
        << thermalObjectiveType << "'." << abort(FatalError);
}

if
(
   !frozenHotRegionEnabled
 && thermalObjectiveType != "legacyMeanTemperature"
)
{
    FatalErrorInFunction
        << "Legacy mode supports thermalObjectiveType "
        << "'legacyMeanTemperature' only. Received '"
        << thermalObjectiveType << "'." << abort(FatalError);
}

volScalarField hotRegionMask
(
    IOobject
    (
        "hotRegionMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hotRegionMask", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField wallRegionMask
(
    IOobject
    (
        "wallRegionMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wallRegionMask", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField coldRegionMask
(
    IOobject
    (
        "coldRegionMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("coldRegionMask", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

labelList hotRegionCells;
labelList wallRegionCells;

if (frozenHotRegionEnabled)
{
    // Legacy fixed-solid zones conduct heat but must not carry the cold-side
    // convective flux or the cold velocity-adjoint coupling.
    if (solid_area)
    {
        forAll(cells_solid, i)
        {
            coldRegionMask[cells_solid[i]] = 0.0;
        }
    }

    const label hotZoneID = mesh.cellZones().findZoneID(hotZoneName);
    const label wallZoneID = mesh.cellZones().findZoneID(wallZoneName);

    if (hotZoneID < 0)
    {
        FatalErrorInFunction
            << "Frozen-hot mode requires cellZone '" << hotZoneName
            << "'." << abort(FatalError);
    }
    if (wallZoneID < 0)
    {
        FatalErrorInFunction
            << "Frozen-hot mode requires cellZone '" << wallZoneName
            << "'." << abort(FatalError);
    }

    hotRegionCells = mesh.cellZones()[hotZoneID];
    wallRegionCells = mesh.cellZones()[wallZoneID];

    label nHotCells = hotRegionCells.size();
    label nWallCells = wallRegionCells.size();
    reduce(nHotCells, sumOp<label>());
    reduce(nWallCells, sumOp<label>());

    if (nHotCells == 0 || nWallCells == 0)
    {
        FatalErrorInFunction
            << "Hot-fluid and separator-wall zones must both be non-empty."
            << abort(FatalError);
    }

    forAll(hotRegionCells, i)
    {
        const label celli = hotRegionCells[i];
        hotRegionMask[celli] = 1.0;
        coldRegionMask[celli] = 0.0;
        designMask[celli] = 0.0;
        x[celli] = 0.999;
        xp[celli] = 0.999;
        xh[celli] = 0.999;
    }

    forAll(wallRegionCells, i)
    {
        const label celli = wallRegionCells[i];
        if (hotRegionMask[celli] > 0.5)
        {
            FatalErrorInFunction
                << "Cell " << celli
                << " belongs to both hot-fluid and separator-wall zones."
                << abort(FatalError);
        }

        wallRegionMask[celli] = 1.0;
        coldRegionMask[celli] = 0.0;
        designMask[celli] = 0.0;
        x[celli] = 0.0;
        xp[celli] = 0.0;
        xh[celli] = 0.0;
    }
}

hotRegionMask.correctBoundaryConditions();
wallRegionMask.correctBoundaryConditions();
coldRegionMask.correctBoundaryConditions();
designMask.correctBoundaryConditions();
x.correctBoundaryConditions();
xp.correctBoundaryConditions();
xh.correctBoundaryConditions();

designVolume = fvc::domainIntegrate(designMask).value();
if (designVolume <= SMALL)
{
    FatalErrorInFunction
        << "Design volume is zero or negative after fixed-region masking."
        << abort(FatalError);
}

surfaceScalarField designFilterFaceMask
(
    IOobject
    (
        "designFilterFaceMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroDesignFilterFaceMask", dimless, 0.0)
);

scalarField& designFilterInternal =
    designFilterFaceMask.primitiveFieldRef();
const labelUList& filterOwner = mesh.owner();
const labelUList& filterNeighbour = mesh.neighbour();

forAll(designFilterInternal, facei)
{
    designFilterInternal[facei] =
        designMask[filterOwner[facei]]
       *designMask[filterNeighbour[facei]];
}

forAll(designFilterFaceMask.boundaryField(), patchi)
{
    designFilterFaceMask.boundaryFieldRef()[patchi] = 0.0;
}

surfaceScalarField hotFaceMask
(
    IOobject
    (
        "hotFaceMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hotFaceMask", dimless, 0.0)
);

surfaceScalarField coldFaceMask
(
    IOobject
    (
        "coldFaceMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("coldFaceMask", dimless, 0.0)
);

scalarField& hotFaceMaskInternal = hotFaceMask.primitiveFieldRef();
scalarField& coldFaceMaskInternal = coldFaceMask.primitiveFieldRef();
const labelUList& owner = mesh.owner();
const labelUList& neighbour = mesh.neighbour();
label directHotColdFaces = 0;

forAll(hotFaceMaskInternal, facei)
{
    const label own = owner[facei];
    const label nei = neighbour[facei];
    const bool ownHot = hotRegionMask[own] > 0.5;
    const bool neiHot = hotRegionMask[nei] > 0.5;
    const bool ownCold = coldRegionMask[own] > 0.5;
    const bool neiCold = coldRegionMask[nei] > 0.5;

    hotFaceMaskInternal[facei] = ownHot && neiHot ? 1.0 : 0.0;
    coldFaceMaskInternal[facei] = ownCold && neiCold ? 1.0 : 0.0;

    if ((ownHot && neiCold) || (ownCold && neiHot))
    {
        ++directHotColdFaces;
    }
}

forAll(mesh.boundary(), patchi)
{
    fvsPatchScalarField& hotPatch = hotFaceMask.boundaryFieldRef()[patchi];
    fvsPatchScalarField& coldPatch = coldFaceMask.boundaryFieldRef()[patchi];
    const fvPatchScalarField& hotMaskPatch =
        hotRegionMask.boundaryField()[patchi];
    const fvPatchScalarField& coldMaskPatch =
        coldRegionMask.boundaryField()[patchi];
    const scalarField hotOwn(hotMaskPatch.patchInternalField());
    const scalarField coldOwn(coldMaskPatch.patchInternalField());

    if (mesh.boundary()[patchi].coupled())
    {
        const scalarField hotNei(hotMaskPatch.patchNeighbourField());
        const scalarField coldNei(coldMaskPatch.patchNeighbourField());

        forAll(hotPatch, facei)
        {
            const bool ownHot = hotOwn[facei] > 0.5;
            const bool neiHot = hotNei[facei] > 0.5;
            const bool ownCold = coldOwn[facei] > 0.5;
            const bool neiCold = coldNei[facei] > 0.5;

            hotPatch[facei] = ownHot && neiHot ? 1.0 : 0.0;
            coldPatch[facei] = ownCold && neiCold ? 1.0 : 0.0;

            if ((ownHot && neiCold) || (ownCold && neiHot))
            {
                ++directHotColdFaces;
            }
        }
    }
    else
    {
        forAll(hotPatch, facei)
        {
            hotPatch[facei] = hotOwn[facei] > 0.5 ? 1.0 : 0.0;
            coldPatch[facei] = coldOwn[facei] > 0.5 ? 1.0 : 0.0;
        }
    }
}

reduce(directHotColdFaces, sumOp<label>());
if (frozenHotRegionEnabled && directHotColdFaces > 0)
{
    FatalErrorInFunction
        << "Detected " << directHotColdFaces
        << " direct hot-fluid/cold-fluid faces. The separator wall is not "
        << "continuous." << abort(FatalError);
}

label hotInletPatchID = -1;
label hotOutletPatchID = -1;
label coldInletPatchID = -1;
label coldOutletPatchID = -1;

if (frozenHotRegionEnabled)
{
    hotInletPatchID = mesh.boundaryMesh().findPatchID(hotInletPatchName);
    hotOutletPatchID = mesh.boundaryMesh().findPatchID(hotOutletPatchName);
    coldInletPatchID = mesh.boundaryMesh().findPatchID(coldInletPatchName);
    coldOutletPatchID = mesh.boundaryMesh().findPatchID(coldOutletPatchName);

    if
    (
        hotInletPatchID < 0 || hotOutletPatchID < 0
     || coldInletPatchID < 0 || coldOutletPatchID < 0
    )
    {
        FatalErrorInFunction
            << "Frozen-hot mode requires valid hot/cold inlet and outlet "
            << "patch names." << abort(FatalError);
    }
}

IOobject UHotIO
(
    hotVelocityFieldName,
    runTime.timeName(),
    mesh,
    frozenHotRegionEnabled
  ? IOobject::MUST_READ
  : IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> UHotFrozenPtr;
if (frozenHotRegionEnabled)
{
    UHotFrozenPtr.reset(new volVectorField(UHotIO, mesh));
}
else
{
    UHotFrozenPtr.reset
    (
        new volVectorField
        (
            UHotIO,
            mesh,
            dimensionedVector("UHotFrozen", U.dimensions(), vector::zero)
        )
    );
}
volVectorField& UHotFrozen = UHotFrozenPtr();

forAll(UHotFrozen.primitiveField(), celli)
{
    if (hotRegionMask[celli] < 0.5)
    {
        UHotFrozen.primitiveFieldRef()[celli] = vector::zero;
    }
}
UHotFrozen.correctBoundaryConditions();

surfaceScalarField phiHotFrozen
(
    IOobject
    (
        "phiHotFrozen",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    hotFaceMask*fvc::flux(UHotFrozen)
);

surfaceScalarField phiThermal
(
    IOobject
    (
        "phiThermal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

volScalarField thermalObjectiveMask
(
    IOobject
    (
        "thermalObjectiveMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalObjectiveMask", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

if (frozenHotRegionEnabled)
{
    thermalObjectiveMask = coldRegionMask;
}
else
{
    thermalObjectiveMask = scalar(1.0);
}
thermalObjectiveMask.correctBoundaryConditions();

scalar thermalObjectiveVolume =
    fvc::domainIntegrate(thermalObjectiveMask).value();
if (thermalObjectiveVolume <= SMALL)
{
    FatalErrorInFunction
        << "Thermal objective region is empty." << abort(FatalError);
}

scalar coldOutletArea = 0.0;
if (frozenHotRegionEnabled)
{
    coldOutletArea =
        sum(mesh.magSf().boundaryField()[coldOutletPatchID]);
    reduce(coldOutletArea, sumOp<scalar>());

    if (coldOutletArea <= SMALL)
    {
        FatalErrorInFunction
            << "Cold outlet patch has zero area." << abort(FatalError);
    }
}

scalarField thermalObjectiveDerivative(mesh.nCells(), 0.0);
scalarField pressureConstraintDerivative(mesh.nCells(), 0.0);

if (thermalObjectiveType == "maximizeColdOutletTemperature")
{
    const labelUList& outletCells =
        mesh.boundary()[coldOutletPatchID].faceCells();
    const scalarField& outletAreas =
        mesh.magSf().boundaryField()[coldOutletPatchID];

    forAll(outletCells, facei)
    {
        thermalObjectiveDerivative[outletCells[facei]] -=
            outletAreas[facei]
           /Foam::max
            (
                coldOutletArea*thermalReferenceTemperature,
                SMALL
            );
    }
}

if (nObjPatch >= 2)
{
    const label inletPatchI = conPatchList[0];
    const label outletPatchI = conPatchList[1];
    scalar inletArea = sum(mesh.magSf().boundaryField()[inletPatchI]);
    scalar outletArea = sum(mesh.magSf().boundaryField()[outletPatchI]);
    reduce(inletArea, sumOp<scalar>());
    reduce(outletArea, sumOp<scalar>());

    if (!p.boundaryField()[inletPatchI].fixesValue())
    {
        const labelUList& inletCells =
            mesh.boundary()[inletPatchI].faceCells();
        const scalarField& inletAreas =
            mesh.magSf().boundaryField()[inletPatchI];

        forAll(inletCells, facei)
        {
            pressureConstraintDerivative[inletCells[facei]] +=
                rhoFluid*inletAreas[facei]
               /Foam::max(pressureDropMaxPa*inletArea, SMALL);
        }
    }

    if (!p.boundaryField()[outletPatchI].fixesValue())
    {
        const labelUList& outletCells =
            mesh.boundary()[outletPatchI].faceCells();
        const scalarField& outletAreas =
            mesh.magSf().boundaryField()[outletPatchI];

        forAll(outletCells, facei)
        {
            pressureConstraintDerivative[outletCells[facei]] -=
                rhoFluid*outletAreas[facei]
               /Foam::max(pressureDropMaxPa*outletArea, SMALL);
        }
    }
}

const scalar hotRhoCpRatio =
    hotRhoCp.value()/Foam::max(rhoc.value(), SMALL);

if (frozenHotRegionEnabled)
{
    const scalar hotVolume = fvc::domainIntegrate(hotRegionMask).value();
    const scalar wallVolume = fvc::domainIntegrate(wallRegionMask).value();
    tmp<volVectorField> tOutsideHotVelocity =
        (scalar(1.0) - hotRegionMask)*UHotFrozen;
    const scalar maxOutsideHot =
        gMax(mag(tOutsideHotVelocity().primitiveField()));

    scalar hotInFlux =
        -sum(phiHotFrozen.boundaryField()[hotInletPatchID]);
    scalar hotOutFlux =
        sum(phiHotFrozen.boundaryField()[hotOutletPatchID]);
    reduce(hotInFlux, sumOp<scalar>());
    reduce(hotOutFlux, sumOp<scalar>());

    const scalar hotMassBalanceError =
        mag(hotOutFlux - hotInFlux)
       /Foam::max(max(mag(hotInFlux), mag(hotOutFlux)), SMALL);

    tmp<volScalarField> tDivHot = fvc::div(phiHotFrozen);
    const scalar maxHotDivergence =
        gMax(mag(tDivHot().primitiveField()));

    Info<< "Frozen-hot mode enabled:" << nl
        << "  hotZone=" << hotZoneName << ", volume=" << hotVolume << nl
        << "  wallZone=" << wallZoneName << ", volume=" << wallVolume << nl
        << "  hotVelocityField=" << hotVelocityFieldName << nl
        << "  hotRhoCpRatio=" << hotRhoCpRatio << nl
        << "  max|UHotFrozen| outside hot zone=" << maxOutsideHot << nl
        << "  hot inlet flux=" << hotInFlux << nl
        << "  hot outlet flux=" << hotOutFlux << nl
        << "  hot mass-balance error=" << hotMassBalanceError << nl
        << "  max|div(phiHotFrozen)|=" << maxHotDivergence << nl
        << "  thermalObjectiveType=" << thermalObjectiveType << nl
        << "  thermalReferenceTemperature="
        << thermalReferenceTemperature << nl
        << "  mmaUpdateEnabled=" << mmaUpdateEnabled << endl;

    if
    (
        hotMassBalanceError > hotMassBalanceTolerance
     || maxHotDivergence > hotDivergenceTolerance
    )
    {
        if (strictFrozenHotValidation)
        {
            FatalErrorInFunction
                << "Frozen hot-side velocity is not sufficiently conservative. "
                << "massBalanceError=" << hotMassBalanceError
                << ", maxDivergence=" << maxHotDivergence
                << abort(FatalError);
        }
        else
        {
            WarningInFunction
                << "Frozen hot-side velocity conservation tolerances were "
                << "exceeded." << endl;
        }
    }
}
else
{
    Info<< "Frozen-hot mode disabled; legacy single-stream model active."
        << endl;
}
