const bool ransFlowModel = (flowModel == "incompressibleRANSFrozen");
bool turbulenceConverged = !ransFlowModel;

volScalarField fluidIndicator
(
    IOobject
    (
        "fluidIndicator",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroFluidIndicator", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField chiTurbulence
(
    IOobject
    (
        "chiTurbulence",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroChiTurbulence", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField nutFrozen
(
    IOobject
    (
        "nutFrozen",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroNutFrozen", nu.dimensions(), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField nuEffFrozen
(
    IOobject
    (
        "nuEffFrozen",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    nu,
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alphaTurbulent
(
    IOobject
    (
        "alphaTurbulent",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroAlphaTurbulent", DTMolecular.dimensions(), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volSymmTensorField DTEffective
(
    IOobject
    (
        "DTEffective",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    DTMolecular
);

#include "updateFrozenTurbulenceFields.H"