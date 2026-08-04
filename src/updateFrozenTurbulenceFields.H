fluidIndicator = dimensionedScalar
(
    "zeroFluidIndicator",
    dimless,
    0.0
);
chiTurbulence = dimensionedScalar
(
    "zeroChiTurbulence",
    dimless,
    0.0
);
if (!freezeTurbulenceForValidation)
{
    nutFrozen = dimensionedScalar
    (
        "zeroNutFrozen",
        nutFrozen.dimensions(),
        0.0
    );
}

if (ransFlowModel)
{
    volScalarField& kRANS =
        mesh.lookupObjectRef<volScalarField>("k");
    volScalarField& omegaRANS =
        mesh.lookupObjectRef<volScalarField>("omega");
    volScalarField& nutRANS =
        mesh.lookupObjectRef<volScalarField>("nut");

    forAll(xh.primitiveField(), celli)
    {
        if (coldRegionMask[celli] <= 0.5)
        {
            kRANS[celli] = 1e-10;
            omegaRANS[celli] = 7.38e4;
            nutRANS[celli] = 0.0;
            nutFrozen[celli] = 0.0;
        }

        const scalar indicator = Foam::min
        (
            Foam::max(xh[celli]/0.999, scalar(0.0)),
            scalar(1.0)
        );
        fluidIndicator[celli] = indicator;
        chiTurbulence[celli] =
            coldRegionMask[celli]
           *Foam::pow(indicator, turbulenceDampingExponent);
        if (!freezeTurbulenceForValidation)
        {
            nutFrozen[celli] = chiTurbulence[celli]
               *Foam::max(nutRANS[celli], scalar(0.0));
        }
    }

}

nuEffFrozen = nu + nutFrozen;
alphaTurbulent = nutFrozen/turbulentPrandtl;

// Identity tensor for extending scalar alphaTurbulent to tensor addition.
// static const ensures single initialization across multiple includes.
static const dimensionedSymmTensor I_dim
(
    "I",
    dimless,
    symmTensor(1, 0, 0, 1, 0, 1)
);
DTEffective = DTMolecular + alphaTurbulent * I_dim;

fluidIndicator.correctBoundaryConditions();
chiTurbulence.correctBoundaryConditions();
nutFrozen.correctBoundaryConditions();
nuEffFrozen.correctBoundaryConditions();
alphaTurbulent.correctBoundaryConditions();
DTEffective.correctBoundaryConditions();