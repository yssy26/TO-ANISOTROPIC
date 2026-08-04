Info<< "Reading thermal properties\n" << endl;
IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoc(thermalProperties.lookup("rhoc"));
dimensionedScalar kf(thermalProperties.lookup("kf"));

// --- Feature flag: anisotropic conductivity ---
const Switch useAnisotropicConductivity =
    optProperties.lookupOrDefault<Switch>("useAnisotropicConductivity", false);

// Solid conductivity tensor KsTensor and fluid identity tensor kfIdentity.
// When useAnisotropicConductivity=false, KsTensor = ks*I (isotropic fallback).
// When useAnisotropicConductivity=true,  KsTensor = solidK (full anisotropic).
symmTensor KsTensor;
if (useAnisotropicConductivity)
{
    dimensionedSymmTensor solidK(thermalProperties.lookup("solidK"));
    const symmTensor& K = solidK.value();
    // Check positive-definiteness via Sylvester's criterion (leading
    // principal minors > 0).  This avoids OpenFOAM 7's eigenValues()
    // which returns wrong results for degenerate (repeated-eigenvalue)
    // tensors such as isotropic k*I — it reports (0, 0, k) instead of
    // (k, k, k) due to numerical issues in the Jacobi iteration.
    const scalar minor1 = K.xx();
    const scalar minor2 = K.xx()*K.yy() - K.xy()*K.xy();
    const scalar minor3 =
        K.xx()*(K.yy()*K.zz() - K.yz()*K.yz())
      - K.xy()*(K.xy()*K.zz() - K.yz()*K.xz())
      + K.xz()*(K.xy()*K.yz() - K.yy()*K.xz());
    if (minor1 <= 0 || minor2 <= 0 || minor3 <= 0)
    {
        FatalErrorInFunction
            << "solidK must be positive definite (Sylvester's criterion)."
            << " solidK=" << K
            << " leading minors=(" << minor1 << ' ' << minor2 << ' '
            << minor3 << ')' << abort(FatalError);
    }
    KsTensor = solidK.value();
    Info<< "Anisotropic conductivity: solidK=" << KsTensor << endl;
}
else
{
    dimensionedScalar ks(thermalProperties.lookup("ks"));
    KsTensor = symmTensor
    (
        ks.value(), 0, 0,
        ks.value(), 0,
        ks.value()
    );
    Info<< "Isotropic conductivity (tensor path): ks=" << ks.value() << endl;
}

const symmTensor kfIdentity
(
    kf.value(), 0, 0,
    kf.value(), 0,
    kf.value()
);

// --- Molecular thermal diffusivity (tensor) ---
// DTMolecular replaces the old scalar DT.  Type is always volSymmTensorField
// so that all code paths are type-checked at compile time.
volSymmTensorField DTMolecular
(
    IOobject
    (
        "DTMolecular",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "DTMolecularInit",
        kf.dimensions()/rhoc.dimensions(),
        symmTensor::zero
    ),
    zeroGradientFvPatchScalarField::typeName
);

// Initialize DTMolecular from current xh (inline, not via
// updateMaterialProperties.H, because frozen-hot-region data is not yet
// available at this include point).
forAll(DTMolecular.primitiveFieldRef(), celli)
{
    const scalar g = xh[celli]*(1.0 + qu)
                   / (qu + xh[celli] + SMALL);
    DTMolecular.primitiveFieldRef()[celli] =
        (KsTensor + (kfIdentity - KsTensor)*g) / rhoc.value();
}
DTMolecular.correctBoundaryConditions();

// --- Thermal diffusivity derivative dD/dxh (tensor) ---
// Declared here (not in sensitivity.H) because AdjHeatTransfer.H needs it
// before sensitivity.H runs.  Computed in updateMaterialProperties.H.
volSymmTensorField dDTDxh
(
    IOobject
    (
        "dDTDxh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "zeroDDTDxh",
        DTMolecular.dimensions(),
        symmTensor::zero
    ),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
mesh.setFluxRequired(T.name());

Info<< "Reading field Tb\n" << endl;
volScalarField Tb
(
    IOobject
    (
        "Tb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField fsensMeanT
(
    IOobject
    (
        "fsensMeanT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   -alphaMax*(1.0 + qu)*qu
   /sqr(qu + xh + SMALL)*(U & Ub),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField fsenshMeanT
(
    IOobject
    (
        "fsenshMeanT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fsensMeanT
);
