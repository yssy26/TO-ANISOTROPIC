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

// mmaUpdateEnabled is assigned from optProperties in
// createFrozenHotRegionFields.H; declare it here (top-level) so the
// anisotropic MMA gate in this file can reference the case value.
Switch mmaUpdateEnabled =
    optProperties.lookupOrDefault<Switch>("mmaUpdateEnabled", true);

// --- Feature flag and validation gates for anisotropic conductivity ---
const Switch useAnisotropicConductivity =
    optProperties.lookupOrDefault<Switch>("useAnisotropicConductivity", false);
const Switch anisotropicConductivityValidated =
    optProperties.lookupOrDefault<Switch>
    (
        "anisotropicConductivityValidated",
        false
    );
const Switch allowOffDiagonalAnisotropy =
    optProperties.lookupOrDefault<Switch>
    (
        "allowOffDiagonalAnisotropy",
        false
    );
const scalar offDiagonalConductivityTolerance =
    optProperties.lookupOrDefault<scalar>
    (
        "offDiagonalConductivityTolerance",
        1e-12
    );

if (offDiagonalConductivityTolerance < 0.0)
{
    FatalErrorInFunction
        << "offDiagonalConductivityTolerance must be non-negative."
        << abort(FatalError);
}

// Solid conductivity tensor KsTensor and fluid identity tensor kfIdentity.
// When useAnisotropicConductivity=false, KsTensor = ks*I (isotropic fallback).
// When useAnisotropicConductivity=true,  KsTensor = solidK (full anisotropic).
symmTensor KsTensor;
if (useAnisotropicConductivity)
{
    dimensionedSymmTensor solidK(thermalProperties.lookup("solidK"));
    const symmTensor& K = solidK.value();

    // Positive-definiteness via Sylvester's criterion. OpenFOAM 7's
    // eigenValues() can be unreliable for repeated eigenvalues such as k*I.
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

    const scalar diagonalScale = Foam::max
    (
        Foam::max(mag(K.xx()), mag(K.yy())),
        Foam::max(mag(K.zz()), scalar(SMALL))
    );
    const scalar maxOffDiagonal = Foam::max
    (
        Foam::max(mag(K.xy()), mag(K.xz())),
        mag(K.yz())
    );
    const scalar relativeOffDiagonal = maxOffDiagonal/diagonalScale;

    if
    (
       !allowOffDiagonalAnisotropy
     && relativeOffDiagonal > offDiagonalConductivityTolerance
    )
    {
        FatalErrorInFunction
            << "Off-diagonal solidK components are disabled until the "
            << "general-tensor flux and adjoint tests are passed. solidK="
            << K << ", relativeOffDiagonal=" << relativeOffDiagonal
            << ", tolerance=" << offDiagonalConductivityTolerance << ". "
            << "Use an axis-aligned diagonal tensor, or explicitly set "
            << "allowOffDiagonalAnisotropy=true for validation-only runs."
            << abort(FatalError);
    }
    if (allowOffDiagonalAnisotropy && relativeOffDiagonal > SMALL)
    {
        WarningInFunction
            << "General off-diagonal anisotropy is enabled. The current "
            << "Stage A acceptance suite must include rotated-tensor flux "
            << "and adjoint checks before production optimization." << endl;
    }

    KsTensor = K;
    Info<< "Anisotropic conductivity: solidK=" << KsTensor << endl;
    Info<< "Principal mesh-axis conductivities: Kxx=" << K.xx()
        << ", Kyy=" << K.yy() << ", Kzz=" << K.zz()
        << ", relativeOffDiagonal=" << relativeOffDiagonal << endl;

    if (!anisotropicConductivityValidated)
    {
        WarningInFunction
            << "anisotropicConductivityValidated=false. Forward validation "
            << "is permitted, but topology updates are blocked until the "
            << "Stage A suite passes." << endl;
    }
    if (mmaUpdateEnabled && !anisotropicConductivityValidated)
    {
        FatalErrorInFunction
            << "MMA is blocked for anisotropic conductivity until Stage A "
            << "validation passes. Set anisotropicConductivityValidated=true "
            << "only after tests/test_anisotropic_material_model.py and the "
            << "case-level acceptance suite report PASS."
            << abort(FatalError);
    }
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
// DTMolecular replaces the old scalar DT. Type is always volSymmTensorField
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

// Initialize DTMolecular from current xh. The frozen-hot-region data is not
// available yet at this include point, so updateMaterialProperties.H cannot
// be used here.
forAll(DTMolecular.primitiveFieldRef(), celli)
{
    const scalar g = xh[celli]*(1.0 + qu)
                   / (qu + xh[celli] + SMALL);
    DTMolecular.primitiveFieldRef()[celli] =
        (KsTensor + (kfIdentity - KsTensor)*g) / rhoc.value();
}
DTMolecular.correctBoundaryConditions();

// --- Thermal diffusivity derivative dD/dxh (tensor) ---
// Declared here because AdjHeatTransfer.H needs it before sensitivity.H runs.
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
