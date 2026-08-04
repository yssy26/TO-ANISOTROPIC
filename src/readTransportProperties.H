Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu(transportProperties.lookup("nu"));
dimensionedScalar alphaMax(transportProperties.lookup("alphaMax"));
dimensionedScalar alphaMAX(transportProperties.lookup("alphaMAX"));
dimensionedScalar alphamax(transportProperties.lookup("alphaMax"));
dimensionedScalar alphamin(transportProperties.lookup("alphamin"));
scalar qu(readScalar(optProperties.lookup("qu")));

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphaMax*qu*(1.0 - xh)/(qu + xh + SMALL)
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

scalarField primalPressureMobility(mesh.nCells(), 0.0);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
turbulence->validate();

Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phia
(
    IOobject
    (
        "phia",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ua) & mesh.Sf()
);

label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell(pa, simple.dict(), paRefCell, paRefValue);
mesh.setFluxRequired(pa.name());

Info<< "Reading field pb\n" << endl;
volScalarField pb
(
    IOobject
    (
        "pb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ub\n" << endl;
volVectorField Ub
(
    IOobject
    (
        "Ub",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phib
(
    IOobject
    (
        "phib",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ub) & mesh.Sf()
);

label pbRefCell = 0;
scalar pbRefValue = 0.0;
setRefCell(pb, simple.dict(), pbRefCell, pbRefValue);
mesh.setFluxRequired(pb.name());

volScalarField gsensPowerDiss
(
    IOobject
    (
        "gsensPowerDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   -alphaMax*(1.0 + qu)*qu
   /sqr(qu + xh + SMALL)*(U & Ua),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gsenshPowerDiss
(
    IOobject
    (
        "gsenshPowerDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gsensPowerDiss
);

Info<< "Reading field pc\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Uc\n" << endl;
volVectorField Uc
(
    IOobject
    (
        "Uc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phic
(
    IOobject
    (
        "phic",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

label pcRefCell = 0;
scalar pcRefValue = 0.0;
setRefCell(pc, simple.dict(), pcRefCell, pcRefValue);
mesh.setFluxRequired(pc.name());

volScalarField gsensPressureDrop
(
    IOobject
    (
        "gsensPressureDrop",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   -alphaMax*(1.0 + qu)*qu
   /sqr(qu + xh + SMALL)*(U & Uc),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gsenshPressureDrop
(
    IOobject
    (
        "gsenshPressureDrop",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gsensPressureDrop
);
