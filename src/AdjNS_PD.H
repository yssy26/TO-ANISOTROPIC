if (adjointMode == "discrete")
{
    vectorField discreteFlowVelocityRhs(mesh.nCells(), vector::zero);
    scalarField discreteFlowPressureRhs(pressureConstraintDerivative);
    scalarField discreteExternalFaceFluxAdjoint(mesh.nFaces(), 0.0);
    volVectorField& discreteAdjointVelocity = Uc;
    volScalarField& discreteAdjointPressure = pc;
    const word discreteAdjointLabel("pressureDrop");

    const Switch useStageB4DiagnosticSolver =
        optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false);
    scalar selectedDiscreteResidual = GREAT;
    if (useStageB4DiagnosticSolver)
    {
        {
            #include "solveDiscreteFlowAdjoint.H"
            selectedDiscreteResidual = discreteRelativeResidual;
        }
    }
    else
    {
        {
            #include "solveDiscreteFlowAdjointProduction.H"
            selectedDiscreteResidual = discreteRelativeResidual;
        }
    }

    adjointResidualU = selectedDiscreteResidual;
    adjointResidualP = 0.0;
}
else
{
    volScalarField pressureAdjointContinuitySource
    (
        IOobject
        (
            "pressureAdjointContinuitySource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(phic)*scalar(0.0)
    );

    forAll(pressureConstraintDerivative, celli)
    {
        pressureAdjointContinuitySource[celli] =
            pressureConstraintDerivative[celli]
           /Foam::max(mesh.V()[celli], SMALL);
    }
    pressureAdjointContinuitySource.correctBoundaryConditions();

    for (i = 0; i < adjointInnerIter; ++i)
    {
        volVectorField adjointTransposeConvectionc((Uc & fvc::grad(U)));

        tmp<fvVectorMatrix> tUcEqn
        (
            fvm::div(-phi, Uc)
          + adjointTransposeConvectionc
          - fvm::laplacian(nuEffFrozen, Uc)
          + fvm::Sp(alpha, Uc)
         ==
            fvOptions(Uc)
        );
        fvVectorMatrix& UcEqn = tUcEqn.ref();
        UcEqn.relax();
        fvOptions.constrain(UcEqn);

        Info<< "AdjNS_PD.H - solving UcEqn, iter " << i << endl;
        {
            SolverPerformance<vector> perfUc =
                solve(UcEqn == -fvc::grad(pc));
            adjointResidualU = mag(perfUc.finalResidual());
        }
        fvOptions.correct(Uc);

        volScalarField rAUc(1.0/UcEqn.A());
        if (gMin(rAUc) <= 0)
        {
            FatalErrorInFunction
                << "rAUc contains non-positive values."
                << abort(FatalError);
        }

        volVectorField HbyAc("HbyAc", Uc);
        HbyAc = rAUc*UcEqn.H();
        surfaceScalarField phiHbyAc("phiHbyAc", fvc::flux(HbyAc));
        adjustPhi(phiHbyAc, Uc, pc);

        tmp<volScalarField> rAtUc(rAUc);
        if (simple.consistent())
        {
            rAtUc = 1.0/(1.0/rAUc - UcEqn.H1());
            phiHbyAc +=
                fvc::interpolate(rAtUc() - rAUc)
               *fvc::snGrad(pc)*mesh.magSf();
            HbyAc -= (rAUc - rAtUc())*fvc::grad(pc);
        }

        tUcEqn.clear();
        constrainPressure(pc, Uc, phiHbyAc, rAtUc());

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pcEqn
            (
                fvm::laplacian(rAtUc(), pc)
             == fvc::div(phiHbyAc) + pressureAdjointContinuitySource
            );
            pcEqn.setReference(pcRefCell, pcRefValue);
            SolverPerformance<scalar> perfPc = pcEqn.solve();
            adjointResidualP = perfPc.finalResidual();

            if (simple.finalNonOrthogonalIter())
            {
                phic = phiHbyAc - pcEqn.flux();
            }
        }

        pc.relax();
        Uc = HbyAc - rAtUc()*fvc::grad(pc);
        Uc.correctBoundaryConditions();
        fvOptions.correct(Uc);

        Uc.storePrevIter();
        pc.storePrevIter();
        phic.storePrevIter();
    }
}
