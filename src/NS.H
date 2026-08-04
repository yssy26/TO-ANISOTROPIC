const label flowCorrectorLimit =
    ransFlowModel ? turbulenceMaxCorrectors : flowInnerIter;
label turbulenceConvergedCount = 0;
turbulenceConverged = !ransFlowModel;

vectorField previousU;
scalarField previousP;
scalarField previousK;
scalarField previousOmega;
scalarField previousNut;

if (ransFlowModel)
{
    previousU = U.primitiveField();
    previousP = p.primitiveField();
    previousK =
        mesh.lookupObject<volScalarField>("k").primitiveField();
    previousOmega =
        mesh.lookupObject<volScalarField>("omega").primitiveField();
    previousNut = nutFrozen.primitiveField();
}

for (i = 0; i < flowCorrectorLimit; ++i)
{
    forwardPrimalOuterIterations = i + 1;
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nuEffFrozen, U)
      + fvm::Sp(alpha, U)
     ==
       -fvc::grad(p)
      + fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    if (ransFlowModel)
    {
        UEqn -= fvc::div
        (
            nuEffFrozen*dev2(Foam::T(fvc::grad(U)))
        );
    }

    UEqn.relax();
    Info<< "NS.H solving, corrector " << i + 1 << endl;
    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        SolverPerformance<vector> perfU =
            solve(UEqn == -fvc::grad(p));
        forwardInitialResidualU = mag(perfU.initialResidual());
        forwardFinalResidualU = mag(perfU.finalResidual());
        fvOptions.correct(U);
    }

    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)
           *fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    primalPressureMobility = rAtU().primitiveField();
    tUEqn.clear();

    constrainPressure(p, U, phiHbyA, rAtU());

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
        SolverPerformance<scalar> perfP = pEqn.solve();
        forwardInitialResidualP = perfP.initialResidual();
        forwardFinalResidualP = perfP.finalResidual();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    p.relax();
    U = HbyA - rAtU()*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    if (ransFlowModel)
    {
        if (!freezeTurbulenceForValidation)
        {
            turbulence->correct();
        }
        #include "updateFrozenTurbulenceFields.H"

        const volScalarField& kField =
            mesh.lookupObject<volScalarField>("k");
        const volScalarField& omegaField =
            mesh.lookupObject<volScalarField>("omega");

        auto coldRelativeScalarChange =
        [&](const scalarField& current, const scalarField& previous)
        {
            scalar numerator = 0.0;
            scalar denominator = 0.0;
            forAll(current, celli)
            {
                if (coldRegionMask[celli] > 0.5)
                {
                    numerator += mesh.V()[celli]
                       *sqr(current[celli] - previous[celli]);
                    denominator += mesh.V()[celli]*sqr(current[celli]);
                }
            }
            reduce(numerator, sumOp<scalar>());
            reduce(denominator, sumOp<scalar>());
            return Foam::sqrt
            (
                numerator/Foam::max(denominator, scalar(VSMALL))
            );
        };
        auto activeTurbulenceRelativeScalarChange =
        [&](const scalarField& current, const scalarField& previous)
        {
            scalar numerator = 0.0;
            scalar denominator = 0.0;
            forAll(current, celli)
            {
                const scalar weight =
                    mesh.V()[celli]*chiTurbulence[celli];
                numerator += weight*sqr(current[celli] - previous[celli]);
                denominator += weight*sqr(current[celli]);
            }
            reduce(numerator, sumOp<scalar>());
            reduce(denominator, sumOp<scalar>());
            return Foam::sqrt
            (
                numerator/Foam::max(denominator, scalar(VSMALL))
            );
        };

        scalar uNumerator = 0.0;
        scalar uDenominator = 0.0;
        forAll(U.primitiveField(), celli)
        {
            if (coldRegionMask[celli] > 0.5)
            {
                uNumerator += mesh.V()[celli]
                   *magSqr(U[celli] - previousU[celli]);
                uDenominator += mesh.V()[celli]*magSqr(U[celli]);
            }
        }
        reduce(uNumerator, sumOp<scalar>());
        reduce(uDenominator, sumOp<scalar>());

        const scalar uChange = Foam::sqrt
        (
            uNumerator/Foam::max(uDenominator, scalar(VSMALL))
        );
        const scalar pChange = coldRelativeScalarChange
        (
            p.primitiveField(),
            previousP
        );
        const scalar kChange = activeTurbulenceRelativeScalarChange
        (
            kField.primitiveField(),
            previousK
        );
        const scalar omegaChange = activeTurbulenceRelativeScalarChange
        (
            omegaField.primitiveField(),
            previousOmega
        );
        const scalar nutChange = activeTurbulenceRelativeScalarChange
        (
            nutFrozen.primitiveField(),
            previousNut
        );
        // The frozen primal/adjoint operators depend on nutFrozen, not on
        // omega directly. Keep omega as a diagnostic while gating the state
        // and the effective frozen transport coefficient that are consumed.
        const scalar maximumRansChange = Foam::max
        (
            Foam::max(uChange, pChange),
            Foam::max(kChange, nutChange)
        );

        // Phase K: expose for CSV logging
        lastRelativeChangeK = kChange;
        lastRelativeChangeOmega = omegaChange;

        if
        (
            i + 1 >= turbulenceMinCorrectors
         && maximumRansChange < turbulenceConvergenceTolerance
        )
        {
            ++turbulenceConvergedCount;
        }
        else
        {
            turbulenceConvergedCount = 0;
        }

        if ((i + 1) % 10 == 0 || turbulenceConvergedCount > 0)
        {
            Info<< "Frozen-RANS relative changes: U=" << uChange
                << ", p=" << pChange
                << ", k=" << kChange
                << ", omega=" << omegaChange
                << ", nut=" << nutChange
                << ", consecutive=" << turbulenceConvergedCount << endl;
        }

        previousU = U.primitiveField();
        previousP = p.primitiveField();
        previousK = kField.primitiveField();
        previousOmega = omegaField.primitiveField();
        previousNut = nutFrozen.primitiveField();

        if (turbulenceConvergedCount >= 5)
        {
            turbulenceConverged = true;
            Info<< "Frozen-RANS flow converged after " << i + 1
                << " correctors." << endl;
            U.storePrevIter();
            p.storePrevIter();
            phi.storePrevIter();
            break;
        }
    }

    U.storePrevIter();
    p.storePrevIter();
    phi.storePrevIter();
}

if (ransFlowModel && !turbulenceConverged)
{
    if (solveFlowAdjoints || mmaUpdateEnabled)
    {
        FatalErrorInFunction
            << "Frozen-RANS flow did not satisfy five consecutive relative-"
            << "change checks below " << turbulenceConvergenceTolerance
            << " after " << turbulenceMaxCorrectors << " correctors. "
            << "Adjoint/design update is blocked." << abort(FatalError);
    }
    else
    {
        WarningInFunction
            << "Forward-only RANS warm-up reached "
            << turbulenceMaxCorrectors << " correctors without satisfying "
            << "the final convergence gate. No adjoint or design update will "
            << "be evaluated; continue the same design at the next outer "
            << "iteration." << endl;
    }
}
