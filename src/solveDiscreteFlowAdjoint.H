if (Pstream::parRun())
{
    FatalErrorInFunction
        << "The reduced discrete flow adjoint currently supports serial "
        << "validation only; processor-patch transpose terms are not yet "
        << "implemented."
        << abort(FatalError);
}

const label discreteNCells = mesh.nCells();
const label discreteNVelocity = 3*discreteNCells;
const label discreteNUnknowns = 4*discreteNCells;
const labelUList& discreteOwner = mesh.owner();
const labelUList& discreteNeighbour = mesh.neighbour();
const scalarField& discreteWeights = mesh.weights().primitiveField();
const vectorField& discreteSf = mesh.Sf().primitiveField();
const scalarField& discretePhi = phi.primitiveField();

// The scalar Laplacian is assembled exactly as in NS.H. In frozen-RANS mode
// the explicit deviatoric stress is applied separately with a matched
// reverse-order transpose below.
tmp<fvVectorMatrix> tDiscretePrimalMomentum
(
    fvm::div(phi, U)
  - fvm::laplacian(nuEffFrozen, U)
  + fvm::Sp(alpha, U)
);
fvVectorMatrix& discretePrimalMomentum = tDiscretePrimalMomentum.ref();

scalarField discreteDiagX(discretePrimalMomentum.diag());
scalarField discreteDiagY(discretePrimalMomentum.diag());
scalarField discreteDiagZ(discretePrimalMomentum.diag());

forAll(discretePrimalMomentum.internalCoeffs(), patchi)
{
    const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
    const vectorField& patchDiag =
        discretePrimalMomentum.internalCoeffs()[patchi];

    forAll(patchCells, facei)
    {
        const label celli = patchCells[facei];
        discreteDiagX[celli] += patchDiag[facei].x();
        discreteDiagY[celli] += patchDiag[facei].y();
        discreteDiagZ[celli] += patchDiag[facei].z();
    }
}

const scalarField& discreteLower = discretePrimalMomentum.lower();
const scalarField& discreteUpper = discretePrimalMomentum.upper();

auto discreteUIndex = [&](const label celli, const direction cmpt)
{
    return 3*celli + cmpt;
};
auto discretePIndex = [&](const label celli)
{
    return discreteNVelocity + celli;
};
auto applyFrozenDeviatoricJT =
[&](const scalarField& input, scalarField& output)
{
    if (!ransFlowModel)
    {
        return;
    }

    tensorField stressAdjoint(discreteNCells, tensor::zero);
    tensorField gradientAdjoint(discreteNCells, tensor::zero);

    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const vector lambdaOwn
        (
            input[discreteUIndex(own, 0)],
            input[discreteUIndex(own, 1)],
            input[discreteUIndex(own, 2)]
        );
        const vector lambdaNei
        (
            input[discreteUIndex(nei, 0)],
            input[discreteUIndex(nei, 1)],
            input[discreteUIndex(nei, 2)]
        );
        const tensor faceStressAdjoint =
            discreteSf[facei]*(lambdaNei - lambdaOwn);

        stressAdjoint[own] +=
            discreteWeights[facei]*faceStressAdjoint;
        stressAdjoint[nei] +=
            (1.0 - discreteWeights[facei])*faceStressAdjoint;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& patchMagSf =
            mesh.magSf().boundaryField()[patchi];
        const scalarField& patchDelta =
            mesh.deltaCoeffs().boundaryField()[patchi];
        const scalarField& patchNuEff =
            nuEffFrozen.boundaryField()[patchi];
        const bool velocityFixed = U.boundaryField()[patchi].fixesValue();

        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            const vector lambdaCell
            (
                input[discreteUIndex(celli, 0)],
                input[discreteUIndex(celli, 1)],
                input[discreteUIndex(celli, 2)]
            );
            const tensor boundaryStressAdjoint =
               -patchSf[facei]*lambdaCell;
            const tensor boundaryGradientAdjoint =
                patchNuEff[facei]
               *dev2(Foam::T(boundaryStressAdjoint));
            const vector normal =
                patchSf[facei]/Foam::max(patchMagSf[facei], SMALL);

            gradientAdjoint[celli] +=
                boundaryGradientAdjoint
              - normal*(normal & boundaryGradientAdjoint);

            if (velocityFixed)
            {
                const vector normalGradientAdjoint =
                   -patchDelta[facei]
                   *(normal & boundaryGradientAdjoint);
                for (direction cmpt = 0; cmpt < 3; ++cmpt)
                {
                    output[discreteUIndex(celli, cmpt)] +=
                        normalGradientAdjoint.component(cmpt);
                }
            }
        }
    }

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        gradientAdjoint[celli] +=
            nuEffFrozen[celli]
           *dev2(Foam::T(stressAdjoint[celli]));
    }

    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const vector faceVelocityAdjoint =
            discreteSf[facei]
          &
            (
                gradientAdjoint[own]/Foam::max(mesh.V()[own], SMALL)
              - gradientAdjoint[nei]/Foam::max(mesh.V()[nei], SMALL)
            );

        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(own, cmpt)] +=
                discreteWeights[facei]
               *faceVelocityAdjoint.component(cmpt);
            output[discreteUIndex(nei, cmpt)] +=
                (1.0 - discreteWeights[facei])
               *faceVelocityAdjoint.component(cmpt);
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            continue;
        }

        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            const vector boundaryVelocityAdjoint =
                patchSf[facei]
              &
                (
                    gradientAdjoint[celli]
                   /Foam::max(mesh.V()[celli], SMALL)
                );
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                output[discreteUIndex(celli, cmpt)] +=
                    boundaryVelocityAdjoint.component(cmpt);
            }
        }
    }
};

auto applyPressureFluxCorrection =
[&](const scalarField& input, scalarField& faceFluxCorrection)
{
    faceFluxCorrection = scalar(0.0);
    vectorField pressureGradientSum(discreteNCells, vector::zero);

    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar pressureFace =
            discreteWeights[facei]*input[discretePIndex(own)]
          + (1.0 - discreteWeights[facei])
           *input[discretePIndex(nei)];
        pressureGradientSum[own] += discreteSf[facei]*pressureFace;
        pressureGradientSum[nei] -= discreteSf[facei]*pressureFace;
    }

    forAll(mesh.boundary(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            pressureGradientSum[celli] +=
                patchSf[facei]*input[discretePIndex(celli)];
        }
    }

    vectorField mobilityGradient(discreteNCells, vector::zero);
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        mobilityGradient[celli] =
            primalPressureMobility[celli]
           *pressureGradientSum[celli]
           /Foam::max(mesh.V()[celli], SMALL);
    }

    const scalarField& deltaCoeffs = mesh.deltaCoeffs().primitiveField();
    const scalarField& faceAreas = mesh.magSf().primitiveField();
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar weight = discreteWeights[facei];
        const scalar mobilityFace =
            weight*primalPressureMobility[own]
          + (1.0 - weight)*primalPressureMobility[nei];
        const scalar normalPressureFlux =
            mobilityFace*deltaCoeffs[facei]*faceAreas[facei]
           *(
                input[discretePIndex(nei)]
              - input[discretePIndex(own)]
            );
        faceFluxCorrection[facei] =
            discreteSf[facei]
          & (
                weight*mobilityGradient[own]
              + (1.0 - weight)*mobilityGradient[nei]
            );
        faceFluxCorrection[facei] -= normalPressureFlux;
    }

    forAll(mesh.boundary(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& patchDelta =
            mesh.deltaCoeffs().boundaryField()[patchi];
        const scalarField& patchArea =
            mesh.magSf().boundaryField()[patchi];
        const label patchStart = mesh.boundary()[patchi].start();
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();

        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            const label globalFace = patchStart + facei;
            faceFluxCorrection[globalFace] =
                patchSf[facei] & mobilityGradient[celli];
            if (pressureFixed)
            {
                faceFluxCorrection[globalFace] +=
                    primalPressureMobility[celli]
                   *patchDelta[facei]*patchArea[facei]
                   *input[discretePIndex(celli)];
            }
        }
    }
};

auto applyPressureFluxCorrectionTranspose =
[&](const scalarField& faceAdjoint, scalarField& output)
{
    vectorField mobilityGradientAdjoint(discreteNCells, vector::zero);
    const scalarField& deltaCoeffs = mesh.deltaCoeffs().primitiveField();
    const scalarField& faceAreas = mesh.magSf().primitiveField();

    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar weight = discreteWeights[facei];
        const scalar multiplier = faceAdjoint[facei];
        mobilityGradientAdjoint[own] +=
            weight*discreteSf[facei]*multiplier;
        mobilityGradientAdjoint[nei] +=
            (1.0 - weight)*discreteSf[facei]*multiplier;

        const scalar mobilityFace =
            weight*primalPressureMobility[own]
          + (1.0 - weight)*primalPressureMobility[nei];
        const scalar normalCoefficient =
            mobilityFace*deltaCoeffs[facei]*faceAreas[facei];
        output[discretePIndex(own)] += normalCoefficient*multiplier;
        output[discretePIndex(nei)] -= normalCoefficient*multiplier;
    }

    forAll(mesh.boundary(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& patchDelta =
            mesh.deltaCoeffs().boundaryField()[patchi];
        const scalarField& patchArea =
            mesh.magSf().boundaryField()[patchi];
        const label patchStart = mesh.boundary()[patchi].start();
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();

        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            const scalar multiplier = faceAdjoint[patchStart + facei];
            mobilityGradientAdjoint[celli] += patchSf[facei]*multiplier;
            if (pressureFixed)
            {
                output[discretePIndex(celli)] +=
                    primalPressureMobility[celli]
                   *patchDelta[facei]*patchArea[facei]*multiplier;
            }
        }
    }

    vectorField pressureGradientAdjoint(discreteNCells, vector::zero);
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        pressureGradientAdjoint[celli] =
            primalPressureMobility[celli]
           *mobilityGradientAdjoint[celli]
           /Foam::max(mesh.V()[celli], SMALL);
    }

    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar faceMultiplier =
            discreteSf[facei]
          & (
                pressureGradientAdjoint[own]
              - pressureGradientAdjoint[nei]
            );
        output[discretePIndex(own)] +=
            discreteWeights[facei]*faceMultiplier;
        output[discretePIndex(nei)] +=
            (1.0 - discreteWeights[facei])*faceMultiplier;
    }

    forAll(mesh.boundary(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            output[discretePIndex(celli)] +=
                patchSf[facei] & pressureGradientAdjoint[celli];
        }
    }
};
scalarField discreteExternalFluxTranspose(discreteNUnknowns, 0.0);
forAll(discreteUpper, facei)
{
    const label own = discreteOwner[facei];
    const label nei = discreteNeighbour[facei];
    const vector faceTranspose =
        discreteSf[facei]*discreteExternalFaceFluxAdjoint[facei];
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        discreteExternalFluxTranspose[discreteUIndex(own, cmpt)] +=
            discreteWeights[facei]*faceTranspose.component(cmpt);
        discreteExternalFluxTranspose[discreteUIndex(nei, cmpt)] +=
            (1.0 - discreteWeights[facei])*faceTranspose.component(cmpt);
    }
}
forAll(mesh.boundary(), patchi)
{
    if (U.boundaryField()[patchi].fixesValue())
    {
        continue;
    }
    const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
    const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
    const label patchStart = mesh.boundary()[patchi].start();
    forAll(patchCells, facei)
    {
        const label celli = patchCells[facei];
        const vector faceTranspose =
            patchSf[facei]
           *discreteExternalFaceFluxAdjoint[patchStart + facei];
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            discreteExternalFluxTranspose[discreteUIndex(celli, cmpt)] +=
                faceTranspose.component(cmpt);
        }
    }
}
applyPressureFluxCorrectionTranspose
(
    discreteExternalFaceFluxAdjoint,
    discreteExternalFluxTranspose
);
for (label celli = 0; celli < discreteNCells; ++celli)
{
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        discreteFlowVelocityRhs[celli].component(cmpt) +=
            discreteExternalFluxTranspose[discreteUIndex(celli, cmpt)];
    }
    discreteFlowPressureRhs[celli] +=
        discreteExternalFluxTranspose[discretePIndex(celli)];
}
auto applyDiscreteFlowJ =
[&](const scalarField& input, scalarField& output)
{
    output = scalar(0.0);
    scalarField pressureFluxCorrection(mesh.nFaces(), 0.0);
    applyPressureFluxCorrection(input, pressureFluxCorrection);
    if (ransFlowModel)
    {
        volVectorField discreteDeltaU
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U
        );

        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            discreteDeltaU[celli] = vector
            (
                input[discreteUIndex(celli, 0)],
                input[discreteUIndex(celli, 1)],
                input[discreteUIndex(celli, 2)]
            );
        }
        forAll(mesh.boundary(), patchi)
        {
            if (U.boundaryField()[patchi].fixesValue())
            {
                discreteDeltaU.boundaryFieldRef()[patchi] == vector::zero;
            }
        }
        discreteDeltaU.correctBoundaryConditions();

        tmp<volVectorField> tDeviatoricResidual = -fvc::div
        (
            nuEffFrozen
           *dev2
            (
                Foam::T(fvc::grad(discreteDeltaU, "grad(U)"))
            ),
            "div((nuEffFrozen*dev2(T(grad(U)))))"
        );
        const vectorField& deviatoricResidual =
            tDeviatoricResidual().primitiveField();
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                output[discreteUIndex(celli, cmpt)] +=
                    mesh.V()[celli]
                   *deviatoricResidual[celli].component(cmpt);
            }
        }
    }

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        output[discreteUIndex(celli, 0)] +=
            discreteDiagX[celli]*input[discreteUIndex(celli, 0)];
        output[discreteUIndex(celli, 1)] +=
            discreteDiagY[celli]*input[discreteUIndex(celli, 1)];
        output[discreteUIndex(celli, 2)] +=
            discreteDiagZ[celli]*input[discreteUIndex(celli, 2)];
    }

    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];

        vector deltaUOwn
        (
            input[discreteUIndex(own, 0)],
            input[discreteUIndex(own, 1)],
            input[discreteUIndex(own, 2)]
        );
        vector deltaUNei
        (
            input[discreteUIndex(nei, 0)],
            input[discreteUIndex(nei, 1)],
            input[discreteUIndex(nei, 2)]
        );

        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(own, cmpt)] +=
                discreteUpper[facei]
               *input[discreteUIndex(nei, cmpt)];
            output[discreteUIndex(nei, cmpt)] +=
                discreteLower[facei]
               *input[discreteUIndex(own, cmpt)];
        }

        const vector deltaUFace =
            discreteWeights[facei]*deltaUOwn
          + (1.0 - discreteWeights[facei])*deltaUNei;
        const scalar deltaPhiFace =
            (discreteSf[facei] & deltaUFace)
          + pressureFluxCorrection[facei];
        const vector velocityJump = U[nei] - U[own];
        const label downwind = discretePhi[facei] >= 0.0 ? nei : own;

        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(downwind, cmpt)] +=
                deltaPhiFace*velocityJump.component(cmpt);
        }

        output[discretePIndex(own)] += deltaPhiFace;
        output[discretePIndex(nei)] -= deltaPhiFace;

        const scalar deltaPFace =
            discreteWeights[facei]*input[discretePIndex(own)]
          + (1.0 - discreteWeights[facei])
           *input[discretePIndex(nei)];
        const vector pressureContribution =
            discreteSf[facei]*deltaPFace;

        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(own, cmpt)] +=
                pressureContribution.component(cmpt);
            output[discreteUIndex(nei, cmpt)] -=
                pressureContribution.component(cmpt);
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        const label patchStart = mesh.boundary()[patchi].start();
        const bool velocityFixed = U.boundaryField()[patchi].fixesValue();
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();

        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];

            if (!velocityFixed)
            {
                vector deltaUCell
                (
                    input[discreteUIndex(celli, 0)],
                    input[discreteUIndex(celli, 1)],
                    input[discreteUIndex(celli, 2)]
                );
                output[discretePIndex(celli)] +=
                    (patchSf[facei] & deltaUCell)
                  + pressureFluxCorrection[patchStart + facei];
            }

            if (!pressureFixed)
            {
                const vector pressureContribution =
                    patchSf[facei]*input[discretePIndex(celli)];
                for (direction cmpt = 0; cmpt < 3; ++cmpt)
                {
                    output[discreteUIndex(celli, cmpt)] +=
                        pressureContribution.component(cmpt);
                }
            }
        }
    }
};

auto applyDiscreteFlowJT =
[&](const scalarField& input, scalarField& output)
{
    output = scalar(0.0);
    scalarField fluxAdjoint(mesh.nFaces(), 0.0);

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        output[discreteUIndex(celli, 0)] +=
            discreteDiagX[celli]*input[discreteUIndex(celli, 0)];
        output[discreteUIndex(celli, 1)] +=
            discreteDiagY[celli]*input[discreteUIndex(celli, 1)];
        output[discreteUIndex(celli, 2)] +=
            discreteDiagZ[celli]*input[discreteUIndex(celli, 2)];
    }

    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];

        vector lambdaOwn
        (
            input[discreteUIndex(own, 0)],
            input[discreteUIndex(own, 1)],
            input[discreteUIndex(own, 2)]
        );
        vector lambdaNei
        (
            input[discreteUIndex(nei, 0)],
            input[discreteUIndex(nei, 1)],
            input[discreteUIndex(nei, 2)]
        );

        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(own, cmpt)] +=
                discreteLower[facei]
               *input[discreteUIndex(nei, cmpt)];
            output[discreteUIndex(nei, cmpt)] +=
                discreteUpper[facei]
               *input[discreteUIndex(own, cmpt)];
        }

        const vector velocityJump = U[nei] - U[own];
        const scalar convectionMultiplier =
            discretePhi[facei] >= 0.0
          ? (lambdaNei & velocityJump)
          : (lambdaOwn & velocityJump);
        const vector fluxTranspose =
            discreteSf[facei]*convectionMultiplier;

        const scalar pressureDifference =
            input[discretePIndex(own)] - input[discretePIndex(nei)];
        const vector continuityTranspose =
            discreteSf[facei]*pressureDifference;
        fluxAdjoint[facei] =
            convectionMultiplier + pressureDifference;

        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(own, cmpt)] +=
                discreteWeights[facei]
               *(
                    fluxTranspose.component(cmpt)
                  + continuityTranspose.component(cmpt)
                );
            output[discreteUIndex(nei, cmpt)] +=
                (1.0 - discreteWeights[facei])
               *(
                    fluxTranspose.component(cmpt)
                  + continuityTranspose.component(cmpt)
                );
        }

        const scalar pressureTranspose =
            discreteSf[facei] & (lambdaOwn - lambdaNei);
        output[discretePIndex(own)] +=
            discreteWeights[facei]*pressureTranspose;
        output[discretePIndex(nei)] +=
            (1.0 - discreteWeights[facei])*pressureTranspose;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        const label patchStart = mesh.boundary()[patchi].start();
        const bool velocityFixed = U.boundaryField()[patchi].fixesValue();
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();

        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];

            if (!velocityFixed)
            {
                fluxAdjoint[patchStart + facei] =
                    input[discretePIndex(celli)];
                const vector continuityTranspose =
                    patchSf[facei]*input[discretePIndex(celli)];
                for (direction cmpt = 0; cmpt < 3; ++cmpt)
                {
                    output[discreteUIndex(celli, cmpt)] +=
                        continuityTranspose.component(cmpt);
                }
            }

            if (!pressureFixed)
            {
                vector lambdaCell
                (
                    input[discreteUIndex(celli, 0)],
                    input[discreteUIndex(celli, 1)],
                    input[discreteUIndex(celli, 2)]
                );
                output[discretePIndex(celli)] +=
                    patchSf[facei] & lambdaCell;
            }
        }
    }
    applyPressureFluxCorrectionTranspose(fluxAdjoint, output);
    applyFrozenDeviatoricJT(input, output);
};

auto discreteGlobalDot = [&](const scalarField& a, const scalarField& bField)
{
    scalar result = sum(a*bField);
    reduce(result, sumOp<scalar>());
    return result;
};

scalar discreteMaxDotError = 0.0;
for (label testI = 0; testI < 3; ++testI)
{
    scalarField delta(discreteNUnknowns, 0.0);
    scalarField lambda(discreteNUnknowns, 0.0);
    scalarField Jdelta(discreteNUnknowns, 0.0);
    scalarField JTlambda(discreteNUnknowns, 0.0);

    forAll(delta, dofi)
    {
        delta[dofi] = Foam::sin((0.013 + 0.002*testI)*(dofi + 1));
        lambda[dofi] = Foam::cos((0.017 + 0.003*testI)*(dofi + 1));
    }

    applyDiscreteFlowJ(delta, Jdelta);
    applyDiscreteFlowJT(lambda, JTlambda);
    const scalar left = discreteGlobalDot(lambda, Jdelta);
    const scalar right = discreteGlobalDot(delta, JTlambda);
    const scalar error =
        mag(left - right)
       /Foam::max(Foam::max(mag(left), mag(right)), scalar(1e-30));
    discreteMaxDotError = Foam::max(discreteMaxDotError, error);
}

Info<< "Reduced cold-flow operator transpose dot-test ("
    << discreteAdjointLabel << ") max relative error="
    << discreteMaxDotError << endl;

if (discreteMaxDotError > 1e-7)
{
    FatalErrorInFunction
        << "Reduced cold-flow transpose dot test failed for "
        << discreteAdjointLabel << " with error " << discreteMaxDotError
        << abort(FatalError);
}

scalar discreteMomentumScale =
    (
        sum(mag(discreteDiagX))
      + sum(mag(discreteDiagY))
      + sum(mag(discreteDiagZ))
    )/Foam::max(3.0*discreteNCells, scalar(1.0));
reduce(discreteMomentumScale, sumOp<scalar>());
discreteMomentumScale = Foam::max(discreteMomentumScale, scalar(SMALL));

scalar discreteAreaScale = sum(mesh.magSf().primitiveField());
label discreteFaceCount = mesh.nInternalFaces();
forAll(mesh.boundary(), patchi)
{
    discreteAreaScale += sum(mesh.magSf().boundaryField()[patchi]);
    discreteFaceCount += mesh.boundary()[patchi].size();
}
reduce(discreteAreaScale, sumOp<scalar>());
reduce(discreteFaceCount, sumOp<label>());
discreteAreaScale /= Foam::max(scalar(discreteFaceCount), scalar(1.0));
discreteAreaScale = Foam::max(discreteAreaScale, scalar(SMALL));

scalarField discreteRhs(discreteNUnknowns, 0.0);
for (label celli = 0; celli < discreteNCells; ++celli)
{
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        discreteRhs[discreteUIndex(celli, cmpt)] =
            discreteFlowVelocityRhs[celli].component(cmpt)
           /discreteMomentumScale;
    }
    discreteRhs[discretePIndex(celli)] =
        discreteFlowPressureRhs[celli]/discreteAreaScale;
}

auto applyScaledDiscreteFlowJT =
[&](const scalarField& input, scalarField& output)
{
    scalarField physicalInput(input);
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        physicalInput[discretePIndex(celli)] *=
            discreteMomentumScale/discreteAreaScale;
    }

    scalarField physicalOutput(discreteNUnknowns, 0.0);
    applyDiscreteFlowJT(physicalInput, physicalOutput);

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(celli, cmpt)] =
                physicalOutput[discreteUIndex(celli, cmpt)]
               /discreteMomentumScale;
        }
        output[discretePIndex(celli)] =
            physicalOutput[discretePIndex(celli)]/discreteAreaScale;
    }
};

const label discreteMaxIterations =
    optProperties.lookupOrDefault<label>("discreteFlowAdjointMaxIter", 2000);
const scalar discreteTolerance =
    optProperties.lookupOrDefault<scalar>("discreteFlowAdjointTolerance", 1e-9);

scalarField discreteSolution(discreteNUnknowns, 0.0);
scalarField discreteResidual(discreteNUnknowns, 0.0);
scalarField discreteOperatorSolution(discreteNUnknowns, 0.0);

const scalar discreteRhsNorm = Foam::sqrt
(
    Foam::max(discreteGlobalDot(discreteRhs, discreteRhs), scalar(1e-300))
);

scalarField discreteSchurDiag(discreteNCells, 1e-8);
forAll(discreteUpper, facei)
{
    const label own = discreteOwner[facei];
    const label nei = discreteNeighbour[facei];
    const scalar weight = discreteWeights[facei];
    const vector sfScaled = discreteSf[facei]/discreteAreaScale;
    const scalar ownerContribution = sqr(weight)
       *(
            sqr(sfScaled.x())
           /Foam::max
            (
                mag(discreteDiagX[own]/discreteMomentumScale),
                scalar(1e-12)
            )
          + sqr(sfScaled.y())
           /Foam::max
            (
                mag(discreteDiagY[own]/discreteMomentumScale),
                scalar(1e-12)
            )
          + sqr(sfScaled.z())
           /Foam::max
            (
                mag(discreteDiagZ[own]/discreteMomentumScale),
                scalar(1e-12)
            )
        );
    const scalar neighbourContribution = sqr(1.0 - weight)
       *(
            sqr(sfScaled.x())
           /Foam::max
            (
                mag(discreteDiagX[nei]/discreteMomentumScale),
                scalar(1e-12)
            )
          + sqr(sfScaled.y())
           /Foam::max
            (
                mag(discreteDiagY[nei]/discreteMomentumScale),
                scalar(1e-12)
            )
          + sqr(sfScaled.z())
           /Foam::max
            (
                mag(discreteDiagZ[nei]/discreteMomentumScale),
                scalar(1e-12)
            )
        );
    const scalar faceContribution =
        Foam::max(ownerContribution + neighbourContribution, scalar(1e-12));
    discreteSchurDiag[own] += faceContribution;
    discreteSchurDiag[nei] += faceContribution;
}

auto applyDiscreteBlockPreconditioner =
[&](const scalarField& input, scalarField& output)
{
    scalarField inverseDiagX(discreteNCells, 0.0);
    scalarField inverseDiagY(discreteNCells, 0.0);
    scalarField inverseDiagZ(discreteNCells, 0.0);

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        inverseDiagX[celli] = 1.0/Foam::max
        (
            mag(discreteDiagX[celli]/discreteMomentumScale),
            scalar(1e-12)
        );
        inverseDiagY[celli] = 1.0/Foam::max
        (
            mag(discreteDiagY[celli]/discreteMomentumScale),
            scalar(1e-12)
        );
        inverseDiagZ[celli] = 1.0/Foam::max
        (
            mag(discreteDiagZ[celli]/discreteMomentumScale),
            scalar(1e-12)
        );
        output[discreteUIndex(celli, 0)] =
            inverseDiagX[celli]*input[discreteUIndex(celli, 0)];
        output[discreteUIndex(celli, 1)] =
            inverseDiagY[celli]*input[discreteUIndex(celli, 1)];
        output[discreteUIndex(celli, 2)] =
            inverseDiagZ[celli]*input[discreteUIndex(celli, 2)];
    }

    scalarField pressureResidual(discreteNCells, 0.0);
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        vector tentativeOwn
        (
            output[discreteUIndex(own, 0)],
            output[discreteUIndex(own, 1)],
            output[discreteUIndex(own, 2)]
        );
        vector tentativeNei
        (
            output[discreteUIndex(nei, 0)],
            output[discreteUIndex(nei, 1)],
            output[discreteUIndex(nei, 2)]
        );
        const scalar pressureTranspose =
            (discreteSf[facei]/discreteAreaScale)
          & (tentativeOwn - tentativeNei);
        pressureResidual[own] +=
            discreteWeights[facei]*pressureTranspose;
        pressureResidual[nei] +=
            (1.0 - discreteWeights[facei])*pressureTranspose;
    }
    forAll(mesh.boundary(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            vector tentativeCell
            (
                output[discreteUIndex(celli, 0)],
                output[discreteUIndex(celli, 1)],
                output[discreteUIndex(celli, 2)]
            );
            pressureResidual[celli] +=
                (patchSf[facei]/discreteAreaScale) & tentativeCell;
        }
    }

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        output[discretePIndex(celli)] =
            (input[discretePIndex(celli)] - pressureResidual[celli])
           /discreteSchurDiag[celli];
    }

    vectorField pressureCoupling(discreteNCells, vector::zero);
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar pressureDifference =
            output[discretePIndex(own)]
          - output[discretePIndex(nei)];
        const vector faceCoupling =
            (discreteSf[facei]/discreteAreaScale)*pressureDifference;
        pressureCoupling[own] +=
            discreteWeights[facei]*faceCoupling;
        pressureCoupling[nei] +=
            (1.0 - discreteWeights[facei])*faceCoupling;
    }
    forAll(mesh.boundary(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            pressureCoupling[celli] +=
                (patchSf[facei]/discreteAreaScale)
               *output[discretePIndex(celli)];
        }
    }

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        output[discreteUIndex(celli, 0)] = inverseDiagX[celli]
           *(
                input[discreteUIndex(celli, 0)]
              - pressureCoupling[celli].x()
            );
        output[discreteUIndex(celli, 1)] = inverseDiagY[celli]
           *(
                input[discreteUIndex(celli, 1)]
              - pressureCoupling[celli].y()
            );
        output[discreteUIndex(celli, 2)] = inverseDiagZ[celli]
           *(
                input[discreteUIndex(celli, 2)]
              - pressureCoupling[celli].z()
            );
    }
};

const label discreteGmresRestart =
    optProperties.lookupOrDefault<label>("discreteFlowGmresRestart", 40);
if (discreteGmresRestart < 2)
{
    FatalErrorInFunction
        << "discreteFlowGmresRestart must be at least 2."
        << abort(FatalError);
}

std::vector<scalarField> discreteBasis
(
    discreteGmresRestart + 1,
    scalarField(discreteNUnknowns, 0.0)
);
std::vector<scalarField> discretePreconditionedBasis
(
    discreteGmresRestart,
    scalarField(discreteNUnknowns, 0.0)
);
std::vector<std::vector<scalar>> discreteHessenberg
(
    discreteGmresRestart + 1,
    std::vector<scalar>(discreteGmresRestart, 0.0)
);
std::vector<scalar> discreteCosine(discreteGmresRestart, 0.0);
std::vector<scalar> discreteSine(discreteGmresRestart, 0.0);
std::vector<scalar> discreteG(discreteGmresRestart + 1, 0.0);
std::vector<scalar> discreteY(discreteGmresRestart, 0.0);

label discreteIteration = 0;
scalar discreteRelativeResidual = GREAT;

while
(
    discreteIteration < discreteMaxIterations
 && discreteRelativeResidual > discreteTolerance
)
{
    applyScaledDiscreteFlowJT
    (
        discreteSolution,
        discreteOperatorSolution
    );
    discreteResidual = discreteRhs - discreteOperatorSolution;
    const scalar discreteResidualNorm = Foam::sqrt
    (
        discreteGlobalDot(discreteResidual, discreteResidual)
    );
    discreteRelativeResidual = discreteResidualNorm/discreteRhsNorm;
    if (discreteRelativeResidual <= discreteTolerance)
    {
        break;
    }

    discreteBasis[0] = discreteResidual/discreteResidualNorm;
    for (label rowi = 0; rowi <= discreteGmresRestart; ++rowi)
    {
        discreteG[rowi] = 0.0;
        for (label coli = 0; coli < discreteGmresRestart; ++coli)
        {
            discreteHessenberg[rowi][coli] = 0.0;
        }
    }
    discreteG[0] = discreteResidualNorm;

    label discreteInnerCount = 0;
    for
    (
        label innerI = 0;
        innerI < discreteGmresRestart
     && discreteIteration < discreteMaxIterations;
        ++innerI, ++discreteIteration
    )
    {
        discreteInnerCount = innerI + 1;
        applyDiscreteBlockPreconditioner
        (
            discreteBasis[innerI],
            discretePreconditionedBasis[innerI]
        );

        scalarField arnoldiVector(discreteNUnknowns, 0.0);
        applyScaledDiscreteFlowJT
        (
            discretePreconditionedBasis[innerI],
            arnoldiVector
        );

        for (label basisI = 0; basisI <= innerI; ++basisI)
        {
            discreteHessenberg[basisI][innerI] =
                discreteGlobalDot(arnoldiVector, discreteBasis[basisI]);
            arnoldiVector -=
                discreteHessenberg[basisI][innerI]
               *discreteBasis[basisI];
        }

        discreteHessenberg[innerI + 1][innerI] = Foam::sqrt
        (
            discreteGlobalDot(arnoldiVector, arnoldiVector)
        );
        if (discreteHessenberg[innerI + 1][innerI] > VSMALL)
        {
            discreteBasis[innerI + 1] =
                arnoldiVector/discreteHessenberg[innerI + 1][innerI];
        }

        for (label rotationI = 0; rotationI < innerI; ++rotationI)
        {
            const scalar upperEntry =
                discreteHessenberg[rotationI][innerI];
            const scalar lowerEntry =
                discreteHessenberg[rotationI + 1][innerI];
            discreteHessenberg[rotationI][innerI] =
                discreteCosine[rotationI]*upperEntry
              + discreteSine[rotationI]*lowerEntry;
            discreteHessenberg[rotationI + 1][innerI] =
               -discreteSine[rotationI]*upperEntry
              + discreteCosine[rotationI]*lowerEntry;
        }

        const scalar rotationNorm = Foam::sqrt
        (
            sqr(discreteHessenberg[innerI][innerI])
          + sqr(discreteHessenberg[innerI + 1][innerI])
        );
        if (rotationNorm <= VSMALL)
        {
            FatalErrorInFunction
                << "GMRES Arnoldi breakdown in " << discreteAdjointLabel
                << " at iteration " << discreteIteration
                << abort(FatalError);
        }

        discreteCosine[innerI] =
            discreteHessenberg[innerI][innerI]/rotationNorm;
        discreteSine[innerI] =
            discreteHessenberg[innerI + 1][innerI]/rotationNorm;
        discreteHessenberg[innerI][innerI] = rotationNorm;
        discreteHessenberg[innerI + 1][innerI] = 0.0;

        const scalar oldG = discreteG[innerI];
        discreteG[innerI] = discreteCosine[innerI]*oldG;
        discreteG[innerI + 1] = -discreteSine[innerI]*oldG;
        discreteRelativeResidual =
            mag(discreteG[innerI + 1])/discreteRhsNorm;

        if (discreteRelativeResidual <= discreteTolerance)
        {
            ++discreteIteration;
            break;
        }
    }

    for (label rowi = discreteInnerCount - 1; rowi >= 0; --rowi)
    {
        scalar rhsEntry = discreteG[rowi];
        for (label coli = rowi + 1; coli < discreteInnerCount; ++coli)
        {
            rhsEntry -=
                discreteHessenberg[rowi][coli]*discreteY[coli];
        }
        discreteY[rowi] =
            rhsEntry/discreteHessenberg[rowi][rowi];
    }

    for (label basisI = 0; basisI < discreteInnerCount; ++basisI)
    {
        discreteSolution +=
            discreteY[basisI]*discretePreconditionedBasis[basisI];
        discreteY[basisI] = 0.0;
    }
}

applyScaledDiscreteFlowJT(discreteSolution, discreteOperatorSolution);
discreteResidual = discreteRhs - discreteOperatorSolution;
discreteRelativeResidual = Foam::sqrt
(
    discreteGlobalDot(discreteResidual, discreteResidual)
)/discreteRhsNorm;

Info<< "Reduced discrete flow adjoint " << discreteAdjointLabel
    << ": GMRES iterations=" << discreteIteration
    << ", relative residual=" << discreteRelativeResidual << endl;

if (discreteRelativeResidual > discreteTolerance)
{
    FatalErrorInFunction
        << "Reduced discrete flow adjoint did not converge for "
        << discreteAdjointLabel << ": residual="
        << discreteRelativeResidual << ", tolerance=" << discreteTolerance
        << abort(FatalError);
}
for (label celli = 0; celli < discreteNCells; ++celli)
{
    discreteAdjointVelocity[celli] = vector
    (
        discreteSolution[discreteUIndex(celli, 0)],
        discreteSolution[discreteUIndex(celli, 1)],
        discreteSolution[discreteUIndex(celli, 2)]
    );
    discreteAdjointPressure[celli] =
        discreteSolution[discretePIndex(celli)]
       *discreteMomentumScale/discreteAreaScale;
}
discreteAdjointVelocity.correctBoundaryConditions();
discreteAdjointPressure.correctBoundaryConditions();