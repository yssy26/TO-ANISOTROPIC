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
{
    scalar sD = 0.0, mnD = GREAT, mxD = -GREAT;
    scalar sA = 0.0, mnA = GREAT, mxA = -GREAT;
    scalar sF = 0.0, mnF = GREAT, mxF = -GREAT;
    forAll(discreteDiagX, celli)
    {
        const scalar d = discreteDiagX[celli];
        const scalar a = alpha[celli];
        const scalar f = phi[celli];
        sD += d; mnD = Foam::min(mnD, d); mxD = Foam::max(mxD, d);
        sA += a; mnA = Foam::min(mnA, a); mxA = Foam::max(mxA, a);
        sF += f; mnF = Foam::min(mnF, f); mxF = Foam::max(mxF, f);
    }
    reduce(sD, sumOp<scalar>()); reduce(mnD, minOp<scalar>()); reduce(mxD, maxOp<scalar>());
    reduce(sA, sumOp<scalar>()); reduce(mnA, minOp<scalar>()); reduce(mxA, maxOp<scalar>());
    reduce(sF, sumOp<scalar>()); reduce(mnF, minOp<scalar>()); reduce(mxF, maxOp<scalar>());
    const scalar nC = Foam::max(scalar(mesh.nCells()), 1.0);
    Info<< "AdjDiag: diagX avg=" << sD/nC << " [min=" << mnD << ",max=" << mxD << "]"
        << " ; alpha avg=" << sA/nC << " [min=" << mnA << ",max=" << mxA << "]"
        << " ; phi avg=" << sF/nC << " [min=" << mnF << ",max=" << mxF << "]"
        << " (" << discreteAdjointLabel << ")" << endl;
}

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

// rAU = 1/A() (same basis as the adjoint assembly); shared by the forward
// (reduced deltaPhi) and transpose (hAdjoint) operators.
scalarField rAUAdj(discreteNCells, 0.0);
for (label celli = 0; celli < discreteNCells; ++celli)
{
    rAUAdj[celli] =
        1.0/Foam::max(discretePrimalMomentum.A().ref()[celli], SMALL);
}

// BFINAL-003 (P<-U relaxation patch): U-equation relaxation factor read from
// fvSolution relaxationFactors/equations/U exactly as fvMatrix::relax() does
// (guarded by relaxEquation, default 1.0 = no relaxation).  The production
// SIMPLE P-row flux tangent uses the RELAXED momentum mapping:
//   HbyA = rAU_rel*H_rel  with  rAU_rel = alphaRel*rAU_u  and the relax-source
//   derivative (D_rel-D0)*dU = (1/alphaRel-1)*A_u*dU entering H, which yields
//   dHbyA/dU = alphaRel*rAU_u*dH_u + (1-alphaRel)*dU.
// No hard-coded 0.4/0.6/2.5: the coefficients follow the actual relaxation
// factor (0.4 for the current case -> 0.4 mobility + 0.6 relax-source).
const scalar alphaRel =
    mesh.relaxEquation("U") ? mesh.equationRelaxationFactor("U") : scalar(1.0);
{
    scalar sU = 0.0, mnU = GREAT, mxU = -GREAT;
    scalar sL = 0.0, mnL = GREAT, mxL = -GREAT;
    scalar sUV = 0.0, mnUV = GREAT, mxUV = -GREAT;
    forAll(discreteUpper, fi)
    {
        const label own = discreteOwner[fi];
        const label nei = discreteNeighbour[fi];
        sU += discreteUpper[fi]; mnU = Foam::min(mnU, discreteUpper[fi]); mxU = Foam::max(mxU, discreteUpper[fi]);
        sL += discreteLower[fi]; mnL = Foam::min(mnL, discreteLower[fi]); mxL = Foam::max(mxL, discreteLower[fi]);
        sUV += discreteUpper[fi]/mesh.V()[own]; mnUV = Foam::min(mnUV, discreteUpper[fi]/mesh.V()[own]); mxUV = Foam::max(mxUV, discreteUpper[fi]/mesh.V()[own]);
    }
    const scalar nF = Foam::max(scalar(mesh.nInternalFaces()), 1.0);
    Info<< "AdjCoeff: upper[avg=" << sU/nF << " min=" << mnU << " max=" << mxU
        << "] lower[avg=" << sL/nF << " min=" << mnL << " max=" << mxL
        << "] upper/V[avg=" << sUV/nF << " min=" << mnUV << " max=" << mxUV
        << "]" << endl;
}

auto discreteUIndex = [&](const label celli, const direction cmpt)
{
    return 3*celli + cmpt;
};
auto discretePIndex = [&](const label celli)
{
    return discreteNVelocity + celli;
};

// Match fvMatrix::setReference(): a pressure gauge exists only when no
// pressure boundary fixes the value.  The current fixedValue-outlet case is
// nonsingular and must retain every physical pressure row.
const bool discretePressureNeedsReference = p.needReference();
const auto pinDiscretePressureReference = [&](scalarField& values)
{
    if (discretePressureNeedsReference)
    {
        values[discretePIndex(pRefCell)] = 0.0;
    }
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
    // SIMPLE (consistent=false) phi reconstruction:
    //   phi = phiHbyA - flux(pEqn), phiHbyA independent of p,
    //   so dphi/dp = -d(flux(pEqn))/dp = +mobility*snGrad(dp)*magSf.
    // Tangential mobility terms do NOT appear in the continuity row.
    faceFluxCorrection = scalar(0.0);
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
        faceFluxCorrection[facei] = normalPressureFlux;
    }

    forAll(mesh.boundary(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
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

// BFINAL-018 (fix 3): the previous pressure-flux-correction transpose
// helper (interior kf with the OPPOSITE sign to the operator's own J_PP
// block) was removed together with the old plain-interpolation U-row
// routing; the source folding below mirrors the operator's P-row flux
// tangent slot by slot instead.
// BFINAL-018 (fix 3, source folding): mirror of the production routing in
// solveDiscreteFlowAdjointProduction.H - the face functional
//   g_f = dJ/dphi_f - (dR_T/dphi_f)^T . Tb
// is routed through the transpose of the SAME relaxed-SIMPLE flux tangent
// the J operator's P rows implement (BFINAL-003):
//   dphi_f = Sf & ( alphaRel*( w*rAU_o*dH_o + (1-w)*rAU_n*dH_n )
//                 + (1-alphaRel)*( w*dU_o + (1-w)*dU_n ) )
//          - kf_f*( dp_n - dp_o )
// The old plain-interpolation U-routing plus the opposite-sign kf helper
// transposed a different, unphysical
// map and is retired (BFINAL-016 B1: 86%/70%/sign-flip vs flow truth).
scalarField discreteExternalFluxTranspose(discreteNUnknowns, 0.0);
vectorField discreteExternalHA(discreteNCells, vector::zero);
forAll(discreteUpper, facei)
{
    const label own = discreteOwner[facei];
    const label nei = discreteNeighbour[facei];
    const scalar faceFunctional =
        discreteExternalFaceFluxAdjoint[facei];
    const vector faceTranspose = discreteSf[facei]*faceFunctional;
    const scalar weight = discreteWeights[facei];
    discreteExternalHA[own] +=
        alphaRel*rAUAdj[own]*weight*faceTranspose;
    discreteExternalHA[nei] +=
        alphaRel*rAUAdj[nei]*(1.0 - weight)*faceTranspose;
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        discreteExternalFluxTranspose[discreteUIndex(own, cmpt)] +=
            (1.0 - alphaRel)*weight*faceTranspose.component(cmpt);
        discreteExternalFluxTranspose[discreteUIndex(nei, cmpt)] +=
            (1.0 - alphaRel)*(1.0 - weight)*faceTranspose.component(cmpt);
    }
    const scalar mobilityFace =
        weight*primalPressureMobility[own]
      + (1.0 - weight)*primalPressureMobility[nei];
    const scalar kf =
        mobilityFace
       *mesh.deltaCoeffs().primitiveField()[facei]
       *mesh.magSf().primitiveField()[facei];
    discreteExternalFluxTranspose[discretePIndex(own)] += kf*faceFunctional;
    discreteExternalFluxTranspose[discretePIndex(nei)] -= kf*faceFunctional;
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
    const bool pressureFixed = p.boundaryField()[patchi].fixesValue();
    forAll(patchCells, facei)
    {
        const label celli = patchCells[facei];
        const scalar faceFunctional =
            discreteExternalFaceFluxAdjoint[patchStart + facei];
        const vector faceTranspose = patchSf[facei]*faceFunctional;
        discreteExternalHA[celli] +=
            alphaRel*rAUAdj[celli]*faceTranspose;
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            discreteExternalFluxTranspose[discreteUIndex(celli, cmpt)] +=
                (1.0 - alphaRel)*faceTranspose.component(cmpt);
        }
        if (pressureFixed)
        {
            discreteExternalFluxTranspose[discretePIndex(celli)] +=
                primalPressureMobility[celli]
               *mesh.deltaCoeffs().boundaryField()[patchi][facei]
               *mesh.magSf().boundaryField()[patchi][facei]
               *faceFunctional;
        }
    }
}
// deltaH^T: transpose of the momentum-H tangent on the hA field (mirrors
// the J^T deltaH block of this file / the production operator).
forAll(discreteUpper, facei)
{
    const label own = discreteOwner[facei];
    const label nei = discreteNeighbour[facei];
    const scalar Vinvo = 1.0/Foam::max(mesh.V()[own], SMALL);
    const scalar Vinn = 1.0/Foam::max(mesh.V()[nei], SMALL);
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        discreteExternalFluxTranspose[discreteUIndex(nei, cmpt)] -=
            discreteUpper[facei]*Vinvo
           *discreteExternalHA[own].component(cmpt);
        discreteExternalFluxTranspose[discreteUIndex(own, cmpt)] -=
            discreteLower[facei]*Vinn
           *discreteExternalHA[nei].component(cmpt);
    }
}
forAll(mesh.boundary(), patchi)
{
    const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
    const vectorField internalCoeffs =
        discretePrimalMomentum.internalCoeffs()[patchi];
    forAll(patchCells, facei)
    {
        const label celli = patchCells[facei];
        const scalar cav =
            (
                internalCoeffs[facei].x()
              + internalCoeffs[facei].y()
              + internalCoeffs[facei].z()
            )/3.0;
        const scalar Vinvc = 1.0/Foam::max(mesh.V()[celli], SMALL);
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            discreteExternalFluxTranspose[discreteUIndex(celli, cmpt)] +=
                (
                    -internalCoeffs[facei].component(cmpt) + cav
                )*Vinvc*discreteExternalHA[celli].component(cmpt);
        }
    }
}
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
scalarField collectedDPhiJ(mesh.nInternalFaces(), 0.0);
bool collectDPhiJ = false;

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

    {
        scalar sR = 0.0, mnR = GREAT, mxR = -GREAT, sA = 0.0;
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            const scalar r = rAUAdj[celli];
            const scalar a = discretePrimalMomentum.A().ref()[celli];
            sR += r; mnR = Foam::min(mnR, r); mxR = Foam::max(mxR, r);
            sA += a;
        }
        reduce(sR, sumOp<scalar>()); reduce(mnR, minOp<scalar>()); reduce(mxR, maxOp<scalar>());
        reduce(sA, sumOp<scalar>());
        Info<< "AdjRAU: rAUAdj[avg=" << sR/discreteNCells
            << " min=" << mnR << " max=" << mxR
            << "] A()[avg=" << sA/discreteNCells << "]" << endl;
    }

    scalarField dHx(discreteNCells, 0.0);
    scalarField dHy(discreteNCells, 0.0);
    scalarField dHz(discreteNCells, 0.0);
    scalarField dHxOff(discreteNCells, 0.0);
    scalarField dHyOff(discreteNCells, 0.0);
    scalarField dHzOff(discreteNCells, 0.0);
    // NOTE: fvMatrix::H() divides by the cell volume at the end
    // (Hphi /= mesh.V()).  The off-diagonal action is accumulated here and
    // divided by V below to match the H() per-volume semantics.
    // FIX(3x): the previous outer "for (direction c = 0; c < 3; ++c)" was
    // redundant (inner statements are already indexed per component), which
    // accumulated every face contribution three times (Gate H4: 0.0100 vs
    // 0.00335 = 3.0x).  Removed.
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        dHx[own] -= discreteUpper[facei]*input[discreteUIndex(nei, 0)];
        dHy[own] -= discreteUpper[facei]*input[discreteUIndex(nei, 1)];
        dHz[own] -= discreteUpper[facei]*input[discreteUIndex(nei, 2)];
        dHx[nei] -= discreteLower[facei]*input[discreteUIndex(own, 0)];
        dHy[nei] -= discreteLower[facei]*input[discreteUIndex(own, 1)];
        dHz[nei] -= discreteLower[facei]*input[discreteUIndex(own, 2)];
        dHxOff[own] -= discreteUpper[facei]*input[discreteUIndex(nei, 0)];
        dHyOff[own] -= discreteUpper[facei]*input[discreteUIndex(nei, 1)];
        dHzOff[own] -= discreteUpper[facei]*input[discreteUIndex(nei, 2)];
        dHxOff[nei] -= discreteLower[facei]*input[discreteUIndex(own, 0)];
        dHyOff[nei] -= discreteLower[facei]*input[discreteUIndex(own, 1)];
        dHzOff[nei] -= discreteLower[facei]*input[discreteUIndex(own, 2)];
    }
    {
        scalar nO = 0.0;
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            nO += magSqr(vector(dHxOff[celli], dHyOff[celli], dHzOff[celli]));
        }
        scalar nIn = 0.0, mxIn = 0.0;
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                const scalar v = input[discreteUIndex(celli, c)];
                nIn += sqr(v);
                mxIn = Foam::max(mxIn, mag(v));
            }
        }
        Info<< "AdjDeltaH offdiag(noV): |dH_off|L2=" << Foam::sqrt(nO)
            << " inputL2=" << Foam::sqrt(nIn) << " inputMax=" << mxIn
            << " in0=" << input[discreteUIndex(0, 0)]
            << " in1=" << input[discreteUIndex(0, 1)]
            << " in2=" << input[discreteUIndex(0, 2)] << endl;
    }
    // boundary-diagonal tangent replicating fvMatrix::H():
    //   Bdiag_c = -internalCoeffs.component(c) + cmptAv(internalCoeffs)
    //   dH[cell].c += Bdiag_c * dU_internal[cell].c      (owner cell value)
    // boundarySource (boundaryCoeffs*psi_b) is a boundary-field-only
    // contribution; for non-coupled patches its frozen-coefficient
    // derivative w.r.t. U is zero, so bcV*dU is NOT added here.
    // No patch is skipped: H()'s boundaryDiag acts on the internal cell
    // value adjacent to fixedValue patches as well.
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const vectorField icV(discretePrimalMomentum.internalCoeffs()[patchi]);
        forAll(fc, f)
        {
            const label celli = fc[f];
            const scalar cav = (icV[f].x() + icV[f].y() + icV[f].z())/3.0;
            dHx[celli] += (-icV[f].x() + cav)*input[discreteUIndex(celli, 0)];
            dHy[celli] += (-icV[f].y() + cav)*input[discreteUIndex(celli, 1)];
            dHz[celli] += (-icV[f].z() + cav)*input[discreteUIndex(celli, 2)];
        }
    }

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        const scalar Vinv = 1.0/Foam::max(mesh.V()[celli], SMALL);
        dHx[celli] *= Vinv;
        dHy[celli] *= Vinv;
        dHz[celli] *= Vinv;
    }
    {
        scalar nH = 0.0, mxH = 0.0, nHa = 0.0;
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            const scalar h = mag(vector(dHx[celli], dHy[celli], dHz[celli]));
            nH += sqr(h);
            mxH = Foam::max(mxH, h);
            nHa += 1.0;
        }
        Info<< "AdjDeltaH: |dH|L2=" << Foam::sqrt(nH)
            << " avg=" << Foam::sqrt(nH/Foam::max(nHa, 1.0))
            << " max=" << mxH << endl;
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

        // ---- reduced one-pass phi tangent (frozen-phi SIMPLE path) ----
        // dH field (precomputed above): dH(own) -= upper*dU_nei etc.
        const scalar wf = discreteWeights[facei];
        const scalar rAUo = rAUAdj[own];
        const scalar rAUn = rAUAdj[nei];
        const scalar mobF =
            wf*primalPressureMobility[own]
          + (1.0 - wf)*primalPressureMobility[nei];
        const scalar dcfF = mesh.deltaCoeffs().primitiveField()[facei];
        const scalar mafF = mesh.magSf().primitiveField()[facei];
        const vector dHown(dHx[own], dHy[own], dHz[own]);
        const vector dHnei(dHx[nei], dHy[nei], dHz[nei]);
        const vector rAUdH = wf*rAUo*dHown + (1.0 - wf)*rAUn*dHnei;
        const scalar deltaPhiFace =
            (discreteSf[facei] & rAUdH)
          + mobF*(input[discretePIndex(nei)] - input[discretePIndex(own)])
           *dcfF*mafF;
        // BFINAL-003 (P<-U relaxation patch): the P-row flux tangent must use
        // the PRODUCTION relaxed-SIMPLE mapping (BFINAL-002 Candidate B):
        //   dphi/dU = Sf & (alphaRel*rAU_u*dH_u + (1-alphaRel)*dUbar)
        // where alphaRel is the actual fvSolution U-equation relaxation factor
        // and the (1-alphaRel)*dUbar term is the relaxation-source derivative
        // (D_rel-D0)*dU entering H.  The relaxation U-part is unchanged.
        // BFINAL-008 (J_PP actual-primal closure): the P-P part below now
        // carries the ACTUAL pressure-residual sign, dR_P/dp = -L_prod:
        //   - mobF*(p_nei-p_own)*dcfF*mafF   (exported +L_int -> -L_int),
        // with the fixedValue-p boundary internalCoeffs diagonal added in the
        // boundary loop below.  The OLD deltaPhiFace is kept ONLY for the
        // momentum velocityJump (convective) term below, which is part of
        // the J_UU block and must stay byte-for-byte unchanged.
        const vector interpDU = wf*deltaUOwn + (1.0 - wf)*deltaUNei;
        const scalar deltaPhiFacePU =
            (discreteSf[facei] & (alphaRel*rAUdH + (1.0 - alphaRel)*interpDU))
          - mobF*(input[discretePIndex(nei)] - input[discretePIndex(own)])
           *dcfF*mafF;
        if (collectDPhiJ)
        {
            collectedDPhiJ[facei] = deltaPhiFacePU;
        }
        const vector velocityJump = U[nei] - U[own];
        const label downwind = discretePhi[facei] >= 0.0 ? nei : own;

        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(downwind, cmpt)] +=
                deltaPhiFace*velocityJump.component(cmpt);
        }
        (void)velocityJump; (void)downwind;

        output[discretePIndex(own)] += deltaPhiFacePU;
        output[discretePIndex(nei)] -= deltaPhiFacePU;

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
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();
        // BFINAL-003 (P<-U relaxation patch): production constrainHbyA pins
        // HbyA_b = U_b on non-assignable (fixedValue/noSlip) U patches so the
        // boundary flux tangent is exactly 0 there.  On assignable
        // (zeroGradient) U patches the boundary face value is extrapolated
        // from the cell: dphi_b = Sf_b & (alphaRel*rAU*dH + (1-alphaRel)*dU_c).
        const bool uAssignable = U.boundaryField()[patchi].assignable();

        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];

            if (uAssignable)
            {
                // reduced boundary dphi (matches fvc::flux(HbyA) boundary,
                // extrapolated-calculated => face value = cell value):
                //   dphi_b = Sf_b & (alphaRel*rAU*dH + (1-alphaRel)*dU),
                // with dH already /V-scaled (the patch loop runs AFTER the
                // /V step on dHx/dHy/dHz).
                const scalar rAUb = rAUAdj[celli];
                const vector dHb
                (
                    dHx[celli],
                    dHy[celli],
                    dHz[celli]
                );
                const vector dUb
                (
                    input[discreteUIndex(celli, 0)],
                    input[discreteUIndex(celli, 1)],
                    input[discreteUIndex(celli, 2)]
                );
                output[discretePIndex(celli)] +=
                    (patchSf[facei]
                     & (alphaRel*rAUb*dHb + (1.0 - alphaRel)*dUb));
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

            // BFINAL-008 (J_PP actual-primal closure): fixedValue-p boundary
            // P-P diagonal from the production fvm::laplacian(rAtU,p)
            // addBoundaryDiag internalCoeffs, residual convention
            // dR_P/dp:  diag(cell) += -internalCoeffs_b
            //                    = +rAtU_c*deltaCoeffs_b*|Sf|_b.
            // boundaryCoeffs source is 0 (p_b = 0 on fixedValue p), so there
            // is NO source term; zeroGradient p patches contribute 0; no
            // cell-center identity row is inserted (the term is a
            // boundary-FACE flux linearization, not a Dirichlet row).
            if (pressureFixed)
            {
                output[discretePIndex(celli)] +=
                    primalPressureMobility[celli]
                   *mesh.deltaCoeffs().boundaryField()[patchi][facei]
                   *mesh.magSf().boundaryField()[patchi][facei]
                   *input[discretePIndex(celli)];
            }
        }
    }
};

auto applyDiscreteFlowJT =
[&](const scalarField& input, scalarField& output)
{
    output = scalar(0.0);

    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        output[discreteUIndex(celli, 0)] +=
            discreteDiagX[celli]*input[discreteUIndex(celli, 0)];
        output[discreteUIndex(celli, 1)] +=
            discreteDiagY[celli]*input[discreteUIndex(celli, 1)];
        output[discreteUIndex(celli, 2)] +=
            discreteDiagZ[celli]*input[discreteUIndex(celli, 2)];
    }

    // Stage 2 accumulator: phi-adjoint -> deltaH^T
    scalarField hAx(mesh.nCells(), 0.0);
    scalarField hAy(mesh.nCells(), 0.0);
    scalarField hAz(mesh.nCells(), 0.0);

    // Stage 1: internal faces
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
        const scalar lambdaPown = input[discretePIndex(own)];
        const scalar lambdaPnei = input[discretePIndex(nei)];

        // U-U offdiag transpose (A^T)
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(own, cmpt)] +=
                discreteLower[facei]
               *input[discreteUIndex(nei, cmpt)];
            output[discreteUIndex(nei, cmpt)] +=
                discreteUpper[facei]
               *input[discreteUIndex(own, cmpt)];
        }

        // direct momentum pressure-gradient transpose (U-P block), keep
        const scalar pressureTranspose =
            discreteSf[facei] & (lambdaOwn - lambdaNei);
        output[discretePIndex(own)] +=
            discreteWeights[facei]*pressureTranspose;
        output[discretePIndex(nei)] +=
            (1.0 - discreteWeights[facei])*pressureTranspose;

        // reduced-flux transpose (BFINAL-003 P<-U relaxation patch):
        //   forward: the P row is routed through the corrected P-U flux
        //     dphi_PU = Sf&(alphaRel*rAU*dH + (1-alphaRel)*dUbar) + kf*(p_n-p_o)
        //     while the velocityJump conv term keeps the OLD unrelaxed dphi,
        //   so the adjoint scalar splits:
        //     phiAdjointRel = convTerm + alphaRel*lambdaPdiff  (-> hA path)
        //     + direct (1-alphaRel)*lambdaPdiff*(wf*Sf)        (-> U rows)
        //   BFINAL-008 (J_PP actual-primal closure): the kf path is the exact
        //   transpose of the corrected forward P-P block (dR_P/dp = -L_int:
        //   diag +kf, off-diag -kf).  Only the lambdaPdiff coefficient flips
        //   sign (P(own) += kf*lambdaPdiff, P(nei) -= kf*lambdaPdiff); the
        //   convTerm (J_UP^T via the OLD velocityJump dphi) keeps its
        //   unchanged -kf/+kf routing; phiAdjointRel/hA/directRel (J_PU) are
        //   byte-identical.
        const vector velocityJump = U[nei] - U[own];
        const label downwind = discretePhi[facei] >= 0.0 ? nei : own;
        const vector lambdaDown = (downwind == nei) ? lambdaNei : lambdaOwn;
        const scalar convTerm = lambdaDown & velocityJump;
        const scalar lambdaPdiff = lambdaPown - lambdaPnei;
        const scalar phiAdjointPP = lambdaPdiff - convTerm;

        const scalar wf = discreteWeights[facei];
        const scalar mobF =
            wf*primalPressureMobility[own]
          + (1.0 - wf)*primalPressureMobility[nei];
        const scalar kf =
            mobF*mesh.deltaCoeffs().primitiveField()[facei]
               *mesh.magSf().primitiveField()[facei];
        output[discretePIndex(own)] += kf*phiAdjointPP;
        output[discretePIndex(nei)] -= kf*phiAdjointPP;

        const scalar phiAdjointRel = convTerm + alphaRel*lambdaPdiff;
        const vector hOwn =
            phiAdjointRel*wf*rAUAdj[own]*discreteSf[facei];
        const vector hNei =
            phiAdjointRel*(1.0 - wf)*rAUAdj[nei]*discreteSf[facei];
        hAx[own] += hOwn.x(); hAy[own] += hOwn.y(); hAz[own] += hOwn.z();
        hAx[nei] += hNei.x(); hAy[nei] += hNei.y(); hAz[nei] += hNei.z();

        // direct relax-source transpose: U(own/nei) <- (1-alphaRel)*lambdaPdiff*Sf
        const vector directRel =
            (1.0 - alphaRel)*lambdaPdiff*discreteSf[facei];
        output[discreteUIndex(own, 0)] += wf*directRel.x();
        output[discreteUIndex(own, 1)] += wf*directRel.y();
        output[discreteUIndex(own, 2)] += wf*directRel.z();
        output[discreteUIndex(nei, 0)] += (1.0 - wf)*directRel.x();
        output[discreteUIndex(nei, 1)] += (1.0 - wf)*directRel.y();
        output[discreteUIndex(nei, 2)] += (1.0 - wf)*directRel.z();
    }

    // boundary: reduced boundary-flux transpose (lambdaP -> hA) and
    // direct U-P transpose (keep)
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        // BFINAL-003 (P<-U relaxation patch): forward boundary dphi is 0 on
        // non-assignable (fixedValue/noSlip) U patches and
        // Sf_b&(alphaRel*rAU*dH + (1-alphaRel)*dU_cell) on assignable ones;
        // the transpose mirrors that: alphaRel-scaled hB + direct
        // (1-alphaRel)*lambdaPc*patchSf term, assignable-only.
        const bool uAssignable = U.boundaryField()[patchi].assignable();
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();

        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];

            if (uAssignable)
            {
                const scalar lambdaPc = input[discretePIndex(celli)];
                const vector hB =
                    lambdaPc*(alphaRel*rAUAdj[celli])*patchSf[facei];
                hAx[celli] += hB.x(); hAy[celli] += hB.y(); hAz[celli] += hB.z();

                // direct relax-source transpose on the boundary face:
                // U(cell) <- (1-alphaRel)*lambdaPc*patchSf
                const vector dB =
                    (1.0 - alphaRel)*lambdaPc*patchSf[facei];
                output[discreteUIndex(celli, 0)] += dB.x();
                output[discreteUIndex(celli, 1)] += dB.y();
                output[discreteUIndex(celli, 2)] += dB.z();
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

            // BFINAL-008 (J_PP actual-primal closure): the fixedValue-p
            // boundary P-P diagonal is symmetric (diagonal-only), so J^T
            // carries the identical +rAtU_c*deltaCoeffs_b*|Sf|_b entry,
            // matching the corrected forward boundary loop above.
            if (pressureFixed)
            {
                output[discretePIndex(celli)] +=
                    primalPressureMobility[celli]
                   *mesh.deltaCoeffs().boundaryField()[patchi][facei]
                   *mesh.magSf().boundaryField()[patchi][facei]
                   *input[discretePIndex(celli)];
            }
        }
    }

    // Stage 2: deltaH^T (hAdjoint -> U)
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar Vinvo = 1.0/Foam::max(mesh.V()[own], SMALL);
        const scalar Vinn = 1.0/Foam::max(mesh.V()[nei], SMALL);
        // dH_own depends on dU_nei via -upper/V_own ; dH_nei on dU_own
        // via -lower/V_nei (transpose of the /V-scaled off-diagonal).
        output[discreteUIndex(nei, 0)] -= discreteUpper[facei]*Vinvo*hAx[own];
        output[discreteUIndex(nei, 1)] -= discreteUpper[facei]*Vinvo*hAy[own];
        output[discreteUIndex(nei, 2)] -= discreteUpper[facei]*Vinvo*hAz[own];
        output[discreteUIndex(own, 0)] -= discreteLower[facei]*Vinn*hAx[nei];
        output[discreteUIndex(own, 1)] -= discreteLower[facei]*Vinn*hAy[nei];
        output[discreteUIndex(own, 2)] -= discreteLower[facei]*Vinn*hAz[nei];
    }
    // boundaryDiag transpose: Bdiag_c/V * hA[cell]
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const vectorField icV(discretePrimalMomentum.internalCoeffs()[patchi]);
        forAll(fc, f)
        {
            const label celli = fc[f];
            const scalar cav = (icV[f].x() + icV[f].y() + icV[f].z())/3.0;
            const scalar Vinvc = 1.0/Foam::max(mesh.V()[celli], SMALL);
            output[discreteUIndex(celli, 0)] +=
                (-icV[f].x() + cav)*Vinvc*hAx[celli];
            output[discreteUIndex(celli, 1)] +=
                (-icV[f].y() + cav)*Vinvc*hAy[celli];
            output[discreteUIndex(celli, 2)] +=
                (-icV[f].z() + cav)*Vinvc*hAz[celli];
        }
    }

    applyFrozenDeviatoricJT(input, output);
};

auto discreteGlobalDot = [&](const scalarField& a, const scalarField& bField)
{
    scalar result = sum(a*bField);
    reduce(result, sumOp<scalar>());
    return result;
};

scalar discreteMaxDotError = 0.0;
scalarField dotErrors(12, 0.0);
for (label testI = 0; testI < 12; ++testI)
{
    scalarField delta(discreteNUnknowns, 0.0);
    scalarField lambda(discreteNUnknowns, 0.0);
    scalarField Jdelta(discreteNUnknowns, 0.0);
    scalarField JTlambda(discreteNUnknowns, 0.0);

    forAll(delta, dofi)
    {
        delta[dofi] = Foam::sin((0.013 + 0.0017*testI)*(dofi + 1));
        lambda[dofi] = Foam::cos((0.017 + 0.0023*testI)*(dofi + 1));
    }

    applyDiscreteFlowJ(delta, Jdelta);
    applyDiscreteFlowJT(lambda, JTlambda);
    const scalar left = discreteGlobalDot(lambda, Jdelta);
    const scalar right = discreteGlobalDot(delta, JTlambda);
    const scalar error =
        mag(left - right)
       /Foam::max(Foam::max(mag(left), mag(right)), scalar(1e-30));
    discreteMaxDotError = Foam::max(discreteMaxDotError, error);
    dotErrors[testI] = error;
}

Info<< "Reduced cold-flow operator transpose dot-test ("
    << discreteAdjointLabel << ") max relative error="
    << discreteMaxDotError
    << " all=" << dotErrors << endl;

// Block-wise dot tests to localize any residual transpose mismatch
{
    const label nB = 12;
    scalar errPU = 0.0, errPP = 0.0, errUU = 0.0, errUP = 0.0, errB = 0.0;
    boolList bndC(mesh.nCells(), false);
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        forAll(fc, f) { bndC[fc[f]] = true; }
    }
    for (label k = 0; k < nB; ++k)
    {
        scalarField delta(discreteNUnknowns, 0.0);
        scalarField lambda(discreteNUnknowns, 0.0);
        scalarField Jd(discreteNUnknowns, 0.0);
        scalarField JtL(discreteNUnknowns, 0.0);
        auto errf = [&]()
        {
            applyDiscreteFlowJ(delta, Jd);
            applyDiscreteFlowJT(lambda, JtL);
            const scalar L = discreteGlobalDot(lambda, Jd);
            const scalar R = discreteGlobalDot(delta, JtL);
            return mag(L - R)/Foam::max(Foam::max(mag(L), mag(R)), scalar(1e-30));
        };
        // PU: delta U-only, lambda P-only
        {
            for (label celli = 0; celli < discreteNCells; ++celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    delta[discreteUIndex(celli, c)] =
                        Foam::sin((0.011 + 0.001*k)*(discreteUIndex(celli, c) + 1));
                }
                lambda[discretePIndex(celli)] =
                    Foam::cos((0.017 + 0.002*k)*(discretePIndex(celli) + 1));
            }
            errPU = Foam::max(errPU, errf());
        }
        // PP: delta p-only, lambda p-only
        {
            delta = scalar(0.0); lambda = scalar(0.0);
            for (label celli = 0; celli < discreteNCells; ++celli)
            {
                delta[discretePIndex(celli)] =
                    Foam::sin((0.013 + 0.001*k)*(discretePIndex(celli) + 1));
                lambda[discretePIndex(celli)] =
                    Foam::cos((0.019 + 0.002*k)*(discretePIndex(celli) + 1));
            }
            errPP = Foam::max(errPP, errf());
        }
        // UU: delta U-only, lambda U-only (offdiag + conv + diag)
        {
            delta = scalar(0.0); lambda = scalar(0.0);
            for (label celli = 0; celli < discreteNCells; ++celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    delta[discreteUIndex(celli, c)] =
                        Foam::sin((0.023 + 0.001*k)*(discreteUIndex(celli, c) + 1));
                    lambda[discreteUIndex(celli, c)] =
                        Foam::cos((0.029 + 0.002*k)*(discreteUIndex(celli, c) + 1));
                }
            }
            errUU = Foam::max(errUU, errf());
        }
        // UP: delta p-only, lambda U-only (direct grad-p + conv dphi(dp) path)
        {
            delta = scalar(0.0); lambda = scalar(0.0);
            for (label celli = 0; celli < discreteNCells; ++celli)
            {
                delta[discretePIndex(celli)] =
                    Foam::sin((0.041 + 0.001*k)*(discretePIndex(celli) + 1));
                for (direction c = 0; c < 3; ++c)
                {
                    lambda[discreteUIndex(celli, c)] =
                        Foam::cos((0.043 + 0.002*k)*(discreteUIndex(celli, c) + 1));
                }
            }
            errUP = Foam::max(errUP, errf());
        }
        // boundary-cell U-only: delta/lambda only on boundary cells
        {
            delta = scalar(0.0); lambda = scalar(0.0);
            for (label celli = 0; celli < discreteNCells; ++celli)
            {
                if (!bndC[celli]) continue;
                for (direction c = 0; c < 3; ++c)
                {
                    delta[discreteUIndex(celli, c)] =
                        Foam::sin((0.031 + 0.001*k)*(discreteUIndex(celli, c) + 1));
                    lambda[discreteUIndex(celli, c)] =
                        Foam::cos((0.037 + 0.002*k)*(discreteUIndex(celli, c) + 1));
                }
            }
            errB = Foam::max(errB, errf());
        }
    }
    Info<< "BlockDot (" << discreteAdjointLabel << "): PU=" << errPU
        << " PP=" << errPP << " UU=" << errUU << " UP=" << errUP
        << " boundaryU=" << errB << endl;
}

if (discreteMaxDotError > 1e-7)
{
    if (optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false))
    {
        WarningInFunction
            << "Reduced cold-flow transpose dot test EXCEEDED tolerance for "
            << discreteAdjointLabel << " (error=" << discreteMaxDotError
            << ") - probe diagnostics only, continuing" << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Reduced cold-flow transpose dot test failed for "
            << discreteAdjointLabel << " with error " << discreteMaxDotError
            << abort(FatalError);
    }
}

// ---------------------------------------------------------------------------
// Stage B4: flow-Jacobian consistency probe (gated by stageB4JacobianProbe).
// Verifies J (applyDiscreteFlowJ) against a finite-difference of the frozen
// momentum residual, component by component:
//   T1: momentum U-U part (dU != 0, dp = 0, dphi = fvc::flux(dU))
//       exercises the fvm matrix + deviatoric + upwind convection
//       non-linearity (velocityJump term).
// Does NOT modify any gradient formula or threshold.
// ---------------------------------------------------------------------------
if (optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false))
{
    Info<< "--- Stage B4: flow Jacobian consistency probe ---" << endl;

    // dPhi_J export: J's explicit flux linearization on the same direction
    // as the NS.H rebuilt probe (sin(0.271*(celli+1))), face-level.
    {
        scalarField stateJ(discreteNUnknowns, 0.0);
        forAll(U, celli)
        {
            const scalar s = Foam::sin(0.271*(celli + 1));
            stateJ[discreteUIndex(celli, 0)] = s;
            stateJ[discreteUIndex(celli, 1)] = 0.5*s;
            stateJ[discreteUIndex(celli, 2)] = -0.3*s;
        }
        scalarField outJ(discreteNUnknowns, 0.0);
        collectDPhiJ = true;
        applyDiscreteFlowJ(stateJ, outJ);
        collectDPhiJ = false;
        {
            scalar nDF = 0.0, mxDF = 0.0;
            forAll(collectedDPhiJ, fi)
            {
                nDF += sqr(collectedDPhiJ[fi]);
                mxDF = Foam::max(mxDF, mag(collectedDPhiJ[fi]));
            }
            Info<< "StageB4 dPhi_J stats: L2=" << Foam::sqrt(nDF)
                << " max=" << mxDF << " ; P-row outJ L2="
                << Foam::sqrt(Foam::sum(sqr(outJ))) << endl;
        }
        OFstream dJOut("stageB2_dphiDudU_J.mtx");
        dJOut.precision(16);
        forAll(collectedDPhiJ, fi)
        {
            dJOut<< collectedDPhiJ[fi] << nl;
        }
        // forward-J P-row output (div of the reduced deltaPhi)
        OFstream jPOut("stageB2_Jdq_P_forward.mtx");
        jPOut.precision(16);
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            jPOut<< outJ[discretePIndex(celli)] << nl;
        }
        Info<< "StageB4 dPhi_J + forward P-row exported (nFaces="
            << mesh.nInternalFaces() << ")" << endl;

        // ---- Gate H1: same-basis dH (reassembled FD) vs my analytic dHx ----
        {
            const scalar hG = 0.1;
            // dUG is the unit direction (NOT scaled by hG); the FD below
            // divides by 2*hG, so the analytic action on dUG is directly
            // comparable.  Previously dUG was pre-multiplied by hG, giving
            // a spurious factor-10 (1/hG) mismatch.
            volVectorField dUG(U);
            dUG.primitiveFieldRef() = vector::zero;
            forAll(dUG, celli)
            {
                dUG[celli] = vector
                (
                    Foam::sin(0.271*(celli + 1)),
                    0.5*Foam::sin(0.271*(celli + 1)),
                    -0.3*Foam::sin(0.271*(celli + 1))
                );
            }
            dUG.correctBoundaryConditions();

            vectorField Hp(mesh.nCells(), vector::zero);
            vectorField Hm(mesh.nCells(), vector::zero);
            for (label sgn4 = 0; sgn4 < 2; ++sgn4)
            {
                const scalar sg4 = (sgn4 == 0) ? 1.0 : -1.0;
                volVectorField UmodG(U);
                forAll(UmodG.primitiveFieldRef(), celli)
                {
                    UmodG.primitiveFieldRef()[celli] += sg4*hG*dUG[celli];
                }
                UmodG.correctBoundaryConditions();
                tmp<fvVectorMatrix> tEqG
                (
                    fvm::div(phi, UmodG)
                  - fvm::laplacian(nuEffFrozen, UmodG)
                  + fvm::Sp(alpha, UmodG)
                );
                const vectorField HG(tEqG.ref().H());
                if (sgn4 == 0) { Hp = HG; }
                else           { Hm = HG; }
            }
            // analytic dH on the same basis (offdiag + boundary, /V)
            scalarField dHax(mesh.nCells(), 0.0);
            scalarField dHay(mesh.nCells(), 0.0);
            scalarField dHaz(mesh.nCells(), 0.0);
            forAll(discreteUpper, fi)
            {
                const label own = discreteOwner[fi];
                const label nei = discreteNeighbour[fi];
                dHax[own] -= discreteUpper[fi]*dUG[nei].x();
                dHay[own] -= discreteUpper[fi]*dUG[nei].y();
                dHaz[own] -= discreteUpper[fi]*dUG[nei].z();
                dHax[nei] -= discreteLower[fi]*dUG[own].x();
                dHay[nei] -= discreteLower[fi]*dUG[own].y();
                dHaz[nei] -= discreteLower[fi]*dUG[own].z();
            }
            forAll(mesh.boundary(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                if (U.boundaryField()[patchi].fixesValue()) continue;
                const vectorField icV(discretePrimalMomentum.internalCoeffs()[patchi]);
                forAll(fc, f)
                {
                    dHax[fc[f]] -= icV[f].x()*dUG[fc[f]].x();
                    dHay[fc[f]] -= icV[f].y()*dUG[fc[f]].y();
                    dHaz[fc[f]] -= icV[f].z()*dUG[fc[f]].z();
                }
            }
            scalar nA4 = 0.0, nD4 = 0.0, dot4 = 0.0;
            forAll(dHax, celli)
            {
                const vector fd
                (
                    (Hp[celli].x() - Hm[celli].x())/(2.0*hG),
                    (Hp[celli].y() - Hm[celli].y())/(2.0*hG),
                    (Hp[celli].z() - Hm[celli].z())/(2.0*hG)
                );
                const vector an
                (
                    dHax[celli]/mesh.V()[celli],
                    dHay[celli]/mesh.V()[celli],
                    dHaz[celli]/mesh.V()[celli]
                );
                nA4 += magSqr(an);
                nD4 += magSqr(fd - an);
                dot4 += (fd & an);
            }
            Info<< "GateH1 same-basis dH: |an|L2=" << Foam::sqrt(nA4)
                << " |fd|L2=" << Foam::sqrt(nA4 + nD4)
                << " relL2=" << Foam::sqrt(nD4/Foam::max(nA4, SMALL))
                << " cos=" << dot4/Foam::max(Foam::sqrt(nA4*nA4), SMALL)
                << endl;

            // ---- 4-term H() decomposition on the same basis (full UEqn) ----
            {
                volVectorField UmodP(U);
                forAll(UmodP.primitiveFieldRef(), celli)
                {
                    UmodP.primitiveFieldRef()[celli] += hG*dUG[celli];
                }
                UmodP.correctBoundaryConditions();
                tmp<fvVectorMatrix> tEqP4
                (
                    fvm::div(phi, UmodP)
                  - fvm::laplacian(nuEffFrozen, UmodP)
                  + fvm::Sp(alpha, UmodP)
                );
                fvVectorMatrix& EqP4 = tEqP4.ref();
                const vectorField HFDP4(EqP4.H());
                vectorField dHoff(mesh.nCells(), vector::zero);
                vectorField dHbnd(mesh.nCells(), vector::zero);
                forAll(EqP4.upper(), fi)
                {
                    const label own = discreteOwner[fi];
                    const label nei = discreteNeighbour[fi];
                    for (direction c = 0; c < 3; ++c)
                    {
                        dHoff[own].component(c) -= EqP4.upper()[fi]*dUG[nei].component(c);
                        dHoff[nei].component(c) -= EqP4.lower()[fi]*dUG[own].component(c);
                    }
                }
                forAll(mesh.boundary(), patchi)
                {
                    const labelUList& fc = mesh.boundary()[patchi].faceCells();
                    const vectorField icV4(EqP4.internalCoeffs()[patchi]);
                    forAll(fc, f)
                    {
                        for (direction c = 0; c < 3; ++c)
                        {
                            dHbnd[fc[f]].component(c) -= icV4[f].component(c)*dUG[fc[f]].component(c);
                        }
                    }
                }
                scalar nOff = 0.0, nBnd = 0.0, nS = 0.0, nBS = 0.0;
                scalar nFD2 = 0.0;
                forAll(dHoff, celli)
                {
                    for (direction c = 0; c < 3; ++c)
                    {
                        const scalar fd = (Hp[celli].component(c) - Hm[celli].component(c))/(2.0*hG);
                        const scalar o = dHoff[celli].component(c)/mesh.V()[celli];
                        const scalar b = dHbnd[celli].component(c)/mesh.V()[celli];
                        nOff += sqr(o);
                        nBnd += sqr(b);
                        nS += sqr(fd - o - b);
                        nFD2 += sqr(fd);
                    }
                }
                Info<< "GateH1 4term: |FD|=" << Foam::sqrt(nFD2)
                    << " |offdiag/V|=" << Foam::sqrt(nOff)
                    << " |boundary/V|=" << Foam::sqrt(nBnd)
                    << " |residual(FD-off-bnd)|=" << Foam::sqrt(nS)
                    << " (source+boundarySource contribution)"
                    << endl;

                // ---- per-operator source() FD (same converged basis) ----
                {
                    vectorField sDP(mesh.nCells(), vector::zero);
                    vectorField sDM(mesh.nCells(), vector::zero);
                    vectorField sLP(mesh.nCells(), vector::zero);
                    vectorField sLM(mesh.nCells(), vector::zero);
                    vectorField sSP(mesh.nCells(), vector::zero);
                    vectorField sSM(mesh.nCells(), vector::zero);
                    for (label sgn5 = 0; sgn5 < 2; ++sgn5)
                    {
                        const scalar sg5 = (sgn5 == 0) ? 1.0 : -1.0;
                        volVectorField UmodS(U);
                        forAll(UmodS.primitiveFieldRef(), celli)
                        {
                            UmodS.primitiveFieldRef()[celli] += sg5*hG*dUG[celli];
                        }
                        UmodS.correctBoundaryConditions();
                        tmp<fvVectorMatrix> tDv(fvm::div(phi, UmodS));
                        tmp<fvVectorMatrix> tLp(fvm::laplacian(nuEffFrozen, UmodS));
                        tmp<fvVectorMatrix> tSp(fvm::Sp(alpha, UmodS));
                        const vectorField sD(tDv.ref().source());
                        const vectorField sL(tLp.ref().source());
                        const vectorField sS(tSp.ref().source());
                        if (sgn5 == 0) { sDP = sD; sLP = sL; sSP = sS; }
                        else           { sDM = sD; sLM = sL; sSM = sS; }
                    }
                    auto l2s5 = [&](const vectorField& a, const vectorField& b)
                    {
                        scalar s = 0.0;
                        forAll(a, i) { s += magSqr(a[i] - b[i]); }
                        return Foam::sqrt(s)/(2.0*hG);
                    };
                    Info<< "GateH1 sourceFD: div=" << l2s5(sDP, sDM)
                        << " laplacian=" << l2s5(sLP, sLM)
                        << " Sp=" << l2s5(sSP, sSM)
                        << " (total source-FD ~2.41e7)"
                        << endl;
                }
            }
        }

        // ---- Gate H2: single-cell impulse, per-face offdiag audit ----
        {
            tmp<fvVectorMatrix> tH2
            (
                fvm::div(phi, U)
              - fvm::laplacian(nuEffFrozen, U)
              + fvm::Sp(alpha, U)
            );
            fvVectorMatrix& EqH2 = tH2.ref();
            const scalarField& upG = EqH2.upper();
            const scalarField& loG = EqH2.lower();

            const label cell0 = 0;
            label nAdj = 0;
            forAll(discreteUpper, fi)
            {
                if (discreteOwner[fi] == cell0 || discreteNeighbour[fi] == cell0)
                {
                    ++nAdj;
                }
            }
            Info<< "H2 cell0=" << cell0 << " adjacent internal faces=" << nAdj << endl;

            for (direction cdir = 0; cdir < 3; ++cdir)
            {
                volVectorField dUImp(U);
                dUImp.primitiveFieldRef() = vector::zero;
                dUImp.primitiveFieldRef()[cell0].component(cdir) = 1.0;

                scalar nGate = 0.0, nDisc = 0.0, nDiff = 0.0, dot = 0.0;
                label nPr = 0;
                forAll(discreteUpper, fi)
                {
                    const label own = discreteOwner[fi];
                    const label nei = discreteNeighbour[fi];
                    if (own != cell0 && nei != cell0) continue;
                    const scalar duO = dUImp[own].component(cdir);
                    const scalar duN = dUImp[nei].component(cdir);
                    // H()-offdiag raw contributions (before /V):
                    // own row -= up*dU_nei ; nei row -= lo*dU_own
                    const scalar gOwn = -upG[fi]*duN;
                    const scalar gNei = -loG[fi]*duO;
                    const scalar dOwn = -discreteUpper[fi]*duN;
                    const scalar dNei = -discreteLower[fi]*duO;
                    if (nPr < 12)
                    {
                        Info<< "H2 c=" << label(cdir)
                            << " face=" << fi << " own=" << own << " nei=" << nei
                            << " upG=" << upG[fi] << " loG=" << loG[fi]
                            << " upD=" << discreteUpper[fi] << " loD=" << discreteLower[fi]
                            << " dUo=" << duO << " dUn=" << duN
                            << " gateOwn=" << gOwn << " gateNei=" << gNei
                            << " discOwn=" << dOwn << " discNei=" << dNei
                            << endl;
                        ++nPr;
                    }
                    nGate += sqr(gOwn) + sqr(gNei);
                    nDisc += sqr(dOwn) + sqr(dNei);
                    nDiff += sqr(gOwn - dOwn) + sqr(gNei - dNei);
                    dot += gOwn*dOwn + gNei*dNei;
                }
                Info<< "H2 c=" << label(cdir) << " nFaces=" << nAdj
                    << " |gate|=" << Foam::sqrt(nGate)
                    << " |disc|=" << Foam::sqrt(nDisc)
                    << " relL2=" << Foam::sqrt(nDiff/Foam::max(nDisc, SMALL))
                    << " cos=" << dot/Foam::max(Foam::sqrt(nDisc*nDisc), SMALL)
                    << endl;
            }

            // ---- Gate H3: explicit applyDiscreteFlowJ(stateJ) ----
            // same direction as Gate H1 dUG (sin 0.271), then compare the
            // internal offdiag print with Gate H1 (expect ~3.2e-3 noV).
            {
                scalarField stateJ3(discreteNUnknowns, 0.0);
                forAll(U, celli)
                {
                    const scalar s = Foam::sin(0.271*(celli + 1));
                    stateJ3[discreteUIndex(celli, 0)] = s;
                    stateJ3[discreteUIndex(celli, 1)] = 0.5*s;
                    stateJ3[discreteUIndex(celli, 2)] = -0.3*s;
                }
                scalarField outJ3(discreteNUnknowns, 0.0);
                applyDiscreteFlowJ(stateJ3, outJ3);
                {
                    scalar nP = 0.0, mxP = 0.0;
                    for (label celli = 0; celli < discreteNCells; ++celli)
                    {
                        const scalar v = stateJ3[discretePIndex(celli)];
                        nP += sqr(v);
                        mxP = Foam::max(mxP, mag(v));
                    }
                    Info<< "GateH3 stateJ3_p: max=" << mxP
                        << " L2=" << Foam::sqrt(nP) << endl;
                }
                Info<< "GateH3 applied applyDiscreteFlowJ(stateJ3); see AdjDeltaH"
                    << " offdiag(noV) above; outL2="
                    << Foam::sqrt(Foam::sum(sqr(outJ3))) << endl;

                // ---- Gate H4: full-field per-face offdiag comparison ----
                // same input (sin 0.271 field), three coefficient sources:
                //  (a) discreteUpper/Lower  (applyDiscreteFlowJ's own coeffs)
                //  (b) tH2 re-assembled      (Gate H1/H2 coeffs)
                // and compare their /V-before L2 against the applyDiscreteFlowJ
                // internal print (expect 0.0100397... for stateJ3).
                {
                    tmp<fvVectorMatrix> tH4
                    (
                        fvm::div(phi, U)
                      - fvm::laplacian(nuEffFrozen, U)
                      + fvm::Sp(alpha, U)
                    );
                    fvVectorMatrix& EqH4 = tH4.ref();
                    const scalarField& up4 = EqH4.upper();
                    const scalarField& lo4 = EqH4.lower();

                    scalarField oDx(mesh.nCells(), 0.0);
                    scalarField oDy(mesh.nCells(), 0.0);
                    scalarField oDz(mesh.nCells(), 0.0);
                    scalarField oRx(mesh.nCells(), 0.0);
                    scalarField oRy(mesh.nCells(), 0.0);
                    scalarField oRz(mesh.nCells(), 0.0);
                    scalar nDiff = 0.0, nDisc = 0.0, nRas = 0.0;
                    forAll(discreteUpper, fi)
                    {
                        const label own = discreteOwner[fi];
                        const label nei = discreteNeighbour[fi];
                        const scalar sN = Foam::sin(0.271*(nei + 1));
                        const scalar sO = Foam::sin(0.271*(own + 1));
                        // (a) discrete coeffs, same formula as applyDiscreteFlowJ
                        oDx[own] -= discreteUpper[fi]*sN;
                        oDy[own] -= discreteUpper[fi]*(0.5*sN);
                        oDz[own] -= discreteUpper[fi]*(-0.3*sN);
                        oDx[nei] -= discreteLower[fi]*sO;
                        oDy[nei] -= discreteLower[fi]*(0.5*sO);
                        oDz[nei] -= discreteLower[fi]*(-0.3*sO);
                        // (b) re-assembled coeffs
                        oRx[own] -= up4[fi]*sN;
                        oRy[own] -= up4[fi]*(0.5*sN);
                        oRz[own] -= up4[fi]*(-0.3*sN);
                        oRx[nei] -= lo4[fi]*sO;
                        oRy[nei] -= lo4[fi]*(0.5*sO);
                        oRz[nei] -= lo4[fi]*(-0.3*sO);
                        nDiff += sqr(discreteUpper[fi] - up4[fi])
                              + sqr(discreteLower[fi] - lo4[fi]);
                        nDisc += sqr(discreteUpper[fi]) + sqr(discreteLower[fi]);
                        nRas += sqr(up4[fi]) + sqr(lo4[fi]);
                    }
                    scalar nA = 0.0, nB = 0.0, nC = 0.0, dotAB = 0.0;
                    forAll(oDx, celli)
                    {
                        const vector a(oDx[celli], oDy[celli], oDz[celli]);
                        const vector b(oRx[celli], oRy[celli], oRz[celli]);
                        nA += magSqr(a);
                        nB += magSqr(b);
                        nC += magSqr(a - b);
                        dotAB += (a & b);
                    }
                    Info<< "GateH4 full-field: |off_discrete|L2=" << Foam::sqrt(nA)
                        << " |off_reassemble|L2=" << Foam::sqrt(nB)
                        << " relL2=" << Foam::sqrt(nC/Foam::max(nB, SMALL))
                        << " cos=" << dotAB/Foam::max(Foam::sqrt(nB*nB), SMALL)
                        << " ; coeffDiffL2=" << Foam::sqrt(nDiff)
                        << " (applyDiscreteFlowJ printed 0.0100397...)"
                        << endl;
                }

                // ---- Gate PR: same-basis reduced P-row FD (multi-h) ----
                // phi reduced path: flux( rAU * H ) with frozen phi (dA=0),
                // same basis as the adjoint; compare against the forward-J
                // P-row (outJ3 P components) for the same direction dUG.
                {
                    const scalar hs[4] = {1e-3, 3e-4, 1e-4, 3e-5};
                    volVectorField dirR(U);
                    dirR.primitiveFieldRef() = vector::zero;
                    forAll(dirR, celli)
                    {
                        const scalar s = Foam::sin(0.271*(celli + 1));
                        dirR[celli] = vector(s, 0.5*s, -0.3*s);
                    }
                    dirR.correctBoundaryConditions();
                    for (label hi2 = 0; hi2 < 4; ++hi2)
                    {
                        const scalar hh = hs[hi2];
                        scalarField rpP(mesh.nCells(), 0.0);
                        scalarField rpM(mesh.nCells(), 0.0);
                        surfaceScalarField phiFp(phi);
                        surfaceScalarField phiFm(phi);
                        vectorField HpR(mesh.nCells(), vector::zero);
                        vectorField HmR(mesh.nCells(), vector::zero);
                        volScalarField rAUc(1.0/discretePrimalMomentum.A().ref());
                        for (label sgn6 = 0; sgn6 < 2; ++sgn6)
                        {
                            const scalar sg6 = (sgn6 == 0) ? 1.0 : -1.0;
                            volVectorField UmodR(U);
                            forAll(UmodR.primitiveFieldRef(), celli)
                            {
                                UmodR.primitiveFieldRef()[celli] += sg6*hh*dirR[celli];
                            }
                            UmodR.correctBoundaryConditions();
                            tmp<fvVectorMatrix> tEqR
                            (
                                fvm::div(phi, UmodR)
                              - fvm::laplacian(nuEffFrozen, UmodR)
                              + fvm::Sp(alpha, UmodR)
                            );
                            const vectorField HR(tEqR.ref().H());
                            if (sgn6 == 0) { HpR = HR; }
                            else           { HmR = HR; }
                            volScalarField rAUf(1.0/tEqR.ref().A());
                            volVectorField HbyA(rAUf*tEqR.ref().H());
                            const surfaceScalarField phiF(fvc::flux(HbyA));
                            if (sgn6 == 0) { phiFp = phiF; }
                            else           { phiFm = phiF; }
                            scalarField rp(mesh.nCells(), 0.0);
                            forAll(mesh.owner(), fi)
                            {
                                rp[mesh.owner()[fi]] += phiF[fi];
                                rp[mesh.neighbour()[fi]] -= phiF[fi];
                            }
                            forAll(phiF.boundaryField(), patchi)
                            {
                                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                                const scalarField& phB = phiF.boundaryField()[patchi];
                                forAll(fc, f) { rp[fc[f]] += phB[f]; }
                            }
                            if (sgn6 == 0) { rpP = rp; }
                            else           { rpM = rp; }
                        }
                        // export the internal-face dphi FD (same direction dirR)
                        {
                            OFstream dFdOut("stageB2_dphiDudU_FD.mtx");
                            dFdOut.precision(16);
                            forAll(mesh.owner(), fi)
                            {
                                dFdOut<< (phiFp[fi] - phiFm[fi])/(2.0*hh) << nl;
                            }
                        }
                        // ---- FD-H oracle: dPhi_from_HFD vs dPhi_FD vs dPhi_J ----
                        {
                            scalarField dHxR(mesh.nCells(), 0.0);
                            scalarField dHyR(mesh.nCells(), 0.0);
                            scalarField dHzR(mesh.nCells(), 0.0);
                            forAll(dHxR, celli)
                            {
                                dHxR[celli] = (HpR[celli].x() - HmR[celli].x())/(2.0*hh);
                                dHyR[celli] = (HpR[celli].y() - HmR[celli].y())/(2.0*hh);
                                dHzR[celli] = (HpR[celli].z() - HmR[celli].z())/(2.0*hh);
                            }
                            scalar nA7 = 0.0, nB7 = 0.0, nC7 = 0.0, dot7 = 0.0;
                            scalar nAA = 0.0, nAB = 0.0, dotA = 0.0;
                            forAll(mesh.owner(), fi)
                            {
                                const label own = mesh.owner()[fi];
                                const label nei = mesh.neighbour()[fi];
                                const scalar wf7 = discreteWeights[fi];
                                const vector dHown(dHxR[own], dHyR[own], dHzR[own]);
                                const vector dHnei(dHxR[nei], dHyR[nei], dHzR[nei]);
                                const scalar dPhiH =
                                    (discreteSf[fi]
                                     & (wf7*rAUc[own]*dHown + (1.0 - wf7)*rAUc[nei]*dHnei));
                                const scalar dPhiFD =
                                    (phiFp[fi] - phiFm[fi])/(2.0*hh);
                                const scalar dPhiJ = collectedDPhiJ[fi];
                                nA7 += sqr(dPhiH);
                                nB7 += sqr(dPhiFD - dPhiH);
                                nC7 += sqr(dPhiJ - dPhiH);
                                dot7 += dPhiJ*dPhiH;
                                nAA += sqr(dPhiJ);
                                nAB += sqr(dPhiFD - dPhiJ);
                                dotA += dPhiFD*dPhiJ;
                            }
                            Info<< "GateOracle h=" << hh
                                << " A(dPhiHfd vs dPhiFD): relL2="
                                << Foam::sqrt(nB7/Foam::max(nA7, SMALL))
                                << " B(dPhiJ vs dPhiHfd): relL2="
                                << Foam::sqrt(nC7/Foam::max(nA7, SMALL))
                                << " cos=" << dot7/Foam::max(Foam::sqrt(nA7*nA7), SMALL)
                                << " ; full(dPhiJ vs dPhiFD): relL2="
                                << Foam::sqrt(nAB/Foam::max(nAA, SMALL))
                                << " cos=" << dotA/Foam::max(Foam::sqrt(nAA*nAA), SMALL)
                                << endl;
                        }

                        // ---- deltaH region split: J-analytic dHx vs dH_FD ----
                        // recompute the J-side dH (same formula as
                        // applyDiscreteFlowJ, on dirR) and compare against the
                        // reassembled dH_FD, split internal vs boundary cells.
                        {
                            scalarField dHxJ(mesh.nCells(), 0.0);
                            scalarField dHyJ(mesh.nCells(), 0.0);
                            scalarField dHzJ(mesh.nCells(), 0.0);
                            forAll(discreteUpper, fi)
                            {
                                const label own = discreteOwner[fi];
                                const label nei = discreteNeighbour[fi];
                                dHxJ[own] -= discreteUpper[fi]*dirR[nei].x();
                                dHyJ[own] -= discreteUpper[fi]*dirR[nei].y();
                                dHzJ[own] -= discreteUpper[fi]*dirR[nei].z();
                                dHxJ[nei] -= discreteLower[fi]*dirR[own].x();
                                dHyJ[nei] -= discreteLower[fi]*dirR[own].y();
                                dHzJ[nei] -= discreteLower[fi]*dirR[own].z();
                            }
                            forAll(mesh.boundary(), patchi)
                            {
                                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                                const vectorField icVJ(discretePrimalMomentum.internalCoeffs()[patchi]);
                                forAll(fc, f)
                                {
                                    const label celli = fc[f];
                                    const scalar cavJ =
                                        (icVJ[f].x() + icVJ[f].y() + icVJ[f].z())/3.0;
                                    dHxJ[celli] += (-icVJ[f].x() + cavJ)*dirR[celli].x();
                                    dHyJ[celli] += (-icVJ[f].y() + cavJ)*dirR[celli].y();
                                    dHzJ[celli] += (-icVJ[f].z() + cavJ)*dirR[celli].z();
                                }
                            }
                            boolList bndC(mesh.nCells(), false);
                            forAll(mesh.boundary(), patchi)
                            {
                                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                                forAll(fc, f) { bndC[fc[f]] = true; }
                            }
                            scalar nAi2 = 0.0, nDi2 = 0.0, doti2 = 0.0;
                            scalar nAb = 0.0, nDb = 0.0, dotb = 0.0;
                            label nIn2 = 0, nBd2 = 0;
                            forAll(dHxJ, celli)
                            {
                                const vector an
                                (
                                    dHxJ[celli]/mesh.V()[celli],
                                    dHyJ[celli]/mesh.V()[celli],
                                    dHzJ[celli]/mesh.V()[celli]
                                );
                                const vector fd
                                (
                                    (HpR[celli].x() - HmR[celli].x())/(2.0*hh),
                                    (HpR[celli].y() - HmR[celli].y())/(2.0*hh),
                                    (HpR[celli].z() - HmR[celli].z())/(2.0*hh)
                                );
                                if (bndC[celli])
                                {
                                    nAb += magSqr(an);
                                    nDb += magSqr(fd - an);
                                    dotb += (fd & an);
                                    ++nBd2;
                                }
                                else
                                {
                                    nAi2 += magSqr(an);
                                    nDi2 += magSqr(fd - an);
                                    doti2 += (fd & an);
                                    ++nIn2;
                                }
                            }
                            Info<< "GateDHsplit internal(" << nIn2 << "): relL2="
                                << Foam::sqrt(nDi2/Foam::max(nAi2, SMALL))
                                << " cos=" << doti2/Foam::max(Foam::sqrt(nAi2*nAi2), SMALL)
                                << " |an|=" << Foam::sqrt(nAi2)
                                << " boundary(" << nBd2 << "): relL2="
                                << Foam::sqrt(nDb/Foam::max(nAb, SMALL))
                                << " cos=" << dotb/Foam::max(Foam::sqrt(nAb*nAb), SMALL)
                                << " |an|=" << Foam::sqrt(nAb)
                                << " (J-dHx vs dH_FD)"
                                << endl;

                            // A/B: boundary variants (icV-only, bcV-only) vs FD
                            {
                                scalarField b1x(mesh.nCells(), 0.0);
                                scalarField b1y(mesh.nCells(), 0.0);
                                scalarField b1z(mesh.nCells(), 0.0);
                                scalarField b2x(mesh.nCells(), 0.0);
                                scalarField b2y(mesh.nCells(), 0.0);
                                scalarField b2z(mesh.nCells(), 0.0);
                                forAll(mesh.boundary(), patchi)
                                {
                                    const labelUList& fc = mesh.boundary()[patchi].faceCells();
                                    if (U.boundaryField()[patchi].fixesValue()) continue;
                                    const vectorField icVJ2(discretePrimalMomentum.internalCoeffs()[patchi]);
                                    const vectorField bcVJ2(discretePrimalMomentum.boundaryCoeffs()[patchi]);
                                    forAll(fc, f)
                                    {
                                        const label celli = fc[f];
                                        b1x[celli] -= icVJ2[f].x()*dirR[celli].x();
                                        b1y[celli] -= icVJ2[f].y()*dirR[celli].y();
                                        b1z[celli] -= icVJ2[f].z()*dirR[celli].z();
                                        b2x[celli] += bcVJ2[f].x()*dirR[celli].x();
                                        b2y[celli] += bcVJ2[f].y()*dirR[celli].y();
                                        b2z[celli] += bcVJ2[f].z()*dirR[celli].z();
                                    }
                                }
                                auto bndL2 = [&](const scalarField& ax, const scalarField& ay, const scalarField& az)
                                {
                                    scalar na = 0.0, nd = 0.0, dt = 0.0;
                                    forAll(ax, celli)
                                    {
                                        if (!bndC[celli]) continue;
                                        const vector an
                                        (
                                            ax[celli]/mesh.V()[celli],
                                            ay[celli]/mesh.V()[celli],
                                            az[celli]/mesh.V()[celli]
                                        );
                                        const vector fd
                                        (
                                            (HpR[celli].x() - HmR[celli].x())/(2.0*hh),
                                            (HpR[celli].y() - HmR[celli].y())/(2.0*hh),
                                            (HpR[celli].z() - HmR[celli].z())/(2.0*hh)
                                        );
                                        na += magSqr(an);
                                        nd += magSqr(fd - an);
                                        dt += (fd & an);
                                    }
                                    return Foam::sqrt(nd/Foam::max(na, SMALL));
                                };
                                Info<< "GateBndAB icV-only relL2="
                                    << bndL2(b1x, b1y, b1z)
                                    << " bcV-only relL2="
                                    << bndL2(b2x, b2y, b2z)
                                    << " (boundary cells, vs dH_FD; icV+bcV was "
                                    << "0.482)" << endl;
                            }
                        }
                        // internal-cell subset mask (cells not touching a patch)
                        {
                            boolList bndCell(mesh.nCells(), false);
                            forAll(mesh.boundary(), patchi)
                            {
                                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                                forAll(fc, f) { bndCell[fc[f]] = true; }
                            }
                            scalar nAi = 0.0, nDi = 0.0, doti = 0.0;
                            label nInt = 0;
                            for (label celli = 0; celli < discreteNCells; ++celli)
                            {
                                if (bndCell[celli]) continue;
                                const scalar fd = (rpP[celli] - rpM[celli])/(2.0*hh);
                                const scalar jp = outJ3[discretePIndex(celli)];
                                nAi += sqr(jp);
                                nDi += sqr(fd - jp);
                                doti += fd*jp;
                                ++nInt;
                            }
                            Info<< "GatePRint h=" << hh << " nInt=" << nInt
                                << " |J_P|=" << Foam::sqrt(nAi)
                                << " |FD_P|=" << Foam::sqrt(nAi + nDi)
                                << " relL2=" << Foam::sqrt(nDi/Foam::max(nAi, SMALL))
                                << " cos=" << doti/Foam::max(Foam::sqrt(nAi*nAi), SMALL)
                                << " (internal cells only)"
                                << endl;
                        }
                        scalar nA6 = 0.0, nD6 = 0.0, dot6 = 0.0;
                        for (label celli = 0; celli < discreteNCells; ++celli)
                        {
                            const scalar fd = (rpP[celli] - rpM[celli])/(2.0*hh);
                            const scalar jp = outJ3[discretePIndex(celli)];
                            nA6 += sqr(jp);
                            nD6 += sqr(fd - jp);
                            dot6 += fd*jp;
                        }
                        Info<< "GatePR h=" << hh
                            << " |J_P|=" << Foam::sqrt(nA6)
                            << " |FD_P|=" << Foam::sqrt(nA6 + nD6)
                            << " relL2=" << Foam::sqrt(nD6/Foam::max(nA6, SMALL))
                            << " cos=" << dot6/Foam::max(Foam::sqrt(nA6*nA6), SMALL)
                            << endl;
                    }
                }
            }
        }
    }

    const scalar epsB4 = 1e-6;

    volVectorField dUField(U);
    dUField.rename("dUField");
    dUField.primitiveFieldRef() = vector::zero;
    forAll(dUField, celli)
    {
        const scalar s = Foam::sin(0.037*(celli + 1));
        dUField[celli] = epsB4*vector(s, 0.5*s, -0.3*s);
    }
    forAll(mesh.boundary(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            dUField.boundaryFieldRef()[patchi] == vector::zero;
        }
    }
    dUField.correctBoundaryConditions();

    const volVectorField Ubase = U;
    const surfaceScalarField phibase = phi;

    // J applied to [dU; dp=0]
    scalarField state(discreteNUnknowns, 0.0);
    forAll(dUField, celli)
    {
        state[discreteUIndex(celli, 0)] = dUField[celli].x();
        state[discreteUIndex(celli, 1)] = dUField[celli].y();
        state[discreteUIndex(celli, 2)] = dUField[celli].z();
    }
    scalarField Jdelta(discreteNUnknowns, 0.0);
    applyDiscreteFlowJ(state, Jdelta);

    // frozen momentum residual (no grad(p); dp = 0 so its variation is 0)
    auto momentumResidualB4 =
    [&](const volVectorField& UU, const surfaceScalarField& phiff)
        -> vectorField
    {
        tmp<fvVectorMatrix> tM
        (
            fvm::div(phiff, UU, "div(phi,U)")
          - fvm::laplacian(nuEffFrozen, UU)
          + fvm::Sp(alpha, UU)
        );
        fvVectorMatrix& M = tM.ref();
        if (ransFlowModel)
        {
            M -= fvc::div
            (
                nuEffFrozen*dev2(Foam::T(fvc::grad(UU)))
            );
        }
        // (M & UU) = (D*UU - H_ldu*UU - source)/V  (operator& semantics;
        //  H_ldu = upper/lower contribution only, no diag, no /V)
        // A*UU = D*UU + H_ldu*UU  =>  R_U = A*UU - source
        //       = V*(M&UU) + 2*H_ldu(UU)
        vectorField hldu(mesh.nCells(), vector::zero);
        for (direction c = 0; c < 3; ++c)
        {
            const scalarField psiCmp(UU.primitiveField().component(c));
            const scalarField hCmp = M.lduMatrix::H(psiCmp);
            forAll(hCmp, celli)
            {
                hldu[celli].component(c) = hCmp[celli];
            }
        }
        return
            mesh.V()*(M & UU)().primitiveField()
          + 2.0*hldu;
    };

    const surfaceScalarField dphiField("dphiField", fvc::flux(dUField));
    // Perturbed field with the SAME boundary types as U: copy U, add the
    // internal perturbation, update zero-gradient boundary values from the
    // neighbouring interior (fixed-value boundaries stay at U values).
    volVectorField Upert(U);
    Upert.primitiveFieldRef() += dUField.primitiveField();
    forAll(mesh.boundary(), patchi)
    {
        if (!U.boundaryField()[patchi].fixesValue())
        {
            const labelUList& pc = mesh.boundary()[patchi].faceCells();
            fvPatchVectorField& pb = Upert.boundaryFieldRef()[patchi];
            forAll(pc, facei)
            {
                pb[facei] = Upert[pc[facei]];
            }
        }
    }
    const surfaceScalarField phiper("phiper", phibase + dphiField);
    const vectorField RU0 = momentumResidualB4(Ubase, phibase);
    const vectorField RU1 = momentumResidualB4(Upert, phiper);

    scalar normR0 = Foam::sqrt
    (
        Foam::max(sum(magSqr(RU0)), scalar(1e-300))
    );
    reduce(normR0, sumOp<scalar>());
    scalar normR1 = Foam::sqrt
    (
        Foam::max(sum(magSqr(RU1)), scalar(1e-300))
    );
    reduce(normR1, sumOp<scalar>());

    // component norms of R0 (V*(M&U), 2H, source, dev2)
    tmp<fvVectorMatrix> tMnorm
    (
        fvm::div(phibase, Ubase, "div(phi,U)")
      - fvm::laplacian(nuEffFrozen, Ubase)
      + fvm::Sp(alpha, Ubase)
    );
    fvVectorMatrix& Mnorm = tMnorm.ref();
    if (ransFlowModel)
    {
        Mnorm -= fvc::div
        (
            nuEffFrozen*dev2(Foam::T(fvc::grad(Ubase))),
            "div((nuEffFrozen*dev2(T(grad(U)))))"
        );
    }
    const vectorField VAMU = mesh.V()*(Mnorm & Ubase)().primitiveField();
    const vectorField Hval = Mnorm.H()().primitiveField();
    const vectorField Sval = Mnorm.source();
    const vectorField VAmuNoSrc =
        (mesh.V()*(Mnorm & Ubase)().primitiveField());
    scalar nVAMU = Foam::sqrt(Foam::max(sum(magSqr(VAMU)), scalar(1e-300)));
    scalar nH = Foam::sqrt(Foam::max(sum(magSqr(Hval)), scalar(1e-300)));
    scalar nS = Foam::sqrt(Foam::max(sum(magSqr(Sval)), scalar(1e-300)));
    reduce(nVAMU, sumOp<scalar>());
    reduce(nH, sumOp<scalar>());
    reduce(nS, sumOp<scalar>());
    Info<< "Stage B4 R0 components: |V*(M&U)|=" << nVAMU
        << " |2H|=" << 2.0*nH << " |source|=" << nS
        << " |diag*U|=" << Foam::sqrt
        (
            Foam::max(sum(magSqr(Mnorm.diag()*Ubase.primitiveField())), scalar(1e-300))
        ) << endl;
    scalar nUpperB4 = Foam::sqrt(Foam::max(sum(sqr(Mnorm.upper())), scalar(1e-300)));
    scalar nLowerB4 = Foam::sqrt(Foam::max(sum(sqr(Mnorm.lower())), scalar(1e-300)));
    scalar nDiagB4 = Foam::sqrt(Foam::max(sum(sqr(Mnorm.diag())), scalar(1e-300)));
    reduce(nUpperB4, sumOp<scalar>());
    reduce(nLowerB4, sumOp<scalar>());
    reduce(nDiagB4, sumOp<scalar>());
    Info<< "Stage B4 matrix norms: |upper|=" << nUpperB4
        << " |lower|=" << nLowerB4 << " |diag|=" << nDiagB4 << endl;
    Info<< "Stage B4 fields: max|phi|=" << gMax(mag(phi.primitiveField()))
        << " max|nuEffFrozen|=" << gMax(nuEffFrozen.primitiveField())
        << " max|alpha|=" << gMax(alpha.primitiveField())
        << " max|U|=" << gMax(mag(U.primitiveField())) << endl;
    scalar normJdelta = 0.0;
    forAll(RU0, celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            normJdelta += sqr(Jdelta[discreteUIndex(celli, c)]);
        }
    }
    reduce(normJdelta, sumOp<scalar>());
    Info<< "Stage B4 norms: |R(U0)|=" << normR0
        << " |R(U1)|=" << normR1
        << " |Jdelta_mom|=" << Foam::sqrt(normJdelta) << endl;

    scalar numT1 = 0.0, denT1 = 0.0, numT1Conv = 0.0, denT1Conv = 0.0;
    scalar maxCellErr = 0.0, maxCellDen = 0.0;
    forAll(RU0, celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            const scalar jval = Jdelta[discreteUIndex(celli, c)];
            const scalar fval = (RU1[celli].component(c)
                               - RU0[celli].component(c));
            numT1 += sqr(jval - fval);
            denT1 += sqr(fval);
            maxCellErr = Foam::max(maxCellErr, mag(jval - fval));
            maxCellDen = Foam::max(maxCellDen, mag(fval));
        }
    }
    reduce(numT1, sumOp<scalar>());
    reduce(denT1, sumOp<scalar>());
    reduce(maxCellErr, maxOp<scalar>());
    reduce(maxCellDen, maxOp<scalar>());
    const scalar relT1 = Foam::sqrt(numT1/Foam::max(denT1, SMALL));
    scalar dotT1 = 0.0, normJ = 0.0, normR = 0.0;
    forAll(RU0, celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            const scalar jval = Jdelta[discreteUIndex(celli, c)];
            const scalar fval = (RU1[celli].component(c)
                               - RU0[celli].component(c));
            dotT1 += jval*fval;
            normJ += sqr(jval);
            normR += sqr(fval);
        }
    }
    reduce(dotT1, sumOp<scalar>());
    reduce(normJ, sumOp<scalar>());
    reduce(normR, sumOp<scalar>());
    const scalar cosT1 = dotT1/Foam::max(Foam::sqrt(normJ*normR), SMALL);
    const scalar ratioT1 = Foam::sqrt(normJ/Foam::max(normR, SMALL));
    Info<< "Stage B4 T1 (momentum U-U, dp=0): relL2=" << relT1
        << " cos=" << cosT1 << " |J|/|FD|=" << ratioT1
        << " maxCellAbsErr=" << maxCellErr
        << " maxCellAbsFD=" << maxCellDen << endl;

    // T-cont: continuity row of J vs fvc::div(dphiField)
    scalar numCont = 0.0, denCont = 0.0;
    {
        const volScalarField dphiDiv = fvc::div(dphiField);
        forAll(dphiDiv, celli)
        {
            const scalar jp = Jdelta[discretePIndex(celli)];
            const scalar fp = mesh.V()[celli]*dphiDiv[celli];
            numCont += sqr(jp - fp);
            denCont += sqr(fp);
        }
    }
    reduce(numCont, sumOp<scalar>());
    reduce(denCont, sumOp<scalar>());
    const scalar relCont = Foam::sqrt(numCont/Foam::max(denCont, SMALL));
    Info<< "Stage B4 T-cont (continuity row): relL2 = " << relCont << endl;

    // T2: pressure coupling (dp != 0, dU = 0)
    // SIMPLE (consistent=false) phi reconstruction:
    //   phi = phiHbyA - flux(pEqn),  phiHbyA = flux(rAU*H(U)) independent of p
    //   => dphi/dp = + rAtU*snGrad(dp)*magSf   (rAtU = rAU here)
    // Compare the J continuity row (deltaPhiFace part, dU=0) with
    // div( rAtU*snGrad(dp)*magSf ).
    volScalarField dpField(p);
    dpField.rename("dpField");
    dpField.primitiveFieldRef() = scalar(0.0);
    forAll(dpField, celli)
    {
        dpField[celli] = epsB4*Foam::cos(0.043*(celli + 1));
    }
    forAll(mesh.boundary(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            dpField.boundaryFieldRef()[patchi] == scalar(0.0);
        }
    }
    dpField.correctBoundaryConditions();

    scalarField state2(discreteNUnknowns, 0.0);
    forAll(dpField, celli)
    {
        state2[discretePIndex(celli)] = dpField[celli];
    }
    scalarField J2(discreteNUnknowns, 0.0);
    applyDiscreteFlowJ(state2, J2);

    volScalarField mobilityField
    (
        IOobject("mobilityField", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("mobilityField", dimTime, scalar(0.0))
    );
    forAll(mobilityField, celli)
    {
        mobilityField[celli] = primalPressureMobility[celli];
    }
    mobilityField.correctBoundaryConditions();
    const surfaceScalarField mobFace("mobFace", fvc::interpolate(mobilityField));
    const surfaceScalarField dphiSIMPLE
    (
        "dphiSIMPLE",
        mobFace*fvc::snGrad(dpField)*mesh.magSf()
    );
    const volScalarField dphiSIMPLEDiv = fvc::div(dphiSIMPLE);

    scalar numT2 = 0.0, denT2 = 0.0, numT2mom = 0.0, denT2mom = 0.0;
    forAll(dphiSIMPLEDiv, celli)
    {
        const scalar jp = J2[discretePIndex(celli)];
        const scalar fp = mesh.V()[celli]*dphiSIMPLEDiv[celli];
        numT2 += sqr(jp - fp);
        denT2 += sqr(fp);
        for (direction c = 0; c < 3; ++c)
        {
            numT2mom += sqr(J2[discreteUIndex(celli, c)]);
        }
    }
    reduce(numT2, sumOp<scalar>());
    reduce(denT2, sumOp<scalar>());
    reduce(numT2mom, sumOp<scalar>());
    reduce(denT2mom, sumOp<scalar>());
    const scalar relT2 = Foam::sqrt(numT2/Foam::max(denT2, SMALL));
    Info<< "Stage B4 T2 (dp!=0 continuity): relL2(J_P, div(dphi_SIMPLE))="
        << relT2 << " |J_mom(dp)|=" << Foam::sqrt(numT2mom) << endl;

    // T-fvm: J diag/upper/lower part vs A(phi0)*dU (frozen-phi fvm matrix)
    scalarField Jfvm(discreteNUnknowns, 0.0);
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        Jfvm[discreteUIndex(celli, 0)] +=
            discreteDiagX[celli]*state[discreteUIndex(celli, 0)];
        Jfvm[discreteUIndex(celli, 1)] +=
            discreteDiagY[celli]*state[discreteUIndex(celli, 1)];
        Jfvm[discreteUIndex(celli, 2)] +=
            discreteDiagZ[celli]*state[discreteUIndex(celli, 2)];
    }
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        for (direction c = 0; c < 3; ++c)
        {
            Jfvm[discreteUIndex(own, c)] +=
                discreteUpper[facei]*state[discreteUIndex(nei, c)];
            Jfvm[discreteUIndex(nei, c)] +=
                discreteLower[facei]*state[discreteUIndex(own, c)];
        }
    }
    vectorField hlduFvm(mesh.nCells(), vector::zero);
    for (direction c = 0; c < 3; ++c)
    {
        const scalarField psiCmp(dUField.primitiveField().component(c));
        const scalarField hCmp = discretePrimalMomentum.lduMatrix::H(psiCmp);
        forAll(hCmp, celli)
        {
            hlduFvm[celli].component(c) = hCmp[celli];
        }
    }
    const vectorField Afd =
        mesh.V()*(discretePrimalMomentum & dUField)().primitiveField()
      + 2.0*hlduFvm;
    scalar numFvm = 0.0, denFvm = 0.0;
    forAll(Afd, celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            numFvm += sqr
            (
                Jfvm[discreteUIndex(celli, c)] - Afd[celli].component(c)
            );
            denFvm += sqr(Afd[celli].component(c));
        }
    }
    reduce(numFvm, sumOp<scalar>());
    reduce(denFvm, sumOp<scalar>());
    const scalar relFvm = Foam::sqrt(numFvm/Foam::max(denFvm, SMALL));
    Info<< "Stage B4 T-fvm (frozen-phi fvm matrix part): relL2 = " << relFvm
        << endl;

    // T-dev: J deviatoric part vs residual deviatoric variation
    // residual deviatoric = fvc::div(nuEff*dev2(T(grad(U))))
    const volVectorField Udev = Ubase + dUField;
    tmp<volVectorField> tDev0 = fvc::div
    (
        nuEffFrozen*dev2(Foam::T(fvc::grad(Ubase))),
        "div((nuEffFrozen*dev2(T(grad(U)))))"
    );
    tmp<volVectorField> tDev1 = fvc::div
    (
        nuEffFrozen*dev2(Foam::T(fvc::grad(Udev))),
        "div((nuEffFrozen*dev2(T(grad(U)))))"
    );
    // J dev part: recompute the deviatoric contribution of applyDiscreteFlowJ
    scalarField Jdev(discreteNUnknowns, 0.0);
    {
        tensorField stressA(discreteNCells, tensor::zero);
        tensorField gradA(discreteNCells, tensor::zero);
        forAll(discreteUpper, facei)
        {
            const label own = discreteOwner[facei];
            const label nei = discreteNeighbour[facei];
            const vector lOwn(state[discreteUIndex(own,0)], state[discreteUIndex(own,1)], state[discreteUIndex(own,2)]);
            const vector lNei(state[discreteUIndex(nei,0)], state[discreteUIndex(nei,1)], state[discreteUIndex(nei,2)]);
            const tensor faceS = discreteSf[facei]*(lNei - lOwn);
            stressA[own] += discreteWeights[facei]*faceS;
            stressA[nei] += (1.0 - discreteWeights[facei])*faceS;
        }
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            gradA[celli] = nuEffFrozen[celli]*dev2(Foam::T(stressA[celli]));
        }
        forAll(discreteUpper, facei)
        {
            const label own = discreteOwner[facei];
            const label nei = discreteNeighbour[facei];
            const vector faceV =
                discreteSf[facei]
              & (
                    gradA[own]/Foam::max(mesh.V()[own], SMALL)
                  - gradA[nei]/Foam::max(mesh.V()[nei], SMALL)
                );
            for (direction c = 0; c < 3; ++c)
            {
                Jdev[discreteUIndex(own, c)] += discreteWeights[facei]*faceV.component(c);
                Jdev[discreteUIndex(nei, c)] += (1.0 - discreteWeights[facei])*faceV.component(c);
            }
        }
    }
    scalar numDev = 0.0, denDev = 0.0;
    {
        const vectorField dRdev =
            tDev1().primitiveField() - tDev0().primitiveField();
        forAll(dRdev, celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                numDev += sqr
                (
                    Jdev[discreteUIndex(celli, c)]
                  - (-mesh.V()[celli]*dRdev[celli].component(c))
                );
                denDev += sqr(mesh.V()[celli]*dRdev[celli].component(c));
            }
        }
    }
    reduce(numDev, sumOp<scalar>());
    reduce(denDev, sumOp<scalar>());
    const scalar relDev = Foam::sqrt(numDev/Foam::max(denDev, SMALL));
    Info<< "Stage B4 T-dev (deviatoric part): relL2 = " << relDev << endl;

    // T-conv: J convection non-linearity (velocityJump term) vs residual
    // phi-variation: [div(phi',U') - div(phi,U')] - [div(phi',U) - div(phi,U)]
    // = -fvc::div(dphi, U') + fvc::div(dphi, U)  (upwind, frozen phi change)
    // Simplified one-sided check: compare J conv term against
    // -[R(phi',U) - R(phi,U)] restricted to the convective operator.
    tmp<fvVectorMatrix> tM0(fvm::div(phibase, Ubase, "div(phi,U)"));
    tmp<fvVectorMatrix> tM1(fvm::div(phiper, Ubase, "div(phi,U)"));
    auto resWithHldu = [&](fvVectorMatrix& MM, const volVectorField& UU)
    {
        vectorField hldu(mesh.nCells(), vector::zero);
        for (direction c = 0; c < 3; ++c)
        {
            const scalarField psiCmp(UU.primitiveField().component(c));
            const scalarField hCmp = MM.lduMatrix::H(psiCmp);
            forAll(hCmp, celli)
            {
                hldu[celli].component(c) = hCmp[celli];
            }
        }
        return mesh.V()*(MM & UU)().primitiveField() + 2.0*hldu;
    };
    const vectorField convDelta =
        resWithHldu(tM1.ref(), Ubase) - resWithHldu(tM0.ref(), Ubase);
    scalar numConv = 0.0, denConv = 0.0;
    forAll(convDelta, celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            const scalar jConvCell = Jdelta[discreteUIndex(celli, c)]
                - Jfvm[discreteUIndex(celli, c)]
                - Jdev[discreteUIndex(celli, c)];
            numConv += sqr(jConvCell - convDelta[celli].component(c));
            denConv += sqr(convDelta[celli].component(c));
        }
    }
    reduce(numConv, sumOp<scalar>());
    reduce(denConv, sumOp<scalar>());
    const scalar relConv = Foam::sqrt(numConv/Foam::max(denConv, SMALL));
    Info<< "Stage B4 T-conv (convection phi-variation part): relL2 = "
        << relConv << endl;
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
pinDiscretePressureReference(discreteRhs);

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
    if (discretePressureNeedsReference)
    {
        output[discretePIndex(pRefCell)] = input[discretePIndex(pRefCell)];
    }
};

// =========================================================================
// Explicit assembly of J^T (physical coefficients, COO -> CSR) as an
// independent oracle.  Every term replicates applyDiscreteFlowJT() exactly.
// =========================================================================
std::vector<label> cooRow, cooCol;
std::vector<scalar> cooVal;
cooRow.reserve(discreteNUnknowns*30);
cooCol.reserve(discreteNUnknowns*30);
cooVal.reserve(discreteNUnknowns*30);

auto addCoo = [&](label r, label c, scalar v)
{
    if (mag(v) < 1e-300)
    {
        return;
    }
    cooRow.push_back(r);
    cooCol.push_back(c);
    cooVal.push_back(v);
};

// ---- U-U diagonal (fvm matrix, physical) ----
for (label celli = 0; celli < discreteNCells; ++celli)
{
    addCoo(discreteUIndex(celli, 0), discreteUIndex(celli, 0), discreteDiagX[celli]);
    addCoo(discreteUIndex(celli, 1), discreteUIndex(celli, 1), discreteDiagY[celli]);
    addCoo(discreteUIndex(celli, 2), discreteUIndex(celli, 2), discreteDiagZ[celli]);
}

const scalarField& cooDeltaCoeffs = mesh.deltaCoeffs().primitiveField();
const scalarField& cooFaceAreas = mesh.magSf().primitiveField();

// ---- internal faces (two-stage reduced-J^T structure) ----
forAll(discreteUpper, facei)
{
    const label own = discreteOwner[facei];
    const label nei = discreteNeighbour[facei];
    const scalar weight = discreteWeights[facei];
    const vector sf = discreteSf[facei];
    const label downwind = discretePhi[facei] >= 0.0 ? nei : own;
    const vector velocityJump = U[nei] - U[own];
    const scalar mobF =
        weight*primalPressureMobility[own]
      + (1.0 - weight)*primalPressureMobility[nei];
    const scalar kf = mobF*cooDeltaCoeffs[facei]*cooFaceAreas[facei];
    const scalar rAUo = rAUAdj[own];
    const scalar rAUn = rAUAdj[nei];
    const scalar Vinvo = 1.0/Foam::max(mesh.V()[own], SMALL);
    const scalar Vinn = 1.0/Foam::max(mesh.V()[nei], SMALL);

    // U-U offdiag transpose (A^T)
    for (direction c = 0; c < 3; ++c)
    {
        addCoo(discreteUIndex(own, c), discreteUIndex(nei, c), discreteLower[facei]);
        addCoo(discreteUIndex(nei, c), discreteUIndex(own, c), discreteUpper[facei]);
    }

    // direct pressure transpose (P row <- U): keep
    for (direction c = 0; c < 3; ++c)
    {
        addCoo(discretePIndex(own), discreteUIndex(own, c),  weight*sf.component(c));
        addCoo(discretePIndex(own), discreteUIndex(nei, c), -weight*sf.component(c));
        addCoo(discretePIndex(nei), discreteUIndex(own, c),  (1.0 - weight)*sf.component(c));
        addCoo(discretePIndex(nei), discreteUIndex(nei, c), -(1.0 - weight)*sf.component(c));
    }

    // kf path (BFINAL-008 J_PP actual-primal closure): the corrected P-P
    // block is dR_P/dp = -L_int, i.e. diag +kf, off-diag -kf (the exact
    // transpose of the corrected matrix-free P-P: P(own) += kf*lambdaPdiff,
    // P(nei) -= kf*lambdaPdiff).  The convTerm entries below (velocityJump
    // via the OLD dphi, J_UP^T) keep the unchanged -kf/+kf routing.
    addCoo(discretePIndex(own), discretePIndex(own),  kf);
    addCoo(discretePIndex(own), discretePIndex(nei), -kf);
    addCoo(discretePIndex(nei), discretePIndex(own), -kf);
    addCoo(discretePIndex(nei), discretePIndex(nei),  kf);
    for (direction c2 = 0; c2 < 3; ++c2)
    {
        addCoo(discretePIndex(own), discreteUIndex(downwind, c2), -kf*velocityJump.component(c2));
        addCoo(discretePIndex(nei), discreteUIndex(downwind, c2),  kf*velocityJump.component(c2));
    }

    // BFINAL-003 (P<-U relaxation patch): direct relax-source transpose
    //   U(own/nei) <- (1-alphaRel)*lambdaPdiff*(weight*Sf)  (internal faces),
    // mirroring the matrix-free applyDiscreteFlowJT direct term.
    for (direction c = 0; c < 3; ++c)
    {
        addCoo(discreteUIndex(own, c), discretePIndex(own),  (1.0 - alphaRel)*weight*sf.component(c));
        addCoo(discreteUIndex(own, c), discretePIndex(nei), -(1.0 - alphaRel)*weight*sf.component(c));
        addCoo(discreteUIndex(nei, c), discretePIndex(own),  (1.0 - alphaRel)*(1.0 - weight)*sf.component(c));
        addCoo(discreteUIndex(nei, c), discretePIndex(nei), -(1.0 - alphaRel)*(1.0 - weight)*sf.component(c));
    }
}

// hAdjoint path (Stage 1 -> Stage 2) expanded over cross-face sums:
//   hA[cell] accumulates over ALL faces of the cell, then
//   Stage2: U(nei_f) -= upper_f/V_own_f*hA[own_f] ; U(own_f) -= lower_f/V_nei_f*hA[nei_f]
{
    std::vector<std::vector<label>> cellFacesOf(mesh.nCells());
    forAll(discreteUpper, fi)
    {
        cellFacesOf[discreteOwner[fi]].push_back(fi);
        cellFacesOf[discreteNeighbour[fi]].push_back(fi);
    }
    // hB coefficient per cell: hA[cell] += lambdaP(cell)*rAU*Sf_b summed over
    // all boundary faces of the cell (Stage-1 boundary path).  BFINAL-003
    // (P<-U relaxation patch): assignable-only (constrainHbyA pins the
    // non-assignable patches to dphi_b = 0) and alphaRel-scaled (relaxed
    // mobility), mirroring the matrix-free applyDiscreteFlowJT boundary path.
    scalarField hBx(mesh.nCells(), 0.0);
    scalarField hBy(mesh.nCells(), 0.0);
    scalarField hBz(mesh.nCells(), 0.0);
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const bool uAssignable = U.boundaryField()[patchi].assignable();
        forAll(fc, f)
        {
            const label celli = fc[f];
            if (uAssignable)
            {
                hBx[celli] += alphaRel*rAUAdj[celli]*pSf[f].x();
                hBy[celli] += alphaRel*rAUAdj[celli]*pSf[f].y();
                hBz[celli] += alphaRel*rAUAdj[celli]*pSf[f].z();
            }
        }
    }
    forAll(discreteUpper, fi)
    {
        const label own = discreteOwner[fi];
        const label nei = discreteNeighbour[fi];
        const scalar Vinvo = 1.0/Foam::max(mesh.V()[own], SMALL);
        const scalar Vinn = 1.0/Foam::max(mesh.V()[nei], SMALL);
        const scalar upperF = discreteUpper[fi];
        const scalar lowerF = discreteLower[fi];

        // hB part of hA[own]/hA[nei] fed through Stage 2:
        //   U(nei,c) -= upper_f/V_own * hBcoef[own] * lambdaP(own)
        //   U(own,c) -= lower_f/V_nei * hBcoef[nei] * lambdaP(nei)
        addCoo(discreteUIndex(nei, 0), discretePIndex(own), -upperF*Vinvo*hBx[own]);
        addCoo(discreteUIndex(nei, 1), discretePIndex(own), -upperF*Vinvo*hBy[own]);
        addCoo(discreteUIndex(nei, 2), discretePIndex(own), -upperF*Vinvo*hBz[own]);
        addCoo(discreteUIndex(own, 0), discretePIndex(nei), -lowerF*Vinn*hBx[nei]);
        addCoo(discreteUIndex(own, 1), discretePIndex(nei), -lowerF*Vinn*hBy[nei]);
        addCoo(discreteUIndex(own, 2), discretePIndex(nei), -lowerF*Vinn*hBz[nei]);

        // U(nei,c) -= upper_f/V_own * hA[own]
        for (const label fp : cellFacesOf[own])
        {
            const label o2 = discreteOwner[fp];
            const label n2 = discreteNeighbour[fp];
            const scalar w2 = discreteWeights[fp];
            const vector sf2 = discreteSf[fp];
            const scalar rAUo2 = rAUAdj[own];
            const label down2 = discretePhi[fp] >= 0.0 ? n2 : o2;
            const vector vj2 = U[n2] - U[o2];
            // hA[cell] lambdaP(own_f)=+coef, lambdaP(nei_f)=-coef regardless
            // of which end the cell is (phi = P_own - P_nei + conv)
            const scalar coefP = (o2 == own) ? w2 : (1.0 - w2);
            // BFINAL-003 (P<-U relaxation patch): alphaRel-scaled P-part
            // (relaxed mobility); convCo (convTerm path) stays unchanged.
            const scalar coP_o2 = alphaRel*coefP*rAUo2;
            const scalar coP_n2 = -alphaRel*coefP*rAUo2;
            const scalar convCo = (o2 == own) ? w2*rAUo2 : (1.0 - w2)*rAUo2;
            for (direction c = 0; c < 3; ++c)
            {
                addCoo(discreteUIndex(nei, c), discretePIndex(o2), -upperF*Vinvo*coP_o2*sf2.component(c));
                addCoo(discreteUIndex(nei, c), discretePIndex(n2), -upperF*Vinvo*coP_n2*sf2.component(c));
                for (direction c2d = 0; c2d < 3; ++c2d)
                {
                    addCoo(discreteUIndex(nei, c), discreteUIndex(down2, c2d),
                        -upperF*Vinvo*convCo*sf2.component(c)*vj2.component(c2d));
                }
            }
        }

        // U(own,c) -= lower_f/V_nei * hA[nei]
        for (const label fp : cellFacesOf[nei])
        {
            const label o2 = discreteOwner[fp];
            const label n2 = discreteNeighbour[fp];
            const scalar w2 = discreteWeights[fp];
            const vector sf2 = discreteSf[fp];
            const scalar rAUn2 = rAUAdj[nei];
            const label down2 = discretePhi[fp] >= 0.0 ? n2 : o2;
            const vector vj2 = U[n2] - U[o2];
            const scalar coefP = (o2 == nei) ? w2 : (1.0 - w2);
            // BFINAL-003 (P<-U relaxation patch): alphaRel-scaled P-part
            // (relaxed mobility); convCo (convTerm path) stays unchanged.
            const scalar coP_o2 = alphaRel*coefP*rAUn2;
            const scalar coP_n2 = -alphaRel*coefP*rAUn2;
            const scalar convCo = (o2 == nei) ? w2*rAUn2 : (1.0 - w2)*rAUn2;
            for (direction c = 0; c < 3; ++c)
            {
                addCoo(discreteUIndex(own, c), discretePIndex(o2), -lowerF*Vinn*coP_o2*sf2.component(c));
                addCoo(discreteUIndex(own, c), discretePIndex(n2), -lowerF*Vinn*coP_n2*sf2.component(c));
                for (direction c2d = 0; c2d < 3; ++c2d)
                {
                    addCoo(discreteUIndex(own, c), discreteUIndex(down2, c2d),
                        -lowerF*Vinn*convCo*sf2.component(c)*vj2.component(c2d));
                }
            }
        }
    }
}

// ---- boundary (reduced boundary-flux transpose + direct U-P) ----
forAll(mesh.boundary(), patchi)
{
    const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
    const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
    const bool pressureFixed = p.boundaryField()[patchi].fixesValue();
    // BFINAL-003 (P<-U relaxation patch): boundary hB assignable-only and
    // alphaRel-scaled; direct (1-alphaRel) relax-source transpose term on
    // assignable faces only (mirrors applyDiscreteFlowJT boundary path).
    const bool uAssignable = U.boundaryField()[patchi].assignable();
    const vectorField icV(discretePrimalMomentum.internalCoeffs()[patchi]);

    forAll(patchCells, facei)
    {
        const label celli = patchCells[facei];

        if (uAssignable)
        {
            const scalar cav = (icV[facei].x() + icV[facei].y() + icV[facei].z())/3.0;
            const scalar Vinvc = 1.0/Foam::max(mesh.V()[celli], SMALL);
            const scalar rAUb = rAUAdj[celli];
            // hB = lambdaP_cell*alphaRel*rAU*Sf_b ; Stage2:
            // U(cell,c) += Bdiag_c/V*hB[c]
            for (direction c = 0; c < 3; ++c)
            {
                const scalar Bdiag = (-icV[facei].component(c) + cav);
                addCoo
                (
                    discreteUIndex(celli, c),
                    discretePIndex(celli),
                    Bdiag*Vinvc*(alphaRel*rAUb)*patchSf[facei].component(c)
                );
            }
            // direct (1-alphaRel)*lambdaPc*patchSf -> U(cell,c)
            for (direction c = 0; c < 3; ++c)
            {
                addCoo
                (
                    discreteUIndex(celli, c),
                    discretePIndex(celli),
                    (1.0 - alphaRel)*patchSf[facei].component(c)
                );
            }
        }

        // direct U-P (P row <- U): keep
        if (!pressureFixed)
        {
            for (direction c = 0; c < 3; ++c)
            {
                addCoo(discretePIndex(celli), discreteUIndex(celli, c), patchSf[facei].component(c));
            }
        }

        // BFINAL-008 (J_PP actual-primal closure): fixedValue-p boundary P-P
        // diagonal = +rAtU_c*deltaCoeffs_b*|Sf|_b (production
        // fvm::laplacian(rAtU,p) addBoundaryDiag internalCoeffs, residual
        // convention dR_P/dp; boundaryCoeffs source = 0 since p_b = 0, so no
        // source entry; zeroGradient p patches contribute 0).  Diagonal-only
        // => symmetric, so the exported J^T carries the corrected boundary
        // term, matching the matrix-free forward/transpose boundary loops.
        if (pressureFixed)
        {
            addCoo
            (
                discretePIndex(celli),
                discretePIndex(celli),
                primalPressureMobility[celli]
               *mesh.deltaCoeffs().boundaryField()[patchi][facei]
               *mesh.magSf().boundaryField()[patchi][facei]
            );
        }
    }
}

// =========================================================================
// Explicit frozen deviatoric-stress transpose (U-U block) for the oracle.
// Replicates applyFrozenDeviatoricJT() as explicit COO entries so that the
// exported explicitJT.mtx is the COMPLETE frozen-RANS J^T.  The ExplicitJToracle
// check further below (explicit vs matrix-free) then verifies this assembly to
// ~1e-16, so an error here cannot silently pass.
// =========================================================================
if (ransFlowModel)
{
    // Per-cell gradient coefficient map: gradCoef[c] = list of
    //   (cell m, coeff[27]), coeff[(i*3+j)*3+k] = d grad[c][i][j]/d lambda[m][k]
    // built from the interior stress and boundary projection (radius-1).
    typedef std::pair<label, std::vector<scalar>> GradTerm;
    std::vector<std::vector<GradTerm>> gradCoef(mesh.nCells());

    auto gradFindOrAdd = [&](label c, label m) -> std::vector<scalar>&
    {
        auto& lst = gradCoef[c];
        for (auto& t : lst)
        {
            if (t.first == m)
            {
                return t.second;
            }
        }
        lst.emplace_back(m, std::vector<scalar>(27, 0.0));
        return lst.back().second;
    };

    auto gradAdd = [&](label c, label m, label i, label j, label k, scalar v)
    {
        if (mag(v) < 1e-300)
        {
            return;
        }
        auto& a = gradFindOrAdd(c, m);
        a[(i*3 + j)*3 + k] += v;
    };

    // ---- interior stress -> grad ----
    // grad[c][i][j] = nu_c * (sigma_c[j][i] - (1/3) delta_ij tr(sigma_c)).
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar w = discreteWeights[facei];
        const vector sf = discreteSf[facei];
        const scalar nuOwn = nuEffFrozen[own];
        const scalar nuNei = nuEffFrozen[nei];

        for (label i = 0; i < 3; ++i)
        {
            for (label j = 0; j < 3; ++j)
            {
                const scalar twoThirdDelta = (i == j) ? (2.0/3.0) : 0.0;

                // owner side: sigma_own += w * sf * (lambda_nei - lambda_own)
                gradAdd(own, nei, i, j, i,  nuOwn*w*sf[j]);
                gradAdd(own, own, i, j, i, -nuOwn*w*sf[j]);
                // neighbour side: sigma_nei += (1-w) * sf * (lambda_nei - lambda_own)
                gradAdd(nei, own, i, j, i, -nuNei*(1.0 - w)*sf[j]);
                gradAdd(nei, nei, i, j, i,  nuNei*(1.0 - w)*sf[j]);
                for (label k = 0; k < 3; ++k)
                {
                    // trace term -(2/3)delta_ij * (sf . (lambda_nei - lambda_own))
                    gradAdd(own, nei, i, j, k, -twoThirdDelta*nuOwn*w*sf[k]);
                    gradAdd(own, own, i, j, k,  twoThirdDelta*nuOwn*w*sf[k]);
                    gradAdd(nei, own, i, j, k,  twoThirdDelta*nuNei*(1.0 - w)*sf[k]);
                    gradAdd(nei, nei, i, j, k, -twoThirdDelta*nuNei*(1.0 - w)*sf[k]);
                }
            }
        }
    }

    // ---- boundary projection -> grad (and velocityFixed direct output) ----
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& patchMagSf = mesh.magSf().boundaryField()[patchi];
        const scalarField& patchDelta = mesh.deltaCoeffs().boundaryField()[patchi];
        const scalarField& patchNuEff = nuEffFrozen.boundaryField()[patchi];
        const bool velocityFixed = U.boundaryField()[patchi].fixesValue();

        forAll(patchCells, facei)
        {
            const label c = patchCells[facei];
            const vector sf = patchSf[facei];
            const scalar mag = patchMagSf[facei];
            const vector n = sf/Foam::max(mag, SMALL);
            const scalar nu = patchNuEff[facei];
            const scalar delta = patchDelta[facei];

            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    const scalar twoThirdDelta = (i == j) ? (2.0/3.0) : 0.0;
                    for (label k = 0; k < 3; ++k)
                    {
                        const scalar dki = (k == i) ? 1.0 : 0.0;
                        const scalar v = nu
                           *(
                                -sf[j]*(dki - n[i]*n[k])
                              + twoThirdDelta*sf[k]
                              - (2.0/3.0)*n[i]*n[j]*sf[k]
                            );
                        gradAdd(c, c, i, j, k, v);
                    }
                }
            }

            if (velocityFixed)
            {
                for (label j = 0; j < 3; ++j)
                {
                    for (label k = 0; k < 3; ++k)
                    {
                        const scalar v =
                            -delta*nu*( -sf[j]*n[k] + (2.0/3.0)*n[j]*sf[k] );
                        addCoo(discreteUIndex(c, j), discreteUIndex(c, k), v);
                    }
                }
            }
        }
    }

    // ---- divergence: grad -> output (interior) ----
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar w = discreteWeights[facei];
        const vector sf = discreteSf[facei];
        const scalar Vinvo = 1.0/Foam::max(mesh.V()[own], SMALL);
        const scalar Vinn = 1.0/Foam::max(mesh.V()[nei], SMALL);

        for (label j = 0; j < 3; ++j)
        {
            for (const auto& term : gradCoef[own])
            {
                const label m = term.first;
                const auto& A = term.second;
                for (label k = 0; k < 3; ++k)
                {
                    scalar cSum = 0.0;
                    for (label i = 0; i < 3; ++i)
                    {
                        cSum += sf[i]*A[(i*3 + j)*3 + k];
                    }
                    if (mag(cSum) < 1e-300)
                    {
                        continue;
                    }
                    addCoo(discreteUIndex(own, j), discreteUIndex(m, k),  w*cSum*Vinvo);
                    addCoo(discreteUIndex(nei, j), discreteUIndex(m, k),  (1.0 - w)*cSum*Vinvo);
                }
            }
            for (const auto& term : gradCoef[nei])
            {
                const label m = term.first;
                const auto& A = term.second;
                for (label k = 0; k < 3; ++k)
                {
                    scalar cSum = 0.0;
                    for (label i = 0; i < 3; ++i)
                    {
                        cSum += sf[i]*A[(i*3 + j)*3 + k];
                    }
                    if (mag(cSum) < 1e-300)
                    {
                        continue;
                    }
                    addCoo(discreteUIndex(own, j), discreteUIndex(m, k), -w*cSum*Vinn);
                    addCoo(discreteUIndex(nei, j), discreteUIndex(m, k), -(1.0 - w)*cSum*Vinn);
                }
            }
        }
    }

    // ---- divergence: grad -> output (boundary, non-fixedValue U) ----
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
            const label c = patchCells[facei];
            const vector sf = patchSf[facei];
            const scalar Vinvc = 1.0/Foam::max(mesh.V()[c], SMALL);
            for (label j = 0; j < 3; ++j)
            {
                for (const auto& term : gradCoef[c])
                {
                    const label m = term.first;
                    const auto& A = term.second;
                    for (label k = 0; k < 3; ++k)
                    {
                        scalar cSum = 0.0;
                        for (label i = 0; i < 3; ++i)
                        {
                            cSum += sf[i]*A[(i*3 + j)*3 + k];
                        }
                        if (mag(cSum) < 1e-300)
                        {
                            continue;
                        }
                        addCoo(discreteUIndex(c, j), discreteUIndex(m, k), cSum*Vinvc);
                    }
                }
            }
        }
    }
}

// ---- COO -> CSR (sort by row then col, merge duplicates) ----
const label cooNnz = cooRow.size();
std::vector<label> perm(cooNnz);
for (label i = 0; i < cooNnz; ++i)
{
    perm[i] = i;
}
std::sort
(
    perm.begin(),
    perm.end(),
    [&](label a, label b)
    {
        if (cooRow[a] != cooRow[b])
        {
            return cooRow[a] < cooRow[b];
        }
        return cooCol[a] < cooCol[b];
    }
);
std::vector<label> csrRowPtr(discreteNUnknowns + 1, 0);
std::vector<label> csrCol;
std::vector<scalar> csrVal;
label prevRow = -1;
for (label pi = 0; pi < cooNnz; ++pi)
{
    const label i = perm[pi];
    const label r = cooRow[i];
    if (r != prevRow)
    {
        for (label rr = prevRow + 1; rr <= r; ++rr)
        {
            csrRowPtr[rr] = csrCol.size();
        }
        prevRow = r;
    }
    if (!csrCol.empty() && csrCol.back() == cooCol[i]
     && csrRowPtr[r + 1] == csrCol.size())
    {
        csrVal.back() += cooVal[i];
    }
    else
    {
        csrCol.push_back(cooCol[i]);
        csrVal.push_back(cooVal[i]);
    }
}
for (label rr = prevRow + 1; rr <= discreteNUnknowns; ++rr)
{
    csrRowPtr[rr] = csrCol.size();
}

// Only all-Neumann pressure systems need the legacy gauge row.  With a
// fixedValue pressure boundary the physical P-P block is already nonsingular.
if (discretePressureNeedsReference)
{
    const label pRefRow = discretePIndex(pRefCell);
    for (label j = csrRowPtr[pRefRow]; j < csrRowPtr[pRefRow + 1]; ++j)
    {
        csrVal[j] = 0.0;
    }
    // locate/add the diagonal entry
    bool found = false;
    for (label j = csrRowPtr[pRefRow]; j < csrRowPtr[pRefRow + 1]; ++j)
    {
        if (csrCol[j] == pRefRow) { csrVal[j] = 1.0; found = true; break; }
    }
    if (!found)
    {
        csrCol.push_back(pRefRow);
        csrVal.push_back(1.0);
        // maintain order: insert into the row by resorting
        std::sort(csrCol.begin() + csrRowPtr[pRefRow], csrCol.begin() + csrRowPtr[pRefRow + 1]);
        std::sort(csrVal.begin() + csrRowPtr[pRefRow], csrVal.begin() + csrRowPtr[pRefRow + 1]);
    }
}

auto applyExplicitMatVec = [&](const scalarField& x, scalarField& y)
{
    y = scalar(0.0);
    for (label r = 0; r < discreteNUnknowns; ++r)
    {
        scalar acc = 0.0;
        for (label j = csrRowPtr[r]; j < csrRowPtr[r + 1]; ++j)
        {
            acc += csrVal[j]*x[csrCol[j]];
        }
        y[r] = acc;
    }
};

// scaled application replicating applyScaledDiscreteFlowJT()
auto applyExplicitScaled = [&](const scalarField& input, scalarField& output)
{
    scalarField physInput(input);
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        physInput[discretePIndex(celli)] *=
            discreteMomentumScale/discreteAreaScale;
    }
    scalarField physOut(discreteNUnknowns, 0.0);
    // The deviatoric transpose is now assembled explicitly into the CSR, so
    // applyExplicitMatVec already carries the FULL operator.
    applyExplicitMatVec(physInput, physOut);
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[discreteUIndex(celli, cmpt)] =
                physOut[discreteUIndex(celli, cmpt)]/discreteMomentumScale;
        }
        output[discretePIndex(celli)] =
            physOut[discretePIndex(celli)]/discreteAreaScale;
    }
    if (discretePressureNeedsReference)
    {
        output[discretePIndex(pRefCell)] = input[discretePIndex(pRefCell)];
    }
};

// ---- oracle verification: explicit vs matrix-free ----
{
    scalar maxRelL2 = 0.0, maxCosErr = 0.0, maxRelLU = 0.0, maxRelLP = 0.0;
    for (label testI = 0; testI < 5; ++testI)
    {
        scalarField vTest(discreteNUnknowns, 0.0);
        forAll(vTest, i)
        {
            vTest[i] = Foam::sin((0.011 + 0.004*testI)*(i + 1))
                     + 0.3*Foam::cos((0.007 + 0.005*testI)*(i + 2));
        }
        scalarField yRef(discreteNUnknowns, 0.0);
        scalarField yExp(discreteNUnknowns, 0.0);
        applyScaledDiscreteFlowJT(vTest, yRef);
        applyExplicitScaled(vTest, yExp);
        scalar nRef = 0.0, nDiff = 0.0, dot = 0.0;
        scalar nRefU = 0.0, nDiffU = 0.0, nRefP = 0.0, nDiffP = 0.0;
        forAll(vTest, i)
        {
            nRef += sqr(yRef[i]);
            nDiff += sqr(yExp[i] - yRef[i]);
            dot += yExp[i]*yRef[i];
        }
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                nRefU += sqr(yRef[discreteUIndex(celli, c)]);
                nDiffU += sqr(yExp[discreteUIndex(celli, c)] - yRef[discreteUIndex(celli, c)]);
            }
            nRefP += sqr(yRef[discretePIndex(celli)]);
            nDiffP += sqr(yExp[discretePIndex(celli)] - yRef[discretePIndex(celli)]);
        }
        reduce(nRef, sumOp<scalar>());
        reduce(nDiff, sumOp<scalar>());
        reduce(dot, sumOp<scalar>());
        reduce(nRefU, sumOp<scalar>());
        reduce(nDiffU, sumOp<scalar>());
        reduce(nRefP, sumOp<scalar>());
        reduce(nDiffP, sumOp<scalar>());
        maxRelL2 = Foam::max
        (
            maxRelL2,
            Foam::sqrt(nDiff/Foam::max(nRef, SMALL))
        );
        maxCosErr = Foam::max
        (
            maxCosErr,
            mag(1.0 - dot/Foam::max(Foam::sqrt(nRef*nRef), SMALL))
        );
        maxRelLU = Foam::max
        (
            maxRelLU,
            Foam::sqrt(nDiffU/Foam::max(nRefU, SMALL))
        );
        maxRelLP = Foam::max
        (
            maxRelLP,
            Foam::sqrt(nDiffP/Foam::max(nRefP, SMALL))
        );
    }
    Info<< "ExplicitJToracle: maxRelL2=" << maxRelL2
        << " maxCosErr=" << maxCosErr
        << " maxRelL_U=" << maxRelLU
        << " maxRelL_P=" << maxRelLP << endl;
}

// ---- export MatrixMarket (physical matrix + physical rhs, tagged) ----
{
    const word mmPath =
        optProperties.lookupOrDefault<fileName>("discreteExplicitMatrixFile", fileName("explicitJT.mtx"));
    const word mmRhsPath =
        optProperties.lookupOrDefault<fileName>("discreteExplicitRhsFile", fileName("explicitRhs"));
    {
        OFstream mmOut(mmPath);
        mmOut.precision(16);
        mmOut<< "%%MatrixMarket matrix coordinate real general" << nl;
        mmOut<< discreteNUnknowns << ' ' << discreteNUnknowns << ' ' << csrVal.size() << nl;
        for (label r = 0; r < discreteNUnknowns; ++r)
        {
            for (label j = csrRowPtr[r]; j < csrRowPtr[r + 1]; ++j)
            {
                mmOut<< r + 1 << ' ' << csrCol[j] + 1 << ' ' << csrVal[j] << nl;
            }
        }
    }
    OFstream rhsOut(mmRhsPath + "_" + discreteAdjointLabel + ".mtx");
    rhsOut.precision(16);
    rhsOut<< "%%MatrixMarket matrix array real general" << nl;
    rhsOut<< discreteNUnknowns << " 1" << nl;
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            rhsOut<< discreteRhs[discreteUIndex(celli, c)]*discreteMomentumScale << nl;
        }
        rhsOut<< discreteRhs[discretePIndex(celli)]*discreteAreaScale << nl;
    }
}

const label discreteMaxIterations =
    optProperties.lookupOrDefault<label>("discreteFlowAdjointMaxIter", 2000);
const scalar discreteTolerance =
    optProperties.lookupOrDefault<scalar>("discreteFlowAdjointTolerance", 1e-12);

scalarField discreteSolution(discreteNUnknowns, 0.0);
scalarField discreteResidual(discreteNUnknowns, 0.0);
scalarField discreteOperatorSolution(discreteNUnknowns, 0.0);
pinDiscretePressureReference(discreteSolution);

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

{
    scalar schurAvg = 0.0, schurMin = GREAT, schurMax = 0.0;
    forAll(discreteSchurDiag, celli)
    {
        schurAvg += discreteSchurDiag[celli];
        schurMin = Foam::min(schurMin, discreteSchurDiag[celli]);
        schurMax = Foam::max(schurMax, discreteSchurDiag[celli]);
    }
    reduce(schurAvg, sumOp<scalar>());
    reduce(schurMin, minOp<scalar>());
    reduce(schurMax, maxOp<scalar>());
    Info<< "SchurDiag: avg=" << schurAvg/Foam::max(scalar(discreteNCells), 1)
        << " min=" << schurMin << " max=" << schurMax << endl;
}
// ---------------------------------------------------------------------------
// Block preconditioner for J^T  (physical space):
//   [A^T  G^T] [zU]   [rU]      A   = fvm momentum matrix (div/laplacian/Sp)
//   [D^T  S^T] [zp] = [rp]      S   = SIMPLE pressure Poisson laplacian(rAtU)
// M:  zU* = (A^T)^-1 rU                 (PBiCGStab, transpose)
//     rp* = rp - D^T zU*                (P-row U coupling of J^T)
//     zp  = S^-1 rp*                    (GAMG)
//     zU  = zU* - (A^T)^-1 (C^T zp)     (C^T = U-row p coupling of J^T)
// Signs follow applyDiscreteFlowJ/applyDiscreteFlowJT exactly.
// ---------------------------------------------------------------------------
volScalarField discretePcMobility(p);
discretePcMobility.rename("discretePcMobility");
scalar pcMobAvg = 0.0, pcMobMin = GREAT, pcMobMax = 0.0;
forAll(discretePcMobility, celli)
{
    const scalar m = primalPressureMobility[celli];
    discretePcMobility[celli] = m;
    pcMobAvg += m;
    pcMobMin = Foam::min(pcMobMin, m);
    pcMobMax = Foam::max(pcMobMax, m);
}
reduce(pcMobAvg, sumOp<scalar>());
reduce(pcMobMin, minOp<scalar>());
reduce(pcMobMax, maxOp<scalar>());
pcMobAvg /= Foam::max(scalar(discreteNCells), 1.0);
Info<< "PC mobility: avg=" << pcMobAvg
    << " min=" << pcMobMin << " max=" << pcMobMax << endl;
// Jacobi-scale the Schur: mobility/avg(mobility) so that diag(S) ~ O(1)
forAll(discretePcMobility, celli)
{
    discretePcMobility[celli] =
        discretePcMobility[celli]/Foam::max(pcMobAvg, SMALL);
}
// clamp to positive definite (avoid near-singular Schur for PCG)
forAll(discretePcMobility, celli)
{
    discretePcMobility[celli] =
        Foam::max(discretePcMobility[celli], 1e-6);
}
discretePcMobility.correctBoundaryConditions();

// diagonal momentum scale from the physical fvm matrix (kept alive)
tmp<fvVectorMatrix> tDiscretePcA
(
    fvm::div(phi, U, "div(phi,U)")
  - fvm::laplacian(nuEffFrozen, U)
  + fvm::Sp(alpha, U)
);
fvVectorMatrix& discretePcA = tDiscretePcA.ref();
scalar pcAvgDiag = 0.0;
forAll(discretePcA.diag(), celli)
{
    pcAvgDiag += mag(discretePcA.diag()[celli]);
}
reduce(pcAvgDiag, sumOp<scalar>());
pcAvgDiag /= Foam::max(scalar(discreteNCells), 1.0);
pcAvgDiag = Foam::max(pcAvgDiag, SMALL);

dictionary discretePcSdict;
discretePcSdict.add("solver", "PCG");
discretePcSdict.add("preconditioner", "DIC");
discretePcSdict.add("tolerance", 1e-8);
discretePcSdict.add("relTol", 0.0);
discretePcSdict.add("maxIter", 500);

auto applyDiscreteBlockPreconditioner =
[&](const scalarField& input, scalarField& output)
{
    // scaled RHS -> physical RHS  (J^T_scaled = D_out*J^T_phys*D_in)
    volVectorField rUField("pcRU", Uc);
    rUField.primitiveFieldRef() = vector::zero;
    volScalarField rpField("pcRP", pc);
    rpField.primitiveFieldRef() = scalar(0.0);
    forAll(rUField, celli)
    {
        rUField[celli] = discreteMomentumScale*vector
        (
            input[discreteUIndex(celli, 0)],
            input[discreteUIndex(celli, 1)],
            input[discreteUIndex(celli, 2)]
        );
        rpField[celli] = discreteAreaScale*input[discretePIndex(celli)];
    }

    output = input;
    return;
    volVectorField zUstar("pcZU", Uc);
    zUstar.primitiveFieldRef() = vector::zero;
    forAll(mesh.boundary(), patchi)
    {
        zUstar.boundaryFieldRef()[patchi] == vector::zero;
    }
    forAll(zUstar, celli)
    {
        const scalar invD =
            1.0/Foam::max(mag(discretePcA.diag()[celli]), 1e-2*pcAvgDiag);
        zUstar[celli] = invD*rUField[celli];
    }

    // 2. rp* = rp - D^T zU* - correction^T(conv(zU*))  (P-row couplings)
    const scalarField& pcDeltaCoeffs = mesh.deltaCoeffs().primitiveField();
    const scalarField& pcFaceAreas = mesh.magSf().primitiveField();
    volScalarField rpStar("pcRPS", pc);
    rpStar.primitiveFieldRef() = scalar(0.0);
    forAll(rpStar, celli)
    {
        rpStar[celli] = rpField[celli];
    }
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar weight = discreteWeights[facei];
        const scalar pT =
            discreteSf[facei] & (zUstar[own] - zUstar[nei]);
        rpStar[own] -= weight*pT;
        rpStar[nei] -= (1.0 - weight)*pT;

        // correction^T(conv(zU*)): J^T P-row lambda_U coupling
        const label downwind =
            discretePhi[facei] >= 0.0 ? nei : own;
        const vector lambdaDown
        (
            zUstar[downwind].x(),
            zUstar[downwind].y(),
            zUstar[downwind].z()
        );
        const vector velocityJump = U[nei] - U[own];
        const scalar convFace = lambdaDown & velocityJump;
        const scalar mobilityFace =
            weight*primalPressureMobility[own]
          + (1.0 - weight)*primalPressureMobility[nei];
        const scalar normCoeff =
            mobilityFace*pcDeltaCoeffs[facei]*pcFaceAreas[facei];
        rpStar[own] += normCoeff*convFace;
        rpStar[nei] -= normCoeff*convFace;
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
            rpStar[celli] -= patchSf[facei] & zUstar[celli];
        }
    }
    // Remove the constant component only for an all-Neumann Schur system.
    // A fixedValue pressure boundary makes this projection non-physical.
    if (discretePressureNeedsReference)
    {
        scalar rpMean = sum(rpStar.primitiveField());
        reduce(rpMean, sumOp<scalar>());
        rpMean /= Foam::max(scalar(discreteNCells), 1.0);
        rpStar.primitiveFieldRef() -= rpMean;
    }

    // 3. zp = S^-1 rp*   (GAMG, reference-pinned; psi = zpField)
    volScalarField zpField("pcZP", pc);
    zpField.primitiveFieldRef() = scalar(0.0);
    forAll(mesh.boundary(), patchi)
    {
        zpField.boundaryFieldRef()[patchi] == scalar(0.0);
    }
    {
        tmp<fvScalarMatrix> tS
        (
            fvm::laplacian(discretePcMobility, zpField)
        );
        fvScalarMatrix& S = tS.ref();
        if (discretePressureNeedsReference)
        {
            S.setReference(pRefCell, 0.0);
        }
        S.source() = rpStar.primitiveField();
        S.solve(discretePcSdict);
        // Enforce the gauge only when the Schur system actually needs one.
        if (discretePressureNeedsReference)
        {
            const scalar zpRef = zpField[pRefCell];
            zpField.primitiveFieldRef() -= zpRef;
        }
        // undo the Jacobi scaling (mobility was divided by its average),
        // so zp is the physical Schur solution again
        zpField.primitiveFieldRef() /= pcMobAvg;
        // gauge check (once): true residual  S*zp - src  and cos(S*zp, rp*)
        {
            static bool sGaugePrinted = false;
            if
            (
                !sGaugePrinted
             && optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false)
            )
            {
                sGaugePrinted = true;
                const scalarField sPsi(zpField.primitiveField());
                const scalarField sSrc(S.source());
                // manual interior residual:  S*zp - src  (boundary ignored)
                scalarField sRes = sSrc;
                const scalarField& sDiag = S.diag();
                const scalarField& sLower = S.lower();
                const scalarField& sUpper = S.upper();
                forAll(sRes, celli)
                {
                    sRes[celli] -= sDiag[celli]*sPsi[celli];
                }
                forAll(sLower, facei)
                {
                    const label own = discreteOwner[facei];
                    const label nei = discreteNeighbour[facei];
                    sRes[own] -= sUpper[facei]*sPsi[nei];
                    sRes[nei] -= sLower[facei]*sPsi[own];
                }
                scalar nRes = 0.0, nSrc = 0.0, nSzp = 0.0, dotSzpRp = 0.0;
                forAll(sRes, celli)
                {
                    nRes += sqr(sRes[celli]);
                    nSrc += sqr(sSrc[celli]);
                    nSzp += sqr(sSrc[celli] - sRes[celli]);
                    dotSzpRp += (sSrc[celli] - sRes[celli])
                               *rpStar[celli];
                }
                reduce(nRes, sumOp<scalar>());
                reduce(nSrc, sumOp<scalar>());
                reduce(nSzp, sumOp<scalar>());
                reduce(dotSzpRp, sumOp<scalar>());
                Info<< "S gauge: relRes=" << Foam::sqrt(nRes/Foam::max(nSrc, SMALL))
                    << " cos(S*zp, rp*)="
                    << dotSzpRp/Foam::max(Foam::sqrt(nSzp*nSrc), SMALL)
                    << " |zp|=" << Foam::sqrt(Foam::sum(sqr(sPsi)))
                    << " pcAvg=" << gAverage(pc)
                    << " pcMax=" << gMax(pc)
                    << " zp[pRef]=" << sPsi[pRefCell]
                    << " zpMax=" << gMax(sPsi) << endl;
            }
        }
    }

    // 4. zUcorr = invDiag * (C^T zp)   (U-row p coupling, physical)
    vectorField rhsCorr(mesh.nCells(), vector::zero);
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        const scalar weight = discreteWeights[facei];
        const vector cT =
            discreteSf[facei]*(zpField[own] - zpField[nei]);
        rhsCorr[own] += weight*cT;
        rhsCorr[nei] += (1.0 - weight)*cT;
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
            rhsCorr[celli] += patchSf[facei]*zpField[celli];
        }
    }
    volVectorField zUcorr("pcZUC", Uc);
    zUcorr.primitiveFieldRef() = vector::zero;
    forAll(mesh.boundary(), patchi)
    {
        zUcorr.boundaryFieldRef()[patchi] == vector::zero;
    }
    forAll(zUcorr, celli)
    {
        const scalar invD =
            1.0/Foam::max(mag(discretePcA.diag()[celli]), 1e-2*pcAvgDiag);
        zUcorr[celli] = invD*rhsCorr[celli];
    }

    // one-shot internal diagnostics (always print while probing)
    {
        if (optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false))
        {
            scalar nIn = 0.0, nZU = 0.0, nRP = 0.0, nZP = 0.0, nZC = 0.0;
            scalar nOutP = 0.0, zpMaxA = 0.0;
            forAll(zUstar, celli)
            {
                nIn += sqr(input[discreteUIndex(celli, 0)])
                     + sqr(input[discreteUIndex(celli, 1)])
                     + sqr(input[discreteUIndex(celli, 2)])
                     + sqr(input[discretePIndex(celli)]);
                nZU += magSqr(zUstar[celli]);
                nRP += sqr(rpStar[celli]);
                nZP += sqr(zpField[celli]);
                nZC += magSqr(zUcorr[celli]);
                const scalar outP =
                    zpField[celli]*(discreteAreaScale/discreteMomentumScale);
                nOutP += sqr(outP);
                zpMaxA = Foam::max(zpMaxA, mag(zpField[celli]));
            }
            reduce(nIn, sumOp<scalar>());
            reduce(nZU, sumOp<scalar>());
            reduce(nRP, sumOp<scalar>());
            reduce(nZP, sumOp<scalar>());
            reduce(nZC, sumOp<scalar>());
            reduce(nOutP, sumOp<scalar>());
            Info<< "PC diag: |in|=" << Foam::sqrt(nIn)
                << " |zU*|=" << Foam::sqrt(nZU)
                << " |rp*|=" << Foam::sqrt(nRP)
                << " |zp|=" << Foam::sqrt(nZP)
                << " zpMax=" << zpMaxA
                << " |outP|=" << Foam::sqrt(nOutP)
                << " |zUcorr|=" << Foam::sqrt(nZC) << endl;
        }
    }

    // 5. output (physical -> scaled):  x_U = y_U,  x_P = y_P*(A/M)
    forAll(zUstar, celli)
    {
        const vector zU = zUstar[celli] - zUcorr[celli];
        output[discreteUIndex(celli, 0)] = zU.x();
        output[discreteUIndex(celli, 1)] = zU.y();
        output[discreteUIndex(celli, 2)] = zU.z();
        output[discretePIndex(celli)] =
            zpField[celli]*(discreteAreaScale/discreteMomentumScale);
    }
    if (discretePressureNeedsReference)
    {
        output[discretePIndex(pRefCell)] = input[discretePIndex(pRefCell)];
    }
};


// PC correctness probe: M*(J^T_scaled*v) should approximate v
if (optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false))
{
    scalarField pcTestV(discreteNUnknowns, 0.0);
    forAll(pcTestV, i)
    {
        pcTestV[i] = Foam::sin(0.021*(i + 1));
    }
    scalarField pcJTv(discreteNUnknowns, 0.0);
    scalarField pcMv(discreteNUnknowns, 0.0);
    applyScaledDiscreteFlowJT(pcTestV, pcJTv);
    applyDiscreteBlockPreconditioner(pcJTv, pcMv);
    scalar cosPC = 0.0, nV = 0.0, nM = 0.0, nJT = 0.0;
    forAll(pcTestV, i)
    {
        cosPC += pcMv[i]*pcTestV[i];
        nV += sqr(pcTestV[i]);
        nM += sqr(pcMv[i]);
        nJT += sqr(pcJTv[i]);
    }
    reduce(cosPC, sumOp<scalar>());
    reduce(nV, sumOp<scalar>());
    reduce(nM, sumOp<scalar>());
    reduce(nJT, sumOp<scalar>());
    Info<< "PC probe: cos(M*JT*v, v)="
        << cosPC/Foam::max(Foam::sqrt(nV*nM), SMALL)
        << " |M*JT*v|/|v|=" << Foam::sqrt(nM/Foam::max(nV, SMALL))
        << " |JT*v|/|v|=" << Foam::sqrt(nJT/Foam::max(nV, SMALL))
        << endl;

    // sub-probe A: pure-U test vector
    {
        scalarField vU(discreteNUnknowns, 0.0);
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            const scalar s = Foam::sin(0.023*(celli + 1));
            vU[discreteUIndex(celli, 0)] = s;
            vU[discreteUIndex(celli, 1)] = 0.5*s;
            vU[discreteUIndex(celli, 2)] = -0.4*s;
        }
        scalarField JTvU(discreteNUnknowns, 0.0);
        scalarField MvU(discreteNUnknowns, 0.0);
        applyScaledDiscreteFlowJT(vU, JTvU);
        applyDiscreteBlockPreconditioner(JTvU, MvU);
        scalar cosU = 0.0;
        scalar nVU = 0.0, nMvU = 0.0;
        forAll(vU, i)
        {
            cosU += MvU[i]*vU[i];
            nVU += sqr(vU[i]);
            nMvU += sqr(MvU[i]);
        }
        reduce(cosU, sumOp<scalar>());
        reduce(nVU, sumOp<scalar>());
        reduce(nMvU, sumOp<scalar>());
        Info<< "PC sub-probe U: cos="
            << cosU/Foam::max(Foam::sqrt(nVU*nMvU), SMALL)
            << " ratio=" << Foam::sqrt(nMvU/Foam::max(nVU, SMALL))
            << endl;
    }

    // sub-probe B: pure-P test vector
    {
        scalarField vP(discreteNUnknowns, 0.0);
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            vP[discretePIndex(celli)] = Foam::cos(0.029*(celli + 1));
        }
        scalarField JTvP(discreteNUnknowns, 0.0);
        scalarField MvP(discreteNUnknowns, 0.0);
        applyScaledDiscreteFlowJT(vP, JTvP);
        applyDiscreteBlockPreconditioner(JTvP, MvP);
        scalar cosP = 0.0, nVP = 0.0, nMvP = 0.0;
        scalar nJTU = 0.0, nJTP = 0.0, nMvPU = 0.0, nMvPP = 0.0;
        forAll(vP, i)
        {
            cosP += MvP[i]*vP[i];
            nVP += sqr(vP[i]);
            nMvP += sqr(MvP[i]);
        }
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                nJTU += sqr(JTvP[discreteUIndex(celli, c)]);
                nMvPU += sqr(MvP[discreteUIndex(celli, c)]);
            }
            nJTP += sqr(JTvP[discretePIndex(celli)]);
            nMvPP += sqr(MvP[discretePIndex(celli)]);
        }
        reduce(cosP, sumOp<scalar>());
        reduce(nVP, sumOp<scalar>());
        reduce(nMvP, sumOp<scalar>());
        reduce(nJTU, sumOp<scalar>());
        reduce(nJTP, sumOp<scalar>());
        reduce(nMvPU, sumOp<scalar>());
        reduce(nMvPP, sumOp<scalar>());
        Info<< "PC sub-probe P: cos="
            << cosP/Foam::max(Foam::sqrt(nVP*nMvP), SMALL)
            << " ratio=" << Foam::sqrt(nMvP/Foam::max(nVP, SMALL))
            << " |JT_U|=" << Foam::sqrt(nJTU)
            << " |JT_P|=" << Foam::sqrt(nJTP)
            << " |M_U|=" << Foam::sqrt(nMvPU)
            << " |M_P|=" << Foam::sqrt(nMvPP)
            << endl;
        // locate the largest |JT_U| entries
        {
            scalar best1 = 0.0, best2 = 0.0, best3 = 0.0;
            label bc1 = -1, bc2 = -1, bc3 = -1;
            for (label celli = 0; celli < discreteNCells; ++celli)
            {
                scalar v = 0.0;
                for (direction c = 0; c < 3; ++c)
                {
                    v += sqr(JTvP[discreteUIndex(celli, c)]);
                }
                v = Foam::sqrt(v);
                if (v > best1) { best3 = best2; bc3 = bc2; best2 = best1; bc2 = bc1; best1 = v; bc1 = celli; }
                else if (v > best2) { best3 = best2; bc3 = bc2; best2 = v; bc2 = celli; }
                else if (v > best3) { best3 = v; bc3 = celli; }
            }
            Info<< "  JT_U top cells: c=" << bc1 << " |U|=" << best1
                << " |vP|=" << mag(vP[discretePIndex(bc1)])
                << " alpha=" << alpha[bc1] << " U=" << U[bc1]
                << " mob=" << primalPressureMobility[bc1] << endl;
            Info<< "  JT_U top cells: c=" << bc2 << " |U|=" << best2
                << " |vP|=" << mag(vP[discretePIndex(bc2)])
                << " alpha=" << alpha[bc2] << " U=" << U[bc2] << endl;
            Info<< "  JT_U top cells: c=" << bc3 << " |U|=" << best3
                << " |vP|=" << mag(vP[discretePIndex(bc3)])
                << " alpha=" << alpha[bc3] << " U=" << U[bc3] << endl;
        Info<< "  scales: M=" << discreteMomentumScale
            << " A=" << discreteAreaScale
            << " M/A=" << discreteMomentumScale/discreteAreaScale
            << " |Sf|avg=" << gAverage(mag(discreteSf))
            << " avgDiag=" << pcAvgDiag << endl;
        Info<< "  P-cell0: vP=" << vP[discretePIndex(0)]
            << " JTvP[P]=" << JTvP[discretePIndex(0)]
            << " MvP[P]=" << MvP[discretePIndex(0)]
            << " MvP[Ux]=" << MvP[discreteUIndex(0, 0)] << endl;
        }
    }
}

const label discreteGmresRestart =
    optProperties.lookupOrDefault<label>("discreteFlowGmresRestart", 40);
if (discreteGmresRestart < 2)
{
    FatalErrorInFunction
        << "discreteFlowGmresRestart must be at least 2."
        << abort(FatalError);
}

// =========================================================================
// Explicit-matrix GMRES + ILU0 (right preconditioned).
// Operator: applyExplicitScaled (CSR + dev2); preconditioner: ILU0(CSR).
// True residual is re-checked with the matrix-free operator at each restart.
// =========================================================================
// ---- ILU0 on the physical CSR (no fill-in) ----
std::vector<scalar> iluVal(csrVal);
scalar iluEps = 1e-10;
{
    scalar diagSum = 0.0;
    for (label i = 0; i < discreteNUnknowns; ++i)
    {
        for (label p = csrRowPtr[i]; p < csrRowPtr[i + 1]; ++p)
        {
            if (csrCol[p] == i) { diagSum += mag(csrVal[p]); break; }
        }
    }
    reduce(diagSum, sumOp<scalar>());
    iluEps =
        1e-6*diagSum/Foam::max(scalar(discreteNUnknowns), 1.0);
    for (label i = 0; i < discreteNUnknowns; ++i)
    {
        for (label p = csrRowPtr[i]; p < csrRowPtr[i + 1]; ++p)
        {
            if (csrCol[p] == i)
            {
                if (mag(iluVal[p]) < iluEps)
                {
                    iluVal[p] = (iluVal[p] >= 0.0 ? 1.0 : -1.0)*iluEps;
                }
                break;
            }
        }
    }
}
{
    label nanCount = 0;
    for (auto v : iluVal) { if (std::isnan(v) || std::isinf(v)) ++nanCount; }
    Info<< "ILU: pre-decomp nnz=" << iluVal.size() << " nan/inf=" << nanCount << endl;
}
for (label i = 0; i < discreteNUnknowns; ++i)
{
    for (label p = csrRowPtr[i]; p < csrRowPtr[i + 1]; ++p)
    {
        const label j = csrCol[p];
        if (j >= i)
        {
            break;
        }
        scalar sum = iluVal[p];
        for (label q = csrRowPtr[j]; q < csrRowPtr[j + 1]; ++q)
        {
            const label k = csrCol[q];
            if (k >= j)
            {
                break;
            }
            label posIK = -1, posKJ = -1;
            for (label r = csrRowPtr[i]; r < csrRowPtr[i + 1]; ++r)
            {
                if (csrCol[r] == k) { posIK = r; break; }
            }
            for (label r = csrRowPtr[k]; r < csrRowPtr[k + 1]; ++r)
            {
                if (csrCol[r] == j) { posKJ = r; break; }
            }
            if (posIK >= 0 && posKJ >= 0)
            {
                sum -= iluVal[posIK]*iluVal[posKJ];
            }
        }
        iluVal[p] = sum;
        label posJI = -1, posII = -1;
        for (label r = csrRowPtr[j]; r < csrRowPtr[j + 1]; ++r)
        {
            if (csrCol[r] == i) { posJI = r; break; }
        }
        for (label r = csrRowPtr[i]; r < csrRowPtr[i + 1]; ++r)
        {
            if (csrCol[r] == i) { posII = r; break; }
        }
        if (posJI >= 0 && posII >= 0)
        {
            iluVal[posII] -= sum*iluVal[posJI];
            if (mag(iluVal[posII]) < iluEps)
            {
                iluVal[posII] =
                    (iluVal[posII] >= 0.0 ? 1.0 : -1.0)*iluEps;
            }
        }
    }
}
{
    label nanCount2 = 0;
    scalar maxAbs = 0.0;
    for (auto v : iluVal)
    {
        if (std::isnan(v) || std::isinf(v)) { ++nanCount2; }
        maxAbs = Foam::max(maxAbs, mag(v));
    }
    Info<< "ILU post-decomp nan=" << nanCount2 << " maxAbs=" << maxAbs << endl;
}
auto applyIluPreconditioner =
[&](const scalarField& r, scalarField& z)
{
    scalarField y(discreteNUnknowns, 0.0);
    for (label i = 0; i < discreteNUnknowns; ++i)
    {
        scalar acc = r[i];
        for (label p = csrRowPtr[i]; p < csrRowPtr[i + 1]; ++p)
        {
            const label j = csrCol[p];
            if (j >= i) { break; }
            acc -= iluVal[p]*y[j];
        }
        y[i] = acc;
    }
    for (label i = discreteNUnknowns - 1; i >= 0; --i)
    {
        scalar acc = y[i];
        scalar diagI = SMALL;
        for (label p = csrRowPtr[i]; p < csrRowPtr[i + 1]; ++p)
        {
            const label j = csrCol[p];
            if (j > i) { acc -= iluVal[p]*z[j]; }
            else if (j == i) { diagI = iluVal[p]; }
        }
        z[i] = acc/Foam::max(mag(diagI), iluEps);
    }
};

// ---- GMRES(restart) + ILU0 ----
const label discreteRestart =
    optProperties.lookupOrDefault<label>("discreteFlowGmresRestart", 50);
std::vector<scalarField> gmresV(discreteRestart + 1);
std::vector<scalarField> gmresZ(discreteRestart);
for (label i = 0; i <= discreteRestart; ++i)
{
    gmresV[i].setSize(discreteNUnknowns, 0.0);
}
for (label i = 0; i < discreteRestart; ++i)
{
    gmresZ[i].setSize(discreteNUnknowns, 0.0);
}
std::vector<std::vector<scalar>> gmresH
(
    discreteRestart + 1,
    std::vector<scalar>(discreteRestart, 0.0)
);
std::vector<scalar> gmresC(discreteRestart, 0.0);
std::vector<scalar> gmresS(discreteRestart, 0.0);
std::vector<scalar> gmresG(discreteRestart + 1, 0.0);
std::vector<scalar> gmresY(discreteRestart, 0.0);
scalarField gmresW(discreteNUnknowns, 0.0);

label discreteIteration = 0;
scalar discreteRelativeResidual = GREAT;

// solve in the PHYSICAL system:  A_phys y = b_phys
scalarField physRhs(discreteNUnknowns, 0.0);
for (label celli = 0; celli < discreteNCells; ++celli)
{
    for (direction c = 0; c < 3; ++c)
    {
        physRhs[discreteUIndex(celli, c)] =
            discreteRhs[discreteUIndex(celli, c)]*discreteMomentumScale;
    }
    physRhs[discretePIndex(celli)] =
        discreteRhs[discretePIndex(celli)]*discreteAreaScale;
}
pinDiscretePressureReference(physRhs);

auto applyExplicitPhys = [&](const scalarField& x, scalarField& y)
{
    // deviatoric transpose is assembled into the CSR; this is now the full op.
    applyExplicitMatVec(x, y);
};

scalarField physSolution(discreteNUnknowns, 0.0);
// optional: load a precomputed physical solution (offline SciPy oracle)
{
    const word solBase =
        optProperties.lookupOrDefault<fileName>("discreteExplicitSolutionFile", fileName("explicitSol"));
    const fileName solFile(solBase + "_" + discreteAdjointLabel + ".mtx");
    if
    (
        optProperties.lookupOrDefault<Switch>("discreteUseExplicitSolution", false)
     && isFile(solFile)
    )
    {
        IFstream solIn(solFile);
        label cnt = 0;
        scalar v = 0.0;
        while (solIn.good() && cnt < discreteNUnknowns)
        {
            solIn>> v;
            physSolution[cnt] = v;
            ++cnt;
        }
        if (cnt >= discreteNUnknowns)
        {
            for (label celli = 0; celli < discreteNCells; ++celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    discreteSolution[discreteUIndex(celli, c)] =
                        physSolution[discreteUIndex(celli, c)];
                }
                discreteSolution[discretePIndex(celli)] =
                    physSolution[discretePIndex(celli)]
                   *(discreteAreaScale/discreteMomentumScale);
            }
            pinDiscretePressureReference(discreteSolution);
            applyScaledDiscreteFlowJT(discreteSolution, discreteOperatorSolution);
            discreteResidual = discreteRhs - discreteOperatorSolution;
            discreteRelativeResidual = Foam::sqrt
            (
                discreteGlobalDot(discreteResidual, discreteResidual)
            )/discreteRhsNorm;
            Info<< "Loaded explicit solution " << solFile
                << " relRes=" << discreteRelativeResidual << endl;
        }
    }
}
scalarField physResidual(discreteNUnknowns, 0.0);
const scalar physRhsNorm = Foam::sqrt
(
    Foam::max(discreteGlobalDot(physRhs, physRhs), scalar(1e-300))
);

applyExplicitPhys(physSolution, physResidual);
physResidual = physRhs - physResidual;
discreteRelativeResidual = Foam::sqrt
(
    discreteGlobalDot(physResidual, physResidual)
)/physRhsNorm;

while
(
    discreteIteration < discreteMaxIterations
 && discreteRelativeResidual > discreteTolerance
 && !optProperties.lookupOrDefault<Switch>("discreteExportOnly", false)
)
{
    applyExplicitPhys(physSolution, physResidual);
    physResidual = physRhs - physResidual;
    const scalar beta = Foam::sqrt
    (
        discreteGlobalDot(physResidual, physResidual)
    );
    discreteRelativeResidual = beta/physRhsNorm;
    if (discreteRelativeResidual <= discreteTolerance)
    {
        break;
    }
    gmresV[0] = physResidual/beta;
    gmresG[0] = beta;
    for (label j = 1; j <= discreteRestart; ++j) { gmresG[j] = 0.0; }
    for (label j = 0; j < discreteRestart; ++j)
    {
        gmresC[j] = 0.0;
        gmresS[j] = 0.0;
        gmresY[j] = 0.0;
    }
    for (label i = 0; i <= discreteRestart; ++i)
    {
        for (label j = 0; j < discreteRestart; ++j)
        {
            gmresH[i][j] = 0.0;
        }
    }

    label innerDone = 0;
    for (label k = 0; k < discreteRestart; ++k)
    {
        applyIluPreconditioner(gmresV[k], gmresZ[k]);
        applyExplicitPhys(gmresZ[k], gmresW);
        for (label i = 0; i <= k; ++i)
        {
            gmresH[i][k] = discreteGlobalDot(gmresW, gmresV[i]);
            forAll(gmresW, ii)
            {
                gmresW[ii] -= gmresH[i][k]*gmresV[i][ii];
            }
        }
        const scalar hk1 = Foam::sqrt(discreteGlobalDot(gmresW, gmresW));
        if (k < 4)
        {
            scalar wMax = 0.0, zMax = 0.0;
            forAll(gmresW, ii)
            {
                wMax = Foam::max(wMax, mag(gmresW[ii]));
                zMax = Foam::max(zMax, mag(gmresZ[k][ii]));
            }
            Info<< "  [GMRES-diag " << discreteAdjointLabel << "] k=" << k
                << " hk1=" << hk1 << " |W|max=" << wMax << " |Z|max=" << zMax
                << " h00=" << gmresH[0][k] << endl;
        }
        if (hk1 < 1e-12)
        {
            innerDone = k + 1;
            break;
        }
        gmresH[k + 1][k] = hk1;
        gmresV[k + 1] = gmresW/hk1;
        for (label i = 0; i < k; ++i)
        {
            const scalar tmp =
                gmresC[i]*gmresH[i][k] + gmresS[i]*gmresH[i + 1][k];
            gmresH[i + 1][k] =
                -gmresS[i]*gmresH[i][k] + gmresC[i]*gmresH[i + 1][k];
            gmresH[i][k] = tmp;
        }
        const scalar rot =
            Foam::sqrt(gmresH[k][k]*gmresH[k][k] + gmresH[k + 1][k]*gmresH[k + 1][k]);
        if (rot < 1e-12)
        {
            innerDone = k + 1;
            break;
        }
        gmresC[k] = gmresH[k][k]/rot;
        gmresS[k] = gmresH[k + 1][k]/rot;
        gmresH[k][k] = rot;
        gmresH[k + 1][k] = 0.0;
        const scalar gk = gmresG[k];
        gmresG[k] = gmresC[k]*gk;
        gmresG[k + 1] = -gmresS[k]*gk;
        innerDone = k + 1;
        discreteRelativeResidual = mag(gmresG[k + 1])/discreteRhsNorm;
        if (discreteRelativeResidual <= discreteTolerance)
        {
            break;
        }
    }

    for (label i = innerDone - 1; i >= 0; --i)
    {
        scalar sum = gmresG[i];
        for (label j = i + 1; j < innerDone; ++j)
        {
            sum -= gmresH[i][j]*gmresY[j];
        }
        gmresY[i] = sum/Foam::max(mag(gmresH[i][i]), SMALL);
    }
    for (label i = 0; i < innerDone; ++i)
    {
        forAll(physSolution, ii)
        {
            physSolution[ii] += gmresY[i]*gmresZ[i][ii];
        }
    }
    pinDiscretePressureReference(physSolution);
    discreteIteration += innerDone;
    Info<< "  [GMRES-ILU " << discreteAdjointLabel << "] iter="
        << discreteIteration << " relRes=" << discreteRelativeResidual
        << endl;
}

// map the physical solution back to the scaled system
for (label celli = 0; celli < discreteNCells; ++celli)
{
    for (direction c = 0; c < 3; ++c)
    {
        discreteSolution[discreteUIndex(celli, c)] =
            physSolution[discreteUIndex(celli, c)];
    }
    discreteSolution[discretePIndex(celli)] =
        physSolution[discretePIndex(celli)]
       *(discreteAreaScale/discreteMomentumScale);
}
pinDiscretePressureReference(discreteSolution);

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
    // DIAGNOSTIC: warn and continue with the partially converged solution
    // (the block Schur preconditioner does not yet match the corrected J).
    WarningInFunction
        << "Reduced discrete flow adjoint " << discreteAdjointLabel
        << " residual=" << discreteRelativeResidual
        << " > tolerance=" << discreteTolerance << " (continuing, diagnostic)"
        << endl;
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

// Uc fingerprint: verify the oracle solution survives to sensitivity.H
if (optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false))
{
    scalar sUc = 0.0, nUc = 0.0, mUc = 0.0;
    forAll(discreteAdjointVelocity, celli)
    {
        const vector ucv = discreteAdjointVelocity[celli];
        sUc += ucv.x() + ucv.y() + ucv.z();
        nUc += magSqr(ucv);
        mUc = Foam::max(mUc, mag(ucv));
    }
    reduce(sUc, sumOp<scalar>());
    reduce(nUc, sumOp<scalar>());
    reduce(mUc, maxOp<scalar>());
    Info<< "UcFingerprint[afterSolve " << discreteAdjointLabel << "]:"
        << " sum=" << sUc << " L2=" << Foam::sqrt(nUc)
        << " max=" << mUc
        << " c0=" << discreteAdjointVelocity[0]
        << " c100=" << discreteAdjointVelocity[100]
        << " c5000=" << discreteAdjointVelocity[5000] << endl;
}

// ---------------------------------------------------------------------------
// Stage B5 (BFINAL-002, DIAGNOSTIC_ONLY): boundary-consistent P-row oracle +
// isolated dphi/dU relaxation oracle.  Guarded by optProperties switch
// "stageB5BoundaryRelaxOracle" (default false).  When disabled this block is
// a no-op (one dictionary lookup only): no fields created/modified, no files
// written, no effect on production values/residuals.  When enabled it runs
// once (static guard) on the first inclusion of this header (thermalCoupling
// adjoint block) using the converged production state; it never modifies
// production state.
// ---------------------------------------------------------------------------
if (optProperties.lookupOrDefault<Switch>("stageB5BoundaryRelaxOracle", false))
{
    static bool stageB5BoundaryRelaxOracleDone = false;
    if (!stageB5BoundaryRelaxOracleDone)
    {
        stageB5BoundaryRelaxOracleDone = true;
        #include "stageB5BoundaryRelaxOracle.H"
    }
}

// ---------------------------------------------------------------------------
// Stage B6 (BFINAL-004, DIAGNOSTIC_ONLY): fixed-state R_x closure oracle
// (RX-A direct alpha derivative, RX-B xh-to-alpha closure, RX-C raw-design
// directional R_x, §9/§10 weighted contributions, volume-gradient regression
// anchor, production gsensPressureDrop replica).  Guarded by optProperties
// switch "stageB6RxDesignOracle" (default false).  When disabled this block is
// a no-op (one dictionary lookup only): no fields created/modified, no files
// written, no effect on production values/residuals.  When enabled it runs
// once (static guard) in the pressureDrop adjoint block only (where Uc/pc are
// the solved exact pressure-drop adjoint); it never modifies production state
// and never solves UEqn/pEqn (only assembles/reads residual(), A(), H(),
// flux() and solves diagnostic filter systems on local fields).
// ---------------------------------------------------------------------------
if (optProperties.lookupOrDefault<Switch>("stageB6RxDesignOracle", false))
{
    static bool stageB6RxDesignOracleDone = false;
    if (!stageB6RxDesignOracleDone)
    {
        stageB6RxDesignOracleDone = true;
        if (discreteAdjointLabel == "pressureDrop")
        {
            #include "stageB6RxDesignOracle.H"
        }
        else
        {
            Info<< "StageB6: skipped (label=" << discreteAdjointLabel
                << "); runs in the pressureDrop adjoint block only." << endl;
        }
    }
}

// ---------------------------------------------------------------------------
// Stage B8 (BFINAL-008, PATCH-validation): actual frozen-primal residual FD
// probes (P1 null-vector P-block + P3 mixed directions) vs the corrected
// matrix-free forward J (applyDiscreteFlowJ).  Guarded by optProperties switch
// "stageB8JPPActualResidualFD" (default false).  When disabled this block is
// a no-op (one dictionary lookup only): no fields created/modified, no files
// written, no effect on production values/residuals.  When enabled it runs
// once (static guard) on the first inclusion of this header using the
// converged production state; it never modifies production state and never
// solves UEqn/pEqn (only assembles/reads residual(), A(), H(), flux() on
// local perturbed copies).
// ---------------------------------------------------------------------------
if (optProperties.lookupOrDefault<Switch>("stageB8JPPActualResidualFD", false))
{
    static bool stageB8JPPActualResidualFDDone = false;
    if (!stageB8JPPActualResidualFDDone)
    {
        stageB8JPPActualResidualFDDone = true;
        #include "stageB8JPPActualResidualFD.H"
    }
}
