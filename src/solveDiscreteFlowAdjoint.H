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
            discreteExternã^»¶‰žËkºwµçZWKœÚ^™J
NÂŸBœ™YXÙJ\ØÜ™]P\™XTØØ[KÝ[SÜØØ[\Š
JNÂœ™YXÙJ\ØÜ™]Q˜XÙPÛÝ[Ý[SÜX™[Š
JNÂ™\ØÜ™]P\™XTØØ[HÏH›Ø[NŽ›X^
ØØ[\Š\ØÜ™]Q˜XÙPÛÝ[
KØØ[\ŠKŒ
JNÂ™\ØÜ™]P\™XTØØ[HH›Ø[NŽ›X^
\ØÜ™]P\™XTØØ[KØØ[\ŠÓPS
JNÂ‚œØØ[\‘šY[\ØÜ™]TšÊ\ØÜ™]S•[šÛ›ÝÛœËŒ
NÂ™›Üˆ
X™[Ù[HHÈÙ[H\ØÜ™]SÙ[ÎÈ
ÊØÙ[JBžÂˆ›Üˆ
\™XÝ[ÛˆÛ\HÈÛ\ÎÈ
ÊØÛ\
BˆÂˆ\ØÜ™]TšÖÙ\ØÜ™]UR[™^
Ù[KÛ\
WHBˆ\ØÜ™]Q›ÝÕ™[ØÚ]TšÖØÙ[WK˜ÛÛ\Û™[
Û\
BˆÙ\ØÜ™]S[ÛY[[TØØ[NÂˆBˆ\ØÜ™]TšÖÙ\ØÜ™]T[™^
Ù[JWHBˆ\ØÜ™]Q›ÝÔ™\ÜÝ\™TšÖØÙ[WKÙ\ØÜ™]P\™XTØØ[NÂŸB‚˜]]È\TØØ[Y\ØÜ™]Q›ÝÒ•B–É—JÛÛœÝØØ[\‘šY[	ˆ[œ]ØØ[\‘šY[	ˆÝ]]
BžÂˆØØ[\‘šY[\ÚXØ[[œ]
[œ]
NÂˆ›Üˆ
X™[Ù[HHÈÙ[H\ØÜ™]SÙ[ÎÈ
ÊØÙ[JBˆÂˆ\ÚXØ[[œ]Ù\ØÜ™]T[™^
Ù[JWH
Bˆ\ØÜ™]S[ÛY[[TØØ[KÙ\ØÜ™]P\™XTØØ[NÂˆB‚ˆØØ[\‘šY[\ÚXØ[Ý]]
\ØÜ™]S•[šÛ›ÝÛœËŒ
NÂˆ\Q\ØÜ™]Q›ÝÒ•
\ÚXØ[[œ]\ÚXØ[Ý]]
NÂ‚ˆ›Üˆ
X™[Ù[HHÈÙ[H\ØÜ™]SÙ[ÎÈ
ÊØÙ[JBˆÂˆ›Üˆ
\™XÝ[ÛˆÛ\HÈÛ\ÎÈ
ÊØÛ\
BˆÂˆÝ]]Ù\ØÜ™]UR[™^
Ù[KÛ\
WHBˆ\ÚXØ[Ý]]Ù\ØÜ™]UR[™^
Ù[KÛ\
WBˆÙ\ØÜ™]S[ÛY[[TØØ[NÂˆBˆÝ]]Ù\ØÜ™]T[™^
Ù[JWHBˆ\ÚXØ[Ý]]Ù\ØÜ™]T[™^
Ù[JWKÙ\ØÜ™]P\™XTØØ[NÂˆBŸNÂ‚˜ÛÛœÝX™[\ØÜ™]SX^]\˜][ÛœÈBˆÜ›Ü\Y\Ë›ÛÚÝ\Ü‘Y˜][X™[Š™\ØÜ™]Q›ÝÐY›Ú[X^]\ˆ‹Œ
NÂ˜ÛÛœÝØØ[\ˆ\ØÜ™]UÛ\˜[˜ÙHBˆÜ›Ü\Y\Ë›ÛÚÝ\Ü‘Y˜][ØØ[\Š™\ØÜ™]Q›ÝÐY›Ú[Û\˜[˜ÙH‹YKNJNÂ‚œØØ[\‘šY[\ØÜ™]TÛÛ][ÛŠ\ØÜ™]S•[šÛ›ÝÛœËŒ
NÂœØØ[\‘šY[\ØÜ™]T™\ÚYX[
\ØÜ™]S•[šÛ›ÝÛœËŒ
NÂœØØ[\‘šY[\ØÜ™]SÜ\˜]Ü”ÛÛ][ÛŠ\ØÜ™]S•[šÛ›ÝÛœËŒ
NÂ‚˜ÛÛœÝØØ[\ˆ\ØÜ™]TšÓ›Ü›HH›Ø[NŽœÜ\Šˆ›Ø[NŽ›X^
\ØÜ™]QÛØ˜[Ý
\ØÜ™]TšË\ØÜ™]TšÊKØØ[\ŠYKLÌ
JBŠNÂ‚œØØ[\‘šY[\ØÜ™]TØÚ\‘XYÊ\ØÜ™]SÙ[ËYKN
NÂ™›Ü[
\ØÜ™]U\\‹˜XÙZJBžÂˆÛÛœÝX™[ÝÛˆH\ØÜ™]SÝÛ™\–Ù˜XÙZWNÂˆÛÛœÝX™[™ZHH\ØÜ™]S™ZYÚ›Ý\–Ù˜XÙZWNÂˆÛÛœÝØØ[\ˆÙZYÚH\ØÜ™]UÙZYÚÖÙ˜XÙZWNÂˆÛÛœÝ™XÝÜˆÙ”ØØ[YH\ØÜ™]TÙ–Ù˜XÙZWKÙ\ØÜ™]P\™XTØØ[NÂˆÛÛœÝØØ[\ˆÝÛ™\ÛÛšX][ÛˆHÜ\ŠÙZYÚ
Bˆ
ŠˆÜ\ŠÙ”ØØ[Yž

JBˆÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖÛÝÛ—KÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
Bˆ
ÈÜ\ŠÙ”ØØ[YžJ
JBˆÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖVÛÝÛ—KÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
Bˆ
ÈÜ\ŠÙ”ØØ[YžŠ
JBˆÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖ–ÛÝÛ—KÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
Bˆ
NÂˆÛÛœÝØØ[\ˆ™ZYÚ›Ý\ÛÛšX][ÛˆHÜ\ŠKŒHÙZYÚ
Bˆ
ŠˆÜ\ŠÙ”ØØ[Yž

JBˆÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖÛ™ZWKÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
Bˆ
ÈÜ\ŠÙ”ØØ[YžJ
JBˆÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖVÛ™ZWKÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
Bˆ
ÈÜ\ŠÙ”ØØ[YžŠ
JBˆÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖ–Û™ZWKÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
Bˆ
NÂˆÛÛœÝØØ[\ˆ˜XÙPÛÛšX][ÛˆBˆ›Ø[NŽ›X^
ÝÛ™\ÛÛšX][Ûˆ
È™ZYÚ›Ý\ÛÛšX][Û‹ØØ[\ŠYKLLŠJNÂˆ\ØÜ™]TØÚ\‘XYÖÛÝÛ—H
ÏH˜XÙPÛÛšX][ÛŽÂˆ\ØÜ™]TØÚ\‘XYÖÛ™ZWH
ÏH˜XÙPÛÛšX][ÛŽÂŸB‚˜]]È\Q\ØÜ™]P›ØÚÔ™XÛÛ™][Û™\ˆB–É—JÛÛœÝØØ[\‘šY[	ˆ[œ]ØØ[\‘šY[	ˆÝ]]
BžÂˆØØ[\‘šY[[™\œÙQXYÖ
\ØÜ™]SÙ[ËŒ
NÂˆØØ[\‘šY[[™\œÙQXYÖJ\ØÜ™]SÙ[ËŒ
NÂˆØØ[\‘šY[[™\œÙQXYÖŠ\ØÜ™]SÙ[ËŒ
NÂ‚ˆ›Üˆ
X™[Ù[HHÈÙ[H\ØÜ™]SÙ[ÎÈ
ÊØÙ[JBˆÂˆ[™\œÙQXYÖØÙ[WHHKŒÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖØÙ[WKÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
NÂˆ[™\œÙQXYÖVØÙ[WHHKŒÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖVØÙ[WKÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
NÂˆ[™\œÙQXYÖ–ØÙ[WHHKŒÑ›Ø[NŽ›X^ˆ
ˆXYÊ\ØÜ™]QXYÖ–ØÙ[WKÙ\ØÜ™]S[ÛY[[TØØ[JKˆØØ[\ŠYKLLŠBˆ
NÂˆÝ]]Ù\ØÜ™]UR[™^
Ù[K
WHBˆ[™\œÙQXYÖØÙ[WJš[œ]Ù\ØÜ™]UR[™^
Ù[K
WNÂˆÝ]]Ù\ØÜ™]UR[™^
Ù[KJWHBˆ[™\œÙQXYÖVØÙ[WJš[œ]Ù\ØÜ™]UR[™^
Ù[KJWNÂˆÝ]]Ù\ØÜ™]UR[™^
Ù[KŠWHBˆ[™\œÙQXYÖ–ØÙ[WJš[œ]Ù\ØÜ™]UR[™^
Ù[KŠWNÂˆB‚ˆØØ[\‘šY[™\ÜÝ\™T™\ÚYX[
\ØÜ™]SÙ[ËŒ
NÂˆ›Ü[
\ØÜ™]U\\‹˜XÙZJBˆÂˆÛÛœÝX™[ÝÛˆH\ØÜ™]SÝÛ™\–Ù˜XÙZWNÂˆÛÛœÝX™[™ZHH\ØÜ™]S™ZYÚ›Ý\–Ù˜XÙZWNÂˆ™XÝÜˆ[]]™SÝÛ‚ˆ
ˆÝ]]Ù\ØÜ™]UR[™^
ÝÛ‹
WKˆÝ]]Ù\ØÜ™]UR[™^
ÝÛ‹JWKˆÝ]]Ù\ØÜ™]UR[™^
ÝÛ‹ŠWBˆ
NÂˆ™XÝÜˆ[]]™S™ZBˆ
ˆÝ]]Ù\ØÜ™]UR[™^
™ZK
WKˆÝ]]Ù\ØÜ™]UR[™^
™ZKJWKˆÝ]]Ù\ØÜ™]UR[™^
™ZKŠWBˆ
NÂˆÛÛœÝØØ[\ˆ™\ÜÝ\™U˜[œÜÜÙHBˆ
\ØÜ™]TÙ–Ù˜XÙZWKÙ\ØÜ™]P\™XTØØ[JBˆ	ˆ
[]]™SÝÛˆH[]]™S™ZJNÂˆ™\ÜÝ\™T™\ÚYX[ÛÝÛ—H
ÏBˆ\ØÜ™]UÙZYÚÖÙ˜XÙZWJœ™\ÜÝ\™U˜[œÜÜÙNÂˆ™\ÜÝ\™T™\ÚYX[Û™ZWH
ÏBˆ
KŒH\ØÜ™]UÙZYÚÖÙ˜XÙZWJJœ™\ÜÝ\™U˜[œÜÜÙNÂˆBˆ›Ü[
Y\Ú˜›Ý[™\žJ
K]ÚJBˆÂˆYˆ
˜›Ý[™\žQšY[

VÜ]ÚWK™š^\Õ˜[YJ
JBˆÂˆÛÛ[YNÂˆBˆÛÛœÝX™[S\Ý	ˆ]ÚÙ[ÈHY\Ú˜›Ý[™\žJ
VÜ]ÚWK™˜XÙPÙ[Ê
NÂˆÛÛœÝ™XÝÜ‘šY[	ˆ]ÚÙˆHY\Ú”ÙŠ
K˜›Ý[™\žQšY[

VÜ]ÚWNÂˆ›Ü[
]ÚÙ[Ë˜XÙZJBˆÂˆÛÛœÝX™[Ù[HH]ÚÙ[ÖÙ˜XÙZWNÂˆ™XÝÜˆ[]]™PÙ[ˆ
ˆÝ]]Ù\ØÜ™]UR[™^
Ù[K
WKˆÝ]]Ù\ØÜ™]UR[™^
Ù[KJWKˆÝ]]Ù\ØÜ™]UR[™^
Ù[KŠWBˆ
NÂˆ™\ÜÝ\™T™\ÚYX[ØÙ[WH
ÏBˆ
]ÚÙ–Ù˜XÙZWKÙ\ØÜ™]P\™XTØØ[JH	ˆ[]]™PÙ[ÂˆBˆB‚ˆ›Üˆ
X™[Ù[HHÈÙ[H\ØÜ™]SÙ[ÎÈ
ÊØÙ[JBˆÂˆÝ]]Ù\ØÜ™]T[™^
Ù[JWHBˆ
[œ]Ù\ØÜ™]T[™^
Ù[JWHH™\ÜÝ\™T™\ÚYX[ØÙ[WJBˆÙ\ØÜ™]TØÚ\‘XYÖØÙ[WNÂˆB‚ˆ™XÝÜ‘šY[™\ÜÝ\™PÛÝ\[™Ê\ØÜ™]SÙ[Ë™XÝÜŽŽž™\›ÊNÂˆ›Ü[
\ØÜ™]U\\‹˜XÙZJBˆÂˆÛÛœÝX™[ÝÛˆH\ØÜ™]SÝÛ™\–Ù˜XÙZWNÂˆÛÛœÝX™[™ZHH\ØÜ™]S™ZYÚ›Ý\–Ù˜XÙZWNÂˆÛÛœÝØØ[\ˆ™\ÜÝ\™QY™™\™[˜ÙHBˆÝ]]Ù\ØÜ™]T[™^
ÝÛŠWBˆHÝ]]Ù\ØÜ™]T[™^
™ZJWNÂˆÛÛœÝ™XÝÜˆ˜XÙPÛÝ\[™ÈBˆ
\ØÜ™]TÙ–Ù˜XÙZWKÙ\ØÜ™]P\™XTØØ[JJœ™\ÜÝ\™QY™™\™[˜ÙNÂˆ™\ÜÝ\™PÛÝ\[™ÖÛÝÛ—H
ÏBˆ\ØÜ™]UÙZYÚÖÙ˜XÙZWJ™˜XÙPÛÝ\[™ÎÂˆ™\ÜÝ\™PÛÝ\[™ÖÛ™ZWH
ÏBˆ
KŒH\ØÜ™]UÙZYÚÖÙ˜XÙZWJJ™˜XÙPÛÝ\[™ÎÂˆBˆ›Ü[
Y\Ú˜›Ý[™\žJ
K]ÚJBˆÂˆYˆ
K˜›Ý[™\žQšY[

VÜ]ÚWK™š^\Õ˜[YJ
JBˆÂˆÛÛ[YNÂˆBˆÛÛœÝX™[S\Ý	ˆ]ÚÙ[ÈHY\Ú˜›Ý[™\žJ
VÜ]ÚWK™˜XÙPÙ[Ê
NÂˆÛÛœÝ™XÝÜ‘šY[	ˆ]ÚÙˆHY\Ú”ÙŠ
K˜›Ý[™\žQšY[

VÜ]ÚWNÂˆ›Ü[
]ÚÙ[Ë˜XÙZJBˆÂˆÛÛœÝX™[Ù[HH]ÚÙ[ÖÙ˜XÙZWNÂˆ™\ÜÝ\™PÛÝ\[™ÖØÙ[WH
ÏBˆ
]ÚÙ–Ù˜XÙZWKÙ\ØÜ™]P\™XTØØ[JBˆ
›Ý]]Ù\ØÜ™]T[™^
Ù[JWNÂˆBˆB‚ˆ›Üˆ
X™[Ù[HHÈÙ[H\ØÜ™]SÙ[ÎÈ
ÊØÙ[JBˆÂˆÝ]]Ù\ØÜ™]UR[™^
Ù[K
WHH[™\œÙQXYÖØÙ[WBˆ
Šˆ[œ]Ù\ØÜ™]UR[™^
Ù[K
WBˆH™\ÜÝ\™PÛÝ\[™ÖØÙ[WKž

Bˆ
NÂˆÝ]]Ù\ØÜ™]UR[™^
Ù[KJWHH[™\œÙQXYÖVØÙ[WBˆ
Šˆ[œ]Ù\ØÜ™]UR[™^
Ù[KJWBˆH™\ÜÝ\™PÛÝ\[™ÖØÙ[WKžJ
Bˆ
NÂˆÝ]]Ù\ØÜ™]UR[™^
Ù[KŠWHH[™\œÙQXYÖ–ØÙ[WBˆ
Šˆ[œ]Ù\ØÜ™]UR[™^
Ù[KŠWBˆH™\ÜÝ\™PÛÝ\[™ÖØÙ[WKžŠ
Bˆ
NÂˆBŸNÂ‚˜ÛÛœÝX™[\ØÜ™]QÛ\™\Ô™\Ý\BˆÜ›Ü\Y\Ë›ÛÚÝ\Ü‘Y˜][X™[Š™\ØÜ™]Q›ÝÑÛ\™\Ô™\Ý\‹
NÂšYˆ
\ØÜ™]QÛ\™\Ô™\Ý\ŠBžÂˆ˜][\œ›Ü’[‘[˜Ý[Û‚ˆ™\ØÜ™]Q›ÝÑÛ\™\Ô™\Ý\]\Ý™H]X\Ý‹ˆ‚ˆX›Ü
˜][\œ›ÜŠNÂŸB‚œÝŽ™XÝÜØØ[\‘šY[ˆ\ØÜ™]P˜\Ú\ÂŠˆ\ØÜ™]QÛ\™\Ô™\Ý\
ÈKˆØØ[\‘šY[
\ØÜ™]S•[šÛ›ÝÛœËŒ
BŠNÂœÝŽ™XÝÜØØ[\‘šY[ˆ\ØÜ™]T™XÛÛ™][Û™Y˜\Ú\ÂŠˆ\ØÜ™]QÛ\™\Ô™\Ý\ˆØØ[\‘šY[
\ØÜ™]S•[šÛ›ÝÛœËŒ
BŠNÂœÝŽ™XÝÜÝŽ™XÝÜØØ[\ˆ\ØÜ™]R\ÜÙ[˜™\™ÂŠˆ\ØÜ™]QÛ\™\Ô™\Ý\
ÈKˆÝŽ™XÝÜØØ[\Š\ØÜ™]QÛ\™\Ô™\Ý\Œ
BŠNÂœÝŽ™XÝÜØØ[\ˆ\ØÜ™]PÛÜÚ[™J\ØÜ™]QÛ\™\Ô™\Ý\Œ
NÂœÝŽ™XÝÜØØ[\ˆ\ØÜ™]TÚ[™J\ØÜ™]QÛ\™\Ô™\Ý\Œ
NÂœÝŽ™XÝÜØØ[\ˆ\ØÜ™]QÊ\ØÜ™]QÛ\™\Ô™\Ý\
ÈKŒ
NÂœÝŽ™XÝÜØØ[\ˆ\ØÜ™]VJ\ØÜ™]QÛ\™\Ô™\Ý\Œ
NÂ‚›X™[\ØÜ™]R]\˜][ÛˆHÂœØØ[\ˆ\ØÜ™]T™[]]™T™\ÚYX[HÔ‘PUÂ‚Ú[BŠˆ\ØÜ™]R]\˜][Ûˆ\ØÜ™]SX^]\˜][ÛœÂˆ	‰ˆ\ØÜ™]T™[]]™T™\ÚYX[ˆ\ØÜ™]UÛ\˜[˜ÙBŠBžÂˆ\TØØ[Y\ØÜ™]Q›ÝÒ•ˆ
ˆ\ØÜ™]TÛÛ][Û‹ˆ\ØÜ™]SÜ\˜]Ü”ÛÛ][Û‚ˆ
NÂˆ\ØÜ™]T™\ÚYX[H\ØÜ™]TšÈH\ØÜ™]SÜ\˜]Ü”ÛÛ][ÛŽÂˆÛÛœÝØØ[\ˆ\ØÜ™]T™\ÚYX[›Ü›HH›Ø[NŽœÜ\ˆ
ˆ\ØÜ™]QÛØ˜[Ý
\ØÜ™]T™\ÚYX[\ØÜ™]T™\ÚYX[
Bˆ
NÂˆ\ØÜ™]T™[]]™T™\ÚYX[H\ØÜ™]T™\ÚYX[›Ü›KÙ\ØÜ™]TšÓ›Ü›NÂˆYˆ
\ØÜ™]T™[]]™T™\ÚYX[H\ØÜ™]UÛ\˜[˜ÙJBˆÂˆœ™XZÎÂˆB‚ˆ\ØÜ™]P˜\Ú\ÖÌHH\ØÜ™]T™\ÚYX[Ù\ØÜ™]T™\ÚYX[›Ü›NÂˆ›Üˆ
X™[›ÝÚHHÈ›ÝÚHH\ØÜ™]QÛ\™\Ô™\Ý\È
ÊÜ›ÝÚJBˆÂˆ\ØÜ™]QÖÜ›ÝÚWHHŒÂˆ›Üˆ
X™[ÛÛHHÈÛÛH\ØÜ™]QÛ\™\Ô™\Ý\È
ÊØÛÛJBˆÂˆ\ØÜ™]R\ÜÙ[˜™\™ÖÜ›ÝÚWVØÛÛWHHŒÂˆBˆBˆ\ØÜ™]QÖÌHH\ØÜ™]T™\ÚYX[›Ü›NÂ‚ˆX™[\ØÜ™]R[›™\ÛÝ[HÂˆ›Ü‚ˆ
ˆX™[[›™\’HHÂˆ[›™\’H\ØÜ™]QÛ\™\Ô™\Ý\ˆ	‰ˆ\ØÜ™]R]\˜][Ûˆ\ØÜ™]SX^]\˜][ÛœÎÂˆ
ÊÚ[›™\’K
ÊÙ\ØÜ™]R]\˜][Û‚ˆ
BˆÂˆ\ØÜ™]R[›™\ÛÝ[H[›™\’H
ÈNÂˆ\Q\ØÜ™]P›ØÚÔ™XÛÛ™][Û™\‚ˆ
ˆ\ØÜ™]P˜\Ú\ÖÚ[›™\’WKˆ\ØÜ™]T™XÛÛ™][Û™Y˜\Ú\ÖÚ[›™\’WBˆ
NÂ‚ˆØØ[\‘šY[\››ÛU™XÝÜŠ\ØÜ™]S•[šÛ›ÝÛœËŒ
NÂˆ\TØØ[Y\ØÜ™]Q›ÝÒ•ˆ
ˆ\ØÜ™]T™XÛÛ™][Û™Y˜\Ú\ÖÚ[›™\’WKˆ\››ÛU™XÝÜ‚ˆ
NÂ‚ˆ›Üˆ
X™[˜\Ú\ÒHHÈ˜\Ú\ÒHH[›™\’NÈ
ÊØ˜\Ú\ÒJBˆÂˆ\ØÜ™]R\ÜÙ[˜™\™ÖØ˜\Ú\ÒWVÚ[›™\’WHBˆ\ØÜ™]QÛØ˜[Ý
\››ÛU™XÝÜ‹\ØÜ™]P˜\Ú\ÖØ˜\Ú\ÒWJNÂˆ\››ÛU™XÝÜˆOBˆ\ØÜ™]R\ÜÙ[˜™\™ÖØ˜\Ú\ÒWVÚ[›™\’WBˆ
™\ØÜ™]P˜\Ú\ÖØ˜\Ú\ÒWNÂˆB‚ˆ\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’H
ÈWVÚ[›™\’WHH›Ø[NŽœÜ\ˆ
ˆ\ØÜ™]QÛØ˜[Ý
\››ÛU™XÝÜ‹\››ÛU™XÝÜŠBˆ
NÂˆYˆ
\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’H
ÈWVÚ[›™\’WHˆ”ÓPS
BˆÂˆ\ØÜ™]P˜\Ú\ÖÚ[›™\’H
ÈWHBˆ\››ÛU™XÝÜ‹Ù\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’H
ÈWVÚ[›™\’WNÂˆB‚ˆ›Üˆ
X™[›Ý][Û’HHÈ›Ý][Û’H[›™\’NÈ
ÊÜ›Ý][Û’JBˆÂˆÛÛœÝØØ[\ˆ\\‘[žHBˆ\ØÜ™]R\ÜÙ[˜™\™ÖÜ›Ý][Û’WVÚ[›™\’WNÂˆÛÛœÝØØ[\ˆÝÙ\‘[žHBˆ\ØÜ™]R\ÜÙ[˜™\™ÖÜ›Ý][Û’H
ÈWVÚ[›™\’WNÂˆ\ØÜ™]R\ÜÙ[˜™\™ÖÜ›Ý][Û’WVÚ[›™\’WHBˆ\ØÜ™]PÛÜÚ[™VÜ›Ý][Û’WJ\\‘[žBˆ
È\ØÜ™]TÚ[™VÜ›Ý][Û’WJ›ÝÙ\‘[žNÂˆ\ØÜ™]R\ÜÙ[˜™\™ÖÜ›Ý][Û’H
ÈWVÚ[›™\’WHBˆY\ØÜ™]TÚ[™VÜ›Ý][Û’WJ\\‘[žBˆ
È\ØÜ™]PÛÜÚ[™VÜ›Ý][Û’WJ›ÝÙ\‘[žNÂˆB‚ˆÛÛœÝØØ[\ˆ›Ý][Û“›Ü›HH›Ø[NŽœÜ\ˆ
ˆÜ\Š\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’WVÚ[›™\’WJBˆ
ÈÜ\Š\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’H
ÈWVÚ[›™\’WJBˆ
NÂˆYˆ
›Ý][Û“›Ü›HH”ÓPS
BˆÂˆ˜][\œ›Ü’[‘[˜Ý[Û‚ˆ‘ÓT‘TÈ\››ÛHœ™XZÙÝÛˆ[ˆˆ\ØÜ™]PY›Ú[X™[ˆˆ]]\˜][Ûˆˆ\ØÜ™]R]\˜][Û‚ˆX›Ü
˜][\œ›ÜŠNÂˆB‚ˆ\ØÜ™]PÛÜÚ[™VÚ[›™\’WHBˆ\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’WVÚ[›™\’WKÜ›Ý][Û“›Ü›NÂˆ\ØÜ™]TÚ[™VÚ[›™\’WHBˆ\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’H
ÈWVÚ[›™\’WKÜ›Ý][Û“›Ü›NÂˆ\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’WVÚ[›™\’WHH›Ý][Û“›Ü›NÂˆ\ØÜ™]R\ÜÙ[˜™\™ÖÚ[›™\’H
ÈWVÚ[›™\’WHHŒÂ‚ˆÛÛœÝØØ[\ˆÛÈH\ØÜ™]QÖÚ[›™\’WNÂˆ\ØÜ™]QÖÚ[›™\’WHH\ØÜ™]PÛÜÚ[™VÚ[›™\’WJ›ÛÎÂˆ\ØÜ™]QÖÚ[›™\’H
ÈWHHY\ØÜ™]TÚ[™VÚ[›™\’WJ›ÛÎÂˆ\ØÜ™]T™[]]™T™\ÚYX[BˆXYÊ\ØÜ™]QÖÚ[›™\’H
ÈWJKÙ\ØÜ™]TšÓ›Ü›NÂ‚ˆYˆ
\ØÜ™]T™[]]™T™\ÚYX[H\ØÜ™]UÛ\˜[˜ÙJBˆÂˆ
ÊÙ\ØÜ™]R]\˜][ÛŽÂˆœ™XZÎÂˆBˆB‚ˆ›Üˆ
X™[›ÝÚHH\ØÜ™]R[›™\ÛÝ[HNÈ›ÝÚHHÈK\›ÝÚJBˆÂˆØØ[\ˆšÑ[žHH\ØÜ™]QÖÜ›ÝÚWNÂˆ›Üˆ
X™[ÛÛHH›ÝÚH
ÈNÈÛÛH\ØÜ™]R[›™\ÛÝ[È
ÊØÛÛJBˆÂˆšÑ[žHOBˆ\ØÜ™]R\ÜÙ[˜™\™ÖÜ›ÝÚWVØÛÛWJ™\ØÜ™]VVØÛÛWNÂˆBˆ\ØÜ™]VVÜ›ÝÚWHBˆšÑ[žKÙ\ØÜ™]R\ÜÙ[˜™\™ÖÜ›ÝÚWVÜ›ÝÚWNÂˆB‚ˆ›Üˆ
X™[˜\Ú\ÒHHÈ˜\Ú\ÒH\ØÜ™]R[›™\ÛÝ[È
ÊØ˜\Ú\ÒJBˆÂˆ\ØÜ™]TÛÛ][Ûˆ
ÏBˆ\ØÜ™]VVØ˜\Ú\ÒWJ™\ØÜ™]T™XÛÛ™][Û™Y˜\Ú\ÖØ˜\Ú\ÒWNÂˆ\ØÜ™]VVØ˜\Ú\ÒWHHŒÂˆBŸB‚˜\TØØ[Y\ØÜ™]Q›ÝÒ•
\ØÜ™]TÛÛ][Û‹\ØÜ™]SÜ\˜]Ü”ÛÛ][ÛŠNÂ™\ØÜ™]T™\ÚYX[H\ØÜ™]TšÈH\ØÜ™]SÜ\˜]Ü”ÛÛ][ÛŽÂ™\ØÜ™]T™[]]™T™\ÚYX[H›Ø[NŽœÜ\Šˆ\ØÜ™]QÛØ˜[Ý
\ØÜ™]T™\ÚYX[\ØÜ™]T™\ÚYX[
BŠKÙ\ØÜ™]TšÓ›Ü›NÂ‚’[™›Ï”™YXÙY\ØÜ™]H›ÝÈY›Ú[ˆ\ØÜ™]PY›Ú[X™[ˆŽˆÓT‘TÈ]\˜][ÛœÏHˆ\ØÜ™]R]\˜][Û‚ˆ‹™[]]™H™\ÚYX[Hˆ\ØÜ™]T™[]]™T™\ÚYX[[™Â‚šYˆ
\ØÜ™]T™[]]™T™\ÚYX[ˆ\ØÜ™]UÛ\˜[˜ÙJBžÂˆ˜][\œ›Ü’[‘[˜Ý[Û‚ˆ”™YXÙY\ØÜ™]H›ÝÈY›Ú[Y›ÝÛÛ™\™ÙH›Üˆ‚ˆ\ØÜ™]PY›Ú[X™[Žˆ™\ÚYX[H‚ˆ\ØÜ™]T™[]]™T™\ÚYX[‹Û\˜[˜ÙOHˆ\ØÜ™]UÛ\˜[˜ÙBˆX›Ü
˜][\œ›ÜŠNÂŸB™›Üˆ
X™[Ù[HHÈÙ[H\ØÜ™]SÙ[ÎÈ
ÊØÙ[JBžÂˆ\ØÜ™]PY›Ú[™[ØÚ]VØÙ[WHH™XÝÜ‚ˆ
ˆ\ØÜ™]TÛÛ][Û–Ù\ØÜ™]UR[™^
Ù[K
WKˆ\ØÜ™]TÛÛ][Û–Ù\ØÜ™]UR[™^
Ù[KJWKˆ\ØÜ™]TÛÛ][Û–Ù\ØÜ™]UR[™^
Ù[KŠWBˆ
NÂˆ\ØÜ™]PY›Ú[™\ÜÝ\™VØÙ[WHBˆ\ØÜ™]TÛÛ][Û–Ù\ØÜ™]T[™^
Ù[JWBˆ
™\ØÜ™]S[ÛY[[TØØ[KÙ\ØÜ™]P\™XTØØ[NÂŸB™\ØÜ™]PY›Ú[™[ØÚ]K˜ÛÜœ™XÝ›Ý[™\žPÛÛ™][ÛœÊ
NÂ™\ØÜ™]PY›Ú[™\ÜÝ\™K˜ÛÜœ™XÝ›Ý[™\žPÛÛ™][ÛœÊ
N