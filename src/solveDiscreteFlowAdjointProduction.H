//-----------------------------------------------------------------------------
// solveDiscreteFlowAdjointProduction.H
// Production reduced discrete cold-flow adjoint for the validated frozen-RANS
// route.  This file intentionally contains only the J^T operator and a robust
// restarted FGMRES solve.  The large solveDiscreteFlowAdjoint.H remains the
// Stage-B4 diagnostic/oracle implementation.
//
// Required symbols from the including scope:
//   mesh, runTime, U, p, phi, alpha, nuEffFrozen, primalPressureMobility,
//   ransFlowModel, pRefCell, optProperties,
//   discreteFlowVelocityRhs, discreteFlowPressureRhs,
//   discreteExternalFaceFluxAdjoint,
//   discreteAdjointVelocity, discreteAdjointPressure, discreteAdjointLabel.
//-----------------------------------------------------------------------------

if (Pstream::parRun())
{
    FatalErrorInFunction
        << "The production reduced discrete flow adjoint is serial-only until "
        << "processor-patch transpose terms are implemented."
        << abort(FatalError);
}

const label prodNCells = mesh.nCells();
const label prodNVelocity = 3*prodNCells;
const label prodNUnknowns = 4*prodNCells;
const labelUList& prodOwner = mesh.owner();
const labelUList& prodNeighbour = mesh.neighbour();
const scalarField& prodWeights = mesh.weights().primitiveField();
const vectorField& prodSf = mesh.Sf().primitiveField();
const scalarField& prodPhi = phi.primitiveField();

const auto prodUIndex = [&](const label celli, const direction cmpt)
{
    return 3*celli + cmpt;
};
const auto prodPIndex = [&](const label celli)
{
    return prodNVelocity + celli;
};

// BFINAL-016 A1 gauge: one-shot mesh/frozen-field export so the frozen
// deviatoric-stress transpose (applyProdFrozenDeviatoricJT below) can be
// replicated offline (numpy/scipy) and added to the exported-matrix
// tangent solves.  Switch-gated ("stageB16MeshExport", default false),
// static-guarded, files only - no numerical effect.
if (optProperties.lookupOrDefault<Switch>("stageB16MeshExport", false))
{
    static bool b16MeshExportDone = false;
    if (!b16MeshExportDone)
    {
        b16MeshExportDone = true;
        {
            OFstream o("b16mesh_owner.mtx");
            for (label f = 0; f < mesh.nInternalFaces(); ++f)
            {
                o << prodOwner[f] << nl;
            }
        }
        {
            OFstream o("b16mesh_neighbour.mtx");
            for (label f = 0; f < mesh.nInternalFaces(); ++f)
            {
                o << prodNeighbour[f] << nl;
            }
        }
        {
            OFstream o("b16mesh_sf.mtx");
            o.precision(16);
            for (label f = 0; f < mesh.nInternalFaces(); ++f)
            {
                o << prodSf[f].x() << ' ' << prodSf[f].y() << ' '
                  << prodSf[f].z() << nl;
            }
        }
        {
            OFstream o("b16mesh_weights.mtx");
            o.precision(16);
            for (label f = 0; f < mesh.nInternalFaces(); ++f)
            {
                o << prodWeights[f] << nl;
            }
        }
        {
            OFstream o("b16mesh_V_nueff.mtx");
            o.precision(16);
            for (label c = 0; c < prodNCells; ++c)
            {
                o << mesh.V()[c] << ' ' << nuEffFrozen[c] << nl;
            }
        }
        {
            // boundary faces: cell, Sf, |Sf|, deltaCoeff, nuEff_b, U-fixed
            OFstream o("b16mesh_boundary.mtx");
            o.precision(16);
            forAll(mesh.boundary(), patchi)
            {
                const labelUList& patchCells =
                    mesh.boundary()[patchi].faceCells();
                const vectorField& patchSf =
                    mesh.Sf().boundaryField()[patchi];
                const scalarField& patchMagSf =
                    mesh.magSf().boundaryField()[patchi];
                const scalarField& patchDelta =
                    mesh.deltaCoeffs().boundaryField()[patchi];
                const scalarField& patchNuEff =
                    nuEffFrozen.boundaryField()[patchi];
                const bool velocityFixed =
                    U.boundaryField()[patchi].fixesValue();
                forAll(patchCells, facei)
                {
                    o << patchCells[facei] << ' '
                      << patchSf[facei].x() << ' '
                      << patchSf[facei].y() << ' '
                      << patchSf[facei].z() << ' '
                      << patchMagSf[facei] << ' '
                      << patchDelta[facei] << ' '
                      << patchNuEff[facei] << ' '
                      << (velocityFixed ? 1 : 0) << nl;
                }
            }
        }
        Info<< "B16MESHEXPORT: written (internalFaces="
            << mesh.nInternalFaces() << " cells=" << prodNCells << ")"
            << endl;
    }
}

// OpenFOAM only inserts a pressure reference when the pressure field has no
// value-fixing boundary.  Keep the reduced adjoint, its RHS, preconditioner,
// residual checks and solution gauge on exactly the same condition.  In the
// current fixedValue-outlet validation case this is false: every physical
// pressure row must remain active.
const bool prodPressureNeedsReference = p.needReference();
const auto pinProdPressureReference = [&](scalarField& values)
{
    if (prodPressureNeedsReference)
    {
        values[prodPIndex(pRefCell)] = 0.0;
    }
};

// Same unrelaxed frozen residual basis as the Stage-B4 validated operator.
tmp<fvVectorMatrix> tProdMomentum
(
    fvm::div(phi, U)
  - fvm::laplacian(nuEffFrozen, U)
  + fvm::Sp(alpha, U)
);
fvVectorMatrix& prodMomentum = tProdMomentum.ref();

scalarField prodDiagX(prodMomentum.diag());
scalarField prodDiagY(prodMomentum.diag());
scalarField prodDiagZ(prodMomentum.diag());
forAll(prodMomentum.internalCoeffs(), patchi)
{
    const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
    const vectorField& patchDiag = prodMomentum.internalCoeffs()[patchi];
    forAll(patchCells, facei)
    {
        const label celli = patchCells[facei];
        prodDiagX[celli] += patchDiag[facei].x();
        prodDiagY[celli] += patchDiag[facei].y();
        prodDiagZ[celli] += patchDiag[facei].z();
    }
}

const scalarField& prodLower = prodMomentum.lower();
const scalarField& prodUpper = prodMomentum.upper();
scalarField prodRAU(prodNCells, 0.0);
for (label celli = 0; celli < prodNCells; ++celli)
{
    prodRAU[celli] = 1.0/Foam::max(prodMomentum.A().ref()[celli], SMALL);
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
// This file is included in the production branch only (the diagnostic branch
// defines its own alphaRel in solveDiscreteFlowAdjoint.H).
const scalar alphaRel =
    mesh.relaxEquation("U") ? mesh.equationRelaxationFactor("U") : scalar(1.0);

// Explicit frozen deviatoric-stress transpose.  Kept identical to the
// Stage-B4 operator so the production solver uses the same validated J^T.
auto applyProdFrozenDeviatoricJT =
[&](const scalarField& input, scalarField& output)
{
    if (!ransFlowModel)
    {
        return;
    }

    tensorField stressAdjoint(prodNCells, tensor::zero);
    tensorField gradientAdjoint(prodNCells, tensor::zero);

    forAll(prodUpper, facei)
    {
        const label own = prodOwner[facei];
        const label nei = prodNeighbour[facei];
        const vector lambdaOwn
        (
            input[prodUIndex(own, 0)],
            input[prodUIndex(own, 1)],
            input[prodUIndex(own, 2)]
        );
        const vector lambdaNei
        (
            input[prodUIndex(nei, 0)],
            input[prodUIndex(nei, 1)],
            input[prodUIndex(nei, 2)]
        );
        const tensor faceStressAdjoint =
            prodSf[facei]*(lambdaNei - lambdaOwn);

        stressAdjoint[own] += prodWeights[facei]*faceStressAdjoint;
        stressAdjoint[nei] += (1.0 - prodWeights[facei])*faceStressAdjoint;
    }

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
            const label celli = patchCells[facei];
            const vector lambdaCell
            (
                input[prodUIndex(celli, 0)],
                input[prodUIndex(celli, 1)],
                input[prodUIndex(celli, 2)]
            );
            const tensor boundaryStressAdjoint = -patchSf[facei]*lambdaCell;
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
                    output[prodUIndex(celli, cmpt)] +=
                        normalGradientAdjoint.component(cmpt);
                }
            }
        }
    }

    for (label celli = 0; celli < prodNCells; ++celli)
    {
        gradientAdjoint[celli] +=
            nuEffFrozen[celli]
           *dev2(Foam::T(stressAdjoint[celli]));
    }

    forAll(prodUpper, facei)
    {
        const label own = prodOwner[facei];
        const label nei = prodNeighbour[facei];
        const vector faceVelocityAdjoint =
            prodSf[facei]
          &
            (
                gradientAdjoint[own]/Foam::max(mesh.V()[own], SMALL)
              - gradientAdjoint[nei]/Foam::max(mesh.V()[nei], SMALL)
            );
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[prodUIndex(own, cmpt)] +=
                prodWeights[facei]*faceVelocityAdjoint.component(cmpt);
            output[prodUIndex(nei, cmpt)] +=
                (1.0 - prodWeights[facei])
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
              & (gradientAdjoint[celli]/Foam::max(mesh.V()[celli], SMALL));
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                output[prodUIndex(celli, cmpt)] +=
                    boundaryVelocityAdjoint.component(cmpt);
            }
        }
    }
};

// Transpose of the pressure-correction contribution used by an explicit
// objective/residual face-flux derivative.
auto applyProdPressureFluxCorrectionTranspose =
[&](const scalarField& faceAdjoint, scalarField& output)
{
    const scalarField& deltaCoeffs = mesh.deltaCoeffs().primitiveField();
    const scalarField& faceAreas = mesh.magSf().primitiveField();

    forAll(prodUpper, facei)
    {
        const label own = prodOwner[facei];
        const label nei = prodNeighbour[facei];
        const scalar weight = prodWeights[facei];
        const scalar multiplier = faceAdjoint[facei];
        const scalar mobilityFace =
            weight*primalPressureMobility[own]
          + (1.0 - weight)*primalPressureMobility[nei];
        const scalar normalCoefficient =
            mobilityFace*deltaCoeffs[facei]*faceAreas[facei];
        output[prodPIndex(own)] -= normalCoefficient*multiplier;
        output[prodPIndex(nei)] += normalCoefficient*multiplier;
    }

    forAll(mesh.boundary(), patchi)
    {
        if (U.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const scalarField& patchDelta = mesh.deltaCoeffs().boundaryField()[patchi];
        const scalarField& patchArea = mesh.magSf().boundaryField()[patchi];
        const label patchStart = mesh.boundary()[patchi].start();
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();

        if (!pressureFixed)
        {
            continue;
        }
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            output[prodPIndex(celli)] +=
                primalPressureMobility[celli]
               *patchDelta[facei]*patchArea[facei]
               *faceAdjoint[patchStart + facei];
        }
    }
};

// Add the explicit face-flux objective/thermal-adjoint source to the physical
// U/p right-hand side before scaling.
scalarField prodExternalTranspose(prodNUnknowns, 0.0);
forAll(prodUpper, facei)
{
    const label own = prodOwner[facei];
    const label nei = prodNeighbour[facei];
    const vector faceTranspose =
        prodSf[facei]*discreteExternalFaceFluxAdjoint[facei];
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        prodExternalTranspose[prodUIndex(own, cmpt)] +=
            prodWeights[facei]*faceTranspose.component(cmpt);
        prodExternalTranspose[prodUIndex(nei, cmpt)] +=
            (1.0 - prodWeights[facei])*faceTranspose.component(cmpt);
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
            prodExternalTranspose[prodUIndex(celli, cmpt)] +=
                faceTranspose.component(cmpt);
        }
    }
}
applyProdPressureFluxCorrectionTranspose
(
    discreteExternalFaceFluxAdjoint,
    prodExternalTranspose
);
for (label celli = 0; celli < prodNCells; ++celli)
{
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        discreteFlowVelocityRhs[celli].component(cmpt) +=
            prodExternalTranspose[prodUIndex(celli, cmpt)];
    }
    discreteFlowPressureRhs[celli] +=
        prodExternalTranspose[prodPIndex(celli)];
}

// Matrix-free physical J^T.  This is the compact production counterpart of
// applyDiscreteFlowJT() in the Stage-B4 diagnostic solver.
auto applyProdFlowJT =
[&](const scalarField& input, scalarField& output)
{
    output = scalar(0.0);

    for (label celli = 0; celli < prodNCells; ++celli)
    {
        output[prodUIndex(celli, 0)] +=
            prodDiagX[celli]*input[prodUIndex(celli, 0)];
        output[prodUIndex(celli, 1)] +=
            prodDiagY[celli]*input[prodUIndex(celli, 1)];
        output[prodUIndex(celli, 2)] +=
            prodDiagZ[celli]*input[prodUIndex(celli, 2)];
    }

    scalarField hAx(prodNCells, 0.0);
    scalarField hAy(prodNCells, 0.0);
    scalarField hAz(prodNCells, 0.0);

    forAll(prodUpper, facei)
    {
        const label own = prodOwner[facei];
        const label nei = prodNeighbour[facei];
        const vector lambdaOwn
        (
            input[prodUIndex(own, 0)],
            input[prodUIndex(own, 1)],
            input[prodUIndex(own, 2)]
        );
        const vector lambdaNei
        (
            input[prodUIndex(nei, 0)],
            input[prodUIndex(nei, 1)],
            input[prodUIndex(nei, 2)]
        );
        const scalar lambdaPown = input[prodPIndex(own)];
        const scalar lambdaPnei = input[prodPIndex(nei)];

        // Frozen fvm momentum A^T off-diagonal.
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[prodUIndex(own, cmpt)] +=
                prodLower[facei]*input[prodUIndex(nei, cmpt)];
            output[prodUIndex(nei, cmpt)] +=
                prodUpper[facei]*input[prodUIndex(own, cmpt)];
        }

        // Direct momentum pressure-gradient transpose.
        const scalar pressureTranspose =
            prodSf[facei] & (lambdaOwn - lambdaNei);
        output[prodPIndex(own)] +=
            prodWeights[facei]*pressureTranspose;
        output[prodPIndex(nei)] +=
            (1.0 - prodWeights[facei])*pressureTranspose;

        // Reduced SIMPLE flux transpose (BFINAL-003 P<-U relaxation patch):
        //   forward: the P row is routed through the corrected P-U flux
        //     dphi_PU = Sf&(alphaRel*rAU*dH + (1-alphaRel)*dUbar) + kf*(p_n-p_o)
        //     while the velocityJump conv term keeps the OLD unrelaxed dphi,
        //   so the adjoint scalar splits:
        //     phiAdjointRel = convTerm + alphaRel*lambdaPdiff  (-> hA path)
        //     + direct (1-alphaRel)*lambdaPdiff*(weight*Sf)    (-> U rows)
        //   BFINAL-008 (J_PP actual-primal closure): the kf path is the exact
        //   transpose of the corrected forward P-P block (dR_P/dp = -L_int:
        //   diag +kf, off-diag -kf).  Only the lambdaPdiff coefficient flips
        //   sign (P(own) += kf*lambdaPdiff, P(nei) -= kf*lambdaPdiff); the
        //   convTerm (J_UP^T via the OLD velocityJump dphi) keeps its
        //   unchanged -kf/+kf routing; phiAdjointRel/hA/directRel (J_PU) are
        //   byte-identical.
        const vector velocityJump = U[nei] - U[own];
        const label downwind = prodPhi[facei] >= 0.0 ? nei : own;
        const vector lambdaDown =
            (downwind == nei) ? lambdaNei : lambdaOwn;
        const scalar convTerm = lambdaDown & velocityJump;
        const scalar lambdaPdiff = lambdaPown - lambdaPnei;
        const scalar phiAdjointPP = lambdaPdiff - convTerm;

        const scalar weight = prodWeights[facei];
        const scalar mobilityFace =
            weight*primalPressureMobility[own]
          + (1.0 - weight)*primalPressureMobility[nei];
        const scalar kf =
            mobilityFace
           *mesh.deltaCoeffs().primitiveField()[facei]
           *mesh.magSf().primitiveField()[facei];

        output[prodPIndex(own)] += kf*phiAdjointPP;
        output[prodPIndex(nei)] -= kf*phiAdjointPP;

        const scalar phiAdjointRel = convTerm + alphaRel*lambdaPdiff;
        const vector hOwn =
            phiAdjointRel*weight*prodRAU[own]*prodSf[facei];
        const vector hNei =
            phiAdjointRel*(1.0 - weight)*prodRAU[nei]*prodSf[facei];
        hAx[own] += hOwn.x(); hAy[own] += hOwn.y(); hAz[own] += hOwn.z();
        hAx[nei] += hNei.x(); hAy[nei] += hNei.y(); hAz[nei] += hNei.z();

        // direct relax-source transpose: U(own/nei) <- (1-alphaRel)*lambdaPdiff*Sf
        const vector directRel =
            (1.0 - alphaRel)*lambdaPdiff*prodSf[facei];
        output[prodUIndex(own, 0)] += weight*directRel.x();
        output[prodUIndex(own, 1)] += weight*directRel.y();
        output[prodUIndex(own, 2)] += weight*directRel.z();
        output[prodUIndex(nei, 0)] += (1.0 - weight)*directRel.x();
        output[prodUIndex(nei, 1)] += (1.0 - weight)*directRel.y();
        output[prodUIndex(nei, 2)] += (1.0 - weight)*directRel.z();
    }

    // Boundary reduced-flux transpose and direct U-p transpose.
    // BFINAL-003 (P<-U relaxation patch): forward boundary dphi is 0 on
    // non-assignable (fixedValue/noSlip) U patches and
    // Sf_b&(alphaRel*rAU*dH + (1-alphaRel)*dU_cell) on assignable ones;
    // the transpose mirrors that: alphaRel-scaled hB + direct
    // (1-alphaRel)*lambdaPc*patchSf term, assignable-only.
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchi];
        const bool uAssignable = U.boundaryField()[patchi].assignable();
        const bool pressureFixed = p.boundaryField()[patchi].fixesValue();
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];

            if (uAssignable)
            {
                const scalar lambdaPc = input[prodPIndex(celli)];
                const vector hB =
                    lambdaPc*(alphaRel*prodRAU[celli])*patchSf[facei];
                hAx[celli] += hB.x(); hAy[celli] += hB.y(); hAz[celli] += hB.z();

                // direct relax-source transpose on the boundary face:
                // U(cell) <- (1-alphaRel)*lambdaPc*patchSf
                const vector dB =
                    (1.0 - alphaRel)*lambdaPc*patchSf[facei];
                output[prodUIndex(celli, 0)] += dB.x();
                output[prodUIndex(celli, 1)] += dB.y();
                output[prodUIndex(celli, 2)] += dB.z();
            }

            if (!pressureFixed)
            {
                const vector lambdaCell
                (
                    input[prodUIndex(celli, 0)],
                    input[prodUIndex(celli, 1)],
                    input[prodUIndex(celli, 2)]
                );
                output[prodPIndex(celli)] += patchSf[facei] & lambdaCell;
            }

            // BFINAL-008 (J_PP actual-primal closure): fixedValue-p boundary
            // P-P diagonal = +rAtU_c*deltaCoeffs_b*|Sf|_b (production
            // fvm::laplacian(rAtU,p) addBoundaryDiag internalCoeffs, residual
            // convention dR_P/dp; boundaryCoeffs source = 0 since p_b = 0;
            // zeroGradient p patches contribute 0; no identity row).  The
            // diagonal-only term is symmetric, so J^T carries it identically.
            if (pressureFixed)
            {
                output[prodPIndex(celli)] +=
                    primalPressureMobility[celli]
                   *mesh.deltaCoeffs().boundaryField()[patchi][facei]
                   *mesh.magSf().boundaryField()[patchi][facei]
                   *input[prodPIndex(celli)];
            }
        }
    }

    // deltaH^T, including fvMatrix::H() /V semantics.
    forAll(prodUpper, facei)
    {
        const label own = prodOwner[facei];
        const label nei = prodNeighbour[facei];
        const scalar Vinvo = 1.0/Foam::max(mesh.V()[own], SMALL);
        const scalar Vinn = 1.0/Foam::max(mesh.V()[nei], SMALL);
        output[prodUIndex(nei, 0)] -= prodUpper[facei]*Vinvo*hAx[own];
        output[prodUIndex(nei, 1)] -= prodUpper[facei]*Vinvo*hAy[own];
        output[prodUIndex(nei, 2)] -= prodUpper[facei]*Vinvo*hAz[own];
        output[prodUIndex(own, 0)] -= prodLower[facei]*Vinn*hAx[nei];
        output[prodUIndex(own, 1)] -= prodLower[facei]*Vinn*hAy[nei];
        output[prodUIndex(own, 2)] -= prodLower[facei]*Vinn*hAz[nei];
    }

    // Boundary H() diagonal transpose.  Do not skip fixedValue patches.
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const vectorField internalCoeffs = prodMomentum.internalCoeffs()[patchi];
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            const scalar cav =
                (internalCoeffs[facei].x()
               + internalCoeffs[facei].y()
               + internalCoeffs[facei].z())/3.0;
            const scalar Vinvc = 1.0/Foam::max(mesh.V()[celli], SMALL);
            output[prodUIndex(celli, 0)] +=
                (-internalCoeffs[facei].x() + cav)*Vinvc*hAx[celli];
            output[prodUIndex(celli, 1)] +=
                (-internalCoeffs[facei].y() + cav)*Vinvc*hAy[celli];
            output[prodUIndex(celli, 2)] +=
                (-internalCoeffs[facei].z() + cav)*Vinvc*hAz[celli];
        }
    }

    applyProdFrozenDeviatoricJT(input, output);
};

auto prodDot = [&](const scalarField& a, const scalarField& b)
{
    scalar result = sum(a*b);
    reduce(result, sumOp<scalar>());
    return result;
};

// Scale U and p blocks to O(1), matching the Stage-B4 convention.
scalar prodMomentumScale =
    (
        sum(mag(prodDiagX))
      + sum(mag(prodDiagY))
      + sum(mag(prodDiagZ))
    )/Foam::max(3.0*prodNCells, scalar(1.0));
reduce(prodMomentumScale, sumOp<scalar>());
prodMomentumScale = Foam::max(prodMomentumScale, scalar(SMALL));

scalar prodAreaScale = sum(mesh.magSf().primitiveField());
label prodFaceCount = mesh.nInternalFaces();
forAll(mesh.boundary(), patchi)
{
    prodAreaScale += sum(mesh.magSf().boundaryField()[patchi]);
    prodFaceCount += mesh.boundary()[patchi].size();
}
reduce(prodAreaScale, sumOp<scalar>());
reduce(prodFaceCount, sumOp<label>());
prodAreaScale /= Foam::max(scalar(prodFaceCount), scalar(1.0));
prodAreaScale = Foam::max(prodAreaScale, scalar(SMALL));

scalarField prodRhs(prodNUnknowns, 0.0);
for (label celli = 0; celli < prodNCells; ++celli)
{
    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        prodRhs[prodUIndex(celli, cmpt)] =
            discreteFlowVelocityRhs[celli].component(cmpt)/prodMomentumScale;
    }
    prodRhs[prodPIndex(celli)] =
        discreteFlowPressureRhs[celli]/prodAreaScale;
}
pinProdPressureReference(prodRhs);

// ---- PROD-SOLVER DIAGNOSTICS (solve layer only) ----
{
    scalar sumRhsU = 0.0, sumRhsP = 0.0;
    scalar nU = 0.0, nP = 0.0;
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            sumRhsU += prodRhs[prodUIndex(celli, c)];
            nU += sqr(prodRhs[prodUIndex(celli, c)]);
        }
        sumRhsP += prodRhs[prodPIndex(celli)];
        nP += sqr(prodRhs[prodPIndex(celli)]);
    }
    reduce(sumRhsU, sumOp<scalar>());
    reduce(sumRhsP, sumOp<scalar>());
    reduce(nU, sumOp<scalar>());
    reduce(nP, sumOp<scalar>());
    Info<< "PRODRHS (" << discreteAdjointLabel << "): sumRhsU="
        << sumRhsU << " sumRhsP=" << sumRhsP
        << " |RhsU|L2=" << Foam::sqrt(nU)
        << " |RhsP|L2=" << Foam::sqrt(nP) << endl;
}

auto applyProdScaledJT =
[&](const scalarField& input, scalarField& output)
{
    scalarField physicalInput(input);
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        physicalInput[prodPIndex(celli)] *=
            prodMomentumScale/prodAreaScale;
    }
    scalarField physicalOutput(prodNUnknowns, 0.0);
    applyProdFlowJT(physicalInput, physicalOutput);
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            output[prodUIndex(celli, cmpt)] =
                physicalOutput[prodUIndex(celli, cmpt)]/prodMomentumScale;
        }
        output[prodPIndex(celli)] =
            physicalOutput[prodPIndex(celli)]/prodAreaScale;
    }
    if (prodPressureNeedsReference)
    {
        // Legacy all-Neumann gauge row.  FixedValue-pressure cases keep the
        // complete physical operator and do not replace any row.
        output[prodPIndex(pRefCell)] = input[prodPIndex(pRefCell)];
    }
};

// Signed Jacobi right-preconditioner.  It is deliberately simple and robust:
// no ILU factorization, no duplicate-CSR assumptions, no pivot sign loss.
scalarField prodJacobiDiag(prodNUnknowns, 0.0);
for (label celli = 0; celli < prodNCells; ++celli)
{
    prodJacobiDiag[prodUIndex(celli, 0)] = prodDiagX[celli]/prodMomentumScale;
    prodJacobiDiag[prodUIndex(celli, 1)] = prodDiagY[celli]/prodMomentumScale;
    prodJacobiDiag[prodUIndex(celli, 2)] = prodDiagZ[celli]/prodMomentumScale;
}
// BFINAL-008 (J_PP actual-primal closure): the production P-P diagonal now
// follows the corrected operator dR_P/dp = -L_prod: interior diag +kf per
// face (was -kf) plus the fixedValue-p boundary internalCoeffs term
// +rAtU_c*deltaCoeffs_b*|Sf|_b on pressure-fixed patches (boundaryCoeffs
// source = 0 since p_b = 0; zeroGradient p patches contribute 0; no identity
// row).  This keeps the Jacobi diagonal consistent with the corrected J^T.
scalarField prodPDiagPhysical(prodNCells, 0.0);
forAll(prodUpper, facei)
{
    const label own = prodOwner[facei];
    const label nei = prodNeighbour[facei];
    const scalar weight = prodWeights[facei];
    const scalar mobilityFace =
        weight*primalPressureMobility[own]
      + (1.0 - weight)*primalPressureMobility[nei];
    const scalar kf =
        mobilityFace
       *mesh.deltaCoeffs().primitiveField()[facei]
       *mesh.magSf().primitiveField()[facei];
    prodPDiagPhysical[own] += kf;
    prodPDiagPhysical[nei] += kf;
}
forAll(mesh.boundary(), patchi)
{
    if (!p.boundaryField()[patchi].fixesValue())
    {
        continue;
    }
    const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
    const scalarField& patchDelta =
        mesh.deltaCoeffs().boundaryField()[patchi];
    const scalarField& patchArea =
        mesh.magSf().boundaryField()[patchi];
    forAll(patchCells, facei)
    {
        const label celli = patchCells[facei];
        prodPDiagPhysical[celli] +=
            primalPressureMobility[celli]*patchDelta[facei]*patchArea[facei];
    }
}
for (label celli = 0; celli < prodNCells; ++celli)
{
    prodJacobiDiag[prodPIndex(celli)] =
        prodPDiagPhysical[celli]
       *prodMomentumScale/sqr(prodAreaScale);
}
if (prodPressureNeedsReference)
{
    prodJacobiDiag[prodPIndex(pRefCell)] = 1.0;
}

scalar prodDiagAverage = 0.0;
forAll(prodJacobiDiag, i)
{
    prodDiagAverage += mag(prodJacobiDiag[i]);
}
reduce(prodDiagAverage, sumOp<scalar>());
prodDiagAverage /= Foam::max(scalar(prodNUnknowns), 1.0);
const scalar prodDiagFloor =
    Foam::max(1e-10, 1e-8*Foam::max(prodDiagAverage, scalar(1e-12)));
forAll(prodJacobiDiag, i)
{
    if (mag(prodJacobiDiag[i]) < prodDiagFloor)
    {
        prodJacobiDiag[i] =
            (prodJacobiDiag[i] < 0.0 ? -prodDiagFloor : prodDiagFloor);
    }
}

// Production setup must remain O(N) or O(nnz).  The former default built an
// exact L1 row sum by applying J^T once for every unknown, which is
// O(N*nnz) and unusable at topology-optimization mesh sizes.  Use a positive
// diagonal scale by default.  The legacy construction is retained only as an
// explicit small-case diagnostic/reference option.
const word prodPreconditionerSetup =
    optProperties.lookupOrDefault<word>
    (
        "discreteProdPreconditionerSetup",
        "pressureGAMG"
    );
scalarField prodL1Diag(prodNUnknowns, 1.0);
if
(
    prodPreconditionerSetup == "diagonal"
 || prodPreconditionerSetup == "pressureGAMG"
)
{
    forAll(prodL1Diag, i)
    {
        prodL1Diag[i] += mag(prodJacobiDiag[i]);
    }
}
else if (prodPreconditionerSetup == "bruteForceL1")
{
    scalarField eCol(prodNUnknowns, 0.0);
    scalarField aCol(prodNUnknowns, 0.0);
    for (label col = 0; col < prodNUnknowns; ++col)
    {
        eCol = scalar(0.0);
        eCol[col] = 1.0;
        applyProdScaledJT(eCol, aCol);
        forAll(aCol, i)
        {
            prodL1Diag[i] += mag(aCol[i]);
        }
    }
}
else
{
    FatalErrorInFunction
        << "Unknown discreteProdPreconditionerSetup='"
        << prodPreconditionerSetup
        << "'. Valid choices are pressureGAMG, diagonal or bruteForceL1."
        << abort(FatalError);
}
{
    scalar lAvg = 0.0, mnL = GREAT, mxL = -GREAT;
    forAll(prodL1Diag, i)
    {
        lAvg += prodL1Diag[i];
        mnL = Foam::min(mnL, prodL1Diag[i]);
        mxL = Foam::max(mxL, prodL1Diag[i]);
    }
    reduce(lAvg, sumOp<scalar>());
    Info<< "PRODPRECSETUP (" << discreteAdjointLabel << "): mode="
        << prodPreconditionerSetup
        << " avg=" << lAvg/prodNUnknowns
        << " min=" << mnL << " max=" << mxL << endl;
}
// Reusable pressure-Poisson approximation to the scaled P-P block.  This
// removes the pressure elliptic condition number instead of merely rescaling
// its diagonal.  Construction is O(nnz); GAMG is used only when the selected
// mode is pressureGAMG.  The physical J^T operator is not changed.
volScalarField prodPressurePrecMobility(p);
prodPressurePrecMobility.rename
(
    "prodPressurePrecMobility_" + discreteAdjointLabel
);
const scalar prodPressurePrecScale =
    prodMomentumScale/sqr(prodAreaScale);
forAll(prodPressurePrecMobility, celli)
{
    prodPressurePrecMobility[celli] =
        primalPressureMobility[celli]*prodPressurePrecScale;
}
forAll(mesh.boundary(), patchi)
{
    const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
    forAll(patchCells, facei)
    {
        prodPressurePrecMobility.boundaryFieldRef()[patchi][facei] =
            primalPressureMobility[patchCells[facei]]
           *prodPressurePrecScale;
    }
}

volScalarField prodPressurePrecPsi(p);
prodPressurePrecPsi.rename
(
    "prodPressurePrecPsi_" + discreteAdjointLabel
);
prodPressurePrecPsi.primitiveFieldRef() = scalar(0.0);
forAll(mesh.boundary(), patchi)
{
    prodPressurePrecPsi.boundaryFieldRef()[patchi] == scalar(0.0);
}

// BFINAL-010 (a) laplacianSchemes guard.  The coefficient face
// interpolation of fvm::laplacian(volScalarField, ...) does NOT read
// interpolationSchemes: the interpolation word comes from the
// laplacianSchemes entry ITstream itself ("Gauss <interp> <snGrad>":
// "Gauss" is consumed by laplacianScheme::New, the next word selects the
// mobility face interpolation - gaussLaplacianScheme.H/laplacianScheme.H).
// The J^T operator kf assumes linear interpolation == mesh.weights(); any
// other word would silently break the matrix/operator equivalence that the
// PRODPRECGAMGCHECK decomposition below enforces.
if (prodPreconditionerSetup == "pressureGAMG")
{
    const dictionary& prodLapSchemes =
        mesh.schemesDict().subDict("laplacianSchemes");
    const wordList prodLapKeys = prodLapSchemes.toc();
    forAll(prodLapKeys, keyi)
    {
        const word& prodLapKey = prodLapKeys[keyi];
        const bool prodLapIsDefault = (prodLapKey == "default");
        const bool prodLapIsPrecEntry =
            (prodLapKey.find("laplacian(prodPressurePrecMobility") == 0);
        if (!prodLapIsDefault && !prodLapIsPrecEntry)
        {
            continue;
        }
        const entry* prodLapEntryPtr =
            prodLapSchemes.lookupEntryPtr(prodLapKey, false, false);
        if (!prodLapEntryPtr)
        {
            continue;
        }
        ITstream& prodLapIs = prodLapEntryPtr->stream();
        const word prodLapFirstWord(prodLapIs);
        word prodLapInterpWord(prodLapFirstWord);
        if (prodLapFirstWord == "Gauss")
        {
            prodLapIs >> prodLapInterpWord;
        }
        if (prodLapInterpWord != "linear")
        {
            FatalErrorInFunction
                << "pressureGAMG preconditioner requires 'linear' as the "
                << "coefficient interpolation word of laplacianSchemes "
                << "entry '" << prodLapKey << "', but found '"
                << prodLapInterpWord << "'. The preconditioner coefficient "
                << "interpolation is taken from the laplacianSchemes entry "
                << "ITstream, while the operator J^T kf assumes linear "
                << "interpolation == mesh.weights(). Change the entry to "
                << "'Gauss linear <snGrad>' or validate the non-linear "
                << "equivalence before enabling pressureGAMG."
                << abort(FatalError);
        }
    }
    Info<< "PRODPRECGAMGSCHEME (" << discreteAdjointLabel
        << "): laplacianSchemes coefficient interpolation 'linear' "
        << "== mesh.weights()" << endl;
}

tmp<fvScalarMatrix> tProdPressurePrecMatrix
(
   -fvm::laplacian(prodPressurePrecMobility, prodPressurePrecPsi)
);
fvScalarMatrix& prodPressurePrecMatrix =
    tProdPressurePrecMatrix.ref();
if (prodPressureNeedsReference)
{
    prodPressurePrecMatrix.setReference(pRefCell, 0.0);
}

const scalar prodPressurePrecTolerance =
    optProperties.lookupOrDefault<scalar>
    (
        "discreteProdPressurePrecTolerance",
        1e-3
    );
const label prodPressurePrecMaxIter =
    optProperties.lookupOrDefault<label>
    (
        "discreteProdPressurePrecMaxIter",
        50
    );
if
(
    prodPressurePrecTolerance <= SMALL
 || prodPressurePrecMaxIter < 1
)
{
    FatalErrorInFunction
        << "Invalid pressureGAMG controls: require "
        << "discreteProdPressurePrecTolerance > 0 and "
        << "discreteProdPressurePrecMaxIter >= 1."
        << abort(FatalError);
}

// BFINAL-011 cycle-1: inner-solver selection for the pressure-Poisson
// preconditioner.  GAMG (unchanged default) reproduces the b010 baseline
// exactly; PCG+DIC is a robustness fallback for rhs classes where the GAMG
// V-cycle diverges (pressureDrop first-apply NaN in the b010 run).  The
// GAMG robustness coefficients are overridable from optProperties so the
// fix round can sweep them without rebuild cycles; defaults reproduce the
// previous behaviour (scaleCorrection unset == OF7 symmetric default true,
// nPreSweeps 0).
const word prodPressurePrecSolver =
    optProperties.lookupOrDefault<word>
    (
        "discreteProdPressurePrecSolver",
        "GAMG"
    );
if
(
    prodPressurePrecSolver != "GAMG"
 && prodPressurePrecSolver != "PCG"
)
{
    FatalErrorInFunction
        << "Invalid discreteProdPressurePrecSolver '"
        << prodPressurePrecSolver << "'. Valid choices are GAMG or PCG."
        << abort(FatalError);
}
dictionary prodPressurePrecSolverDict;
if (prodPressurePrecSolver == "PCG")
{
    prodPressurePrecSolverDict.add("solver", "PCG");
    prodPressurePrecSolverDict.add("preconditioner", "DIC");
    prodPressurePrecSolverDict.add
    (
        "tolerance",
        prodPressurePrecTolerance
    );
    prodPressurePrecSolverDict.add("relTol", 0.0);
    prodPressurePrecSolverDict.add
    (
        "maxIter",
        prodPressurePrecMaxIter
    );
}
else
{
    prodPressurePrecSolverDict.add("solver", "GAMG");
    prodPressurePrecSolverDict.add
    (
        "tolerance",
        prodPressurePrecTolerance
    );
    prodPressurePrecSolverDict.add("relTol", 0.0);
    prodPressurePrecSolverDict.add
    (
        "maxIter",
        prodPressurePrecMaxIter
    );
    prodPressurePrecSolverDict.add("smoother", "GaussSeidel");
    prodPressurePrecSolverDict.add("cacheAgglomeration", true);
    prodPressurePrecSolverDict.add("agglomerator", "faceAreaPair");
    prodPressurePrecSolverDict.add("nCellsInCoarsestLevel", 10);
    prodPressurePrecSolverDict.add("mergeLevels", 1);
    prodPressurePrecSolverDict.add
    (
        "nPreSweeps",
        optProperties.lookupOrDefault<label>
        (
            "discreteProdPressurePrecNPreSweeps",
            0
        )
    );
    prodPressurePrecSolverDict.add("nPostSweeps", 2);
    prodPressurePrecSolverDict.add
    (
        "scaleCorrection",
        optProperties.lookupOrDefault<bool>
        (
            "discreteProdPressurePrecScaleCorrection",
            true
        )
    );
    prodPressurePrecSolverDict.add
    (
        "directSolveCoarsest",
        optProperties.lookupOrDefault<bool>
        (
            "discreteProdPressurePrecDirectSolveCoarsest",
            false
        )
    );
}

// BFINAL-010 (b): four-group matrix/operator equivalence decomposition.
// The preconditioner fvScalarMatrix from -fvm::laplacian(mobility, psi)
// must reproduce the scaled J^T P-P block term by term:
//   interior diag  +kf per face   (negSumDiag diag(), boundary excluded)
//   boundary diag  +mob_c*deltaCoeffs_b*|Sf|_b on fixedValue-p patches
//                   (internalCoeffs accumulation, addBoundaryDiag semantics)
//   off-diagonal  -kf             (upper()/lower())
//   effective diag interior+boundary == prodJacobiDiag[P]
// The old bare prodJacobiDiag-vs-diag() comparison mixed the boundary-less
// fvMatrix diagonal with the boundary-inclusive operator diagonal and is
// replaced by this four-group split.  The four relL2 values are reported
// before the first Krylov iteration; non-finite values or an
// effective/off-diag relL2 above 1e-8 abort the run (hard gate; thresholds
// are safety invariants and deliberately not configurable).
if (prodPreconditionerSetup == "pressureGAMG")
{
    scalarField prodPrecOpInteriorDiag(prodNCells, 0.0);
    scalarField prodPrecOpBoundaryDiag(prodNCells, 0.0);
    scalarField prodPrecOpOffDiagRef(mesh.nInternalFaces(), 0.0);
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = prodOwner[facei];
        const label nei = prodNeighbour[facei];
        const scalar weight = prodWeights[facei];
        const scalar mobilityFace =
            weight*primalPressureMobility[own]
          + (1.0 - weight)*primalPressureMobility[nei];
        const scalar kf =
            mobilityFace
           *prodPressurePrecScale
           *mesh.deltaCoeffs().primitiveField()[facei]
           *mesh.magSf().primitiveField()[facei];
        prodPrecOpInteriorDiag[own] += kf;
        prodPrecOpInteriorDiag[nei] += kf;
        prodPrecOpOffDiagRef[facei] = -kf;
    }
    forAll(mesh.boundary(), patchi)
    {
        if (!p.boundaryField()[patchi].fixesValue())
        {
            continue;
        }
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const scalarField& patchDelta =
            mesh.deltaCoeffs().boundaryField()[patchi];
        const scalarField& patchArea =
            mesh.magSf().boundaryField()[patchi];
        forAll(patchCells, facei)
        {
            const label celli = patchCells[facei];
            prodPrecOpBoundaryDiag[celli] +=
                primalPressureMobility[celli]
               *prodPressurePrecScale
               *patchDelta[facei]*patchArea[facei];
        }
    }

    // Read the matrix through a const lduMatrix& alias: the NON-const
    // lower() materialises the symmetric-storage lower array as a separate
    // copy, silently converting the matrix to asymmetric storage.  That
    // side effect (present in the b010 gate) both rejects PCG ("unknown
    // asymmetric matrix solver") and flips the GAMG scaleCorrection default
    // (true only for matrix.symmetric()).  The const overload returns the
    // upper alias for a symmetric matrix and leaves storage untouched.
    const lduMatrix& prodPrecLduMatrix = prodPressurePrecMatrix;
    const scalarField& prodPrecRawDiag = prodPrecLduMatrix.diag();
    scalarField prodPrecBoundaryDiag(prodNCells, 0.0);
    forAll(mesh.boundary(), patchi)
    {
        const scalarField& patchInternalCoeffs =
            prodPressurePrecMatrix.internalCoeffs()[patchi];
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        forAll(patchCells, facei)
        {
            prodPrecBoundaryDiag[patchCells[facei]] +=
                patchInternalCoeffs[facei];
        }
    }
    const scalarField& prodPrecUpper = prodPrecLduMatrix.upper();
    const scalarField& prodPrecLower = prodPrecLduMatrix.lower();

    scalar prodPrecInteriorDiff2 = 0.0;
    scalar prodPrecInteriorRef2 = 0.0;
    scalar prodPrecBoundaryDiff2 = 0.0;
    scalar prodPrecBoundaryRef2 = 0.0;
    scalar prodPrecEffectiveDiff2 = 0.0;
    scalar prodPrecEffectiveRef2 = 0.0;
    scalar prodPrecUpperDiff2 = 0.0;
    scalar prodPrecLowerDiff2 = 0.0;
    scalar prodPrecOffDiagRef2 = 0.0;
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        if (prodPressureNeedsReference && celli == pRefCell)
        {
            // Legacy gauge row replaced by setReference; excluded from the
            // operator/preconditioner coefficient comparison.
            continue;
        }
        const scalar opInterior = prodPrecOpInteriorDiag[celli];
        const scalar opBoundary = prodPrecOpBoundaryDiag[celli];
        const scalar precBoundary = prodPrecBoundaryDiag[celli];
        const scalar effectiveRef = prodJacobiDiag[prodPIndex(celli)];
        prodPrecInteriorDiff2 += sqr(prodPrecRawDiag[celli] - opInterior);
        prodPrecInteriorRef2 += sqr(opInterior);
        prodPrecBoundaryDiff2 += sqr(precBoundary - opBoundary);
        prodPrecBoundaryRef2 += sqr(opBoundary);
        prodPrecEffectiveDiff2 +=
            sqr(prodPrecRawDiag[celli] + precBoundary - effectiveRef);
        prodPrecEffectiveRef2 += sqr(effectiveRef);
    }
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        prodPrecUpperDiff2 +=
            sqr(prodPrecUpper[facei] - prodPrecOpOffDiagRef[facei]);
        prodPrecLowerDiff2 +=
            sqr(prodPrecLower[facei] - prodPrecOpOffDiagRef[facei]);
        prodPrecOffDiagRef2 += sqr(prodPrecOpOffDiagRef[facei]);
    }
    reduce(prodPrecInteriorDiff2, sumOp<scalar>());
    reduce(prodPrecInteriorRef2, sumOp<scalar>());
    reduce(prodPrecBoundaryDiff2, sumOp<scalar>());
    reduce(prodPrecBoundaryRef2, sumOp<scalar>());
    reduce(prodPrecEffectiveDiff2, sumOp<scalar>());
    reduce(prodPrecEffectiveRef2, sumOp<scalar>());
    reduce(prodPrecUpperDiff2, sumOp<scalar>());
    reduce(prodPrecLowerDiff2, sumOp<scalar>());
    reduce(prodPrecOffDiagRef2, sumOp<scalar>());

    const scalar interiorDiagRelL2 =
        Foam::sqrt
        (
            prodPrecInteriorDiff2
           /Foam::max(prodPrecInteriorRef2, scalar(1e-300))
        );
    const bool prodPrecHasFixedPBoundary =
        (prodPrecBoundaryRef2 > 0.0);
    const scalar boundaryDiagRelL2 =
        prodPrecHasFixedPBoundary
      ? Foam::sqrt
        (
            prodPrecBoundaryDiff2
           /Foam::max(prodPrecBoundaryRef2, scalar(1e-300))
        )
      : scalar(0.0);
    const scalar effectiveDiagRelL2 =
        Foam::sqrt
        (
            prodPrecEffectiveDiff2
           /Foam::max(prodPrecEffectiveRef2, scalar(1e-300))
        );
    const scalar offDiagRelL2 =
        Foam::sqrt
        (
            Foam::max(prodPrecUpperDiff2, prodPrecLowerDiff2)
           /Foam::max(prodPrecOffDiagRef2, scalar(1e-300))
        );

    Info<< "PRODPRECGAMGSETUP (" << discreteAdjointLabel
        << "): mode=pressureGAMG solver=" << prodPressurePrecSolver
        << " innerTol=" << prodPressurePrecTolerance
        << " innerMaxIter=" << prodPressurePrecMaxIter << endl;
    Info<< "PRODPRECGAMGCHECK (" << discreteAdjointLabel
        << "): interiorDiagRelL2=" << interiorDiagRelL2
        << " boundaryDiagRelL2=" << boundaryDiagRelL2
        << (prodPrecHasFixedPBoundary ? "" : " (no fixedValue-p boundary)")
        << " effectiveDiagRelL2=" << effectiveDiagRelL2
        << " offDiagRelL2=" << offDiagRelL2 << endl;

    auto prodPrecIsFinite = [](const scalar v)
    {
        return (v == v) && (Foam::mag(v) < VGREAT);
    };
    if
    (
        !prodPrecIsFinite(interiorDiagRelL2)
     || !prodPrecIsFinite(boundaryDiagRelL2)
     || !prodPrecIsFinite(effectiveDiagRelL2)
     || !prodPrecIsFinite(offDiagRelL2)
     || effectiveDiagRelL2 > 1e-8
     || offDiagRelL2 > 1e-8
    )
    {
        FatalErrorInFunction
            << "pressureGAMG preconditioner matrix is NOT equivalent to "
            << "the scaled J^T P-P block for label '"
            << discreteAdjointLabel << "': interiorDiagRelL2="
            << interiorDiagRelL2 << " boundaryDiagRelL2="
            << boundaryDiagRelL2 << " effectiveDiagRelL2="
            << effectiveDiagRelL2 << " offDiagRelL2=" << offDiagRelL2
            << ". Locate the failing group: interior diag -> "
            << "-fvm::laplacian internal-face kf assembly; boundary diag "
            << "-> fixedValue-p internalCoeffs (addBoundaryDiag) vs "
            << "+mob_c*deltaCoeffs_b*|Sf|_b; effective diag -> "
            << "prodJacobiDiag[P] coverage; off-diag -> upper()/lower() "
            << "vs -kf. Also check that the laplacianSchemes coefficient "
            << "interpolation word is 'linear'."
            << abort(FatalError);
    }
}

// Block lower-triangular right preconditioner (simplified Schur complement):
//   M = [D_U  UP; 0  S],  z = M^{-1} r :
//     zP = S^{-1} rP                       (S ~ L1 row-sum of the P block)
//     zU = D_U^{-1} (rU - UP zP)           (D_U ~ L1 of the U block)
// UP zP is evaluated with the matrix-free operator applied to [0; zP].
// Solve-layer only; applyProdFlowJT is unchanged.
label prodPressurePrecApplyCount = 0;
auto applyProdPreconditioner =
[&](const scalarField& input, scalarField& output)
{
    output.setSize(prodNUnknowns);
    if (prodPreconditionerSetup == "pressureGAMG")
    {
        prodPressurePrecPsi.primitiveFieldRef() = scalar(0.0);
        forAll(mesh.boundary(), patchi)
        {
            prodPressurePrecPsi.boundaryFieldRef()[patchi] == scalar(0.0);
        }
        prodPressurePrecMatrix.source() = scalar(0.0);
        for (label celli = 0; celli < prodNCells; ++celli)
        {
            prodPressurePrecMatrix.source()[celli] =
                input[prodPIndex(celli)];
        }
        if (prodPressureNeedsReference)
        {
            prodPressurePrecMatrix.source()[pRefCell] = 0.0;
        }
        prodPressurePrecMatrix.solve(prodPressurePrecSolverDict);
        // BFINAL-011 cycle-1: fail fast on a non-finite inner solve instead
        // of letting nan propagate silently through FGMRES into the
        // gradient layer (b010 reached the MMA gate three stages later).
        // Pure fast-fail instrumentation; no numerical change.
        ++prodPressurePrecApplyCount;
        label prodPrecPsiNonFinite = 0;
        for (label celli = 0; celli < prodNCells; ++celli)
        {
            const scalar v = prodPressurePrecPsi[celli];
            if ((v != v) || (Foam::mag(v) > VGREAT))
            {
                ++prodPrecPsiNonFinite;
            }
        }
        reduce(prodPrecPsiNonFinite, sumOp<label>());
        if (prodPrecPsiNonFinite > 0)
        {
            FatalErrorInFunction
                << "Inner pressure-Poisson preconditioner solve produced a "
                << "non-finite psi for label '" << discreteAdjointLabel
                << "' at preconditioner apply "
                << prodPressurePrecApplyCount
                << " (inner solver '" << prodPressurePrecSolver
                << "', innerTol=" << prodPressurePrecTolerance
                << ", innerMaxIter=" << prodPressurePrecMaxIter
                << "; non-finite cells: " << prodPrecPsiNonFinite
                << "). The preconditioner matrix passed the "
                << "PRODPRECGAMGCHECK equivalence gate, so this is an inner "
                << "solver robustness failure on this rhs class (compare "
                << "discreteProdPressurePrecSolver PCG, or GAMG "
                << "scaleCorrection/nPreSweeps overrides), not a matrix "
                << "assembly defect."
                << abort(FatalError);
        }
        for (label celli = 0; celli < prodNCells; ++celli)
        {
            output[prodPIndex(celli)] =
                prodPressurePrecPsi[celli];
        }
    }
    else
    {
        for (label celli = 0; celli < prodNCells; ++celli)
        {
            output[prodPIndex(celli)] =
                input[prodPIndex(celli)]
               /prodL1Diag[prodPIndex(celli)];
        }
    }
    scalarField zIn(prodNUnknowns, 0.0);
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        zIn[prodPIndex(celli)] = output[prodPIndex(celli)];
    }
    scalarField upZ(prodNUnknowns, 0.0);
    applyProdScaledJT(zIn, upZ);
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            const label dof = prodUIndex(celli, c);
            output[dof] =
                (input[dof] - upZ[dof])/prodL1Diag[dof];
        }
    }
    if (prodPressureNeedsReference)
    {
        output[prodPIndex(pRefCell)] = input[prodPIndex(pRefCell)];
    }
};

{
    scalar jU = 0.0, jP = 0.0, mnJ = GREAT, mxJ = -GREAT;
    label jZero = 0;
    forAll(prodJacobiDiag, i)
    {
        const scalar v = mag(prodJacobiDiag[i]);
        jU += (i < prodNVelocity ? v : 0.0);
        jP += (i >= prodNVelocity ? v : 0.0);
        mnJ = Foam::min(mnJ, v);
        mxJ = Foam::max(mxJ, v);
        if (v < 1e-14) ++jZero;
    }
    reduce(jU, sumOp<scalar>());
    reduce(jP, sumOp<scalar>());
    Info<< "PRODDIAG (" << discreteAdjointLabel << "): jacobiUavg="
        << jU/Foam::max(3.0*prodNCells, 1.0)
        << " jacobiPavg=" << jP/Foam::max(scalar(prodNCells), 1.0)
        << " jacobiMin=" << mnJ << " jacobiMax=" << mxJ
        << " jacobiNearZero=" << jZero << endl;
}

// Solvability projection is only valid for an all-Neumann pressure field.
// A fixedValue pressure boundary removes the constant mode and its physical
// RHS must not be mean-projected.
if (prodPressureNeedsReference)
{
    scalarField vOnes(prodNUnknowns, 0.0);
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        vOnes[prodPIndex(celli)] = 1.0;
    }
    scalarField zOnes(prodNUnknowns, 0.0);
    applyProdScaledJT(vOnes, zOnes);
    scalar sumZ = 0.0, nZ = 0.0;
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        if (celli == pRefCell) continue;
        sumZ += zOnes[prodPIndex(celli)];
        nZ += sqr(zOnes[prodPIndex(celli)]);
    }
    reduce(sumZ, sumOp<scalar>());
    reduce(nZ, sumOp<scalar>());
    scalar rhsMean = 0.0;
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        if (celli == pRefCell) continue;
        rhsMean += prodRhs[prodPIndex(celli)];
    }
    reduce(rhsMean, sumOp<scalar>());
    rhsMean /= Foam::max(scalar(prodNCells - 1), 1.0);
    Info<< "PRODSOLV (" << discreteAdjointLabel << "): sum(J^T*1)_P="
        << sumZ << " |J^T*1_P|L2=" << Foam::sqrt(nZ)
        << " rhsPmean=" << rhsMean << endl;
    if (mag(sumZ) <= 1e-4*(mag(rhsMean)*prodNCells + SMALL))
    {
        for (label celli = 0; celli < prodNCells; ++celli)
        {
            if (celli == pRefCell) continue;
            prodRhs[prodPIndex(celli)] -= rhsMean;
        }
        pinProdPressureReference(prodRhs);
        Info<< "PRODSOLV (" << discreteAdjointLabel
            << "): removed P-mean from RHS (null-space projection)." << endl;
    }
}

const label prodRestart =
    optProperties.lookupOrDefault<label>("discreteFlowAdjointRestart", 80);
const label prodMaxIterations =
    optProperties.lookupOrDefault<label>("discreteFlowAdjointMaxIter", 4000);
const scalar prodTolerance =
    optProperties.lookupOrDefault<scalar>("discreteFlowAdjointTolerance", 1e-9);
const Switch discreteCheckpointGradient =
    optProperties.lookupOrDefault<Switch>("discreteCheckpointGradient", false);
const Switch prodEnableRitzPilot =
    optProperties.lookupOrDefault<Switch>
    (
        "discreteProdEnableRitzPilot",
        false
    );

if (prodRestart < 2 || prodMaxIterations < 1 || prodTolerance <= SMALL)
{
    FatalErrorInFunction
        << "Invalid production discrete-adjoint controls: require restart>=2, "
        << "maxIter>=1 and tolerance>0."
        << abort(FatalError);
}

scalarField prodSolution(prodNUnknowns, 0.0);
scalarField prodResidual(prodNUnknowns, 0.0);
scalarField prodOperator(prodNUnknowns, 0.0);
pinProdPressureReference(prodSolution);

// Gradient-convergence checkpoint proxy (diagnostic only).  The flow-adjoint
// contribution to dJ/dx and d(gDP)/dx is linear in lambda_U through
//   -dAlphaDxh * (U & lambda_U) * V
// Since dAlphaDxh and V are fixed during the adjoint solve, the functional
//   proxy = sum((U & lambda_U) * Dprobe)
// converges identically to the true flow-adjoint gradient.  Dprobe mirrors
// the Stage-B2 D1 physical direction (normalized, design-masked).
scalarField prodProbeD(prodNCells, 0.0);
{
    const vectorField& C = mesh.C().primitiveField();
    const scalar Lx = gMax(C.component(0)) - gMin(C.component(0));
    const scalar Ly = gMax(C.component(1)) - gMin(C.component(1));
    const scalar Lz = gMax(C.component(2)) - gMin(C.component(2));
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        if (designMask[celli] <= 0.5)
        {
            continue;
        }
        const scalar xx = C[celli].x()/Foam::max(Lx, SMALL);
        const scalar yy = C[celli].y()/Foam::max(Ly, SMALL);
        const scalar zz = C[celli].z()/Foam::max(Lz, SMALL);
        prodProbeD[celli] =
            Foam::sin(2.0*Foam::constant::mathematical::pi*xx)
           *Foam::cos(Foam::constant::mathematical::pi*yy)
          + 0.5*Foam::sin(3.0*Foam::constant::mathematical::pi*zz);
        if (x[celli] < 0.02 || x[celli] > 0.98)
        {
            prodProbeD[celli] = 0.0;
        }
    }
    scalar mx = 0.0;
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        mx = Foam::max(mx, mag(prodProbeD[celli]));
    }
    reduce(mx, maxOp<scalar>());
    if (mx > 1e-30)
    {
        for (label celli = 0; celli < prodNCells; ++celli)
        {
            prodProbeD[celli] /= mx;
        }
    }
}

auto prodGradProxy = [&](const scalarField& sol) -> scalar
{
    scalar proxy = 0.0;
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        const vector lamU
        (
            sol[prodUIndex(celli, 0)],
            sol[prodUIndex(celli, 1)],
            sol[prodUIndex(celli, 2)]
        );
        proxy += (U[celli] & lamU)*prodProbeD[celli];
    }
    reduce(proxy, sumOp<scalar>());
    return proxy;
};

Foam::OFstream ckptStream
(
    runTime.path()/fileName("adjointCheckpoint_" + discreteAdjointLabel + ".tsv")
);
if (discreteCheckpointGradient)
{
    ckptStream.precision(16);
    ckptStream << "iter\ttrueRelRes\tgradProxy\n";
}
// Near-null deflation state: if the residual is dominated by a single
// near-eigenvector (cos(Au,u) ~ 1), deflate that direction explicitly.
scalarField prodDeflationVec(prodNUnknowns, 0.0);
bool prodDeflationActive = false;

const scalar prodRhsNorm = Foam::sqrt
(
    Foam::max(prodDot(prodRhs, prodRhs), scalar(1e-300))
);
scalar discreteRelativeResidual = 0.0;
label prodIteration = 0;

const word prodSolverType =
    optProperties.lookupOrDefault<word>("discreteProdSolverType", "fgmres");
bool prodSolved = false;

if (prodSolverType == "bicgstab" && prodRhsNorm > 1e-145)
{
    // Right-preconditioned BiCGSTAB: no restart information loss, suitable
    // for the strongly non-normal U-p system where restarted FGMRES stalls.
    scalarField r(prodNUnknowns, 0.0);
    scalarField rhat(prodNUnknowns, 0.0);
    scalarField p(prodNUnknowns, 0.0);
    scalarField v(prodNUnknowns, 0.0);
    scalarField s(prodNUnknowns, 0.0);
    scalarField t(prodNUnknowns, 0.0);
    scalarField ph(prodNUnknowns, 0.0);
    scalarField sh(prodNUnknowns, 0.0);
    applyProdScaledJT(prodSolution, prodOperator);
    r = prodRhs - prodOperator;
    pinProdPressureReference(r);
    rhat = r;
    scalar rhoOld = 1.0, alpha = 1.0, omega = 1.0;
    for (label it = 0; it < prodMaxIterations; ++it)
    {
        const scalar rho = prodDot(rhat, r);
        if (mag(rho) < 1e-30)
        {
            Info<< "  [BCGS-PROD " << discreteAdjointLabel
                << "] breakdown rho at iter=" << it << endl;
            break;
        }
        const scalar beta =
            (rho/rhoOld)*(alpha/omega);
        forAll(p, i)
        {
            p[i] = r[i] + beta*(p[i] - omega*v[i]);
        }
        applyProdPreconditioner(p, ph);
        applyProdScaledJT(ph, v);
        const scalar den = prodDot(rhat, v);
        if (mag(den) < 1e-30)
        {
            Info<< "  [BCGS-PROD " << discreteAdjointLabel
                << "] breakdown den at iter=" << it << endl;
            break;
        }
        alpha = rho/den;
        forAll(s, i)
        {
            s[i] = r[i] - alpha*v[i];
        }
        pinProdPressureReference(s);
        const scalar sNorm = Foam::sqrt(prodDot(s, s));
        discreteRelativeResidual = sNorm/prodRhsNorm;
        if (discreteRelativeResidual <= prodTolerance)
        {
            forAll(prodSolution, i)
            {
                prodSolution[i] += alpha*ph[i];
            }
            pinProdPressureReference(prodSolution);
            prodIteration = it + 1;
            prodSolved = true;
            break;
        }
        applyProdPreconditioner(s, sh);
        applyProdScaledJT(sh, t);
        const scalar t2 = prodDot(t, t);
        if (t2 < 1e-30)
        {
            forAll(prodSolution, i)
            {
                prodSolution[i] += alpha*ph[i];
            }
            pinProdPressureReference(prodSolution);
            prodIteration = it + 1;
            prodSolved = true;
            break;
        }
        omega = prodDot(t, s)/t2;
        forAll(prodSolution, i)
        {
            prodSolution[i] += alpha*ph[i] + omega*sh[i];
        }
        pinProdPressureReference(prodSolution);
        forAll(r, i)
        {
            r[i] = s[i] - omega*t[i];
        }
        pinProdPressureReference(r);
        const scalar rNorm = Foam::sqrt(prodDot(r, r));
        discreteRelativeResidual = rNorm/prodRhsNorm;
        rhoOld = rho;
        if (it % 80 == 0)
        {
            Info<< "  [BCGS-PROD " << discreteAdjointLabel << "] iter="
                << it + 1 << " trueRelRes=" << discreteRelativeResidual
                << endl;
        }
        if (discreteRelativeResidual <= prodTolerance)
        {
            prodIteration = it + 1;
            prodSolved = true;
            break;
        }
    }
    if (!prodSolved)
    {
        Info<< "  [BCGS-PROD " << discreteAdjointLabel
            << "] maxIter reached, trueRelRes=" << discreteRelativeResidual
            << endl;
    }
}

if (prodRhsNorm > 1e-145 && !prodSolved)
{
    std::vector<scalarField> V(prodRestart + 1);
    std::vector<scalarField> Z(prodRestart);
    for (label i = 0; i <= prodRestart; ++i)
    {
        V[i].setSize(prodNUnknowns, 0.0);
    }
    for (label i = 0; i < prodRestart; ++i)
    {
        Z[i].setSize(prodNUnknowns, 0.0);
    }
    std::vector<std::vector<scalar>> H
    (
        prodRestart + 1,
        std::vector<scalar>(prodRestart, 0.0)
    );
    std::vector<scalar> C(prodRestart, 0.0);
    std::vector<scalar> S(prodRestart, 0.0);
    std::vector<scalar> G(prodRestart + 1, 0.0);
    std::vector<scalar> Y(prodRestart, 0.0);
    scalarField W(prodNUnknowns, 0.0);

    discreteRelativeResidual = GREAT;
    while
    (
        prodIteration < prodMaxIterations
     && discreteRelativeResidual > prodTolerance
    )
    {
        applyProdScaledJT(prodSolution, prodOperator);
        prodResidual = prodRhs - prodOperator;
        pinProdPressureReference(prodResidual);
        const scalar beta = Foam::sqrt(prodDot(prodResidual, prodResidual));
        discreteRelativeResidual = beta/prodRhsNorm;
        if (discreteRelativeResidual <= prodTolerance)
        {
            break;
        }

        // Near-null deflation: first cycle estimates the dominant residual
        // mode and checks cos(Au,u); later cycles deflate it explicitly.
        if (prodIteration >= prodRestart)
        {
            if (!prodDeflationActive)
            {
                prodDeflationVec = prodResidual/beta;
                scalarField au(prodNUnknowns, 0.0);
                applyProdScaledJT(prodDeflationVec, au);
                const scalar cosAU = prodDot(au, prodDeflationVec)
                   /Foam::max
                    (
                        Foam::sqrt(prodDot(au, au)),
                        scalar(1e-300)
                    );
                Info<< "PRODDEFL (" << discreteAdjointLabel << ") cos(Au,u)="
                    << cosAU << " at iter=" << prodIteration << endl;
                if (mag(cosAU) > 0.9)
                {
                    prodDeflationActive = true;
                    Info<< "PRODDEFL (" << discreteAdjointLabel
                        << "): near-null deflation ENABLED." << endl;
                }
            }
            if (prodDeflationActive)
            {
                scalarField au(prodNUnknowns, 0.0);
                applyProdScaledJT(prodDeflationVec, au);
                const scalar uAu = prodDot(prodDeflationVec, au);
                const scalar uR = prodDot(prodDeflationVec, prodResidual);
                if (mag(uAu) > 1e-30)
                {
                    const scalar cDef = uR/uAu;
                    forAll(prodSolution, dofi)
                    {
                        prodSolution[dofi] +=
                            cDef*prodDeflationVec[dofi];
                    }
                    pinProdPressureReference(prodSolution);
                    // refresh the residual after deflation
                    applyProdScaledJT(prodSolution, prodOperator);
                    prodResidual = prodRhs - prodOperator;
                    pinProdPressureReference(prodResidual);
                }
            }
        }

        const scalar betaDef = Foam::sqrt(prodDot(prodResidual, prodResidual));
        discreteRelativeResidual = betaDef/prodRhsNorm;
        if (discreteRelativeResidual <= prodTolerance)
        {
            break;
        }

        V[0] = prodResidual/betaDef;
        for (label i = 0; i <= prodRestart; ++i)
        {
            G[i] = 0.0;
            if (i < prodRestart)
            {
                C[i] = 0.0;
                S[i] = 0.0;
                Y[i] = 0.0;
            }
        }
        G[0] = beta;
        for (label i = 0; i <= prodRestart; ++i)
        {
            for (label j = 0; j < prodRestart; ++j)
            {
                H[i][j] = 0.0;
            }
        }

        label innerDone = 0;
        for
        (
            label k = 0;
            k < prodRestart && prodIteration + k < prodMaxIterations;
            ++k
        )
        {
            applyProdPreconditioner(V[k], Z[k]);
            applyProdScaledJT(Z[k], W);

            // Two-pass modified Gram-Schmidt improves robustness for the
            // strongly non-normal U-p operator.
            for (label pass = 0; pass < 2; ++pass)
            {
                for (label i = 0; i <= k; ++i)
                {
                    const scalar hij = prodDot(W, V[i]);
                    H[i][k] += hij;
                    forAll(W, dofi)
                    {
                        W[dofi] -= hij*V[i][dofi];
                    }
                }
            }

            const scalar hk1 = Foam::sqrt(prodDot(W, W));
            H[k + 1][k] = hk1;
            if (hk1 > 1e-14)
            {
                V[k + 1] = W/hk1;
            }

            for (label i = 0; i < k; ++i)
            {
                const scalar tmp =
                    C[i]*H[i][k] + S[i]*H[i + 1][k];
                H[i + 1][k] =
                   -S[i]*H[i][k] + C[i]*H[i + 1][k];
                H[i][k] = tmp;
            }

            const scalar rot = Foam::sqrt
            (
                sqr(H[k][k]) + sqr(H[k + 1][k])
            );
            if (rot <= 1e-30)
            {
                innerDone = k + 1;
                break;
            }
            C[k] = H[k][k]/rot;
            S[k] = H[k + 1][k]/rot;
            H[k][k] = rot;
            H[k + 1][k] = 0.0;

            const scalar gk = G[k];
            G[k] = C[k]*gk;
            G[k + 1] = -S[k]*gk;
            innerDone = k + 1;

            const scalar projectedResidual =
                mag(G[k + 1])/prodRhsNorm;
            if
            (
                (prodIteration < 640 || (prodIteration % 1600) < 160)
             && (k <= 12 || projectedResidual < 0.1)
            )
            {
                Info<< "  [FGMRES-TRACE " << discreteAdjointLabel
                    << "] cycleIter=" << prodIteration << " k=" << k
                    << " projRes=" << projectedResidual
                    << " h00=" << H[0][0] << endl;
            }
            if
            (
                projectedResidual <= prodTolerance
             || hk1 <= 1e-14
            )
            {
                break;
            }
        }

        if (innerDone <= 0)
        {
            FatalErrorInFunction
                << "FGMRES failed to build a Krylov basis for "
                << discreteAdjointLabel << abort(FatalError);
        }

        for (label i = innerDone - 1; i >= 0; --i)
        {
            scalar rhsI = G[i];
            for (label j = i + 1; j < innerDone; ++j)
            {
                rhsI -= H[i][j]*Y[j];
            }
            const scalar pivot = H[i][i];
            if (mag(pivot) <= 1e-30)
            {
                FatalErrorInFunction
                    << "FGMRES triangular solve encountered a zero pivot at "
                    << i << " for " << discreteAdjointLabel
                    << abort(FatalError);
            }
            // Preserve the pivot sign.  Never divide by mag(pivot).
            Y[i] = rhsI/pivot;
        }

        for (label i = 0; i < innerDone; ++i)
        {
            forAll(prodSolution, dofi)
            {
                prodSolution[dofi] += Y[i]*Z[i][dofi];
            }
        }
        pinProdPressureReference(prodSolution);
        prodIteration += innerDone;

        // ---- harmonic-Ritz deflation pilot (diagnostic) ----
        if
        (
            prodEnableRitzPilot
         && (prodIteration <= 320
          || prodIteration >= prodMaxIterations - 320)
        )
        {
            const label nH = innerDone;
            // inverse-power iteration on H (upper Hessenberg) for the m
            // smallest-|theta| Ritz pairs; measure residual coverage.
            const label mRitz = 16;
            scalarField ritzTheta(mRitz, 0.0);
            scalarField ritzCover(mRitz, 0.0);
            // copy residual (preconditioned-system residual is in G, but we
            // use the physical residual direction V[0]*beta... use raw r):
            scalarField rPilot(prodNUnknowns, 0.0);
            applyProdScaledJT(prodSolution, prodOperator);
            forAll(rPilot, i)
            {
                rPilot[i] = prodRhs[i] - prodOperator[i];
            }
            pinProdPressureReference(rPilot);
            const scalar rNorm = Foam::max
            (
                Foam::sqrt(prodDot(rPilot, rPilot)),
                scalar(1e-300)
            );
            std::vector<scalarField> ritzVecs;
            for (label mi = 0; mi < mRitz; ++mi)
            {
                scalarField y(nH, 0.0);
                for (label i = 0; i < nH; ++i)
                {
                    y[i] = Foam::sin(0.7*(i + 1) + 0.31*mi);
                }
                scalar theta = 0.0;
                for (label it = 0; it < 200; ++it)
                {
                    // solve H yNew = y (back substitution, Hessenberg)
                    scalarField yNew(nH, 0.0);
                    for (label i = nH - 1; i >= 0; --i)
                    {
                        scalar s = y[i];
                        for (label j = i + 1; j < nH; ++j)
                        {
                            s -= H[i][j]*yNew[j];
                        }
                        const scalar piv = H[i][i];
                        yNew[i] = (mag(piv) > 1e-30) ? s/piv : s*1e30;
                    }
                    // orthogonalize against previous Ritz vectors
                    for (const auto& rv : ritzVecs)
                    {
                        scalar dot = 0.0;
                        for (label i = 0; i < nH; ++i)
                        {
                            dot += yNew[i]*rv[i];
                        }
                        for (label i = 0; i < nH; ++i)
                        {
                            yNew[i] -= dot*rv[i];
                        }
                    }
                    const scalar ny = Foam::sqrt(sum(sqr(yNew)));
                    if (ny < 1e-14) { yNew[0] = 1.0; }
                    else { yNew /= ny; }
                    // Rayleigh quotient
                    scalar num = 0.0, den = 0.0;
                    for (label i = 0; i < nH; ++i)
                    {
                        scalar hi = 0.0;
                        for (label j = 0; j < nH; ++j)
                        {
                            hi += H[i][j]*yNew[j];
                        }
                        num += yNew[i]*hi;
                        den += yNew[i]*yNew[i];
                    }
                    const scalar thNew = num/Foam::max(den, 1e-30);
                    const scalar dth = mag(thNew - theta);
                    theta = thNew;
                    y = yNew;
                    if (dth < 1e-8*mag(theta) + 1e-10) { break; }
                }
                ritzTheta[mi] = theta;
                // Ritz vector u = V*y ; coverage of the physical residual
                scalarField uP(prodNUnknowns, 0.0);
                for (label i = 0; i < nH; ++i)
                {
                    const scalar yi = y[i];
                    const scalarField& vi = V[i];
                    forAll(vi, dofi)
                    {
                        uP[dofi] += yi*vi[dofi];
                    }
                }
                const scalar un = Foam::sqrt(prodDot(uP, uP));
                if (un > 1e-30) { uP /= un; }
                const scalar c = prodDot(uP, rPilot)/rNorm;
                ritzCover[mi] = c;
                ritzVecs.push_back(y);
            }
            Info<< "RITZPILOT (" << discreteAdjointLabel << ") iter="
                << prodIteration << " |r|=" << rNorm
                << " thetas=";
            for (label mi = 0; mi < mRitz; ++mi)
            {
                Info<< ritzTheta[mi] << (mi < mRitz - 1 ? "," : "");
            }
            Info<< " cov=";
            for (label mi = 0; mi < mRitz; ++mi)
            {
                Info<< ritzCover[mi] << (mi < mRitz - 1 ? "," : "");
            }
            Info<< endl;
        }

        applyProdScaledJT(prodSolution, prodOperator);
        prodResidual = prodRhs - prodOperator;
        pinProdPressureReference(prodResidual);
        discreteRelativeResidual =
            Foam::sqrt(prodDot(prodResidual, prodResidual))/prodRhsNorm;

        Info<< "  [FGMRES-PROD " << discreteAdjointLabel << "] iter="
            << prodIteration << " trueRelRes="
            << discreteRelativeResidual << endl;

        if (discreteCheckpointGradient)
        {
            const scalar gproxy = prodGradProxy(prodSolution);
            ckptStream << prodIteration << '\t'
                       << discreteRelativeResidual << '\t'
                       << gproxy << '\n';
            Info<< "  [GRADPROXY " << discreteAdjointLabel << "] iter="
                << prodIteration << " trueRelRes="
                << discreteRelativeResidual << " proxy=" << gproxy << endl;
        }
    }
}

applyProdScaledJT(prodSolution, prodOperator);
prodResidual = prodRhs - prodOperator;
pinProdPressureReference(prodResidual);
discreteRelativeResidual =
    Foam::sqrt(prodDot(prodResidual, prodResidual))/prodRhsNorm;
if (discreteCheckpointGradient)
{
    const scalar gproxy = prodGradProxy(prodSolution);
    ckptStream << prodIteration << '\t'
               << discreteRelativeResidual << '\t'
               << gproxy << '\n';
    Info<< "  [GRADPROXY-FINAL " << discreteAdjointLabel << "] trueRelRes="
        << discreteRelativeResidual << " proxy=" << gproxy << endl;
}
{
    scalar nU = 0.0, nP = 0.0, mxU = 0.0, mxP = 0.0;
    label mxUCell = -1, mxPCell = -1;
    for (label celli = 0; celli < prodNCells; ++celli)
    {
        for (direction c = 0; c < 3; ++c)
        {
            const scalar v = prodResidual[prodUIndex(celli, c)];
            nU += sqr(v);
            if (mag(v) > mxU) { mxU = mag(v); mxUCell = celli; }
        }
        const scalar vp = prodResidual[prodPIndex(celli)];
        nP += sqr(vp);
        if (mag(vp) > mxP) { mxP = mag(vp); mxPCell = celli; }
    }
    reduce(nU, sumOp<scalar>());
    reduce(nP, sumOp<scalar>());
    Info<< "PRODRESID (" << discreteAdjointLabel << "): |rU|L2="
        << Foam::sqrt(nU) << " |rP|L2=" << Foam::sqrt(nP)
        << " max|rU|=" << mxU << " at cell " << mxUCell
        << " max|rP|=" << mxP << " at cell " << mxPCell << endl;
}

Info<< "Production reduced discrete flow adjoint " << discreteAdjointLabel
    << ": FGMRES iterations=" << prodIteration
    << ", true relative residual=" << discreteRelativeResidual << endl;

if
(
   !std::isfinite(discreteRelativeResidual)
 || discreteRelativeResidual > prodTolerance
)
{
    if (optProperties.lookupOrDefault<Switch>("discreteProdConvergeFatal", true))
    {
        FatalErrorInFunction
            << "Production reduced discrete flow adjoint did not converge for "
            << discreteAdjointLabel << ". trueRelRes="
            << discreteRelativeResidual << ", tolerance=" << prodTolerance
            << ". Gradients are blocked rather than using a partial adjoint."
            << abort(FatalError);
    }
    else
    {
        WarningInFunction
            << "Production reduced discrete flow adjoint did not converge for "
            << discreteAdjointLabel << " (DIAGNOSTIC MODE). trueRelRes="
            << discreteRelativeResidual << ", tolerance=" << prodTolerance
            << ". Continuing for cross-objective discrimination."
            << endl;
    }
}

for (label celli = 0; celli < prodNCells; ++celli)
{
    discreteAdjointVelocity[celli] = vector
    (
        prodSolution[prodUIndex(celli, 0)],
        prodSolution[prodUIndex(celli, 1)],
        prodSolution[prodUIndex(celli, 2)]
    );
    discreteAdjointPressure[celli] =
        prodSolution[prodPIndex(celli)]
       *prodMomentumScale/prodAreaScale;
}
discreteAdjointVelocity.correctBoundaryConditions();
discreteAdjointPressure.correctBoundaryConditions();
