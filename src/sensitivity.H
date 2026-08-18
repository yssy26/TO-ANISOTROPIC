Info<< "sensitivity analysis" << endl << flush;

volScalarField dAlphaDxh
(
    IOobject
    (
        "dAlphaDxh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   -alphaMax*(1.0 + qu)*qu
   /sqr(qu + xh + SMALL)
);

// dDTDxh is now declared in readThermalProperties.H and computed in
// updateMaterialProperties.H (tensor form: volSymmTensorField).
// The direction-derivative folding happens in AdjHeatTransfer.H, so
// thermalDiffusionDerivativeDTCell already stores dJ/dxh from thermal
// diffusion — no separate dDTDxh multiplication needed below.

// alpha is clipped by alphamin in update.H; the derivative of the clipped
// branch is zero. Fixed regions are removed by designMask.
forAll(dAlphaDxh.primitiveField(), celli)
{
    if
    (
        designMask[celli] <= 0.5
     || alpha[celli] <= alphamin.value()*(1.0 + 1e-10)
    )
    {
        dAlphaDxh[celli] = 0.0;
    }
}
dAlphaDxh.correctBoundaryConditions();

// The sign and scale of Tb are supplied by the selected thermal objective
// source in AdjHeatTransfer.H.
if (adjointMode == "discrete")
{
    // Uc fingerprint before the gDP sensitivity assembly
    if (optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false))
    {
        scalar sUc = 0.0, nUc = 0.0, mUc = 0.0;
        forAll(Uc, celli)
        {
            const vector ucv = Uc[celli];
            sUc += ucv.x() + ucv.y() + ucv.z();
            nUc += magSqr(ucv);
            mUc = Foam::max(mUc, mag(ucv));
        }
        reduce(sUc, sumOp<scalar>());
        reduce(nUc, sumOp<scalar>());
        reduce(mUc, maxOp<scalar>());
        Info<< "UcFingerprint[beforeSens]:"
            << " sum=" << sUc << " L2=" << Foam::sqrt(nUc)
            << " max=" << mUc
            << " c0=" << Uc[0]
            << " c100=" << Uc[100]
            << " c5000=" << Uc[5000] << endl;
    }
    // Brinkman term: fvm::Sp(alpha, U) assembles alpha*V_cell on the
    // diagonal, so the adjoint sensitivity dJ/dxh includes mesh.V().
    fsenshMeanT = -dAlphaDxh*(U & Ub);
    fsenshMeanT.primitiveFieldRef() *= mesh.V();
    // BFINAL-005 (PATCH_RX_PRESSURE_ROW): transparent decomposition of the
    // pressure-drop design sensitivity into the momentum-row contribution
    // (-lambda_U^T R_U,x, unchanged) and the pressure/continuity-row
    // contribution (-lambda_P^T R_P,x, new):
    //   gsenshPressureDrop_total
    //     = gsenshPressureDrop_momentum + gsenshPressureDrop_pressureRow
    // Momentum row (UNCHANGED formula):
    //   -dAlphaDxh*(U&Uc)*V  (fvm::Sp assembles alpha*V_cell on the diagonal)
    volScalarField gsenshPressureDropMomentum
    (
        IOobject
        (
            "gsenshPressureDropMomentum",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroGPDMomentum", gsenshPressureDrop.dimensions(), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    gsenshPressureDropMomentum = -dAlphaDxh*(U & Uc);
    gsenshPressureDropMomentum.primitiveFieldRef() *= mesh.V();
    // Pressure/continuity row (BFINAL-005, task §5-§10): the exact transpose
    // of the BFINAL-004-validated R_P operator applied to the pressure
    // adjoint, T = J_P^T pc, contracted with the production dAlphaDxh:
    //   gsenshPressureDropPressureRow = -T*dAlphaDxh   (NO extra V)
    // so that sum(gsenshPressureDropPressureRow*z) == -pc^T R_P,x d with
    // z = dxh/dx*d through the EXISTING filter_chainrule.H path.  The helper
    // computes T = J_P^T pc UNCONDITIONALLY (production term active with the
    // stageB6RxDesignOracle switch OFF as well); only its runtime dot-test
    // self-check and diagnostics are gated behind that switch.
    #include "rxPressureRowTranspose.H"
    volScalarField gsenshPressureDropPressureRow
    (
        IOobject
        (
            "gsenshPressureDropPressureRow",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroGPDPressRow", gsenshPressureDrop.dimensions(), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    gsenshPressureDropPressureRow = -rxPressureRowT*dAlphaDxh;
    // Total (single field feeds the existing filter_chainrule.H path).
    gsenshPressureDrop = gsenshPressureDropMomentum + gsenshPressureDropPressureRow;

    if (freezeColdFlowForValidation)
    {
        fsenshMeanT = dimensionedScalar
        (
            "zeroFrozenColdFlowSensitivity",
            fsenshMeanT.dimensions(),
            0.0
        );
        gsenshPressureDropMomentum = dimensionedScalar
        (
            "zeroFrozenPressureSensitivity",
            gsenshPressureDropMomentum.dimensions(),
            0.0
        );
        gsenshPressureDropPressureRow = dimensionedScalar
        (
            "zeroFrozenPressureRowSensitivity",
            gsenshPressureDropPressureRow.dimensions(),
            0.0
        );
        gsenshPressureDrop = dimensionedScalar
        (
            "zeroFrozenPressureSensitivity",
            gsenshPressureDrop.dimensions(),
            0.0
        );
    }

    // BFINAL-005 G1/G2 evidence: export the three decomposed xh-level
    // pressure-drop sensitivities (momentum / pressure-row / total) BEFORE
    // the projection/filter chain mutates the total field.
    if (optProperties.lookupOrDefault<Switch>("stageB6RxDesignOracle", false))
    {
        OFstream rxM("rxpr_prod_gsensh_momentum.mtx");
        OFstream rxP("rxpr_prod_gsensh_pressurerow.mtx");
        OFstream rxT("rxpr_prod_gsensh_total.mtx");
        rxM.precision(16); rxP.precision(16); rxT.precision(16);
        scalar rxNM = 0.0, rxNP = 0.0, rxNT = 0.0;
        forAll(gsenshPressureDropMomentum, celli)
        {
            rxM << gsenshPressureDropMomentum[celli] << nl;
            rxP << gsenshPressureDropPressureRow[celli] << nl;
            rxT << gsenshPressureDrop[celli] << nl;
            rxNM += sqr(gsenshPressureDropMomentum[celli]);
            rxNP += sqr(gsenshPressureDropPressureRow[celli]);
            rxNT += sqr(gsenshPressureDrop[celli]);
        }
        rxM.flush(); rxP.flush(); rxT.flush();
        reduce(rxNM, sumOp<scalar>());
        reduce(rxNP, sumOp<scalar>());
        reduce(rxNT, sumOp<scalar>());
        Info<< "RxPressureRow production xh fields: |momentum|L2="
            << Foam::sqrt(rxNM)
            << " |pressureRow|L2=" << Foam::sqrt(rxNP)
            << " |total|L2=" << Foam::sqrt(rxNT) << endl;
    }

    // thermalDiffusionDerivativeDTCell already contains the full dJ/dxh
    // from thermal diffusion (dDTDxh direction derivative folded in at
    // face level in AdjHeatTransfer.H).  It is a per-cell quantity —
    // do NOT multiply by mesh.V().
    forAll(fsenshMeanT.primitiveField(), celli)
    {
        fsenshMeanT[celli] += thermalDiffusionDerivativeDTCell[celli];
    }
}
else  // adjointMode == "continuous"
{
    // Continuous adjoint sensitivities (derived via Lagrangian dJ/dxh =
    // -<adj, dR/dxh>).  Both terms carry a negative sign:
    //   - Brinkman: dJ/dxh = -dAlphaDxh * (U . Ub) * mesh.V()
    //     (fvm::Sp assembles alpha*V_cell on the diagonal)
    //   - Thermal:  exact face-based Laplacian sensitivity (same as discrete)
    //     thermalDiffusionDerivativeDTCell is a face-sum, already per-cell.
    fsenshMeanT = -dAlphaDxh*(U & Ub);
    fsenshMeanT.primitiveFieldRef() *= mesh.V();
    gsenshPressureDrop = -dAlphaDxh*(U & Uc);
    gsenshPressureDrop.primitiveFieldRef() *= mesh.V();

    if (freezeColdFlowForValidation)
    {
        fsenshMeanT = dimensionedScalar
        (
            "zeroFrozenColdFlowSensitivity",
            fsenshMeanT.dimensions(),
            0.0
        );
        gsenshPressureDrop = dimensionedScalar
        (
            "zeroFrozenPressureSensitivity",
            gsenshPressureDrop.dimensions(),
            0.0
        );
    }

    forAll(fsenshMeanT.primitiveField(), celli)
    {
        fsenshMeanT[celli] += thermalDiffusionDerivativeDTCell[celli];
    }
}

gsenshPowerDiss =
    dAlphaDxh*(U & Ua)/Foam::max(PowerDiss0, SMALL);
gsenshPowerDiss.primitiveFieldRef() *= mesh.V();

#include "filter_chainrule.H"

// BFINAL-005 G2/G3 evidence: export the patched raw-design pressure-drop
// sensitivity (post projection/filter chain, pre MMA scaling).  The
// projection of this field onto a raw-design direction d is D_total_production
// = sum_active(gsensPressureDrop*d).
if (optProperties.lookupOrDefault<Switch>("stageB6RxDesignOracle", false))
{
    OFstream rxG("rxpr_prod_gsens.mtx");
    rxG.precision(16);
    forAll(gsensPressureDrop, celli) { rxG << gsensPressureDrop[celli] << nl; }
    rxG.flush();
}

gx[0] = V;
gx[1] = PressureDropPa/Foam::max(pressureDropMaxPa, SMALL) - 1.0;

forAll(activeDesignCells, activeI)
{
    const label celli = activeDesignCells[activeI];
    xmma[activeI] = x[celli];
    dfdx[activeI] = objectiveGradientScale
       *(
            objectiveMode == "outletTemperature"
          ? fsensMeanT[celli]
          : (1.0 - w1)*fsensMeanT[celli] + w1*gsensPowerDiss[celli]
        );
    dgdx[0][activeI] = gsensVol[celli];
    dgdx[1][activeI] =
        pressureGradientScale*gsensPressureDrop[celli];

    if
    (
       !std::isfinite(xmma[activeI])
     || !std::isfinite(dfdx[activeI])
     || !std::isfinite(dgdx[0][activeI])
     || !std::isfinite(dgdx[1][activeI])
    )
    {
        FatalErrorInFunction
            << "Non-finite MMA data at active design index " << activeI
            << ", cell " << celli << abort(FatalError);
    }
}

if (!std::isfinite(gx[0]) || !std::isfinite(gx[1]))
{
    FatalErrorInFunction
        << "Non-finite constraint value supplied to MMA."
        << abort(FatalError);
}

scalar dfdxMax = -GREAT;
scalar dfdxMin = GREAT;
scalar dg0Max = -GREAT;
scalar dg0Min = GREAT;
scalar dg1Max = -GREAT;
scalar dg1Min = GREAT;

for (label activeI = 0; activeI < nDesignLocal; ++activeI)
{
    dfdxMax = Foam::max(dfdxMax, scalar(dfdx[activeI]));
    dfdxMin = Foam::min(dfdxMin, scalar(dfdx[activeI]));
    dg0Max = Foam::max(dg0Max, scalar(dgdx[0][activeI]));
    dg0Min = Foam::min(dg0Min, scalar(dgdx[0][activeI]));
    dg1Max = Foam::max(dg1Max, scalar(dgdx[1][activeI]));
    dg1Min = Foam::min(dg1Min, scalar(dgdx[1][activeI]));
}

// Every rank participates, including ranks with zero local active cells.
reduce(dfdxMax, maxOp<scalar>());
reduce(dfdxMin, minOp<scalar>());
reduce(dg0Max, maxOp<scalar>());
reduce(dg0Min, minOp<scalar>());
reduce(dg1Max, maxOp<scalar>());
reduce(dg1Min, minOp<scalar>());

if (Pstream::master())
{
    Info<< "Active gradient ranges: dJ=[" << dfdxMin << ", "
        << dfdxMax << "], dV=[" << dg0Min << ", " << dg0Max
        << "], dDP=[" << dg1Min << ", " << dg1Max << "]" << endl << flush;
}

// These fields correspond to the current solved design. Write them before an
// MMA update creates the next, not-yet-solved design.
if (runTime.writeTime() || finalSolvedIteration)
{
    dAlphaDxh.write();
    dDTDxh.write();
    fsenshMeanT.write();
    gsenshPowerDiss.write();
    gsenshPressureDrop.write();
    fsensMeanT.write();
    gsensPowerDiss.write();
    gsensPressureDrop.write();
    gsensVol.write();
    // Stage A diagnostic: write the raw design-objective gradient dJ/dx.
    // Instrumentation only; no equation or optimizer value is changed.
    volScalarField dfdxField
    (
        IOobject
        (
            "dfdx",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroDfdx", dimless, 0.0)
    );
    forAll(activeDesignCells, activeI)
    {
        dfdxField.primitiveFieldRef()[activeDesignCells[activeI]] =
            dfdx[activeI];
    }
    dfdxField.write();

}

if (mmaUpdateEnabled && !finalSolvedIteration)
{
    // --- Phase I: Transaction mechanism ---
    // When enableSSTAcceptanceCheck is true, the MMA candidate is evaluated
    // with full SST re-convergence before being accepted.  Rejected
    // candidates trigger a rollback of design, primal, and MMA state.
    if (enableSSTAcceptanceCheck)
    {
        #include "saveAcceptedState.H"
    }

    bool candidateAccepted = !enableSSTAcceptanceCheck;
    label candidateAttempt = 0;

    do
    {
    if (enableSSTAcceptanceCheck)
    {
        candidateAttempt++;
        candidateAccepted = false;
        mma.movelimit = currentMoveLimit;

        // Re-populate xmma from accepted x (needed after rollback)
        forAll(activeDesignCells, activeI)
        {
            xmma[activeI] = x[activeDesignCells[activeI]];
        }
    }

    // Save MMA history for potential rollback
    MMA::StateSnapshot mmaTxSnapshot;
    if (enableSSTAcceptanceCheck)
    {
        mmaTxSnapshot = mma.saveState();
    }

    Info<< "run MMA on " << nDesignGlobal
        << " active global design variables" << endl << flush;

    const std::vector<double> oldXmma(xmma);
    if (enableSSTAcceptanceCheck)
    {
        mma.propose(xmma, dfdx, gx, dgdx);

        // --- Transaction test: amplify first candidate step to force rejection
        if (transactionTestMode && candidateAttempt == 1)
        {
            Info<< "  TX_TEST: Amplifying candidate step by factor "
                << transactionTestAmplification
                << " to force rejection" << endl << flush;
            forAll(activeDesignCells, activeI)
            {
                const scalar step =
                    xmma[activeI] - oldXmma[activeI];
                xmma[activeI] = oldXmma[activeI]
                  + transactionTestAmplification*step;
                // Clamp to [0, 1]
                xmma[activeI] = Foam::max(0.0, Foam::min(1.0, xmma[activeI]));
            }
        }
    }
    else
    {
        mma.MMAsolver(xmma, dfdx, gx, dgdx);
    }

    word directionOnlyUpdateMode("mma");
    if
    (
        directionOnlyOptimizationApproved
     && !gradientValidated
     && directionOnlyFeasibilityFirst
     && gx[1] > 0.0
    )
    {
        directionOnlyUpdateMode = "pressureFeasibility";
        forAll(activeDesignCells, activeI)
        {
            // Only the direction is used; the fixed-volume projection and
            // movlim normalization below set its admissible magnitude.
            xmma[activeI] = oldXmma[activeI] - dgdx[1][activeI];
        }
    }

    scalar candidateDesignFluidVolume = 0.0;
    forAll(activeDesignCells, activeI)
    {
        const label celli = activeDesignCells[activeI];
        candidateDesignFluidVolume += xmma[activeI]*mesh.V()[celli];
    }
    reduce(candidateDesignFluidVolume, sumOp<scalar>());

    // Keep the prescribed solid fraction exactly and let MMA redistribute
    // material spatially. This gate is intentionally separate from the
    // direction-only gradient approval.
    const scalar targetDesignFluidFraction = 1.0 - solidFractionMin;
    const scalar volumeRestorationShift =
        targetDesignFluidFraction
      - candidateDesignFluidVolume/Foam::max(designVolume, SMALL);

    scalar mmaCandidateMaxRawDesignChange = 0.0;
    scalar directionOnlyBoundScale = GREAT;
    forAll(activeDesignCells, activeI)
    {
        const scalar candidateChange =
            xmma[activeI] + volumeRestorationShift - oldXmma[activeI];
        mmaCandidateMaxRawDesignChange = Foam::max
        (
            mmaCandidateMaxRawDesignChange,
            mag(candidateChange)
        );

        if (candidateChange > SMALL)
        {
            directionOnlyBoundScale = Foam::min
            (
                directionOnlyBoundScale,
                (1.0 - oldXmma[activeI])/candidateChange
            );
        }
        else if (candidateChange < -SMALL)
        {
            directionOnlyBoundScale = Foam::min
            (
                directionOnlyBoundScale,
                -oldXmma[activeI]/candidateChange
            );
        }
    }
    reduce(mmaCandidateMaxRawDesignChange, maxOp<scalar>());
    reduce(directionOnlyBoundScale, minOp<scalar>());

    scalar directionOnlyStepScale = 1.0;
    if
    (
        directionOnlyOptimizationApproved
     && !gradientValidated
     && mmaCandidateMaxRawDesignChange > SMALL
    )
    {
        directionOnlyStepScale = Foam::min
        (
            movlim/mmaCandidateMaxRawDesignChange,
            (1.0 - 1e-6)*directionOnlyBoundScale
        );
    }

    scalar tentativeDesignFluidVolume = 0.0;
    forAll(activeDesignCells, activeI)
    {
        const label celli = activeDesignCells[activeI];
        const scalar candidateChange =
            xmma[activeI] + volumeRestorationShift - oldXmma[activeI];
        const scalar tentativeDesign =
            oldXmma[activeI] + directionOnlyStepScale*candidateChange;
        xmma[activeI] = tentativeDesign;
        tentativeDesignFluidVolume += tentativeDesign*mesh.V()[celli];
    }
    reduce(tentativeDesignFluidVolume, sumOp<scalar>());

    // Remove only the accumulated roundoff from scaling the zero-volume
    // direction. Predictions below use this final accepted design.
    const scalar postScaleVolumeCorrection =
        targetDesignFluidFraction
      - tentativeDesignFluidVolume/Foam::max(designVolume, SMALL);

    scalar restoredDesignFluidVolume = 0.0;
    scalar predictedObjectiveChange = 0.0;
    scalar predictedVolumeChange = 0.0;
    scalar predictedPressureChange = 0.0;
    scalar maxRawDesignChange = 0.0;
    forAll(activeDesignCells, activeI)
    {
        const label celli = activeDesignCells[activeI];
        const scalar restoredDesign =
            xmma[activeI] + postScaleVolumeCorrection;
        if (restoredDesign < 0.0 || restoredDesign > 1.0)
        {
            FatalErrorInFunction
                << "Direction-only step would move active design index "
                << activeI << ", cell " << celli << " outside [0,1]: "
                << restoredDesign << abort(FatalError);
        }

        xmma[activeI] = restoredDesign;
        x[celli] = restoredDesign;
        restoredDesignFluidVolume += restoredDesign*mesh.V()[celli];

        const scalar designChange = restoredDesign - oldXmma[activeI];
        predictedObjectiveChange += dfdx[activeI]*designChange;
        predictedVolumeChange += dgdx[0][activeI]*designChange;
        predictedPressureChange += dgdx[1][activeI]*designChange;
        maxRawDesignChange = Foam::max
        (
            maxRawDesignChange,
            mag(designChange)
        );
    }
    reduce(restoredDesignFluidVolume, sumOp<scalar>());
    reduce(predictedObjectiveChange, sumOp<scalar>());
    reduce(predictedVolumeChange, sumOp<scalar>());
    reduce(predictedPressureChange, sumOp<scalar>());
    reduce(maxRawDesignChange, maxOp<scalar>());

    const scalar restoredFluidFraction =
        restoredDesignFluidVolume/Foam::max(designVolume, SMALL);
    if (mag(restoredFluidFraction - targetDesignFluidFraction) > 1e-12)
    {
        FatalErrorInFunction
            << "Volume restoration failed after roundoff correction: target="
            << targetDesignFluidFraction << ", restored="
            << restoredFluidFraction << abort(FatalError);
    }

    if (Pstream::master())
    {
        static Foam::OFstream mmaStepLog
        (
            runTime.path()/"mma_step_log.dat"
        );
        static bool mmaStepLogHeaderWritten = false;
        mmaStepLog.precision(16);
        if (!mmaStepLogHeaderWritten)
        {
            mmaStepLog
                << "# Iteration\tUpdateMode"
                << "\tMaxRawDesignChange"
                << "\tMmaCandidateMaxRawDesignChange"
                << "\tDirectionOnlyStepScale"
                << "\tPredictedObjectiveChange"
                << "\tPredictedVolumeChange"
                << "\tPredictedPressureChange"
                << "\tVolumeRestorationShift"
                << "\tPostScaleVolumeCorrection"
                << "\tObjectiveGradientScale"
                << "\tPressureGradientScale" << endl << flush;
            mmaStepLogHeaderWritten = true;
        }
        mmaStepLog
            << opt << "\t"
            << directionOnlyUpdateMode << "\t"
            << maxRawDesignChange << "\t"
            << mmaCandidateMaxRawDesignChange << "\t"
            << directionOnlyStepScale << "\t"
            << predictedObjectiveChange << "\t"
            << predictedVolumeChange << "\t"
            << predictedPressureChange << "\t"
            << volumeRestorationShift << "\t"
            << postScaleVolumeCorrection << "\t"
            << objectiveGradientScale << "\t"
            << pressureGradientScale << endl << flush;
    }

    #include "enforceFixedDesignRegions.H"
    #include "filter_x.H"

    // --- Phase I: Candidate evaluation and acceptance ---
    if (enableSSTAcceptanceCheck)
    {
        #include "evaluateCandidate.H"

        if (candidateAccepted)
        {
            // --- Transaction test: verify acceptance on retry
            if (transactionTestMode && candidateAttempt > 1)
            {
                Info<< "  TX_TEST: Candidate ACCEPTED on attempt "
                    << candidateAttempt
                    << " (smaller step after move limit reduction)"
                    << endl << flush;
                Info<< "  TX_TEST: Transaction mechanism test PASSED" << endl
                    << flush;
            }

            mma.commit(xmma);
            consecutiveAccepts++;
            if (consecutiveAccepts >= 3)
            {
                currentMoveLimit = Foam::min
                (
                    currentMoveLimit*moveLimitGrowth,
                    maximumMoveLimit
                );
                consecutiveAccepts = 0;
            }
            // Update accepted-state snapshot to the new design
            #include "saveAcceptedState.H"
        }
        else
        {
            // Rollback to accepted state
            #include "restoreAcceptedState.H"
            mma.restoreState(mmaTxSnapshot);

            // --- Transaction test: verify full rollback after restoration
            if (transactionTestMode && candidateAttempt == 1)
            {
                scalar maxXErr = 0.0, maxUErr = 0.0, maxPErr = 0.0;
                scalar maxTErr = 0.0, maxKErr = 0.0, maxOmegaErr = 0.0;

                // Compare restored fields against the pre-proposal backup
                forAll(xTxBackup, celli)
                {
                    maxXErr = Foam::max(maxXErr,
                        mag(x[celli] - xTxBackup[celli]));
                    maxUErr = Foam::max(maxUErr,
                        mag(U[celli] - UTxBackup[celli]));
                    maxPErr = Foam::max(maxPErr,
                        mag(p[celli] - pTxBackup[celli]));
                    maxTErr = Foam::max(maxTErr,
                        mag(T[celli] - TTxBackup[celli]));
                }
                reduce(maxXErr, maxOp<scalar>());
                reduce(maxUErr, maxOp<scalar>());
                reduce(maxPErr, maxOp<scalar>());
                reduce(maxTErr, maxOp<scalar>());

                // Tensor fields: Frobenius norm of the symmTensor difference
                // (6 components), consistent with the dDTDxh norm used at
                // validateFrozenGradient.H:382.  [v2 R5]
                scalar maxDTErr = 0.0, maxDTEffErr = 0.0;
                forAll(DTMolecular, celli)
                {
                    maxDTErr = Foam::max(maxDTErr,
                        mag(DTMolecular[celli] - DTTxBackup[celli]));
                    maxDTEffErr = Foam::max(maxDTEffErr,
                        mag(DTEffective[celli] - DTEffectiveTxBackup[celli]));
                }
                // [v2 R1] reduce() must stay inside the transactionTestMode
                // block at the same scope level as reduce(maxXErr) above --
                // NEVER inside Pstream::master() (MPI deadlock constraint).
                reduce(maxDTErr, maxOp<scalar>());
                reduce(maxDTEffErr, maxOp<scalar>());

                if (ransFlowModel)
                {
                    const volScalarField& kField =
                        mesh.lookupObject<volScalarField>("k");
                    const volScalarField& omegaField =
                        mesh.lookupObject<volScalarField>("omega");
                    forAll(kTxBackup, celli)
                    {
                        maxKErr = Foam::max(maxKErr,
                            mag(kField[celli] - kTxBackup[celli]));
                        maxOmegaErr = Foam::max(maxOmegaErr,
                            mag(omegaField[celli] - omegaTxBackup[celli]));
                    }
                    reduce(maxKErr, maxOp<scalar>());
                    reduce(maxOmegaErr, maxOp<scalar>());
                }

                // Verify MMA state matches the pre-proposal snapshot
                MMA::StateSnapshot mmaAfterRollback = mma.saveState();
                scalar maxMmaXold1Err = 0.0, maxMmaXold2Err = 0.0;
                scalar maxMmaLowErr = 0.0, maxMmaUppErr = 0.0;
                bool mmaIterMatch =
                    (mmaAfterRollback.iter == mmaTxSnapshot.iter);
                for
                (
                    label i = 0;
                    i < label(mmaTxSnapshot.xold1.size());
                    ++i
                )
                {
                    maxMmaXold1Err = Foam::max(maxMmaXold1Err,
                        mag(mmaAfterRollback.xold1[i]
                           - mmaTxSnapshot.xold1[i]));
                    maxMmaXold2Err = Foam::max(maxMmaXold2Err,
                        mag(mmaAfterRollback.xold2[i]
                           - mmaTxSnapshot.xold2[i]));
                    maxMmaLowErr = Foam::max(maxMmaLowErr,
                        mag(mmaAfterRollback.low[i]
                           - mmaTxSnapshot.low[i]));
                    maxMmaUppErr = Foam::max(maxMmaUppErr,
                        mag(mmaAfterRollback.upp[i]
                           - mmaTxSnapshot.upp[i]));
                }
                reduce(maxMmaXold1Err, maxOp<scalar>());
                reduce(maxMmaXold2Err, maxOp<scalar>());
                reduce(maxMmaLowErr, maxOp<scalar>());
                reduce(maxMmaUppErr, maxOp<scalar>());

                Info<< "  TX_TEST: Rollback verification (vs accepted state):"
                    << endl
                    << "    x:       maxErr=" << maxXErr << endl
                    << "    U:       maxErr=" << maxUErr << endl
                    << "    p:       maxErr=" << maxPErr << endl
                    << "    T:       maxErr=" << maxTErr << endl
                    << "    DT:      maxErr=" << maxDTErr << endl
                    << "    DTEff:   maxErr=" << maxDTEffErr << endl;
                if (ransFlowModel)
                {
                    Info<< "    k:       maxErr=" << maxKErr << endl
                        << "    omega:   maxErr=" << maxOmegaErr << endl;
                }
                Info<< "    MMA:     iterMatch=" << mmaIterMatch
                    << " xold1Err=" << maxMmaXold1Err
                    << " xold2Err=" << maxMmaXold2Err
                    << " lowErr=" << maxMmaLowErr
                    << " uppErr=" << maxMmaUppErr << endl;

                const scalar rollbackTol = 1e-12;
                bool rollbackOK =
                    maxXErr < rollbackTol
                 && maxUErr < rollbackTol
                 && maxPErr < rollbackTol
                 && maxTErr < rollbackTol
                 && maxDTErr < rollbackTol
                 && maxDTEffErr < rollbackTol
                 && mmaIterMatch
                 && maxMmaXold1Err < rollbackTol
                 && maxMmaXold2Err < rollbackTol
                 && maxMmaLowErr < rollbackTol
                 && maxMmaUppErr < rollbackTol;
                if (ransFlowModel)
                {
                    rollbackOK = rollbackOK
                     && maxKErr < rollbackTol
                     && maxOmegaErr < rollbackTol;
                }

                if (rollbackOK)
                {
                    Info<< "  TX_TEST: Rollback verification PASSED"
                        << " (all fields restored to machine precision)"
                        << endl << flush;
                }
                else
                {
                    Info<< "  TX_TEST: Rollback verification FAILED!"
                        << endl << flush;
                }
            }

            currentMoveLimit = Foam::max
            (
                currentMoveLimit*moveLimitReduction,
                minimumMoveLimit
            );
            consecutiveAccepts = 0;

            if (transactionTestMode && candidateAttempt == 1)
            {
                Info<< "  TX_TEST: Move limit reduced to "
                    << currentMoveLimit << endl << flush;
            }

            if (candidateAttempt > maxCandidateRetries)
            {
                Info<< "Transaction: all " << maxCandidateRetries
                    << " retries exhausted. Keeping accepted design."
                    << endl << flush;
                candidateAccepted = true;  // force exit from retry loop
                // Stop the optimization loop
                runTime.setTime
                (
                    runTime.endTime(),
                    runTime.endTime().value()
                );
            }
        }
    }
    } while (enableSSTAcceptanceCheck && !candidateAccepted);
}
else if (mmaUpdateEnabled)
{
    Info<< "Final solved design written without generating an unverified "
        << "MMA candidate." << endl << flush;
}
else
{
    Info<< "mmaUpdateEnabled=false: primal/adjoint fields evaluated without "
        << "advancing the design." << endl << flush;
}
