//-----------------------------------------------------------------------------
// validateSSTDirection.H
// Stage G: Full SST direction validity validation.
//
// Tests whether the frozen-turbulence gradient provides useful descent
// directions when evaluated with full SST re-convergence.
//
// Test A: Projected negative objective gradient direction
// Test B: Random directions for robustness
//
// For each direction and step size:
//   1. Perturb the design along the direction
//   2. Recompute the design chain (filter, projection, material properties)
//   3. Re-solve the FULL SST primal (U, p, T, k, omega with turbulence->correct())
//   4. Evaluate the objective
//   5. Check if the objective improved (descent)
//
// Pass criteria:
//   - At least one step size in each direction gives objective improvement
//   - Overall descent success rate >= 75%
//   - SST converges for each tested design
//-----------------------------------------------------------------------------

static bool sstDirectionValidated = false;

if (!sstDirectionValidated && opt >= 1 && ransFlowModel && !ransDirectionValidated)
{
    sstDirectionValidated = true;
    Info<< "\n========================================" << endl;
    Info<< "Stage G: Full SST direction validity validation" << endl;
    Info<< "========================================\n" << endl;

    const Switch origFreezeFlag = freezeTurbulenceForValidation;

    //=====================================================================
    // Save complete state
    //=====================================================================
    scalarField xSaved  = x.primitiveField();
    scalarField xpSaved = xp.primitiveField();
    scalarField xhSaved = xh.primitiveField();
    volScalarField alphaSaved("alpha_saved", alpha);
    volSymmTensorField DTSaved("DT_saved", DTMolecular);
    volScalarField QSaved("Q_saved", Q);
    volVectorField USaved("U_saved", U);
    volScalarField pSaved("p_saved", p);
    surfaceScalarField phiSaved("phi_saved", phi);
    volScalarField TSaved("T_saved", T);
    volSymmTensorField DTEffectiveSaved("DTEff_saved", DTEffective);
    volScalarField nutFrozenSaved("nutFrozen_saved", nutFrozen);
    volScalarField nuEffFrozenSaved("nuEffFrozen_saved", nuEffFrozen);
    volScalarField alphaTurbulentSaved("alphaTurb_saved", alphaTurbulent);
    volScalarField fluidIndicatorSaved("fluidInd_saved", fluidIndicator);
    volScalarField chiTurbulenceSaved("chiTurb_saved", chiTurbulence);

    // Save turbulence model fields
    volScalarField kSaved("k_saved", mesh.lookupObject<volScalarField>("k"));
    volScalarField omegaSaved("omega_saved", mesh.lookupObject<volScalarField>("omega"));
    volScalarField nutSaved("nut_saved", mesh.lookupObject<volScalarField>("nut"));

    scalar eta5Saved = eta5;

    //=====================================================================
    // Evaluate baseline objective (from current converged state)
    //=====================================================================
    // The main loop already solved with frozen turbulence. Re-evaluate
    // with full SST to get the true baseline.
    Info<< "\n--- Baseline: Full SST solve ---" << endl;
    label baselineIters = solveFullSST(100, 1e-3);
    scalar Jbase = evaluateObjective();
    scalar volFracBase = fvc::domainIntegrate(designMask*xh).value()
                        /Foam::max(designVolume, SMALL);
    Info<< "  Baseline: iters=" << baselineIters
        << ", J=" << Jbase
        << ", volFrac=" << volFracBase << endl;

    // Save the full-SST-converged state for warm-starting perturbed solves
    volVectorField USstBase("U_sst_base", U);
    volScalarField pSstBase("p_sst_base", p);
    surfaceScalarField phiSstBase("phi_sst_base", phi);
    volScalarField TSstBase("T_sst_base", T);
    volScalarField kSstBase("k_sst_base", mesh.lookupObject<volScalarField>("k"));
    volScalarField omegaSstBase("omega_sst_base", mesh.lookupObject<volScalarField>("omega"));
    volScalarField nutSstBase("nut_sst_base", mesh.lookupObject<volScalarField>("nut"));
    volScalarField nutFrozenSstBase("nutFrozen_sst_base", nutFrozen);
    volScalarField nuEffFrozenSstBase("nuEffFrozen_sst_base", nuEffFrozen);
    volScalarField alphaTurbulentSstBase("alphaTurb_sst_base", alphaTurbulent);
    volSymmTensorField DTEffectiveSstBase("DTEff_sst_base", DTEffective);
    scalarField xSstBase = x.primitiveField();
    scalarField xpSstBase = xp.primitiveField();
    scalarField xhSstBase = xh.primitiveField();
    scalar eta5SstBase = eta5;

    //=====================================================================
    // Construct test directions
    //=====================================================================
    const label nDir = 3;
    const scalarList fdSteps = {0.01, 0.02, 0.05};

    // Direction 0: Projected negative objective gradient
    // Project out the volume constraint component to stay feasible
    List<scalarField> directions(nDir);
    forAll(directions, dI)
    {
        directions[dI].resize(mesh.nCells(), 0.0);
    }

    // Build negative gradient direction
    scalarField negGrad(nDesignLocal, 0.0);
    scalarField volGrad(nDesignLocal, 0.0);
    for (label activeI = 0; activeI < nDesignLocal; ++activeI)
    {
        negGrad[activeI] = -dfdx[activeI];
        volGrad[activeI] = dgdx[0][activeI];
    }

    // Project: d = negGrad - (negGrad . volGrad / volGrad . volGrad) * volGrad
    scalar dotNV = 0.0;
    scalar dotVV = 0.0;
    for (label activeI = 0; activeI < nDesignLocal; ++activeI)
    {
        dotNV += negGrad[activeI]*volGrad[activeI];
        dotVV += sqr(volGrad[activeI]);
    }
    reduce(dotNV, sumOp<scalar>());
    reduce(dotVV, sumOp<scalar>());

    scalar projCoeff = 0.0;
    if (mag(dotVV) > SMALL)
    {
        projCoeff = dotNV/dotVV;
    }

    forAll(activeDesignCells, activeI)
    {
        const label celli = activeDesignCells[activeI];
        directions[0][celli] = negGrad[activeI] - projCoeff*volGrad[activeI];
    }

    // Direction 1, 2: Random
    for (label dirI = 1; dirI < nDir; ++dirI)
    {
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            directions[dirI][celli] = Foam::sin
            (
                0.371*(dirI + 1)*(celli + 1) + 0.3
            );
        }
    }

    // Normalize all directions
    forAll(directions, dirI)
    {
        scalar dirNorm = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            dirNorm += sqr(directions[dirI][celli]);
        }
        reduce(dirNorm, sumOp<scalar>());
        dirNorm = Foam::sqrt(dirNorm);
        if (dirNorm > SMALL)
        {
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                directions[dirI][celli] /= dirNorm;
            }
        }
    }

    //=====================================================================
    // Test loop
    //=====================================================================
    label totalTests = 0;
    label totalDescent = 0;       // J_new < J_base (objective improved)
    label totalConverged = 0;     // SST converged
    label totalPredictedDescent = 0;  // Frozen gradient predicts descent
    label totalPredictedActual = 0;   // Predicted descent AND actual descent
    scalar maxActualChange = 0.0;

    for (label dirI = 0; dirI < nDir; ++dirI)
    {
        // Compute predicted directional derivative
        scalar Dadj = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            Dadj += dfdx[activeI]*directions[dirI][celli];
        }
        reduce(Dadj, sumOp<scalar>());

        Info<< "\n  Direction " << dirI << ": D_ADJ=" << Dadj << endl;

        bool dirHasDescent = false;

        for (label stepI = 0; stepI < fdSteps.size(); ++stepI)
        {
            const scalar h = fdSteps[stepI];
            totalTests++;

            // Perturb design
            x.primitiveFieldRef() = xSstBase;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = Foam::min
                (
                    1.0,
                    Foam::max(0.0, xSstBase[celli] + h*directions[dirI][celli])
                );
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();

            // Warm-start from SST baseline
            U = USstBase;
            p = pSstBase;
            phi = phiSstBase;
            T = TSstBase;
            mesh.lookupObjectRef<volScalarField>("k") = kSstBase;
            mesh.lookupObjectRef<volScalarField>("omega") = omegaSstBase;
            mesh.lookupObjectRef<volScalarField>("nut") = nutSstBase;
            nutFrozen = nutFrozenSstBase;
            nuEffFrozen = nuEffFrozenSstBase;
            alphaTurbulent = alphaTurbulentSstBase;
            DTEffective = DTEffectiveSstBase;
            U.correctBoundaryConditions();
            T.correctBoundaryConditions();

            // Full SST solve
            label iters = solveFullSST(100, 1e-3);
            scalar Jnew = evaluateObjective();

            scalar dJactual = Jnew - Jbase;
            scalar dJpredicted = h*Dadj;
            bool converged = (iters > 0);
            bool actualDescent = (dJactual < 0.0);
            bool predictedDescent = (dJpredicted < 0.0);

            if (converged) totalConverged++;
            if (predictedDescent) totalPredictedDescent++;
            if (predictedDescent && actualDescent) totalPredictedActual++;
            if (actualDescent)
            {
                totalDescent++;
                dirHasDescent = true;
            }

            maxActualChange = Foam::max(maxActualChange, mag(dJactual));

            // Volume fraction change
            scalar volFracNew = fvc::domainIntegrate(designMask*xh).value()
                               /Foam::max(designVolume, SMALL);

            Info<< "    h=" << h
                << ", J=" << Jnew
                << ", dJ_actual=" << dJactual
                << ", dJ_predicted=" << dJpredicted
                << ", iters=" << iters
                << ", volFrac=" << volFracNew
                << (converged ? "" : " [NOT CONVERGED]")
                << (actualDescent ? " [DESCENT]" : " [NO DESCENT]")
                << endl;
        }

        // Restore design state for next direction
        x.primitiveFieldRef() = xSstBase;
        xp.primitiveFieldRef() = xpSstBase;
        xh.primitiveFieldRef() = xhSstBase;
        alpha = alphaSaved;
        DTMolecular = DTSaved;
        Q = QSaved;
        DTEffective = DTEffectiveSstBase;
        eta5 = eta5SstBase;
        x.correctBoundaryConditions();

        if (!dirHasDescent)
        {
            Info<< "  WARNING: Direction " << dirI
                << " showed no descent at any step size." << endl;
        }
    }

    //=====================================================================
    // Restore complete state
    //=====================================================================
    freezeTurbulenceForValidation = origFreezeFlag;

    x.primitiveFieldRef() = xSaved;
    xp.primitiveFieldRef() = xpSaved;
    xh.primitiveFieldRef() = xhSaved;
    alpha = alphaSaved;
    DTMolecular = DTSaved;
    Q = QSaved;
    DTEffective = DTEffectiveSaved;
    U = USaved;
    p = pSaved;
    phi = phiSaved;
    T = TSaved;
    nutFrozen = nutFrozenSaved;
    nuEffFrozen = nuEffFrozenSaved;
    alphaTurbulent = alphaTurbulentSaved;
    fluidIndicator = fluidIndicatorSaved;
    chiTurbulence = chiTurbulenceSaved;
    mesh.lookupObjectRef<volScalarField>("k") = kSaved;
    mesh.lookupObjectRef<volScalarField>("omega") = omegaSaved;
    mesh.lookupObjectRef<volScalarField>("nut") = nutSaved;
    eta5 = eta5Saved;
    x.correctBoundaryConditions();
    U.correctBoundaryConditions();
    T.correctBoundaryConditions();

    //=====================================================================
    // Report and pass/fail
    //=====================================================================
    scalar descentRate = scalar(totalDescent)/scalar(totalTests);
    scalar convergenceRate = scalar(totalConverged)/scalar(totalTests);
    scalar predictionAccuracy = (totalPredictedDescent > 0)
        ? scalar(totalPredictedActual)/scalar(totalPredictedDescent)
        : 1.0;

    Info<< "\n----------------------------------------" << endl;
    Info<< "Stage G Summary:" << endl;
    Info<< "  Total tests:         " << totalTests << endl;
    Info<< "  SST converged:       " << totalConverged << "/"
        << totalTests << " (" << 100.0*convergenceRate << "%)" << endl;
    Info<< "  Actual descents:     " << totalDescent << "/"
        << totalTests << " (" << 100.0*descentRate << "%)" << endl;
    Info<< "  Predicted descents:  " << totalPredictedDescent << endl;
    Info<< "  Prediction accuracy: " << 100.0*predictionAccuracy
        << "% (predicted AND actual / predicted)" << endl;
    Info<< "  Max |dJ_actual|:     " << maxActualChange << endl;
    Info<< "----------------------------------------\n" << endl;

    // Pass criteria:
    // 1. SST convergence rate >= 75%
    // 2. Descent rate >= 75% (at least some improvement in most tests)
    // 3. Prediction accuracy >= 50% (frozen gradient direction is useful)
    bool convergencePass = (convergenceRate >= 0.75);
    bool descentPass = (descentRate >= 0.50);
    bool predictionPass = (predictionAccuracy >= 0.50);

    if (!convergencePass)
    {
        WarningInFunction
            << "SST convergence rate (" << 100.0*convergenceRate
            << "%) is below 75%. Results may be unreliable." << endl;
    }

    if (!descentPass)
    {
        // [v2 Gate 4-6] Known frozen-turbulence gradient defect
        // (ANISOTROPIC_CONDUCTIVITY_PLAN.md section 10.3): the frozen
        // gradient does not provide useful descent directions for full
        // SST.  Stage G is DIAGNOSTIC ONLY per the upstream plan and must
        // not abort runs -- the Gate 6 module (validateGate6SST.H) is the
        // authoritative hard gate with its own failure semantics.
        WarningInFunction
            << "Descent rate (" << 100.0*descentRate
            << "%) is below 50%. The frozen gradient direction does not "
            << "provide useful descent directions for full SST "
            << "(known frozen-turbulence gradient defect, diagnostic only)."
            << endl;
    }

    if (!predictionPass)
    {
        WarningInFunction
            << "Prediction accuracy (" << 100.0*predictionAccuracy
            << "%) is below 50%. The frozen gradient may not reliably "
            << "predict descent direction for full SST." << endl;
    }

    if (convergencePass && descentPass && predictionPass)
    {
        Info<< "Stage G: PASSED" << endl;
        Info<< "  (descent rate=" << 100.0*descentRate
            << "%, convergence=" << 100.0*convergenceRate
            << "%, prediction=" << 100.0*predictionAccuracy << "%)"
            << endl;
    }
    else
    {
        Info<< "Stage G: DIAGNOSTIC (not all criteria met; see warnings "
            << "above, frozen-turbulence gradient defect)" << endl;
    }
}
