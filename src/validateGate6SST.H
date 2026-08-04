//-----------------------------------------------------------------------------
// validateGate6SST.H
// Gate 6: Full-SST central-difference hard gate.
//
// Replaces the single-sided descent test (validateSSTDirection.H, which stays
// as a diagnostic) with a central-difference hard gate:
//   - >= 4 volume-neutral design directions
//   - +h / -h perturbations, each re-converged to FULL SST from the same
//     baseline (solveFullSST, 100 iters / 1e-3)
//   - FD plateau confirmation across the {1e-3, 3e-4, 1e-4} ladder
//   - sign consistency for the temperature objective and the pressure-drop
//     constraint, plus mass and thermal balance checks per perturbed state
//
// Pass  -> ransDirectionValidated=true, gradientValidated=true, GATE6: PASSED
// Fail  -> all flags stay false, GATE6: FAILED with details (Warning only).
//
// Relies on the shared lambdas in validateCommon.H (included before this
// file in MTO_HF.C).  [v2 R7] gate6PlateauTolerance default 0.25 is TEMPORARY
// and must be calibrated with the first run's actual spread data.
//-----------------------------------------------------------------------------

static bool gate6Validated = false;
if (!gate6Validated && opt >= 1 && ransFlowModel && !ransDirectionValidated)
{
    gate6Validated = true;
    Info<< "\n========================================" << endl;
    Info<< "GATE6: Full-SST central-difference hard gate" << endl;
    Info<< "========================================\n" << endl;

    const Switch origFreezeG6 = freezeTurbulenceForValidation;

    // [v2 R7] Temporary default; calibrate with first-run spread data.
    const scalar gate6PlateauTolerance = 0.25;
    const scalar massBalanceTolerance = 1e-3;
    const scalar thermalBalanceTolerance = 1e-2;
    const scalarList fdStepsG6 = {1e-3, 3e-4, 1e-4};
    const label nDirG6 = 4;

    //=====================================================================
    // Save complete pre-gate state (same set as Stage G)
    //=====================================================================
    scalarField xG6Saved  = x.primitiveField();
    scalarField xpG6Saved = xp.primitiveField();
    scalarField xhG6Saved = xh.primitiveField();
    volScalarField alphaG6Saved("alpha_g6s", alpha);
    volSymmTensorField DTG6Saved("DT_g6s", DTMolecular);
    volScalarField QG6Saved("Q_g6s", Q);
    volVectorField UG6Saved("U_g6s", U);
    volScalarField pG6Saved("p_g6s", p);
    surfaceScalarField phiG6Saved("phi_g6s", phi);
    volScalarField TG6Saved("T_g6s", T);
    volSymmTensorField DTEffG6Saved("DTEff_g6s", DTEffective);
    volScalarField nutFG6Saved("nutF_g6s", nutFrozen);
    volScalarField nuEffG6Saved("nuEff_g6s", nuEffFrozen);
    volScalarField alphaTG6Saved("alphaT_g6s", alphaTurbulent);
    volScalarField fluidIG6Saved("fluidI_g6s", fluidIndicator);
    volScalarField chiTG6Saved("chiT_g6s", chiTurbulence);
    volScalarField kG6Saved("k_g6s", mesh.lookupObject<volScalarField>("k"));
    volScalarField omegaG6Saved("omega_g6s", mesh.lookupObject<volScalarField>("omega"));
    volScalarField nutG6Saved("nut_g6s", mesh.lookupObject<volScalarField>("nut"));
    scalar eta5G6Saved = eta5;

    // Baseline: full-SST re-convergence.  200-iteration budget: starting
    // from the frozen state, 100 iterations are insufficient (the nutFrozen
    // relative change needs ~130-150 iterations; Stage G's own baseline
    // shows the same).  Perturbed solves warm-start from the converged
    // baseline and use the plan's 100-iteration budget.
    label baselineIters = solveFullSST(200, 1e-3);
    Info<< "GATE6-DBG: baseline solve done, iters=" << baselineIters
        << endl << flush;
    const scalar Jbase = evaluateObjective();
    Info<< "GATE6-DBG: Jbase=" << Jbase << endl << flush;
    const scalar DPbase = evaluatePressureDrop();
    Info<< "GATE6-DBG: DPbase=" << DPbase << endl << flush;
    const scalar massBase = evaluateMassBalance();
    Info<< "GATE6-DBG: massBase=" << massBase << endl << flush;
    const scalar heatBase = evaluateHeatBalance();
    Info<< "GATE6-DBG: heatBase=" << heatBase << endl << flush;
    Info<< "  GATE6 baseline: iters=" << baselineIters
        << ", J=" << Jbase
        << ", DP=" << DPbase
        << ", massImbalance=" << massBase
        << ", thermalBalanceError=" << heatBase << endl;

    // Warm-start SST-converged baseline state
    volVectorField UG6Sst("U_g6sst", U);
    volScalarField pG6Sst("p_g6sst", p);
    surfaceScalarField phiG6Sst("phi_g6sst", phi);
    volScalarField TG6Sst("T_g6sst", T);
    volScalarField kG6Sst("k_g6sst", mesh.lookupObject<volScalarField>("k"));
    volScalarField omegaG6Sst("omega_g6sst", mesh.lookupObject<volScalarField>("omega"));
    volScalarField nutG6Sst("nut_g6sst", mesh.lookupObject<volScalarField>("nut"));
    volScalarField nutFG6Sst("nutF_g6sst", nutFrozen);
    volScalarField nuEffG6Sst("nuEff_g6sst", nuEffFrozen);
    volScalarField alphaTG6Sst("alphaT_g6sst", alphaTurbulent);
    volSymmTensorField DTEffG6Sst("DTEff_g6sst", DTEffective);
    scalarField xG6Sst = x.primitiveField();
    scalarField xpG6Sst = xp.primitiveField();
    scalarField xhG6Sst = xh.primitiveField();
    scalar eta5G6Sst = eta5;

    auto restoreG6SstBase = [&]()
    {
        x.primitiveFieldRef() = xG6Sst;
        xp.primitiveFieldRef() = xpG6Sst;
        xh.primitiveFieldRef() = xhG6Sst;
        alpha = alphaG6Saved;
        DTMolecular = DTG6Saved;
        Q = QG6Saved;
        DTEffective = DTEffG6Sst;
        U = UG6Sst;
        p = pG6Sst;
        phi = phiG6Sst;
        T = TG6Sst;
        nutFrozen = nutFG6Sst;
        nuEffFrozen = nuEffG6Sst;
        alphaTurbulent = alphaTG6Sst;
        fluidIndicator = fluidIG6Saved;
        chiTurbulence = chiTG6Saved;
        mesh.lookupObjectRef<volScalarField>("k") = kG6Sst;
        mesh.lookupObjectRef<volScalarField>("omega") = omegaG6Sst;
        mesh.lookupObjectRef<volScalarField>("nut") = nutG6Sst;
        eta5 = eta5G6Sst;
        x.correctBoundaryConditions();
        U.correctBoundaryConditions();
        T.correctBoundaryConditions();
    };

    //=====================================================================
    // Gate criteria accumulators
    //=====================================================================
    bool gate6Pass = true;
    bool gate6SSTAllConverged = (baselineIters > 0);
    bool gate6SignAllOK = true;
    bool gate6PlateauAllOK = true;
    bool gate6BalanceAllOK = true;
    label gate6SignJBad = 0, gate6SignDPBad = 0;
    label gate6ConvBad = 0, gate6MassBad = 0, gate6HeatBad = 0;

    //=====================================================================
    // Directions
    //=====================================================================
    for (label dirI = 0; dirI < nDirG6; ++dirI)
    {
        // -- direction construction: volume neutral (same as Gate 4 M2) --
        scalarField direction(mesh.nCells(), 0.0);
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            if (dirI == 0)
            {
                // negative objective-gradient direction
                direction[celli] = -dfdx[activeI];
            }
            else if (dirI == 1)
            {
                direction[celli] =
                    Foam::sin(0.271*(dirI + 1)*(celli + 1) + 0.5);
            }
            else
            {
                direction[celli] =
                    Foam::sin(0.173*(dirI + 1)*(celli + 1) + 1.3);
            }
        }

        const scalar hMaxG6 = fdStepsG6[0];
        scalar maxAbsDir = 0.0;
        forAll(activeDesignCells, activeI)
        {
            maxAbsDir = Foam::max
            (
                maxAbsDir,
                mag(direction[activeDesignCells[activeI]])
            );
        }
        reduce(maxAbsDir, maxOp<scalar>());
        maxAbsDir = Foam::max(maxAbsDir, SMALL);

        for (label gsIter = 0; gsIter < 2; ++gsIter)
        {
            scalar dotDV = 0.0, dotVV = 0.0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                dotDV += direction[celli]*gsensVol[celli];
                dotVV += sqr(gsensVol[celli]);
            }
            reduce(dotDV, sumOp<scalar>());
            reduce(dotVV, sumOp<scalar>());
            if (mag(dotVV) >= SMALL)
            {
                const scalar projCoeff = dotDV/dotVV;
                forAll(activeDesignCells, activeI)
                {
                    const label celli = activeDesignCells[activeI];
                    direction[celli] -= projCoeff*gsensVol[celli];
                }
            }
            const scalar margin = 2.0*hMaxG6*maxAbsDir;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                if
                (
                    xG6Saved[celli] < margin
                 || xG6Saved[celli] > 1.0 - margin
                )
                {
                    direction[celli] = 0.0;
                }
            }
        }

        scalar dirNorm = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            dirNorm += sqr(direction[celli]);
        }
        reduce(dirNorm, sumOp<scalar>());
        dirNorm = Foam::sqrt(dirNorm);
        if (dirNorm < SMALL)
        {
            Info<< "  GATE6 WARNING: direction " << dirI
                << " collapsed after volume neutralization; skipped."
                << endl;
            continue;
        }
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            direction[celli] /= dirNorm;
        }

        // adjoint predictions
        scalar AD_J = 0.0, AD_DP = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            AD_J += dfdx[activeI]*direction[celli];
            AD_DP += dgdx[1][activeI]*direction[celli];
        }
        reduce(AD_J, sumOp<scalar>());
        reduce(AD_DP, sumOp<scalar>());
        AD_DP /= Foam::max(pressureDropMaxPa, SMALL);

        Info<< "\n  GATE6 Direction " << dirI
            << ": AD_J=" << AD_J << ", AD_DP=" << AD_DP << endl;

        DynamicList<scalar> fdJSamples(3), fdDPSamples(3);

        for (label stepI = 0; stepI < fdStepsG6.size(); ++stepI)
        {
            const scalar h = fdStepsG6[stepI];

            // --- +h perturbation (from SST baseline) ---
            restoreG6SstBase();
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = Foam::min
                (
                    1.0,
                    Foam::max(0.0, xG6Saved[celli] + h*direction[celli])
                );
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            label itersP = solveFullSST(100, 1e-3);
            const scalar Jplus = evaluateObjective();
            const scalar DPplus = evaluatePressureDrop();
            const scalar massPlus = evaluateMassBalance();
            const scalar heatPlus = evaluateHeatBalance();

            // --- -h perturbation (from SST baseline) ---
            restoreG6SstBase();
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = Foam::min
                (
                    1.0,
                    Foam::max(0.0, xG6Saved[celli] - h*direction[celli])
                );
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            label itersM = solveFullSST(100, 1e-3);
            const scalar Jminus = evaluateObjective();
            const scalar DPminus = evaluatePressureDrop();
            const scalar massMinus = evaluateMassBalance();
            const scalar heatMinus = evaluateHeatBalance();

            const scalar FD_J = (Jplus - Jminus)/(2.0*h);
            const scalar FD_DP =
                (DPplus - DPminus)/(2.0*h)/Foam::max(pressureDropMaxPa, SMALL);

            // sign consistency (near-zero adjoint exempted and reported)
            const bool signJOK =
                   (AD_J*FD_J > 0.0)
                || (mag(AD_J) < 1e-10);
            const bool signDPOK =
                   (AD_DP*FD_DP > 0.0)
                || (mag(AD_DP) < 1e-10);
            if (!signJOK) ++gate6SignJBad;
            if (!signDPOK) ++gate6SignDPBad;
            if (!signJOK || !signDPOK) gate6SignAllOK = false;

            const bool convOK = (itersP > 0 && itersM > 0);
            if (!convOK) ++gate6ConvBad;
            if (!convOK) gate6SSTAllConverged = false;

            const bool massOK =
                   massPlus <= massBalanceTolerance
                && massMinus <= massBalanceTolerance;
            const bool heatOK =
                   heatPlus <= thermalBalanceTolerance
                && heatMinus <= thermalBalanceTolerance;
            if (!massOK) ++gate6MassBad;
            if (!heatOK) ++gate6HeatBad;
            if (!massOK || !heatOK) gate6BalanceAllOK = false;

            fdJSamples.append(FD_J);
            fdDPSamples.append(FD_DP);

            Info<< "    h=" << h
                << ": iters(+/-)=" << itersP << "/" << itersM
                << ", FD_J=" << FD_J << " (AD " << AD_J << ")"
                << ", FD_DP=" << FD_DP << " (AD " << AD_DP << ")"
                << ", sign(J/DP)=" << (signJOK ? "OK" : "BAD")
                << "/" << (signDPOK ? "OK" : "BAD")
                << ", mass=" << massPlus << "/" << massMinus
                << ", heat=" << heatPlus << "/" << heatMinus << endl;
        }

        // FD plateau: same sign across the ladder and bounded spread
        {
            scalar fdJmin = GREAT, fdJmax = -GREAT, fdJmed = 0.0;
            forAll(fdJSamples, i)
            {
                fdJmin = Foam::min(fdJmin, mag(fdJSamples[i]));
                fdJmax = Foam::max(fdJmax, mag(fdJSamples[i]));
                fdJmed += mag(fdJSamples[i]);
            }
            fdJmed /= Foam::max(label(fdJSamples.size()), 1);
            const bool sameSignJ = (fdJSamples[0] > 0)
                ? (fdJSamples[1] > 0 && fdJSamples[2] > 0)
                : (fdJSamples[0] < 0 && fdJSamples[1] < 0 && fdJSamples[2] < 0);
            const scalar spreadJ = (fdJmax - fdJmin)
                / Foam::max(fdJmed, 1e-30);
            const bool plateauJ = sameSignJ && spreadJ <= gate6PlateauTolerance;

            scalar fdDPmin = GREAT, fdDPmax = -GREAT, fdDPmed = 0.0;
            forAll(fdDPSamples, i)
            {
                fdDPmin = Foam::min(fdDPmin, mag(fdDPSamples[i]));
                fdDPmax = Foam::max(fdDPmax, mag(fdDPSamples[i]));
                fdDPmed += mag(fdDPSamples[i]);
            }
            fdDPmed /= Foam::max(label(fdDPSamples.size()), 1);
            const bool sameSignDP = (fdDPSamples[0] > 0)
                ? (fdDPSamples[1] > 0 && fdDPSamples[2] > 0)
                : (fdDPSamples[0] < 0 && fdDPSamples[1] < 0 && fdDPSamples[2] < 0);
            const scalar spreadDP = (fdDPmax - fdDPmin)
                / Foam::max(fdDPmed, 1e-30);
            const bool plateauDP = sameSignDP && spreadDP <= gate6PlateauTolerance;

            if (!plateauJ || !plateauDP) gate6PlateauAllOK = false;

            Info<< "  GATE6 Direction " << dirI
                << " plateau: J spread=" << spreadJ
                << " (tol " << gate6PlateauTolerance << ")"
                << (plateauJ ? " OK" : " BAD")
                << ", DP spread=" << spreadDP
                << (plateauDP ? " OK" : " BAD") << endl;
        }

        // restore pre-gate state at the end of each direction
        x.primitiveFieldRef() = xG6Saved;
        xp.primitiveFieldRef() = xpG6Saved;
        xh.primitiveFieldRef() = xhG6Saved;
        alpha = alphaG6Saved;
        DTMolecular = DTG6Saved;
        Q = QG6Saved;
        DTEffective = DTEffG6Saved;
        U = UG6Saved;
        p = pG6Saved;
        phi = phiG6Saved;
        T = TG6Saved;
        nutFrozen = nutFG6Saved;
        nuEffFrozen = nuEffG6Saved;
        alphaTurbulent = alphaTG6Saved;
        fluidIndicator = fluidIG6Saved;
        chiTurbulence = chiTG6Saved;
        mesh.lookupObjectRef<volScalarField>("k") = kG6Saved;
        mesh.lookupObjectRef<volScalarField>("omega") = omegaG6Saved;
        mesh.lookupObjectRef<volScalarField>("nut") = nutG6Saved;
        eta5 = eta5G6Saved;
        x.correctBoundaryConditions();
        U.correctBoundaryConditions();
        T.correctBoundaryConditions();
    }

    //=====================================================================
    // Final judgement
    //=====================================================================
    freezeTurbulenceForValidation = origFreezeG6;

    gate6Pass =
           gate6SSTAllConverged
        && gate6SignAllOK
        && gate6PlateauAllOK
        && gate6BalanceAllOK;

    Info<< "\n----------------------------------------" << endl;
    Info<< "GATE6 Summary:" << endl;
    Info<< "  SST convergence:      "
        << (gate6SSTAllConverged ? "PASS" : "FAIL")
        << " (baseline iters=" << baselineIters
        << ", bad=" << gate6ConvBad << ")" << endl;
    Info<< "  Sign consistency:     "
        << (gate6SignAllOK ? "PASS" : "FAIL")
        << " (J bad=" << gate6SignJBad
        << ", DP bad=" << gate6SignDPBad << ")" << endl;
    Info<< "  FD plateau:           "
        << (gate6PlateauAllOK ? "PASS" : "FAIL")
        << " (plateau tol=" << gate6PlateauTolerance << ")" << endl;
    Info<< "  Mass balance (1e-3):  "
        << (gate6MassBad == 0 ? "PASS" : "FAIL")
        << " (bad=" << gate6MassBad << ")" << endl;
    Info<< "  Thermal balance (1e-2):"
        << (gate6HeatBad == 0 ? "PASS" : "FAIL")
        << " (bad=" << gate6HeatBad << ")" << endl;
    Info<< "----------------------------------------" << endl;

    if (gate6Pass)
    {
        ransDirectionValidated = true;
        gradientValidated = true;
        Info<< "GATE6: PASSED" << endl;
        Info<< "  (full-SST central differences consistent with the "
            << "adjoint gradient across " << nDirG6 << " directions)" << endl;
    }
    else
    {
        WarningInFunction
            << "GATE6: FAILED -- the frozen gradient does not reproduce the "
            << "full-SST central-difference response on all checks. "
            << "mmaUpdateEnabled/gradientValidated/ransDirectionValidated "
            << "remain false (known frozen-turbulence gradient defect; "
            << "ANISOTROPIC_CONDUCTIVITY_PLAN.md section 10.3)." << endl;
        Info<< "GATE6: FAILED" << endl;
    }
    Info<< "\n========================================\n" << endl;
}
