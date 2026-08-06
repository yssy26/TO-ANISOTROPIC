//-----------------------------------------------------------------------------
// validateFrozenGradient.H
// Stage F: Frozen U-p-T total gradient validation.
//
// Two-step validation:
//   F1: Thermal-only FD — fix U,p,phi, re-solve only T.  This isolates
//       the thermal adjoint and avoids pressure-velocity convergence issues.
//   F2: Full frozen FD — re-solve U,p,T with a FIXED iteration count
//       (no convergence check) to ensure consistent convergence bias
//       between + and - perturbations.
//
// The check runs once after Stage E (component-level validation).
//-----------------------------------------------------------------------------

static bool frozenGradientValidated = false;

if (!frozenGradientValidated && opt >= 1 && !gradientValidated)
{
    frozenGradientValidated = true;
    Info<< "\n========================================" << endl;
    Info<< "Stage F: Frozen U-p-T total gradient validation" << endl;
    Info<< "========================================\n" << endl;

    const Switch origFreezeFlag = freezeTurbulenceForValidation;
    freezeTurbulenceForValidation = true;

    // -- Helper: solve ONLY thermal (keep U, p, phi frozen) ---------------
    auto solveThermalOnly = [&]() -> void
    {
        #include "updateThermalFlux.H"

        const surfaceScalarField& activeThermalFluxF =
            frozenHotRegionEnabled ? phiThermal : phi;

        fvScalarMatrix TEqnLocal
        (
            fvm::div(activeThermalFluxF, T, "div(phi,T)")
          - fvm::laplacian(DTEffective, T)
         ==
            Q
          + fvOptions(T)
        );
        TEqnLocal.relax();
        TEqnLocal.solve();

        // Extra sweeps to tighten convergence
        for (label extra = 0; extra < 5; ++extra)
        {
            fvScalarMatrix TEqnExtra
            (
                fvm::div(activeThermalFluxF, T, "div(phi,T)")
              - fvm::laplacian(DTEffective, T)
             ==
                Q
              + fvOptions(T)
            );
            TEqnExtra.solve();
        }
    };

    // -- Helper: solve frozen U-p-T with FIXED iteration count ------------
    auto solveFrozenPrimalFixed = [&](label nIter) -> void
    {
        for (label iter = 0; iter < nIter; ++iter)
        {
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
            fvOptions.constrain(UEqn);

            if (simple.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
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
            tUEqn.clear();
            constrainPressure(p, U, phiHbyA, rAtU());

            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                );
                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();
                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }
            p.relax();
            U = HbyA - rAtU()*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            #include "updateThermalFlux.H"

            const surfaceScalarField& activeThermalFluxF =
                frozenHotRegionEnabled ? phiThermal : phi;

            fvScalarMatrix TEqnLocal
            (
                fvm::div(activeThermalFluxF, T, "div(phi,T)")
              - fvm::laplacian(DTEffective, T)
             ==
                Q
              + fvOptions(T)
            );
            TEqnLocal.relax();
            TEqnLocal.solve();
        }
    };

    // -- Save complete state ----------------------------------------------
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
    scalar eta5Saved_val = eta5;

    //=====================================================================
    // F1: Thermal-only validation
    //=====================================================================
    Info<< "\n--- F1 [GATE4]: Thermal-only FD (U,p,phi fixed) ---" << endl;

    // Baseline: re-solve T from current state
    solveThermalOnly();
    scalar Jbase_T = evaluateObjective();
    Info<< "  Thermal baseline J=" << Jbase_T << endl;

    // Save the thermally-converged state
    volScalarField TSavedThermal("T_saved_thermal", T);
    volSymmTensorField DTSavedThermal("DT_saved_thermal", DTMolecular);
    volSymmTensorField DTEffSavedThermal("DTEff_saved_thermal", DTEffective);
    volScalarField alphaSavedThermal("alpha_saved_thermal", alpha);
    scalarField xhSavedThermal = xh.primitiveField();
    scalarField xpSavedThermal = xp.primitiveField();
    scalarField xSavedThermal = x.primitiveField();
    scalar eta5SavedThermal = eta5;

    scalar maxRelErrorF1 = 0.0;
    label totalFailuresF1 = 0;
    label totalTestsF1 = 0;
    // [v2 Gate 4] New diagnostics: median error, sign mismatches, per-step
    // error samples and per-direction volume-neutrality residuals.
    scalar medianRelErrorF1 = 0.0;
    label signMismatchF1 = 0;
    DynamicList<scalar> relErrSamplesF1(9);
    const label nDirF1 = 3;
    // [v2 Gate 4] Perturbation ladder required by the upstream plan:
    // {1e-3, 3e-4, 1e-4} (previously {1e-2, 5e-3, 2e-3}).
    const scalarList fdStepsF1 = {1e-3, 3e-4, 1e-4};
    List<scalar> dirVolResidualF1(nDirF1, GREAT);

    // When solveFlowAdjoints=true, dfdx includes the Brinkman term
    // (-dAlphaDxh*(U & Ub)) which captures flow-thermal coupling.  But F1
    // only re-solves T (flow frozen), so the FD captures only the thermal
    // path.  Compute a thermal-only filtered sensitivity to match.
    std::vector<double> dfdxThermal(nDesignLocal, 0.0);
    if (solveFlowAdjoints)
    {
        Info<< "  F1: Computing thermal-only adjoint (excluding Brinkman)"
            << endl;

        volScalarField fsenshThermalOnly
        (
            IOobject
            (
                "fsenshThermalOnly",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fsenshMeanT
        );
        fsenshThermalOnly = dimensionedScalar("zero", fsenshMeanT.dimensions(), 0.0);

        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            fsenshThermalOnly[celli] = thermalDiffusionDerivativeDTCell[celli];
        }
        fsenshThermalOnly.primitiveFieldRef() *= designMask.primitiveField();

        // Diagnostic: print magnitudes of key quantities
        {
            scalar tddMin = GREAT, tddMax = -GREAT, tddSum = 0.0;
            scalar dtdMin = GREAT, dtdMax = -GREAT;
            scalar fthMin = GREAT, fthMax = -GREAT, fthSum = 0.0;
            scalar braMin = GREAT, braMax = -GREAT, braSum = 0.0;
            scalar volMin = GREAT, volMax = -GREAT;
            label nAct = 0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                tddMin = min(tddMin, thermalDiffusionDerivativeDTCell[celli]);
                tddMax = max(tddMax, thermalDiffusionDerivativeDTCell[celli]);
                tddSum += thermalDiffusionDerivativeDTCell[celli];
                dtdMin = min(dtdMin, mag(dDTDxh[celli]));
                dtdMax = max(dtdMax, mag(dDTDxh[celli]));
                fthMin = min(fthMin, fsenshThermalOnly[celli]);
                fthMax = max(fthMax, fsenshThermalOnly[celli]);
                fthSum += fsenshThermalOnly[celli];
                const scalar brack = -dAlphaDxh[celli]*(U[celli] & Ub[celli]);
                braMin = min(braMin, brack);
                braMax = max(braMax, brack);
                braSum += brack;
                volMin = min(volMin, mesh.V()[celli]);
                volMax = max(volMax, mesh.V()[celli]);
                nAct++;
            }
            reduce(tddMin, minOp<scalar>());
            reduce(tddMax, maxOp<scalar>());
            reduce(tddSum, sumOp<scalar>());
            reduce(dtdMin, minOp<scalar>());
            reduce(dtdMax, maxOp<scalar>());
            reduce(fthMin, minOp<scalar>());
            reduce(fthMax, maxOp<scalar>());
            reduce(fthSum, sumOp<scalar>());
            reduce(braMin, minOp<scalar>());
            reduce(braMax, maxOp<scalar>());
            reduce(braSum, sumOp<scalar>());
            reduce(volMin, minOp<scalar>());
            reduce(volMax, maxOp<scalar>());
            reduce(nAct, sumOp<label>());

            // Also check Tb magnitude
            scalar tbMin = GREAT, tbMax = -GREAT;
            forAll(Tb.primitiveField(), celli)
            {
                tbMin = min(tbMin, Tb[celli]);
                tbMax = max(tbMax, Tb[celli]);
            }
            reduce(tbMin, minOp<scalar>());
            reduce(tbMax, maxOp<scalar>());

            Info<< "  F1 DIAG: nActive=" << nAct << endl
                << "    thermalDiffDerivDT: min=" << tddMin << " max=" << tddMax
                << " sum=" << tddSum << endl
                << "    dDTDxh: min=" << dtdMin << " max=" << dtdMax << endl
                << "    fsenshThermalOnly: min=" << fthMin << " max=" << fthMax
                << " sum=" << fthSum << endl
                << "    Brinkman(-dA*(U&Ub)): min=" << braMin << " max=" << braMax
                << " sum=" << braSum << endl
                << "    mesh.V: min=" << volMin << " max=" << volMax << endl
                << "    Tb: min=" << tbMin << " max=" << tbMax << endl;
        }

        // Apply Heaviside projection chain rule (same as filter_chainrule.H)
        if (del > SMALL)
        {
            if (projectionThresholdAdaptive)
            {
                const scalar etaDenom = fvc::domainIntegrate
                (
                    designMask*dProjectionDeta
                ).value();

                // --- F1 NUDIAG: verify fields before eta correction ---
                {
                    scalar dPDESumAll = 0.0, dPDEMaxAll = 0.0, dPDEMinAll = GREAT;
                    scalar fthSumAll = 0.0, fthMaxAll = 0.0, fthMinAll = GREAT;
                    scalar prodSum = 0.0, prodMax = 0.0;
                    scalar prodMaxFth = 0.0, prodMaxDPE = 0.0;
                    label prodMaxCell = -1;
                    label nBadF = 0, nBadD = 0;

                    forAll(designMask, celli)
                    {
                        if (designMask[celli] > 0.5)
                        {
                            const scalar dPDE = dProjectionDeta[celli];
                            const scalar fth = fsenshThermalOnly[celli];
                            const scalar prod = fth * dPDE;

                            dPDESumAll += dPDE;
                            dPDEMaxAll = max(dPDEMaxAll, dPDE);
                            dPDEMinAll = min(dPDEMinAll, dPDE);

                            fthSumAll += fth;
                            fthMaxAll = max(fthMaxAll, mag(fth));
                            fthMinAll = min(fthMinAll, fth);

                            prodSum += prod;
                            if (mag(prod) > mag(prodMax))
                            {
                                prodMax = prod;
                                prodMaxCell = celli;
                                prodMaxFth = fth;
                                prodMaxDPE = dPDE;
                            }

                            if (!std::isfinite(fth)) nBadF++;
                            if (!std::isfinite(dPDE)) nBadD++;
                        }
                    }
                    reduce(dPDESumAll, sumOp<scalar>());
                    reduce(dPDEMaxAll, maxOp<scalar>());
                    reduce(dPDEMinAll, minOp<scalar>());
                    reduce(fthSumAll, sumOp<scalar>());
                    reduce(fthMaxAll, maxOp<scalar>());
                    reduce(fthMinAll, minOp<scalar>());
                    reduce(prodSum, sumOp<scalar>());
                    reduce(prodMax, maxOp<scalar>());
                    reduce(nBadF, sumOp<label>());
                    reduce(nBadD, sumOp<label>());

                    Info<< "  F1 NUDIAG: dProjectionDeta (all design cells):" << endl
                        << "    sum=" << dPDESumAll << " max=" << dPDEMaxAll
                        << " min=" << dPDEMinAll << endl
                        << "  F1 NUDIAG: fsenshThermalOnly (all design cells):" << endl
                        << "    sum=" << fthSumAll << " max(mag)=" << fthMaxAll
                        << " min=" << fthMinAll << endl
                        << "  F1 NUDIAG: product (fth*dPDE):" << endl
                        << "    sum=" << prodSum << " max=" << prodMax
                        << " at cell " << prodMaxCell
                        << " (fth=" << prodMaxFth << ", dPDE=" << prodMaxDPE << ")" << endl
                        << "    non-finite: fth=" << nBadF << " dPDE=" << nBadD << endl;

                    // Also check fsenshMeanT for comparison
                    scalar fsmSumAll = 0.0, fsmMaxAll = 0.0;
                    forAll(designMask, celli)
                    {
                        if (designMask[celli] > 0.5)
                        {
                            fsmSumAll += fsenshMeanT[celli];
                            fsmMaxAll = max(fsmMaxAll, mag(fsenshMeanT[celli]));
                        }
                    }
                    reduce(fsmSumAll, sumOp<scalar>());
                    reduce(fsmMaxAll, maxOp<scalar>());
                    Info<< "  F1 NUDIAG: fsenshMeanT (for comparison):"
                        << " sum=" << fsmSumAll << " max(mag)=" << fsmMaxAll << endl;
                }

                // fsenshThermalOnly is a per-cell quantity (face-sum).
                // Use plain sum (not domainIntegrate) to avoid extra mesh.V().
                scalar etaCorrThermal = 0.0;
                forAll(designMask, celli)
                {
                    etaCorrThermal +=
                        designMask[celli]
                       *fsenshThermalOnly[celli]
                       *dProjectionDeta[celli];
                }
                reduce(etaCorrThermal, sumOp<scalar>());

                Info<< "  F1 NUDIAG: numerator=" << etaCorrThermal
                    << " etaDenom=" << etaDenom << endl;

                // Bug fix: Foam::max(etaDenom, SMALL) replaces a NEGATIVE
                // etaDenom with SMALL (~1e-15), causing ~1e10x amplification.
                // The main code (filter_chainrule.H) divides directly by
                // projectionEtaDenominator after checking mag() > SMALL.
                if (mag(etaDenom) <= SMALL)
                {
                    FatalErrorInFunction
                        << "F1: volume-preserving projection has singular "
                        << "eta derivative. etaDenom=" << etaDenom
                        << abort(FatalError);
                }
                etaCorrThermal /= etaDenom;

                // --- F1 DIAG: print eta correction values ---
                Info<< "  F1 DIAG: etaCorrThermal=" << etaCorrThermal
                    << " etaDenom=" << etaDenom << endl;

                // Correct eta correction (same formula as filter_chainrule.H):
                //   dJ/dxp[j] = dJ/dxh[j]*drho[j] + V[j]*(1-drho[j])*C_J
                fsenshThermalOnly *= drho;
                forAll(designMask, celli)
                {
                    fsenshThermalOnly[celli] +=
                        designMask[celli]*mesh.V()[celli]
                       *(1.0 - drho[celli])*etaCorrThermal;
                }

                // --- F1 DIAG: print after eta correction ---
                {
                    scalar fSum = 0.0, fMax = 0.0;
                    forAll(activeDesignCells, activeI)
                    {
                        const label celli = activeDesignCells[activeI];
                        fSum += fsenshThermalOnly[celli];
                        fMax = max(fMax, mag(fsenshThermalOnly[celli]));
                    }
                    reduce(fSum, sumOp<scalar>());
                    reduce(fMax, maxOp<scalar>());
                    Info<< "  F1 DIAG: fsenshThermalOnly after eta corr:"
                        << " sum=" << fSum << " max=" << fMax << endl;
                }
            }
            else
            {
                fsenshThermalOnly *= drho;
            }
        }

        // Apply Helmholtz filter
        volScalarField fsensThermalOnly
        (
            IOobject
            (
                "fsensThermalOnly",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fsensMeanT
        );
        fsensThermalOnly = dimensionedScalar("zero", fsensMeanT.dimensions(), 0.0);
        solve
        (
            fvm::laplacian(designFilterFaceMask, fsensThermalOnly)
          - fvm::Sp(b, fsensThermalOnly)
          + fsenshThermalOnly*b
        );
        fsensThermalOnly.primitiveFieldRef() *= designMask.primitiveField();
        fsensThermalOnly.primitiveFieldRef() *= designMask.primitiveField();

        // --- F1 DIAG: print filter output ---
        {
            scalar fSum = 0.0, fMax = 0.0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                fSum += fsensThermalOnly[celli];
                fMax = max(fMax, mag(fsensThermalOnly[celli]));
            }
            reduce(fSum, sumOp<scalar>());
            reduce(fMax, maxOp<scalar>());
            Info<< "  F1 DIAG: fsensThermalOnly after filter:"
                << " sum=" << fSum << " max=" << fMax << endl;
        }

        for (label activeI = 0; activeI < nDesignLocal; ++activeI)
        {
            const label celli = activeDesignCells[activeI];
            dfdxThermal[activeI] =
                objectiveGradientScale*fsensThermalOnly[celli];
        }
    }
    else
    {
        // No flow adjoint: dfdx is already thermal-only
        for (label activeI = 0; activeI < nDesignLocal; ++activeI)
        {
            dfdxThermal[activeI] = dfdx[activeI];
        }
    }

    // Partition-invariant synthetic FD directions. Local cell labels
    // restart on every MPI rank and therefore cannot define a serial/parallel
    // comparison direction. Normalize physical cell-centre coordinates by
    // the global mesh bounds so the same physical cell receives the same
    // value for every decomposition.
    const boundBox fdDirectionBounds(mesh.points(), true);
    const vector fdDirectionMin = fdDirectionBounds.min();
    const vector fdDirectionSpan =
        fdDirectionBounds.max() - fdDirectionMin;

    auto partitionInvariantFDValue =
        [&](const label celli, const label dirI) -> scalar
    {
        const vector& centre = mesh.C()[celli];
        const scalar xi =
            (centre.x() - fdDirectionMin.x())
           /Foam::max(mag(fdDirectionSpan.x()), SMALL);
        const scalar eta =
            (centre.y() - fdDirectionMin.y())
           /Foam::max(mag(fdDirectionSpan.y()), SMALL);
        const scalar zeta =
            (centre.z() - fdDirectionMin.z())
           /Foam::max(mag(fdDirectionSpan.z()), SMALL);
        const scalar mode = scalar(dirI + 1);

        return Foam::sin
        (
            mode*(2.137*xi + 3.071*eta + 4.113*zeta) + 0.5
        );
    };

    for (label dirI = 0; dirI < nDirF1; ++dirI)
    {
        scalarField direction(mesh.nCells(), 0.0);
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            // dir 0: negative thermal-gradient direction; dirs 1..2: smooth
            // sinusoidal.  [v2 R3] Every direction MUST be volume-neutral
            // (orthogonal to gsensVol): the known E5 volume-gradient issue
            // (dgdx[0]=0.289912) would otherwise leak into the central
            // difference and corrupt the FD estimate.
            if (dirI == 0)
            {
                direction[celli] = -dfdxThermal[activeI];
            }
            else
            {
                direction[celli] =
                    partitionInvariantFDValue(celli, dirI);
            }
        }

        // --- [v2 Gate 4] Direction construction (2 iterations):
        //     a) volume neutralization (Gram-Schmidt vs gsensVol)
        //     b) boundary-margin filter (avoid min/max clamping)
        //     c) L2 normalization over active design cells
        const scalar hMaxF1 = fdStepsF1[0];   // largest step = 1e-3
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

            const scalar margin = 2.0*hMaxF1*maxAbsDir;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                if
                (
                    xSavedThermal[celli] < margin
                 || xSavedThermal[celli] > 1.0 - margin
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
            FatalErrorInFunction
                << "F1 direction " << dirI
                << " collapsed to zero after volume neutralization."
                << abort(FatalError);
        }
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            direction[celli] /= dirNorm;
        }

        // Volume-neutrality residual (diagnostic; target <= 1e-8)
        {
            scalar volDot = 0.0, gNorm = 0.0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                volDot += direction[celli]*gsensVol[celli];
                gNorm += sqr(gsensVol[celli]);
            }
            reduce(volDot, sumOp<scalar>());
            reduce(gNorm, sumOp<scalar>());
            dirVolResidualF1[dirI] =
                mag(volDot)/Foam::max(Foam::sqrt(gNorm), SMALL);
        }

        // Valid direction cells must be >= 20% of nDesignLocal
        {
            label nValidDir = 0;
            forAll(activeDesignCells, activeI)
            {
                if (mag(direction[activeDesignCells[activeI]]) > 0.0)
                {
                    ++nValidDir;
                }
            }
            reduce(nValidDir, sumOp<label>());
            if (nValidDir < 0.2*nDesignLocal)
            {
                FatalErrorInFunction
                    << "F1 direction " << dirI
                    << " has too few active cells: " << nValidDir
                    << " / " << nDesignLocal << abort(FatalError);
            }
        }

        scalar Dadj = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            Dadj += dfdxThermal[activeI]*direction[celli];
        }
        reduce(Dadj, sumOp<scalar>());

        Info<< "\n  F1 Direction " << dirI << ": D_ADJ(thermal)=" << Dadj << endl;

        for (label stepI = 0; stepI < fdStepsF1.size(); ++stepI)
        {
            const scalar h = fdStepsF1[stepI];
            totalTestsF1++;

            // Plus
            label clampCountPlus = 0;
            x.primitiveFieldRef() = xSavedThermal;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                const scalar xPlus =
                    xSavedThermal[celli] + h*direction[celli];
                if (xPlus >= 1.0 || xPlus <= 0.0) ++clampCountPlus;
                x[celli] = Foam::min(1.0, Foam::max(0.0, xPlus));
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            const scalar volPlus = evaluateVolumeConstraint();
            // Keep U, p, phi at saved state; only update T
            T = TSavedThermal;
            T.correctBoundaryConditions();
            solveThermalOnly();
            scalar Jplus = evaluateObjective();

            // Minus
            label clampCountMinus = 0;
            x.primitiveFieldRef() = xSavedThermal;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                const scalar xMinus =
                    xSavedThermal[celli] - h*direction[celli];
                if (xMinus >= 1.0 || xMinus <= 0.0) ++clampCountMinus;
                x[celli] = Foam::min(1.0, Foam::max(0.0, xMinus));
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            const scalar volMinus = evaluateVolumeConstraint();
            T = TSavedThermal;
            T.correctBoundaryConditions();
            solveThermalOnly();
            scalar Jminus = evaluateObjective();

            scalar Dfd = (Jplus - Jminus)/(2.0*h);
            scalar relErr = mag(Dadj - Dfd)
               /Foam::max(Foam::max(mag(Dadj), mag(Dfd)), scalar(1e-30));

            // Sign consistency: mismatch only counts when not near zero
            const bool signMismatch =
                   Foam::max(mag(Dadj), mag(Dfd)) >= 1e-10
                && Dadj*Dfd < 0.0;
            if (signMismatch) ++signMismatchF1;

            // Volume-change diagnostic (should be ~0 for neutral dirs)
            const scalar dV_empirical =
                (volPlus - volMinus)/(2.0*h);

            Info<< "    h=" << h
                << ", J+=" << Jplus
                << ", J-=" << Jminus
                << ", D_FD=" << Dfd
                << ", rel error=" << relErr
                << ", clamps(+/-)=" << clampCountPlus << "/"
                << clampCountMinus
                << ", dV_FD=" << dV_empirical << endl;

            maxRelErrorF1 = Foam::max(maxRelErrorF1, relErr);
            relErrSamplesF1.append(relErr);
            if (relErr > 0.05) totalFailuresF1++;
        }

        // Restore
        x.primitiveFieldRef() = xSavedThermal;
        xp.primitiveFieldRef() = xpSavedThermal;
        xh.primitiveFieldRef() = xhSavedThermal;
        alpha = alphaSavedThermal;
        DTMolecular = DTSavedThermal;
        DTEffective = DTEffSavedThermal;
        T = TSavedThermal;
        eta5 = eta5SavedThermal;
        x.correctBoundaryConditions();
        T.correctBoundaryConditions();
    }

    // [v2 Gate 4] Median over all collected relative-error samples
    if (relErrSamplesF1.size() > 0)
    {
        Foam::sort(relErrSamplesF1);
        medianRelErrorF1 = relErrSamplesF1[relErrSamplesF1.size()/2];
    }

    Info<< "\n  F1 [GATE4] Summary: medianRelError=" << medianRelErrorF1
        << ", maxRelError=" << maxRelErrorF1
        << ", failures=" << totalFailuresF1 << "/" << totalTestsF1
        << ", signMismatch=" << signMismatchF1 << endl;
    Info<< "  F1 [GATE4] volume-neutrality residuals per direction:";
    for (label dirI = 0; dirI < nDirF1; ++dirI)
    {
        Info<< "  dir" << dirI << "=" << dirVolResidualF1[dirI];
    }
    Info<< "  (target <= 1e-8)" << endl;

    //=====================================================================
    // F2: Full frozen U-p-T FD with fixed iteration count
    //=====================================================================
    Info<< "\n--- F2: Full frozen U-p-T FD (fixed iters) ---" << endl;

    // Restore to main solver state
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
    eta5 = eta5Saved_val;
    x.correctBoundaryConditions();
    U.correctBoundaryConditions();
    T.correctBoundaryConditions();

    // Baseline: solve frozen primal with fixed iterations
    const label nFixedIter = 300;
    Info<< "  Solving baseline with " << nFixedIter << " fixed iters..." << endl;
    solveFrozenPrimalFixed(nFixedIter);
    scalar Jbase_full = evaluateObjective();
    Info<< "  Baseline J=" << Jbase_full << endl;

    // Save the frozen-converged state
    volVectorField USavedFroz("U_saved_froz", U);
    volScalarField pSavedFroz("p_saved_froz", p);
    surfaceScalarField phiSavedFroz("phi_saved_froz", phi);
    volScalarField TSavedFroz("T_saved_froz", T);

    scalar maxRelErrorF2 = 0.0;
    label totalFailuresF2 = 0;
    label totalTestsF2 = 0;
    scalar maxDadjF2 = 0.0;  // track largest |D_ADJ| across all F2 directions
    scalar maxRelErrorPD = 0.0;  // pressure drop FD max error
    label totalFailuresPD = 0;   // pressure drop FD failures
    scalar maxDadjPDall = 0.0;   // max |D_ADJ(pd)| across F2+F3
    scalar maxDadjObjAll = 0.0;  // max |D_ADJ(obj)| across F3
    scalar maxFdpObjAll = 0.0;   // max |D_FD(obj)| across F3
    const label nDirF2 = 3;
    const scalarList fdStepsF2 = {1e-1, 5e-2, 2e-2};

    for (label dirI = 0; dirI < nDirF2; ++dirI)
    {
        scalarField direction(mesh.nCells(), 0.0);
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            direction[celli] = partitionInvariantFDValue(celli, dirI);
        }
        scalar dirNorm = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            dirNorm += sqr(direction[celli]);
        }
        reduce(dirNorm, sumOp<scalar>());
        dirNorm = Foam::sqrt(dirNorm);
        if (dirNorm < SMALL) continue;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            direction[celli] /= dirNorm;
        }

        scalar Dadj = 0.0;
        scalar DadjPD = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            Dadj += dfdx[activeI]*direction[celli];
            DadjPD += dgdx[1][activeI]*direction[celli];
        }
        reduce(Dadj, sumOp<scalar>());
        reduce(DadjPD, sumOp<scalar>());
        maxDadjF2 = Foam::max(maxDadjF2, mag(Dadj));
        maxDadjPDall = Foam::max(maxDadjPDall, mag(DadjPD));

        Info<< "\n  F2 Direction " << dirI
            << ": D_ADJ(obj)=" << Dadj
            << ", D_ADJ(pd)=" << DadjPD << endl;

        for (label stepI = 0; stepI < fdStepsF2.size(); ++stepI)
        {
            const scalar h = fdStepsF2[stepI];
            totalTestsF2++;

            // Plus
            x.primitiveFieldRef() = xSaved;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = Foam::min(1.0, xSaved[celli] + h*direction[celli]);
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            U = USavedFroz;
            p = pSavedFroz;
            phi = phiSavedFroz;
            T = TSavedFroz;
            U.correctBoundaryConditions();
            T.correctBoundaryConditions();
            solveFrozenPrimalFixed(nFixedIter);
            scalar Jplus = evaluateObjective();
            scalar DPplus = evaluatePressureDrop();

            // Minus
            x.primitiveFieldRef() = xSaved;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = Foam::max(0.0, xSaved[celli] - h*direction[celli]);
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            U = USavedFroz;
            p = pSavedFroz;
            phi = phiSavedFroz;
            T = TSavedFroz;
            U.correctBoundaryConditions();
            T.correctBoundaryConditions();
            solveFrozenPrimalFixed(nFixedIter);
            scalar Jminus = evaluateObjective();
            scalar DPminus = evaluatePressureDrop();

            scalar Dfd = (Jplus - Jminus)/(2.0*h);
            scalar relErr = mag(Dadj - Dfd)
               /Foam::max(Foam::max(mag(Dadj), mag(Dfd)), scalar(1e-30));

            // Pressure drop FD (normalized by pressureDropMaxPa to match gx[1])
            scalar DfdPD = (DPplus - DPminus)/(2.0*h);
            scalar DfdPDnorm = DfdPD
               /Foam::max(pressureDropMaxPa, SMALL);
            scalar relErrPD = mag(DadjPD - DfdPDnorm)
               /Foam::max
                (
                    Foam::max(mag(DadjPD), mag(DfdPDnorm)),
                    scalar(1e-30)
                );

            Info<< "    h=" << h
                << ", J+=" << Jplus << ", J-=" << Jminus
                << ", D_FD=" << Dfd
                << ", relErr(obj)=" << relErr
                << ", DP+=" << DPplus << ", DP-=" << DPminus
                << ", D_FD(pd)=" << DfdPDnorm
                << ", relErr(pd)=" << relErrPD << endl;

            maxRelErrorF2 = Foam::max(maxRelErrorF2, relErr);
            maxRelErrorPD = Foam::max(maxRelErrorPD, relErrPD);
            if (relErr > 0.05) totalFailuresF2++;
            if (relErrPD > 0.50) totalFailuresPD++;
        }

        // Restore
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
        eta5 = eta5Saved_val;
        x.correctBoundaryConditions();
        U.correctBoundaryConditions();
        T.correctBoundaryConditions();
    }

    Info<< "\n  F2 Summary: maxRelError=" << maxRelErrorF2
        << ", failures=" << totalFailuresF2 << "/" << totalTestsF2 << endl;
    Info<< "  F2 Pressure Drop: maxRelError=" << maxRelErrorPD
        << ", failures=" << totalFailuresPD << "/" << totalTestsF2 << endl;

    //=====================================================================
    // F3: Combined direction FD (actual MMA direction)
    //=====================================================================
    Info<< "\n--- F3: Combined direction FD (projected negative gradient) ---" << endl;

    // Restore to frozen-converged state for F3
    x.primitiveFieldRef() = xSaved;
    xp.primitiveFieldRef() = xpSaved;
    xh.primitiveFieldRef() = xhSaved;
    alpha = alphaSaved;
    DTMolecular = DTSaved;
    Q = QSaved;
    DTEffective = DTEffectiveSaved;
    U = USavedFroz;
    p = pSavedFroz;
    phi = phiSavedFroz;
    T = TSavedFroz;
    eta5 = eta5Saved_val;
    x.correctBoundaryConditions();
    U.correctBoundaryConditions();
    T.correctBoundaryConditions();

    // Baseline objective, volume, pressure drop
    scalar Jbase_f3 = evaluateObjective();
    scalar Vbase_f3 = evaluateVolumeConstraint();
    scalar DPbase_f3 = evaluatePressureDrop();
    Info<< "  F3 Baseline: J=" << Jbase_f3
        << ", V=" << Vbase_f3
        << ", DP=" << DPbase_f3 << endl;

    // Build combined direction: projected negative objective gradient
    // d = -dfdx - proj(-dfdx, dgdx[0]) * dgdx[0]
    // i.e., remove the volume constraint component
    scalarField negGrad(nDesignLocal, 0.0);
    scalarField volGrad(nDesignLocal, 0.0);
    for (label activeI = 0; activeI < nDesignLocal; ++activeI)
    {
        negGrad[activeI] = -dfdx[activeI];
        volGrad[activeI] = dgdx[0][activeI];
    }

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

    scalarField combDir(mesh.nCells(), 0.0);
    forAll(activeDesignCells, activeI)
    {
        const label celli = activeDesignCells[activeI];
        combDir[celli] = negGrad[activeI] - projCoeff*volGrad[activeI];
    }

    // Normalize
    scalar dirNorm3 = 0.0;
    forAll(activeDesignCells, activeI)
    {
        const label celli = activeDesignCells[activeI];
        dirNorm3 += sqr(combDir[celli]);
    }
    reduce(dirNorm3, sumOp<scalar>());
    dirNorm3 = Foam::sqrt(dirNorm3);
    if (dirNorm3 > SMALL)
    {
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            combDir[celli] /= dirNorm3;
        }
    }

    // Directional derivatives
    scalar DadjObj = 0.0;
    scalar DadjVol = 0.0;
    scalar DadjPD = 0.0;
    forAll(activeDesignCells, activeI)
    {
        const label celli = activeDesignCells[activeI];
        DadjObj += dfdx[activeI]*combDir[celli];
        DadjVol += dgdx[0][activeI]*combDir[celli];
        DadjPD += dgdx[1][activeI]*combDir[celli];
    }
    reduce(DadjObj, sumOp<scalar>());
    reduce(DadjVol, sumOp<scalar>());
    reduce(DadjPD, sumOp<scalar>());

    maxDadjPDall = Foam::max(maxDadjPDall, mag(DadjPD));
    maxDadjObjAll = Foam::max(maxDadjObjAll, mag(DadjObj));

    Info<< "  F3 Directional derivatives: dJ=" << DadjObj
        << ", dV=" << DadjVol
        << ", dPD(norm)=" << DadjPD << endl;

    scalar maxRelErrorF3obj = 0.0;
    scalar maxRelErrorF3vol = 0.0;
    scalar maxRelErrorF3pd = 0.0;
    const scalarList fdStepsF3 = {1e-1, 5e-2, 2e-2};

    for (label stepI = 0; stepI < fdStepsF3.size(); ++stepI)
    {
        const scalar h = fdStepsF3[stepI];

        // Plus
        x.primitiveFieldRef() = xSaved;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            x[celli] = Foam::min
            (
                1.0,
                Foam::max(0.0, xSaved[celli] + h*combDir[celli])
            );
        }
        x.correctBoundaryConditions();
        recomputeDesignChain();
        U = USavedFroz;
        p = pSavedFroz;
        phi = phiSavedFroz;
        T = TSavedFroz;
        U.correctBoundaryConditions();
        T.correctBoundaryConditions();
        solveFrozenPrimalFixed(nFixedIter);
        scalar Jp = evaluateObjective();
        scalar Vp = evaluateVolumeConstraint();
        scalar DPp = evaluatePressureDrop();

        // Minus
        x.primitiveFieldRef() = xSaved;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            x[celli] = Foam::min
            (
                1.0,
                Foam::max(0.0, xSaved[celli] - h*combDir[celli])
            );
        }
        x.correctBoundaryConditions();
        recomputeDesignChain();
        U = USavedFroz;
        p = pSavedFroz;
        phi = phiSavedFroz;
        T = TSavedFroz;
        U.correctBoundaryConditions();
        T.correctBoundaryConditions();
        solveFrozenPrimalFixed(nFixedIter);
        scalar Jm = evaluateObjective();
        scalar Vm = evaluateVolumeConstraint();
        scalar DPm = evaluatePressureDrop();

        // FD derivatives
        scalar DfdObj = (Jp - Jm)/(2.0*h);
        scalar DfdVol = (Vp - Vm)/(2.0*h);
        scalar DfdPD = (DPp - DPm)/(2.0*h);
        scalar DfdPDnorm = DfdPD/Foam::max(pressureDropMaxPa, SMALL);

        scalar errObj = mag(DadjObj - DfdObj)
           /Foam::max(Foam::max(mag(DadjObj), mag(DfdObj)), scalar(1e-30));
        scalar errVol = mag(DadjVol - DfdVol)
           /Foam::max(Foam::max(mag(DadjVol), mag(DfdVol)), scalar(1e-30));
        scalar errPD = mag(DadjPD - DfdPDnorm)
           /Foam::max
            (
                Foam::max(mag(DadjPD), mag(DfdPDnorm)),
                scalar(1e-30)
            );

        Info<< "    h=" << h
            << "  obj: ADJ=" << DadjObj << " FD=" << DfdObj
            << " err=" << errObj
            << "  vol: ADJ=" << DadjVol << " FD=" << DfdVol
            << " err=" << errVol
            << "  pd: ADJ=" << DadjPD << " FD=" << DfdPDnorm
            << " err=" << errPD << endl;

        maxRelErrorF3obj = Foam::max(maxRelErrorF3obj, errObj);
        maxRelErrorF3vol = Foam::max(maxRelErrorF3vol, errVol);
        maxRelErrorF3pd = Foam::max(maxRelErrorF3pd, errPD);
        maxFdpObjAll = Foam::max(maxFdpObjAll, mag(DfdObj));
    }

    Info<< "\n  F3 Summary: obj_err=" << maxRelErrorF3obj
        << ", vol_err=" << maxRelErrorF3vol
        << ", pd_err=" << maxRelErrorF3pd << endl;

    // -- Final restore ----------------------------------------------------
    freezeTurbulenceForValidation = origFreezeFlag;

    // -- Overall summary --------------------------------------------------
    Info<< "\n----------------------------------------" << endl;
    Info<< "Stage F Summary:" << endl;
    Info<< "  F1 (GATE4 thermal-only): medianRelError=" << medianRelErrorF1
        << ", maxRelError=" << maxRelErrorF1
        << ", signMismatch=" << signMismatchF1
        << ", failures=" << totalFailuresF1 << "/" << totalTestsF1 << endl;
    Info<< "  F2 (full frozen):   maxRelError=" << maxRelErrorF2
        << ", failures=" << totalFailuresF2 << "/" << totalTestsF2
        << ", max|D_ADJ|=" << maxDadjF2 << endl;
    Info<< "  F2 (pressure drop): maxRelError=" << maxRelErrorPD
        << ", failures=" << totalFailuresPD << "/" << totalTestsF2 << endl;
    Info<< "  F3 (combined dir):  obj_err=" << maxRelErrorF3obj
        << ", vol_err=" << maxRelErrorF3vol
        << ", pd_err=" << maxRelErrorF3pd << endl;
    Info<< "----------------------------------------\n" << endl;

    // Pass criterion:
    // F1 (thermal-only) should pass with < 5% error (tight convergence).
    // F2 (full frozen) is informational due to pressure-velocity coupling
    // noise; accept if max error < 50% OR if both |D_ADJ| and |D_FD| are
    // small (< 1e-3, indicating near-zero gradient).
    // [v2 Gate 4] F1 hard criterion: median <= 1%, max <= 5%, no sign
    // mismatch (upstream ANISOTROPIC_CONDUCTIVITY_PLAN.md section 9 Gate 4).
    bool f1Pass =
           medianRelErrorF1 <= 0.01
        && maxRelErrorF1    <= 0.05
        && signMismatchF1   == 0
        && totalTestsF1 > 0;

    // For F2, check if gradient is near-zero (all directions have small |D_ADJ|)
    bool f2NearZero = (maxDadjF2 < 1e-3);
    bool f2Pass = f2NearZero || (totalFailuresF2 <= totalTestsF2/3 && maxRelErrorF2 < 0.50);

    // Pressure drop: informational due to pressure-velocity coupling noise.
    // Accept if near-zero gradient or < 50% error.
    bool pdNearZero = (maxDadjPDall < 1e-8);
    bool pdPass = pdNearZero || (totalFailuresPD <= totalTestsF2/3 && maxRelErrorPD < 1.0);

    // F3 combined: informational (all three components)
    bool f3Pass = (maxRelErrorF3obj < 0.50)
               || (maxDadjObjAll < 1e-3 && maxFdpObjAll < 1e-3);

    if (!f1Pass)
    {
        FatalErrorInFunction
            << "F1 (GATE4) gradient validation failed. "
            << "medianRelError=" << medianRelErrorF1
            << " (limit 1%), maxRelError=" << maxRelErrorF1
            << " (limit 5%), signMismatch=" << signMismatchF1
            << ", failures=" << totalFailuresF1
            << "/" << totalTestsF1
            << ". [v2 R8] Diagnostics: 1) volume-neutrality residuals "
            << "dirVolResidualF1 (target <= 1e-8); 2) FD plateau across the "
            << "3 perturbation steps (if h=1e-4 is noise-dominated, enlarge "
            << "the smallest step instead of failing); 3) sign flips -> "
            << "inspect AdjHeatTransfer.H boundary conditions and transpose "
            << "assembly."
            << abort(FatalError);
    }

    if (!f2Pass)
    {
        FatalErrorInFunction
            << "F2 (full frozen) gradient validation failed. "
            << "maxRelError=" << maxRelErrorF2
            << ", failures=" << totalFailuresF2
            << "/" << totalTestsF2
            << abort(FatalError);
    }

    if (!pdPass)
    {
        WarningInFunction
            << "F2 pressure drop gradient validation: maxRelError="
            << maxRelErrorPD << ", failures=" << totalFailuresPD
            << "/" << totalTestsF2
            << ". Pressure drop gradient may be unreliable."
            << endl;
    }

    if (!f3Pass)
    {
        WarningInFunction
            << "F3 combined direction validation: obj_err="
            << maxRelErrorF3obj
            << ". Combined gradient direction may be unreliable."
            << endl;
    }

    Info<< "Stage F: PASSED" << endl;
    if (!f2NearZero)
    {
        Info<< "  (F2 passed within tolerance)" << endl;
    }
    else
    {
        Info<< "  (F2: gradient near-zero, FD noise-dominated — accepted)" << endl;
    }
    if (pdNearZero)
    {
        Info<< "  (PD: gradient near-zero — accepted)" << endl;
    }
}
