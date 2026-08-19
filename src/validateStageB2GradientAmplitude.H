//-----------------------------------------------------------------------------
// validateStageB2GradientAmplitude.H
// Stage B1.5/B2/B3 frozen-turbulence gradient-amplitude validation.
//
// Formal acceptance rules:
//   * D1/D2/D3 are judged independently (no global best-step loophole).
//   * Every +/- finite-difference primal must converge.
//   * J:   same sign, per-direction best relative error <= 5%.
//   * gDP: same sign, per-direction best relative error <= 10%.
//   * gV:  same sign, per-direction best relative error <= 1e-6.
//   * Each metric requires a stable adjacent-step FD plateau (<=30% spread).
//
// This module never advances MMA and never edits optProperties.  A PASS writes
// a small unlock snippet that may be copied into the production case only
// after the report has been reviewed.
//-----------------------------------------------------------------------------

if (stageB2Enabled)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Stage B2/B3 requires serial execution (reduced discrete flow "
            << "adjoint is serial-only)." << abort(FatalError);
    }

    if (mmaUpdateEnabled)
    {
        FatalErrorInFunction
            << "Stage B2/B3 must run with mmaUpdateEnabled=false."
            << abort(FatalError);
    }

    Info<< "\n========================================\n"
        << " Stage B1.5/B2/B3: gradient amplitude validation\n"
        << "========================================\n" << endl;

    const fileName b2OutDir = runTime.path()/"stageB2";
    mkDir(b2OutDir);

    struct B2State
    {
        scalarField x, xp, xh;
        volVectorField U;
        volScalarField p, T, alpha, Q;
        volScalarField nutFrozen, nuEffFrozen, alphaTurbulent;
        volSymmTensorField DTMolecular, DTEffective;
        surfaceScalarField phi;
        volScalarField k, omega, nut;
        scalar eta5;
    };

    const B2State entryState =
    {
        x.primitiveField(), xp.primitiveField(), xh.primitiveField(),
        U, p, T, alpha, Q,
        nutFrozen, nuEffFrozen, alphaTurbulent,
        DTMolecular, DTEffective,
        phi,
        mesh.lookupObject<volScalarField>("k"),
        mesh.lookupObject<volScalarField>("omega"),
        mesh.lookupObject<volScalarField>("nut"),
        eta5
    };

    auto restoreStateInto = [&](const B2State& s)
    {
        x.primitiveFieldRef() = s.x;
        xp.primitiveFieldRef() = s.xp;
        xh.primitiveFieldRef() = s.xh;
        U = s.U;
        p = s.p;
        T = s.T;
        alpha = s.alpha;
        Q = s.Q;
        nutFrozen = s.nutFrozen;
        nuEffFrozen = s.nuEffFrozen;
        alphaTurbulent = s.alphaTurbulent;
        DTMolecular = s.DTMolecular;
        DTEffective = s.DTEffective;
        phi = s.phi;
        mesh.lookupObjectRef<volScalarField>("k") = s.k;
        mesh.lookupObjectRef<volScalarField>("omega") = s.omega;
        mesh.lookupObjectRef<volScalarField>("nut") = s.nut;
        eta5 = s.eta5;
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();
        T.correctBoundaryConditions();
        alpha.correctBoundaryConditions();
        Q.correctBoundaryConditions();
        nutFrozen.correctBoundaryConditions();
        nuEffFrozen.correctBoundaryConditions();
        alphaTurbulent.correctBoundaryConditions();
        DTMolecular.correctBoundaryConditions();
        DTEffective.correctBoundaryConditions();
    };

    // ---------------------------------------------------------------------
    // Phase 1: full-SST baseline, then freeze turbulence.
    // ---------------------------------------------------------------------
    Info<< "Stage B2 phase 1: full-SST baseline solve" << endl;
    freezeTurbulenceForValidation = false;
    const label nFullSST = solveFullSST(1000, 1e-6);
    if (nFullSST < 0)
    {
        FatalErrorInFunction
            << "Stage B2 phase 1: full-SST baseline did not converge."
            << abort(FatalError);
    }
    Info<< "Stage B2 phase 1: full-SST converged in " << nFullSST
        << " iterations." << endl;

    freezeTurbulenceForValidation = true;
    #include "updateFrozenTurbulenceFields.H"

    const B2State baselineState =
    {
        x.primitiveField(), xp.primitiveField(), xh.primitiveField(),
        U, p, T, alpha, Q,
        nutFrozen, nuEffFrozen, alphaTurbulent,
        DTMolecular, DTEffective,
        phi,
        mesh.lookupObject<volScalarField>("k"),
        mesh.lookupObject<volScalarField>("omega"),
        mesh.lookupObject<volScalarField>("nut"),
        eta5
    };

    const scalar Jbase = evaluateObjective();
    const scalar PDbase = evaluatePressureDrop();
    const scalar gDPbase = PDbase/Foam::max(pressureDropMaxPa, SMALL) - 1.0;
    const scalar gVbase = evaluateVolumeConstraint();
    Info<< "Stage B2 baseline: J=" << Jbase << " PD=" << PDbase
        << " gDP=" << gDPbase << " gV=" << gVbase << endl;

    // ---------------------------------------------------------------------
    // Converged frozen U-p-T solver.  A finite-difference point is invalid
    // if this routine returns -1; invalid points never enter the formal gate.
    // ---------------------------------------------------------------------
    auto solveFrozenPrimalConverged = [&](label maxIter) -> label
    {
        label consecutive = 0;
        const scalar J0 = evaluateObjective();
        const scalar gDP0 = evaluatePressureDrop()
                          /Foam::max(pressureDropMaxPa, SMALL) - 1.0;
        // Stage B1 noise-floor based convergence.  The frozen primal is
        // bit-reproducible but runs a deterministic limit cycle: J oscillates
        // at ~1e-7 (relative) and never satisfies a 1e-10 platform window.
        // The reproducible converged value is the sliding-window mean (the
        // Stage B1 "platform value").  We accept convergence when the
        // non-overlapping window mean of both J and gDP stops drifting by
        // more than a small multiple of the measured noise floor.
        const label platformWindow = 20;
        DynamicList<scalar> Jhist;
        DynamicList<scalar> gDPhist;
        Jhist.append(J0);
        gDPhist.append(gDP0);
        const scalar JtolRel = 1e-6;
        const scalar gDPtolRel = 1e-6;

        for (label iter = 0; iter < maxIter; ++iter)
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

            const scalar Jnow = evaluateObjective();
            const scalar gDPnow = evaluatePressureDrop()
                                  /Foam::max(pressureDropMaxPa, SMALL) - 1.0;
            Jhist.append(Jnow);
            gDPhist.append(gDPnow);

            if (Jhist.size() >= 2*platformWindow)
            {
                scalar JprevW = 0.0, JcurW = 0.0;
                scalar gDPprevW = 0.0, gDPcurW = 0.0;
                for (label i = 0; i < platformWindow; ++i)
                {
                    JprevW += Jhist[Jhist.size() - 2*platformWindow + i];
                    JcurW  += Jhist[Jhist.size() - platformWindow + i];
                    gDPprevW += gDPhist[Jhist.size() - 2*platformWindow + i];
                    gDPcurW  += gDPhist[Jhist.size() - platformWindow + i];
                }
                JprevW /= platformWindow;
                JcurW  /= platformWindow;
                gDPprevW /= platformWindow;
                gDPcurW  /= platformWindow;

                const scalar dJrel =
                    mag(JcurW - JprevW)/Foam::max(mag(JprevW), SMALL);
                const scalar dgDPrel =
                    mag(gDPcurW - gDPprevW)/Foam::max(mag(gDPprevW), SMALL);

                if (dJrel <= JtolRel && dgDPrel <= gDPtolRel)
                {
                    ++consecutive;
                }
                else
                {
                    consecutive = 0;
                }
                if (consecutive >= 2)
                {
                    return iter + 1;
                }
            }
        }
        return -1;
    };

    auto runFullB2Chain = [&]()
    {
        {
            #include "updateThermalFlux.H"
        }
        {
            #include "HeatTransfer.H"
        }
        {
            #include "computeObjective.H"
        }
        {
            #include "AdjHeatTransfer.H"
        }
        if (!freezeColdFlowForValidation && solveFlowAdjoints)
        {
            #include "AdjNS_HT.H"
            #include "AdjNS_PD.H"
        }
        {
            #include "sensitivity.H"
        }
    };

    auto gradientProjection = [&](const std::vector<double>& g,
                                  const scalarField& Dk) -> scalar
    {
        scalar projection = 0.0;
        forAll(activeDesignCells, activeI)
        {
            projection += g[activeI]*Dk[activeDesignCells[activeI]];
        }
        reduce(projection, sumOp<scalar>());
        return projection;
    };

    // ---------------------------------------------------------------------
    // Deterministic physical-coordinate directions.
    // ---------------------------------------------------------------------
    scalarField D1(mesh.nCells(), 0.0);
    scalarField D2(mesh.nCells(), 0.0);
    scalarField D3(mesh.nCells(), 0.0);
    {
        const vectorField& C = mesh.C().primitiveField();
        const scalar Lx = gMax(C.component(0)) - gMin(C.component(0));
        const scalar Ly = gMax(C.component(1)) - gMin(C.component(1));
        const scalar Lz = gMax(C.component(2)) - gMin(C.component(2));

        forAll(C, celli)
        {
            if (designMask[celli] <= 0.5)
            {
                continue;
            }
            const scalar xx = C[celli].x()/Foam::max(Lx, SMALL);
            const scalar yy = C[celli].y()/Foam::max(Ly, SMALL);
            const scalar zz = C[celli].z()/Foam::max(Lz, SMALL);
            D1[celli] =
                Foam::sin(2.0*Foam::constant::mathematical::pi*xx)
               *Foam::cos(Foam::constant::mathematical::pi*yy)
               +0.5*Foam::sin(3.0*Foam::constant::mathematical::pi*zz);
            D2[celli] =
                Foam::sin(4.0*Foam::constant::mathematical::pi*xx)
               +0.3*Foam::cos(2.0*Foam::constant::mathematical::pi*xx);
            D3[celli] =
                Foam::cos(2.0*Foam::constant::mathematical::pi*yy)
               *Foam::sin(2.0*Foam::constant::mathematical::pi*zz);

            if (x[celli] < 0.02 || x[celli] > 0.98)
            {
                D1[celli] = 0.0;
                D2[celli] = 0.0;
                D3[celli] = 0.0;
            }
        }

        auto normalizeDirection = [&](scalarField& d)
        {
            scalar mx = 0.0;
            forAll(d, celli)
            {
                mx = Foam::max(mx, mag(d[celli]));
            }
            reduce(mx, maxOp<scalar>());
            if (mx > 1e-30)
            {
                d /= mx;
            }
        };
        normalizeDirection(D1);
        normalizeDirection(D2);
        normalizeDirection(D3);
    }

    // ---------------------------------------------------------------------
    // Baseline adjoint projections. The production discrete solver is used
    // here because validateMmaUnlockGate forbids B2+B4 in the same run.
    // ---------------------------------------------------------------------
    runFullB2Chain();
    const scalar projJ_D1 = gradientProjection(dfdx, D1);
    const scalar projJ_D2 = gradientProjection(dfdx, D2);
    const scalar projJ_D3 = gradientProjection(dfdx, D3);
    const scalar projDP_D1 = gradientProjection(dgdx[1], D1);
    const scalar projDP_D2 = gradientProjection(dgdx[1], D2);
    const scalar projDP_D3 = gradientProjection(dgdx[1], D3);
    const scalar projV_D1 = gradientProjection(dgdx[0], D1);
    const scalar projV_D2 = gradientProjection(dgdx[0], D2);
    const scalar projV_D3 = gradientProjection(dgdx[0], D3);

    Info<< "Stage B2 adjoint projections (J): D1=" << projJ_D1
        << " D2=" << projJ_D2 << " D3=" << projJ_D3 << endl;
    Info<< "Stage B2 adjoint projections (gDP): D1=" << projDP_D1
        << " D2=" << projDP_D2 << " D3=" << projDP_D3 << endl;
    Info<< "Stage B2 adjoint projections (gV): D1=" << projV_D1
        << " D2=" << projV_D2 << " D3=" << projV_D3 << endl;

    Info<< "Stage B1.5 direction self-check: ||D1||inf=" << gMax(mag(D1))
        << " ||D2||inf=" << gMax(mag(D2))
        << " ||D3||inf=" << gMax(mag(D3))
        << " (same fields used for FD and projections)" << endl;

    // BFINAL-013 P2: per-term contraction decomposition at exactly this
    // baseline state (switch-gated, read-only).
    if
    (
        optProperties.lookupOrDefault<Switch>
        (
            "stageB13GradientProbe",
            false
        )
    )
    {
        #include "stageB13ContractionProbe.H"
    }

    struct FDProbe
    {
        scalar J, PD, gDP, gV;
        label nIter;
        bool valid;
    };

    auto probeFD = [&](const scalarField& Dk, const scalar h,
                       const scalar sign) -> FDProbe
    {
        restoreStateInto(baselineState);
        forAll(Dk, celli)
        {
            x.primitiveFieldRef()[celli] += sign*h*Dk[celli];
        }
        recomputeDesignChain();

        FDProbe result;
        result.J = GREAT;
        result.PD = GREAT;
        result.gDP = GREAT;
        result.gV = GREAT;
        result.nIter = solveFrozenPrimalConverged(1000);
        result.valid = result.nIter > 0;

        if (result.valid)
        {
            result.J = evaluateObjective();
            result.PD = evaluatePressureDrop();
            result.gDP =
                result.PD/Foam::max(pressureDropMaxPa, SMALL) - 1.0;
            result.gV = evaluateVolumeConstraint();
        }
        return result;
    };

    auto probeVolumeOnly = [&](const scalarField& Dk, const scalar h,
                               const scalar sign) -> scalar
    {
        restoreStateInto(baselineState);
        forAll(Dk, celli)
        {
            x.primitiveFieldRef()[celli] += sign*h*Dk[celli];
        }
        recomputeDesignChain();
        return evaluateVolumeConstraint();
    };

    struct FDRow
    {
        word dir;
        scalar h;
        scalar FD_J, FD_PD, FD_gDP, FD_gV;
        scalar ADJ_J, ADJ_gDP, ADJ_gV;
        scalar relJ, relDP, relV;
        label nIterPlus, nIterMinus;
        bool valid, signJ, signDP, signV;
    };
    DynamicList<FDRow> fdRows;

    auto scanDir = [&](const word& dirName, const scalarField& Dk,
                       const scalarList& steps,
                       const scalar aJ, const scalar aDP, const scalar aV)
    {
        forAll(steps, si)
        {
            const scalar h = steps[si];
            const FDProbe pp = probeFD(Dk, h, +1.0);
            const FDProbe pm = probeFD(Dk, h, -1.0);

            FDRow row;
            row.dir = dirName;
            row.h = h;
            row.ADJ_J = aJ;
            row.ADJ_gDP = aDP;
            row.ADJ_gV = aV;
            row.nIterPlus = pp.nIter;
            row.nIterMinus = pm.nIter;
            row.valid = pp.valid && pm.valid;
            row.FD_J = GREAT;
            row.FD_PD = GREAT;
            row.FD_gDP = GREAT;
            row.FD_gV = GREAT;
            row.relJ = GREAT;
            row.relDP = GREAT;
            row.relV = GREAT;
            row.signJ = false;
            row.signDP = false;
            row.signV = false;

            if (row.valid)
            {
                row.FD_J = (pp.J - pm.J)/(2.0*h);
                row.FD_PD = (pp.PD - pm.PD)/(2.0*h);
                row.FD_gDP = (pp.gDP - pm.gDP)/(2.0*h);
                row.FD_gV = (pp.gV - pm.gV)/(2.0*h);
                row.relJ =
                    mag(row.FD_J - aJ)/Foam::max(mag(aJ), 1e-30);
                row.relDP =
                    mag(row.FD_gDP - aDP)/Foam::max(mag(aDP), 1e-30);
                row.relV =
                    mag(row.FD_gV - aV)/Foam::max(mag(aV), 1e-30);
                row.signJ =
                    (mag(aJ) <= 1e-12)
                  ? mag(row.FD_J) <= 1e-10
                  : row.FD_J*aJ > 0.0;
                row.signDP =
                    (mag(aDP) <= 1e-12)
                  ? mag(row.FD_gDP) <= 1e-10
                  : row.FD_gDP*aDP > 0.0;
                row.signV =
                    (mag(aV) <= 1e-12)
                  ? mag(row.FD_gV) <= 1e-10
                  : row.FD_gV*aV > 0.0;
            }

            fdRows.append(row);
            if (row.valid)
            {
                Info<< "  " << dirName << " h=" << h
                    << ": FD_J=" << row.FD_J << " ADJ_J=" << aJ
                    << " rel=" << row.relJ << " sign=" << row.signJ
                    << " | FD_gDP=" << row.FD_gDP << " ADJ_gDP=" << aDP
                    << " rel=" << row.relDP << " sign=" << row.signDP
                    << " | FD_gV=" << row.FD_gV << " ADJ_gV=" << aV
                    << " rel=" << row.relV << " sign=" << row.signV
                    << " | nIter(+/-)=" << row.nIterPlus << '/'
                    << row.nIterMinus << endl;
            }
            else
            {
                Info<< "  " << dirName << " h=" << h
                    << ": INVALID finite-difference point; frozen primal "
                    << "iterations(+/-)=" << row.nIterPlus << '/'
                    << row.nIterMinus << endl;
            }
        }
    };

    // BFINAL-012: the D1 diagnostic ladder spans 1e-5 -> 1e-2 (four
    // decades, per the task spec) to expose both the truncation-dominated
    // large-step end and the noise-floor-dominated small-step end.  The
    // formal acceptance steps (accSteps, D1/D2/D3 verdicts) are unchanged
    // so the gate semantics stay identical to the historical B2/B3 runs.
    const scalarList diagSteps = {1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2};
    const scalarList accSteps = {1e-4, 3e-4, 1e-3};

    Info<< "\n--- Stage B2: D1 multi-step diagnostic scan ---" << endl;
    scanDir("D1", D1, diagSteps, projJ_D1, projDP_D1, projV_D1);

    // ---------------------------------------------------------------------
    // Chain decomposition.
    // ---------------------------------------------------------------------
    Info<< "\n--- Stage B2: chain decomposition ---" << endl;
    Info<< "  Volume chain (no flow solve): FD_gV vs <gsensVol,Dk>:" << endl;
    for (label d = 0; d < 3; ++d)
    {
        const scalarField& Dk = d == 0 ? D1 : (d == 1 ? D2 : D3);
        const scalar aV = d == 0 ? projV_D1 : (d == 1 ? projV_D2 : projV_D3);
        const scalar hV = 1e-4;
        const scalar gp = probeVolumeOnly(Dk, hV, +1.0);
        const scalar gm = probeVolumeOnly(Dk, hV, -1.0);
        const scalar FD_gV = (gp - gm)/(2.0*hV);
        const scalar relV =
            mag(FD_gV - aV)/Foam::max(mag(aV), 1e-30);
        Info<< "    D" << d + 1 << ": FD_gV=" << FD_gV
            << " ADJ=" << aV << " rel=" << relV
            << (relV <= 1e-6 ? " PASS" : " FAIL") << endl;
    }

    forAll(fdRows, i)
    {
        if
        (
            fdRows[i].dir == "D1"
         && mag(fdRows[i].h - 1e-3) <= 1e-14
         && fdRows[i].valid
        )
        {
            const scalar lhs = fdRows[i].FD_gDP;
            const scalar rhs =
                fdRows[i].FD_PD/Foam::max(pressureDropMaxPa, SMALL);
            const scalar ident =
                mag(lhs - rhs)/Foam::max(mag(lhs), 1e-30);
            Info<< "  Pressure identity dgDP/dh vs (1/PDmax)dPD/dh: "
                << lhs << " vs " << rhs << " relDiff=" << ident
                << (ident < 1e-12 ? " PASS" : " CHECK") << endl;
        }
    }

    // ---------------------------------------------------------------------
    // B3 formal multi-direction scan.
    // D1 formal rows already exist in diagSteps; add D2/D3 here.
    // ---------------------------------------------------------------------
    Info<< "\n--- Stage B3: multi-direction formal acceptance ---" << endl;
    scanDir("D2", D2, accSteps, projJ_D2, projDP_D2, projV_D2);
    scanDir("D3", D3, accSteps, projJ_D3, projDP_D3, projV_D3);

    auto isFormalStep = [&](const scalar h)
    {
        forAll(accSteps, i)
        {
            if (mag(h - accSteps[i]) <= 1e-14)
            {
                return true;
            }
        }
        return false;
    };

    auto fdSpread = [&](const scalar a, const scalar b)
    {
        return mag(a - b)
          / Foam::max(Foam::max(mag(a), mag(b)), scalar(1e-12));
    };

    struct DirVerdict
    {
        word dir;
        label nRows;
        scalar bestRelJ, bestRelDP, bestRelV;
        bool allValid;
        bool signJ, signDP, signV;
        bool plateauJ, plateauDP, plateauV;
        bool passJ, passDP, passV;
    };

    auto buildVerdict = [&](const word& dirName) -> DirVerdict
    {
        DirVerdict verdict;
        verdict.dir = dirName;
        verdict.nRows = 0;
        verdict.bestRelJ = GREAT;
        verdict.bestRelDP = GREAT;
        verdict.bestRelV = GREAT;
        verdict.allValid = true;
        verdict.signJ = true;
        verdict.signDP = true;
        verdict.signV = true;
        verdict.plateauJ = false;
        verdict.plateauDP = false;
        verdict.plateauV = false;

        forAll(fdRows, i)
        {
            const FDRow& row = fdRows[i];
            if (row.dir != dirName || !isFormalStep(row.h))
            {
                continue;
            }
            ++verdict.nRows;
            verdict.allValid = verdict.allValid && row.valid;
            if (!row.valid)
            {
                continue;
            }
            verdict.bestRelJ = Foam::min(verdict.bestRelJ, row.relJ);
            verdict.bestRelDP = Foam::min(verdict.bestRelDP, row.relDP);
            verdict.bestRelV = Foam::min(verdict.bestRelV, row.relV);
            verdict.signJ = verdict.signJ && row.signJ;
            verdict.signDP = verdict.signDP && row.signDP;
            verdict.signV = verdict.signV && row.signV;
        }
        verdict.allValid = verdict.allValid && verdict.nRows == accSteps.size();

        // Plateau means at least one adjacent formal step pair is stable.
        for (label pairI = 0; pairI < accSteps.size() - 1; ++pairI)
        {
            const FDRow* rowA = nullptr;
            const FDRow* rowB = nullptr;
            forAll(fdRows, i)
            {
                if (fdRows[i].dir != dirName || !fdRows[i].valid)
                {
                    continue;
                }
                if (mag(fdRows[i].h - accSteps[pairI]) <= 1e-14)
                {
                    rowA = &fdRows[i];
                }
                if (mag(fdRows[i].h - accSteps[pairI + 1]) <= 1e-14)
                {
                    rowB = &fdRows[i];
                }
            }
            if (rowA && rowB)
            {
                verdict.plateauJ = verdict.plateauJ
                    || fdSpread(rowA->FD_J, rowB->FD_J) <= 0.30;
                verdict.plateauDP = verdict.plateauDP
                    || fdSpread(rowA->FD_gDP, rowB->FD_gDP) <= 0.30;
                verdict.plateauV = verdict.plateauV
                    || fdSpread(rowA->FD_gV, rowB->FD_gV) <= 0.30;
            }
        }

        verdict.passJ =
               verdict.allValid
            && verdict.signJ
            && verdict.bestRelJ <= 0.05
            && verdict.plateauJ;
        verdict.passDP =
               verdict.allValid
            && verdict.signDP
            && verdict.bestRelDP <= 0.10
            && verdict.plateauDP;
        verdict.passV =
               verdict.allValid
            && verdict.signV
            && verdict.bestRelV <= 1e-6
            && verdict.plateauV;
        return verdict;
    };

    const DirVerdict verdictD1 = buildVerdict("D1");
    const DirVerdict verdictD2 = buildVerdict("D2");
    const DirVerdict verdictD3 = buildVerdict("D3");

    const bool formalConvergenceOK =
           verdictD1.allValid
        && verdictD2.allValid
        && verdictD3.allValid;
    const bool passJ =
           verdictD1.passJ
        && verdictD2.passJ
        && verdictD3.passJ;
    const bool passDP =
           verdictD1.passDP
        && verdictD2.passDP
        && verdictD3.passDP;
    const bool passV =
           verdictD1.passV
        && verdictD2.passV
        && verdictD3.passV;
    const bool plateauOK =
           verdictD1.plateauJ && verdictD1.plateauDP && verdictD1.plateauV
        && verdictD2.plateauJ && verdictD2.plateauDP && verdictD2.plateauV
        && verdictD3.plateauJ && verdictD3.plateauDP && verdictD3.plateauV;
    const bool frozenGradientPass =
        formalConvergenceOK && passJ && passDP && passV && plateauOK;

    auto printVerdict = [&](const DirVerdict& v)
    {
        Info<< "  " << v.dir
            << ": valid=" << v.allValid
            << " bestRel(J/DP/V)=" << v.bestRelJ << '/'
            << v.bestRelDP << '/' << v.bestRelV
            << " sign(J/DP/V)=" << v.signJ << '/'
            << v.signDP << '/' << v.signV
            << " plateau(J/DP/V)=" << v.plateauJ << '/'
            << v.plateauDP << '/' << v.plateauV
            << " pass(J/DP/V)=" << v.passJ << '/'
            << v.passDP << '/' << v.passV << endl;
    };

    Info<< "\nStage B3 per-direction verdicts:" << endl;
    printVerdict(verdictD1);
    printVerdict(verdictD2);
    printVerdict(verdictD3);
    Info<< "  formalConvergenceOK=" << formalConvergenceOK << endl;
    Info<< "  plateauOK=" << plateauOK << endl;
    Info<< "  FROZEN_GRADIENT_STATUS="
        << (frozenGradientPass ? "PASS" : "FAIL") << endl;

    // ---------------------------------------------------------------------
    // Machine-readable outputs.
    // ---------------------------------------------------------------------
    {
        OFstream tsv(b2OutDir/"stage_b2_fd_scan.tsv");
        tsv.precision(16);
        tsv << "Dir\tStep\tValid\tNIterPlus\tNIterMinus\t"
            << "FD_J\tADJ_J\trelJ\tsignJ\t"
            << "FD_gDP\tADJ_gDP\trelDP\tsignDP\t"
            << "FD_gV\tADJ_gV\trelV\tsignV\n";
        forAll(fdRows, i)
        {
            const FDRow& r = fdRows[i];
            tsv << r.dir << '\t' << r.h << '\t' << r.valid << '\t'
                << r.nIterPlus << '\t' << r.nIterMinus << '\t'
                << r.FD_J << '\t' << r.ADJ_J << '\t' << r.relJ << '\t'
                << r.signJ << '\t'
                << r.FD_gDP << '\t' << r.ADJ_gDP << '\t' << r.relDP << '\t'
                << r.signDP << '\t'
                << r.FD_gV << '\t' << r.ADJ_gV << '\t' << r.relV << '\t'
                << r.signV << '\n';
        }

        OFstream sum(b2OutDir/"stage_b2_summary.tsv");
        sum.precision(16);
        sum << "Metric\tValue\n";
        sum << "baselineJ\t" << Jbase << '\n';
        sum << "baselinePD\t" << PDbase << '\n';
        sum << "baselinegDP\t" << gDPbase << '\n';
        sum << "baselinegV\t" << gVbase << '\n';
        sum << "projJ_D1\t" << projJ_D1 << '\n';
        sum << "projJ_D2\t" << projJ_D2 << '\n';
        sum << "projJ_D3\t" << projJ_D3 << '\n';
        sum << "projDP_D1\t" << projDP_D1 << '\n';
        sum << "projDP_D2\t" << projDP_D2 << '\n';
        sum << "projDP_D3\t" << projDP_D3 << '\n';
        sum << "projV_D1\t" << projV_D1 << '\n';
        sum << "projV_D2\t" << projV_D2 << '\n';
        sum << "projV_D3\t" << projV_D3 << '\n';

        auto writeDir = [&](const DirVerdict& v)
        {
            sum << v.dir << "_valid\t" << v.allValid << '\n';
            sum << v.dir << "_bestRelJ\t" << v.bestRelJ << '\n';
            sum << v.dir << "_bestRelDP\t" << v.bestRelDP << '\n';
            sum << v.dir << "_bestRelV\t" << v.bestRelV << '\n';
            sum << v.dir << "_signJ\t" << v.signJ << '\n';
            sum << v.dir << "_signDP\t" << v.signDP << '\n';
            sum << v.dir << "_signV\t" << v.signV << '\n';
            sum << v.dir << "_plateauJ\t" << v.plateauJ << '\n';
            sum << v.dir << "_plateauDP\t" << v.plateauDP << '\n';
            sum << v.dir << "_plateauV\t" << v.plateauV << '\n';
            sum << v.dir << "_passJ\t" << v.passJ << '\n';
            sum << v.dir << "_passDP\t" << v.passDP << '\n';
            sum << v.dir << "_passV\t" << v.passV << '\n';
        };
        writeDir(verdictD1);
        writeDir(verdictD2);
        writeDir(verdictD3);
        sum << "formalConvergenceOK\t" << formalConvergenceOK << '\n';
        sum << "plateauOK\t" << plateauOK << '\n';
        sum << "passJ\t" << passJ << '\n';
        sum << "passDP\t" << passDP << '\n';
        sum << "passV\t" << passV << '\n';
        sum << "status\t" << (frozenGradientPass ? "PASS" : "FAIL") << '\n';
    }

    if (frozenGradientPass)
    {
        OFstream unlock(b2OutDir/"FROZEN_GRADIENT_UNLOCK.txt");
        unlock
            << "# Copy this line to constant/optProperties only after "
            << "reviewing stage_b2_summary.tsv\n"
            << "frozenGradientValidated true;\n";
    }

    Info<< "Stage B2/B3 outputs written to " << b2OutDir << endl;
    restoreStateInto(entryState);
    Info<< "--- Stage B2/B3 complete (state restored) ---\n" << endl;
}
