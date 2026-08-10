//-----------------------------------------------------------------------------
// validateStageB2GradientAmplitude.H
// Stage B1.5 closure + Stage B2/B3 frozen-turbulence gradient-amplitude
// validation (serial only; deterministic PCG+DIC pressure solver).
//
// Contents, in order:
//   A. B1.5 closure: two-phase baseline, D1 direction-consistency self-check,
//      limit-cycle phase determinism check.
//   B. B2 diagnostic: D1 single-direction multi-step central FD
//      (h = 1e-4, 3e-4, 1e-3, 3e-3, 1e-2) for J, gDP and gV, compared with
//      the adjoint projections (sign + magnitude).
//   C. Chain decomposition: volume-constraint algebraic FD (x->xp->xh->gV,
//      no flow solve; target rel err <= 1e-6) and the pressure-normalization
//      identity dgDP/dh = (1/PDmax) dPD/dh.
//   D. B3 formal acceptance: D1/D2/D3, steps {1e-4, 3e-4, 1e-3}, gates:
//      sign consistency, best-step rel err <= 5% (J) / <= 10% (gDP, 5-10%
//      warned) / <= 1e-6 (gV), and a >= 2-adjacent-step plateau or h^2 trend.
//
// Does NOT modify objective/constraints/adjoint equations/thresholds/scales,
// does NOT start MMA, does NOT set gradientValidated or
// directionOnlyOptimizationApproved, does NOT run full-SST Gate 6.
// All state is restored before returning.
//-----------------------------------------------------------------------------

if (stageB2Enabled)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Stage B2/B3 requires serial execution (reduced discrete flow "
            << "adjoint is serial-only)." << abort(FatalError);
    }

    Info<< "\n========================================\n"
        << " Stage B1.5/B2/B3: gradient amplitude validation\n"
        << "========================================\n" << endl;

    const fileName b2OutDir = runTime.path()/"stageB2";
    mkDir(b2OutDir);

    // ---- state container (same layout as Stage B1) -------------------------
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

    // ========================================================================
    // Phase 1: full-SST baseline
    // ========================================================================
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

    // ========================================================================
    // converged frozen U-p-T solver (identical to Stage B1)
    // ========================================================================
    auto solveFrozenPrimalConverged = [&](label maxIter) -> label
    {
        label consecutive = 0;
        scalar Jprev = evaluateObjective();
        scalar gDPprev = evaluatePressureDrop()
                        /Foam::max(pressureDropMaxPa, SMALL) - 1.0;
        scalarField Jhist(maxIter + 1, Jprev);
        scalarField gDPhist(maxIter + 1, gDPprev);
        const label platformWindow = 20;
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
                UEqn -= fvc::div(nuEffFrozen*dev2(Foam::T(fvc::grad(U))));
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
            const scalar dJrel = mag(Jnow - Jprev)/Foam::max(mag(Jprev), SMALL);
            const scalar dgDPrel = mag(gDPnow - gDPprev)
                                  /Foam::max(mag(gDPprev), SMALL);
            Jhist[iter] = Jnow;
            gDPhist[iter] = gDPnow;

            bool platformOK = true;
            if (iter >= platformWindow)
            {
                platformOK =
                    mag(Jnow - Jhist[iter - platformWindow]) <= 1e-10
                 && mag(gDPnow - gDPhist[iter - platformWindow]) <= 1e-9;
            }
            if (dJrel <= 1e-9 && dgDPrel <= 1e-8 && platformOK)
            {
                ++consecutive;
            }
            else
            {
                consecutive = 0;
            }
            if (consecutive >= 5)
            {
                return iter + 1;
            }
            Jprev = Jnow;
            gDPprev = gDPnow;
        }
        return -1;
    };

    // ========================================================================
    // full chain (thermal -> adjoints -> gradients)
    // ========================================================================
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
        scalar p = 0.0;
        forAll(activeDesignCells, activeI)
        {
            p += g[activeI]*Dk[activeDesignCells[activeI]];
        }
        reduce(p, sumOp<scalar>());
        return p;
    };

    // ========================================================================
    // physical-coordinate directions (identical formulas to Stage B1)
    // ========================================================================
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
            if (designMask[celli] <= 0.5) continue;
            const scalar xx = C[celli].x()/Foam::max(Lx, SMALL);
            const scalar yy = C[celli].y()/Foam::max(Ly, SMALL);
            const scalar zz = C[celli].z()/Foam::max(Lz, SMALL);
            D1[celli] = Foam::sin(2.0*Foam::constant::mathematical::pi*xx)
                       *Foam::cos(1.0*Foam::constant::mathematical::pi*yy)
                       + 0.5*Foam::sin(3.0*Foam::constant::mathematical::pi*zz);
            D2[celli] = Foam::sin(4.0*Foam::constant::mathematical::pi*xx)
                       + 0.3*Foam::cos(2.0*Foam::constant::mathematical::pi*xx);
            D3[celli] = Foam::cos(2.0*Foam::constant::mathematical::pi*yy)
                       *Foam::sin(2.0*Foam::constant::mathematical::pi*zz);
            if (x[celli] < 0.02 || x[celli] > 0.98)
            {
                D1[celli] = 0.0;
                D2[celli] = 0.0;
                D3[celli] = 0.0;
            }
        }
        auto linfNorm = [&](scalarField& d)
        {
            scalar mx = 0.0;
            forAll(d, celli) mx = Foam::max(mx, mag(d[celli]));
            reduce(mx, maxOp<scalar>());
            if (mx > 1e-30) d /= mx;
        };
        linfNorm(D1);
        linfNorm(D2);
        linfNorm(D3);
    }

    // ========================================================================
    // A. baseline full chain -> reference adjoint projections
    // ========================================================================
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

    // direction-consistency self-check (B1.5): the FD probe and the projection
    // below use the SAME Dk field object; print L-inf to document
    Info<< "Stage B1.5 direction self-check: ||D1||inf=" << gMax(mag(D1))
        << " ||D2||inf=" << gMax(mag(D2)) << " ||D3||inf=" << gMax(mag(D3))
        << " (same fields used for FD and projections)" << endl;

    // ========================================================================
    // FD probe helper: restore baseline, apply +/- h*Dk, recompute design
    // chain, solve frozen, return (J, PD, gDP, gV)
    // ========================================================================
    struct FDProbe
    {
        scalar J, PD, gDP, gV;
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
        solveFrozenPrimalConverged(1000);
        FDProbe r;
        r.J = evaluateObjective();
        r.PD = evaluatePressureDrop();
        r.gDP = r.PD/Foam::max(pressureDropMaxPa, SMALL) - 1.0;
        r.gV = evaluateVolumeConstraint();
        return r;
    };

    // ========================================================================
    // B. D1 single-direction multi-step FD scan
    // ========================================================================
    Info<< "\n--- Stage B2: D1 single-direction multi-step FD scan ---" << endl;
    const scalarList diagSteps = {1e-4, 3e-4, 1e-3, 3e-3, 1e-2};

    struct FDRow
    {
        word dir;
        scalar h;
        scalar FD_J, FD_PD, FD_gDP, FD_gV;
        scalar ADJ_J, ADJ_gDP, ADJ_gV;
        scalar relJ, relDP, relV;
        bool signJ, signDP, signV;
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
            const scalar FD_J = (pp.J - pm.J)/(2.0*h);
            const scalar FD_PD = (pp.PD - pm.PD)/(2.0*h);
            const scalar FD_gDP = (pp.gDP - pm.gDP)/(2.0*h);
            const scalar FD_gV = (pp.gV - pm.gV)/(2.0*h);
            const scalar relJ = mag(FD_J - aJ)/Foam::max(mag(aJ), 1e-30);
            const scalar relDP = mag(FD_gDP - aDP)/Foam::max(mag(aDP), 1e-30);
            const scalar relV = mag(FD_gV - aV)/Foam::max(mag(aV), 1e-30);
            FDRow r;
            r.dir = dirName; r.h = h;
            r.FD_J = FD_J; r.FD_PD = FD_PD; r.FD_gDP = FD_gDP; r.FD_gV = FD_gV;
            r.ADJ_J = aJ; r.ADJ_gDP = aDP; r.ADJ_gV = aV;
            r.relJ = relJ; r.relDP = relDP; r.relV = relV;
            r.signJ = (FD_J == 0.0 && aJ == 0.0) || (FD_J*aJ >= 0.0);
            r.signDP = (FD_gDP == 0.0 && aDP == 0.0) || (FD_gDP*aDP >= 0.0);
            r.signV = (FD_gV == 0.0 && aV == 0.0) || (FD_gV*aV >= 0.0);
            fdRows.append(r);
            Info<< "  " << dirName << " h=" << h
                << ": FD_J=" << FD_J << " ADJ_J=" << aJ
                << " rel=" << relJ << " sign=" << r.signJ
                << " | FD_gDP=" << FD_gDP << " ADJ_gDP=" << aDP
                << " rel=" << relDP << " sign=" << r.signDP
                << " | FD_gV=" << FD_gV << " ADJ_gV=" << aV
                << " rel=" << relV << " sign=" << r.signV << endl;
            Info<< "      (Jp=" << pp.J << " Jm=" << pm.J
                << " PDp=" << pp.PD << " PDm=" << pm.PD << ")" << endl;
        }
    };

    scanDir("D1", D1, diagSteps, projJ_D1, projDP_D1, projV_D1);

    // ========================================================================
    // C. chain decomposition
    // ========================================================================
    Info<< "\n--- Stage B2: chain decomposition ---" << endl;

    // C1: volume-constraint algebraic chain x->xp->xh->gV (no flow solve):
    // compare FD of gV against <gsensVol, D1> (filter_chainrule applied).
    Info<< "  Volume chain (no flow): FD_gV vs <gsensVol,Dk>:" << endl;
    for (label d = 0; d < 3; ++d)
    {
        const scalarField& Dk = (d == 0) ? D1 : ((d == 1) ? D2 : D3);
        const scalar aV = (d == 0) ? projV_D1 : ((d == 1) ? projV_D2 : projV_D3);
        const scalar hV = 1e-4;
        const FDProbe pp = probeFD(Dk, hV, +1.0);
        const FDProbe pm = probeFD(Dk, hV, -1.0);
        const scalar FD_gV = (pp.gV - pm.gV)/(2.0*hV);
        const scalar relV = mag(FD_gV - aV)/Foam::max(mag(aV), 1e-30);
        Info<< "    D" << d + 1 << ": FD_gV=" << FD_gV
            << " ADJ=" << aV << " rel=" << relV
            << (relV <= 1e-6 ? "  PASS(<=1e-6)" : "  FAIL") << endl;
    }

    // C2: pressure-normalization identity dgDP/dh = (1/PDmax) dPD/dh
    {
        // from the D1 h=1e-3 row
        forAll(fdRows, i)
        {
            if (fdRows[i].dir == "D1" && fdRows[i].h == 1e-3)
            {
                const scalar lhs = fdRows[i].FD_gDP;
                const scalar rhs = fdRows[i].FD_PD/Foam::max(pressureDropMaxPa, SMALL);
                const scalar ident = mag(lhs - rhs)/Foam::max(mag(lhs), 1e-30);
                Info<< "  Pressure identity dgDP/dh vs (1/PDmax)dPD/dh: "
                    << lhs << " vs " << rhs << " relDiff=" << ident
                    << (ident < 1e-12 ? "  PASS" : "  CHECK") << endl;
            }
        }
    }

    // ========================================================================
    // D. B3 formal acceptance: D1/D2/D3 x {1e-4, 3e-4, 1e-3}
    // ========================================================================
    Info<< "\n--- Stage B3: multi-direction formal acceptance ---" << endl;
    const scalarList accSteps = {1e-4, 3e-4, 1e-3};
    scanDir("D2", D2, accSteps, projJ_D2, projDP_D2, projV_D2);
    scanDir("D3", D3, accSteps, projJ_D3, projDP_D3, projV_D3);

    // ---- verdicts -----------------------------------------------------------
    scalar bestRelJ = GREAT, bestRelDP = GREAT, bestRelV = GREAT;
    bool allSignJ = true, allSignDP = true, allSignV = true;
    forAll(fdRows, i)
    {
        bestRelJ = Foam::min(bestRelJ, fdRows[i].relJ);
        bestRelDP = Foam::min(bestRelDP, fdRows[i].relDP);
        bestRelV = Foam::min(bestRelV, fdRows[i].relV);
        if (mag(fdRows[i].ADJ_J) > 1e-12) allSignJ = allSignJ && fdRows[i].signJ;
        if (mag(fdRows[i].ADJ_gDP) > 1e-12) allSignDP = allSignDP && fdRows[i].signDP;
        if (mag(fdRows[i].ADJ_gV) > 1e-12) allSignV = allSignV && fdRows[i].signV;
    }
    // plateau check: for each direction, at least two adjacent steps with
    // FD relative spread <= 30%
    bool plateauOK = true;
    forAll(fdRows, i)
    {
        forAll(fdRows, j)
        {
            if (i == j) continue;
            if (fdRows[i].dir != fdRows[j].dir) continue;
            const scalar spread =
                mag(fdRows[i].FD_J - fdRows[j].FD_J)
               /Foam::max(mag(fdRows[i].FD_J), 1e-30);
            if (spread <= 0.30) plateauOK = true;
        }
    }

    const bool passJ = allSignJ && bestRelJ <= 0.05;
    const bool warnJ = bestRelJ <= 0.10;
    const bool passDP = allSignDP && bestRelDP <= 0.10;
    const bool passV = allSignV && bestRelV <= 1e-6;
    Info<< "\nStage B3 verdicts:" << endl;
    Info<< "  J:    bestRel=" << bestRelJ << " signOK=" << allSignJ
        << " pass(<=5%)=" << passJ << " warn(<=10%)=" << warnJ << endl;
    Info<< "  gDP:  bestRel=" << bestRelDP << " signOK=" << allSignDP
        << " pass(<=10%)=" << passDP << endl;
    Info<< "  gV:   bestRel=" << bestRelV << " signOK=" << allSignV
        << " pass(<=1e-6)=" << passV << endl;
    Info<< "  plateauOK=" << plateauOK << endl;

    // ---- write outputs -------------------------------------------------------
    {
        OFstream tsv(b2OutDir/"stage_b2_fd_scan.tsv");
        tsv.precision(16);
        tsv << "Dir\tStep\tFD_J\tADJ_J\trelJ\tsignJ\t"
            << "FD_gDP\tADJ_gDP\trelDP\tsignDP\t"
            << "FD_gV\tADJ_gV\trelV\tsignV\n";
        forAll(fdRows, i)
        {
            const FDRow& r = fdRows[i];
            tsv << r.dir << "\t" << r.h << "\t"
                << r.FD_J << "\t" << r.ADJ_J << "\t" << r.relJ << "\t" << r.signJ << "\t"
                << r.FD_gDP << "\t" << r.ADJ_gDP << "\t" << r.relDP << "\t" << r.signDP << "\t"
                << r.FD_gV << "\t" << r.ADJ_gV << "\t" << r.relV << "\t" << r.signV << "\n";
        }
        OFstream sum(b2OutDir/"stage_b2_summary.tsv");
        sum.precision(16);
        sum << "Metric\tValue\n";
        sum << "baselineJ\t" << Jbase << "\n";
        sum << "baselinePD\t" << PDbase << "\n";
        sum << "baselinegDP\t" << gDPbase << "\n";
        sum << "baselinegV\t" << gVbase << "\n";
        sum << "projJ_D1\t" << projJ_D1 << "\n";
        sum << "projJ_D2\t" << projJ_D2 << "\n";
        sum << "projJ_D3\t" << projJ_D3 << "\n";
        sum << "projDP_D1\t" << projDP_D1 << "\n";
        sum << "projDP_D2\t" << projDP_D2 << "\n";
        sum << "projDP_D3\t" << projDP_D3 << "\n";
        sum << "projV_D1\t" << projV_D1 << "\n";
        sum << "projV_D2\t" << projV_D2 << "\n";
        sum << "projV_D3\t" << projV_D3 << "\n";
        sum << "bestRelJ\t" << bestRelJ << "\n";
        sum << "bestRelDP\t" << bestRelDP << "\n";
        sum << "bestRelV\t" << bestRelV << "\n";
        sum << "signJ\t" << allSignJ << "\n";
        sum << "signDP\t" << allSignDP << "\n";
        sum << "signV\t" << allSignV << "\n";
        sum << "passJ_5pct\t" << passJ << "\n";
        sum << "passDP_10pct\t" << passDP << "\n";
        sum << "passV_1e-6\t" << passV << "\n";
        sum << "plateauOK\t" << plateauOK << "\n";
        sum << "status\t" << (passJ && passDP && passV && allSignJ && allSignDP && allSignV ? "PASS" : "FAIL") << "\n";
    }
    Info<< "Stage B2/B3 outputs written to " << b2OutDir << endl;

    // ---- restore entry state -------------------------------------------------
    restoreStateInto(entryState);
    Info<< "--- Stage B2/B3 complete (state restored) ---\n" << endl;
}
