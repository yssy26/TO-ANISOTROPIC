//-----------------------------------------------------------------------------
// validateStageBRepeatability.H
// Stage B1: frozen-turbulence baseline repeatability and numerical noise-floor
// validation (serial only; the reduced discrete flow adjoint has no processor
// patch transpose).
//
// Runs only when stageBEnabled is true.  Does NOT change the design, does NOT
// touch MMA, does NOT set gradientValidated, does NOT run the full SST Gate 6,
// and does not modify any objective/constraint/adjoint equation or FD
// threshold.  All state is restored before returning, so the downstream
// validation modules (validateSSTDirection.H, validateGate6SSTStrict.H) see
// the same fields as without this module.
//-----------------------------------------------------------------------------

if (stageBEnabled)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Stage B1 requires serial execution: the reduced discrete "
            << "flow adjoint (solveDiscreteFlowAdjoint.H) has no processor-"
            << "patch transpose." << abort(FatalError);
    }

    Info<< "\n========================================\n"
        << " Stage B1: frozen-turbulence repeatability & noise floor\n"
        << "========================================\n" << endl;

    const fileName b1OutDir = runTime.path()/"stageB1";

    // ---- state container ---------------------------------------------------
    struct B1State
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

    // ---- save the entry state so downstream modules are unaffected --------
    const B1State entryState =
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

    auto restoreStateInto = [&](const B1State& s)
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

    // ---- field statistics (min/max/L1/L2/stable checksum) ------------------
    auto fieldStatsScalar = [&](const char* name, const scalarField& f)
    {
        scalar mn = GREAT, mx = -GREAT;
        scalar l1 = 0.0, l2 = 0.0, nBad = 0.0;
        unsigned long long hsh = 1469598103934665603ULL;
        forAll(f, i)
        {
            if (std::isfinite(f[i]))
            {
                if (f[i] < mn) mn = f[i];
                if (f[i] > mx) mx = f[i];
                const scalar m2 = mag(f[i]);
                l1 += m2;
                l2 += (m2 < 1e150) ? m2*m2 : 1e300;
            }
            else
            {
                nBad += 1.0;
            }
            const double* p = reinterpret_cast<const double*>(&f[i]);
            for (size_t b = 0; b < sizeof(double); ++b)
            {
                hsh ^= static_cast<unsigned char>(p[b]);
                hsh *= 1099511628211ULL;
            }
        }
        reduce(mn, minOp<scalar>());
        reduce(mx, maxOp<scalar>());
        reduce(l1, sumOp<scalar>());
        reduce(l2, sumOp<scalar>());
        reduce(nBad, sumOp<scalar>());
        std::ostringstream ossCks;
        ossCks << std::hex << hsh;
        Info<< "  stat[" << name << "] min=" << mn
            << " max=" << mx << " L1=" << l1
            << " L2=" << Foam::sqrt(l2) << " nBad=" << nBad
            << " cksum=0x" << ossCks.str().c_str() << endl;
    };
    auto fieldStatsVector = [&](const char* name, const vectorField& f)
    {
        scalar mn = GREAT, mx = -GREAT;
        scalar l1 = 0.0, l2 = 0.0, nBad = 0.0;
        unsigned long long hsh = 1469598103934665603ULL;
        forAll(f, i)
        {
            const scalar m2 = mag(f[i]);
            if (std::isfinite(m2))
            {
                if (m2 < mn) mn = m2;
                if (m2 > mx) mx = m2;
                l1 += m2;
                l2 += (m2 < 1e150) ? m2*m2 : 1e300;
            }
            else
            {
                nBad += 1.0;
            }
            const double* p = reinterpret_cast<const double*>(&f[i]);
            for (size_t b = 0; b < sizeof(vector); ++b)
            {
                hsh ^= static_cast<unsigned char>(p[b]);
                hsh *= 1099511628211ULL;
            }
        }
        reduce(mn, minOp<scalar>());
        reduce(mx, maxOp<scalar>());
        reduce(l1, sumOp<scalar>());
        reduce(l2, sumOp<scalar>());
        reduce(nBad, sumOp<scalar>());
        std::ostringstream ossCks;
        ossCks << std::hex << hsh;
        Info<< "  stat[" << name << "] min=" << mn
            << " max=" << mx << " L1=" << l1
            << " L2=" << Foam::sqrt(l2) << " nBad=" << nBad
            << " cksum=0x" << ossCks.str().c_str() << endl;
    };
    auto fieldStatsTensor = [&](const char* name, const symmTensorField& f)
    {
        scalar mn = GREAT, mx = -GREAT;
        scalar l1 = 0.0, l2 = 0.0, nBad = 0.0;
        unsigned long long hsh = 1469598103934665603ULL;
        forAll(f, i)
        {
            const scalar m2 = mag(f[i]);
            if (std::isfinite(m2))
            {
                if (m2 < mn) mn = m2;
                if (m2 > mx) mx = m2;
                l1 += m2;
                l2 += (m2 < 1e150) ? m2*m2 : 1e300;
            }
            else
            {
                nBad += 1.0;
            }
            const double* p = reinterpret_cast<const double*>(&f[i]);
            for (size_t b = 0; b < sizeof(symmTensor); ++b)
            {
                hsh ^= static_cast<unsigned char>(p[b]);
                hsh *= 1099511628211ULL;
            }
        }
        reduce(mn, minOp<scalar>());
        reduce(mx, maxOp<scalar>());
        reduce(l1, sumOp<scalar>());
        reduce(l2, sumOp<scalar>());
        reduce(nBad, sumOp<scalar>());
        std::ostringstream ossCks;
        ossCks << std::hex << hsh;
        Info<< "  stat[" << name << "] min=" << mn
            << " max=" << mx << " L1=" << l1
            << " L2=" << Foam::sqrt(l2) << " nBad=" << nBad
            << " cksum=0x" << ossCks.str().c_str() << endl;
    };

    // ========================================================================
    // Phase 1: full-SST baseline (turbulence active)
    // ========================================================================
    Info<< "Stage B1 phase 1: full-SST baseline solve (freezeTurbulence=false)"
        << endl;
    freezeTurbulenceForValidation = false;
    const label nFullSST = solveFullSST(1000, 1e-6);
    if (nFullSST < 0)
    {
        FatalErrorInFunction
            << "Stage B1 phase 1: full-SST baseline did not converge."
            << abort(FatalError);
    }
    Info<< "Stage B1 phase 1: full-SST converged in " << nFullSST
        << " iterations." << endl;

    // ========================================================================
    // Phase 2: freeze turbulence, save baseline + statistics
    // ========================================================================
    freezeTurbulenceForValidation = true;
    #include "updateFrozenTurbulenceFields.H"

    const B1State baselineState =
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

    Info<< "Stage B1 phase 2: baseline field statistics" << endl;
    fieldStatsVector("U", U.primitiveField());
    fieldStatsScalar("p", p.primitiveField());
    fieldStatsScalar("T", T.primitiveField());
    fieldStatsScalar("k", mesh.lookupObject<volScalarField>("k").primitiveField());
    fieldStatsScalar("omega", mesh.lookupObject<volScalarField>("omega").primitiveField());
    fieldStatsScalar("nut", mesh.lookupObject<volScalarField>("nut").primitiveField());
    fieldStatsScalar("nutFrozen", nutFrozen.primitiveField());
    fieldStatsScalar("nuEffFrozen", nuEffFrozen.primitiveField());
    fieldStatsScalar("alphaTurbulent", alphaTurbulent.primitiveField());
    fieldStatsTensor("DTMolecular", DTMolecular.primitiveField());
    fieldStatsTensor("DTEffective", DTEffective.primitiveField());
    fieldStatsScalar("x", x.primitiveField());
    fieldStatsScalar("xp", xp.primitiveField());
    fieldStatsScalar("xh", xh.primitiveField());
    fieldStatsScalar("alpha", alpha.primitiveField());

    const scalar Jbase = evaluateObjective();
    const scalar DPbase = evaluatePressureDrop();
    const scalar gDPbase = DPbase/Foam::max(pressureDropMaxPa, SMALL) - 1.0;
    const scalar gVbase = evaluateVolumeConstraint();
    Info<< "Stage B1 baseline: J=" << Jbase << " PressureDropPa=" << DPbase
        << " gDP=" << gDPbase << " gV=" << gVbase << endl;

    // ========================================================================
    // converged frozen U-p-T solver (replaces fixed-iteration practice)
    // ========================================================================
    auto solveFrozenPrimalConverged = [&](label maxIter) -> label
    {
        label consecutive = 0;
        scalar Jprev = evaluateObjective();
        scalar gDPprev = evaluatePressureDrop()
                        /Foam::max(pressureDropMaxPa, SMALL) - 1.0;
        // Long-window platform check: the iteration-to-iteration criterion
        // (dJrel <= 1e-9) can trigger while J still drifts on a longer
        // scale.  Keep the last 21 J/gDP values and require both the local
        // criterion and a 20-iteration platform window <= 1e-10 (J) and
        // <= 1e-9 (gDP) before declaring convergence.
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
            SolverPerformance<vector> perfU;
            if (simple.momentumPredictor())
            {
                perfU = solve(UEqn == -fvc::grad(p));
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

            SolverPerformance<scalar> perfP;
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
                );
                pEqn.setReference(pRefCell, pRefValue);
                perfP = pEqn.solve();
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
            SolverPerformance<scalar> perfT = TEqnLocal.solve();

            scalar contErr = 0.0;
            {
                scalar sumPhi = 0.0;
                forAll(phi.primitiveField(), f)
                {
                    sumPhi += phi.primitiveField()[f];
                }
                forAll(mesh.boundary(), patchi)
                {
                    sumPhi += sum(phi.boundaryField()[patchi]);
                }
                reduce(sumPhi, sumOp<scalar>());
                scalar flowScale = 0.0;
                forAll(phi.primitiveField(), f)
                {
                    flowScale = Foam::max(flowScale, mag(phi.primitiveField()[f]));
                }
                reduce(flowScale, maxOp<scalar>());
                contErr = mag(sumPhi)/Foam::max(flowScale, VSMALL);
            }

            const scalar Jnow = evaluateObjective();
            const scalar gDPnow = evaluatePressureDrop()
                                  /Foam::max(pressureDropMaxPa, SMALL) - 1.0;
            const scalar dJrel = mag(Jnow - Jprev)/Foam::max(mag(Jprev), SMALL);
            const scalar dgDPrel = mag(gDPnow - gDPprev)
                                  /Foam::max(mag(gDPprev), SMALL);
            Jhist[iter] = Jnow;
            gDPhist[iter] = gDPnow;

            if ((iter + 1) % 20 == 0 || iter == 0)
            {
                Info<< "  B1 frozen iter " << iter + 1
                    << ": Ures=" << mag(perfU.finalResidual())
                    << " Pres=" << mag(perfP.finalResidual())
                    << " Tres=" << mag(perfT.finalResidual())
                    << " contErr=" << contErr
                    << " dJrel=" << dJrel << " dgDPrel=" << dgDPrel
                    << " J=" << Jnow << " gDP=" << gDPnow << endl;
            }

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
                Info<< "  B1 frozen converged after " << iter + 1
                    << " iterations." << endl;
                U.storePrevIter();
                p.storePrevIter();
                phi.storePrevIter();
                return iter + 1;
            }
            Jprev = Jnow;
            gDPprev = gDPnow;
        }
        WarningInFunction
            << "Stage B1 frozen solve did NOT converge within " << maxIter
            << " iterations. Reporting the reached platform value." << endl;
        return -1;
    };

    // ========================================================================
    // full chain: thermal flux -> T -> objective derivs -> Tb -> Ub -> Uc ->
    // gradients (block-scoped includes; no effect on the main loop scope)
    // ========================================================================
    auto runFullB1Chain = [&]() -> scalar
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
        return evaluateObjective();
    };

    // ========================================================================
    // gradient projection helpers (MMA-consistent, no mesh.V())
    // ========================================================================
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

    struct B1Trial
    {
        word tag;
        label nIter;
        scalar tbRes, ubRes;
        scalar J, DP, gDP, gV;
        scalar df1, df2, dfi;
        scalar dg01, dg02, dg0i;
        scalar dg11, dg12, dg1i;
        scalar projD1, projD2, projD3;
        scalar projDP1, projDP2, projDP3;
        scalar projV1, projV2, projV3;
    };
    DynamicList<B1Trial> trialsA;
    DynamicList<B1Trial> trialsB;

    auto recordTrial = [&](B1Trial& tr, const word& lbl)
    {
        tr.tag = lbl;
        tr.tbRes = adjointResidualT;
        tr.ubRes = adjointResidualU;
        tr.J = evaluateObjective();
        tr.DP = evaluatePressureDrop();
        tr.gDP = tr.DP/Foam::max(pressureDropMaxPa, SMALL) - 1.0;
        tr.gV = evaluateVolumeConstraint();
        tr.df1 = 0.0; tr.df2 = 0.0; tr.dfi = 0.0;
        tr.dg01 = 0.0; tr.dg02 = 0.0; tr.dg0i = 0.0;
        tr.dg11 = 0.0; tr.dg12 = 0.0; tr.dg1i = 0.0;
        forAll(activeDesignCells, activeI)
        {
            tr.df1 += mag(dfdx[activeI]);
            tr.df2 += sqr(dfdx[activeI]);
            tr.dfi = Foam::max(tr.dfi, mag(dfdx[activeI]));
            tr.dg01 += mag(dgdx[0][activeI]);
            tr.dg02 += sqr(dgdx[0][activeI]);
            tr.dg0i = Foam::max(tr.dg0i, mag(dgdx[0][activeI]));
            tr.dg11 += mag(dgdx[1][activeI]);
            tr.dg12 += sqr(dgdx[1][activeI]);
            tr.dg1i = Foam::max(tr.dg1i, mag(dgdx[1][activeI]));
        }
        reduce(tr.df1, sumOp<scalar>());
        reduce(tr.df2, sumOp<scalar>());
        reduce(tr.dfi, maxOp<scalar>());
        reduce(tr.dg01, sumOp<scalar>());
        reduce(tr.dg02, sumOp<scalar>());
        reduce(tr.dg0i, maxOp<scalar>());
        reduce(tr.dg11, sumOp<scalar>());
        reduce(tr.dg12, sumOp<scalar>());
        reduce(tr.dg1i, maxOp<scalar>());
        tr.df2 = Foam::sqrt(tr.df2);
        tr.dg02 = Foam::sqrt(tr.dg02);
        tr.dg12 = Foam::sqrt(tr.dg12);
    };

    // ========================================================================
    // physical-coordinate design directions D1/D2/D3 (L-inf normalized)
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
        Info<< "Stage B1 directions: ||D1||inf=" << gMax(mag(D1))
            << " ||D2||inf=" << gMax(mag(D2))
            << " ||D3||inf=" << gMax(mag(D3)) << endl;
    }

    // ========================================================================
    // Test A: identical initial state, stageBRepeatA full chains
    // ========================================================================
    Info<< "\n--- Stage B1 Test A: identical initial state ---" << endl;
    for (label rep = 0; rep < stageBRepeatA; ++rep)
    {
        restoreStateInto(baselineState);
        Info<< "  Test A rep " << rep + 1 << ": solving frozen primal..." << endl;
        const label nIt = solveFrozenPrimalConverged(1000);
        B1Trial tr;
        tr.nIter = nIt;
        runFullB1Chain();
        recordTrial(tr, "A" + Foam::name(rep + 1));
        trialsA.append(tr);
        Info<< "  Test A rep " << rep + 1 << ": J=" << tr.J
            << " gDP=" << tr.gDP << " gV=" << tr.gV
            << " |dfdx|2=" << tr.df2
            << " |dgdx0|2=" << tr.dg02
            << " |dgdx1|2=" << tr.dg12 << endl;
    }

    // ========================================================================
    // Test B: perturbed U/p/T initial state, stageBRepeatB directions
    // ========================================================================
    Info<< "\n--- Stage B1 Test B: perturbed initial state ---" << endl;
    const vectorField& CC = mesh.C().primitiveField();
    for (label rep = 0; rep < stageBRepeatB; ++rep)
    {
        restoreStateInto(baselineState);
        forAll(U, celli)
        {
            const scalar s1 = Foam::sin(1.7*CC[celli].x() + 0.3*rep);
            const scalar s2 = Foam::cos(1.1*CC[celli].y() + 0.2*rep);
            const scalar s3 = Foam::sin(0.9*CC[celli].z() + 0.4*rep);
            U[celli] += 1e-4*U[celli].x()*vector(s1, s2, s3);
            p[celli] += 1e-4*mag(p[celli])*s2;
            T[celli] += 1e-5*T[celli]*s3;
        }
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();
        T.correctBoundaryConditions();
        Info<< "  Test B rep " << rep + 1 << ": solving frozen primal..." << endl;
        const label nIt = solveFrozenPrimalConverged(1000);
        B1Trial tr;
        tr.nIter = nIt;
        runFullB1Chain();
        recordTrial(tr, "B" + Foam::name(rep + 1));
        trialsB.append(tr);
        Info<< "  Test B rep " << rep + 1 << ": J=" << tr.J
            << " gDP=" << tr.gDP << " |dfdx|2=" << tr.df2 << endl;
    }

    // ========================================================================
    // projections onto D1/D2/D3 (fill into the recorded trials)
    // ========================================================================
    forAll(trialsA, i)
    {
        B1Trial& tr = trialsA[i];
        tr.projD1 = gradientProjection(dfdx, D1);
        tr.projD2 = gradientProjection(dfdx, D2);
        tr.projD3 = gradientProjection(dfdx, D3);
        tr.projV1 = gradientProjection(dgdx[0], D1);
        tr.projV2 = gradientProjection(dgdx[0], D2);
        tr.projV3 = gradientProjection(dgdx[0], D3);
        tr.projDP1 = gradientProjection(dgdx[1], D1);
        tr.projDP2 = gradientProjection(dgdx[1], D2);
        tr.projDP3 = gradientProjection(dgdx[1], D3);
    }
    forAll(trialsB, i)
    {
        B1Trial& tr = trialsB[i];
        tr.projD1 = gradientProjection(dfdx, D1);
        tr.projD2 = gradientProjection(dfdx, D2);
        tr.projD3 = gradientProjection(dfdx, D3);
        tr.projV1 = gradientProjection(dgdx[0], D1);
        tr.projV2 = gradientProjection(dgdx[0], D2);
        tr.projV3 = gradientProjection(dgdx[0], D3);
        tr.projDP1 = gradientProjection(dgdx[1], D1);
        tr.projDP2 = gradientProjection(dgdx[1], D2);
        tr.projDP3 = gradientProjection(dgdx[1], D3);
    }

    // ========================================================================
    // write trial TSV
    // ========================================================================
    mkDir(b1OutDir);
    {
        OFstream tsva(b1OutDir/"stage_b1_repeatability.tsv");
        tsva.precision(16);
        tsva << "Trial\tIter\tTbRes\tUbRes\tJ\tPressureDropPa\tgDP\tgV"
             << "\t|dfdx|1\t|dfdx|2\t|dfdx|inf"
             << "\t|dgdx0|1\t|dgdx0|2\t|dgdx0|inf"
             << "\t|dgdx1|1\t|dgdx1|2\t|dgdx1|inf"
             << "\tprojD1\tprojD2\tprojD3"
             << "\tprojV1\tprojV2\tprojV3"
             << "\tprojDP1\tprojDP2\tprojDP3\n";
        forAll(trialsA, i)
        {
            const B1Trial& tr = trialsA[i];
            tsva << tr.tag << "\t" << tr.nIter << "\t"
                 << tr.tbRes << "\t" << tr.ubRes << "\t" << tr.J << "\t"
                 << tr.DP << "\t" << tr.gDP << "\t" << tr.gV << "\t"
                 << tr.df1 << "\t" << tr.df2 << "\t" << tr.dfi << "\t"
                 << tr.dg01 << "\t" << tr.dg02 << "\t" << tr.dg0i << "\t"
                 << tr.dg11 << "\t" << tr.dg12 << "\t" << tr.dg1i << "\t"
                 << tr.projD1 << "\t" << tr.projD2 << "\t" << tr.projD3 << "\t"
                 << tr.projV1 << "\t" << tr.projV2 << "\t" << tr.projV3 << "\t"
                 << tr.projDP1 << "\t" << tr.projDP2 << "\t" << tr.projDP3
                 << "\n";
        }
        forAll(trialsB, i)
        {
            const B1Trial& tr = trialsB[i];
            tsva << tr.tag << "\t" << tr.nIter << "\t"
                 << tr.tbRes << "\t" << tr.ubRes << "\t" << tr.J << "\t"
                 << tr.DP << "\t" << tr.gDP << "\t" << tr.gV << "\t"
                 << tr.df1 << "\t" << tr.df2 << "\t" << tr.dfi << "\t"
                 << tr.dg01 << "\t" << tr.dg02 << "\t" << tr.dg0i << "\t"
                 << tr.dg11 << "\t" << tr.dg12 << "\t" << tr.dg1i << "\t"
                 << tr.projD1 << "\t" << tr.projD2 << "\t" << tr.projD3 << "\t"
                 << tr.projV1 << "\t" << tr.projV2 << "\t" << tr.projV3 << "\t"
                 << tr.projDP1 << "\t" << tr.projDP2 << "\t" << tr.projDP3
                 << "\n";
        }
    }

    // ========================================================================
    // frozen-field drift check vs baseline
    // ========================================================================
    Info<< "\n--- Stage B1 frozen-field drift check ---" << endl;
    auto driftMetrics = [&](const char* name, const volScalarField& f,
                            const volScalarField& fb)
    {
        scalar maxAbs = 0.0;
        scalar maxRel = 0.0;
        scalar l2num = 0.0;
        scalar l2den = 0.0;
        forAll(f, celli)
        {
            const scalar d = f[celli] - fb[celli];
            maxAbs = Foam::max(maxAbs, mag(d));
            maxRel = Foam::max
            (
                maxRel,
                mag(d)/Foam::max(mag(fb[celli]), VSMALL)
            );
            l2num += sqr(d);
            l2den += sqr(fb[celli]);
        }
        reduce(maxAbs, maxOp<scalar>());
        reduce(maxRel, maxOp<scalar>());
        reduce(l2num, sumOp<scalar>());
        reduce(l2den, sumOp<scalar>());
        Info<< "  drift[" << name << "] maxAbs=" << maxAbs
            << " maxRel=" << maxRel
            << " L2rel=" << Foam::sqrt(l2num/Foam::max(l2den, VSMALL)) << endl;
    };
    driftMetrics("nutFrozen", nutFrozen, baselineState.nutFrozen);
    driftMetrics("nuEffFrozen", nuEffFrozen, baselineState.nuEffFrozen);
    driftMetrics("alphaTurbulent", alphaTurbulent, baselineState.alphaTurbulent);
    driftMetrics("k", mesh.lookupObject<volScalarField>("k"), baselineState.k);
    driftMetrics("omega", mesh.lookupObject<volScalarField>("omega"), baselineState.omega);
    driftMetrics("nut", mesh.lookupObject<volScalarField>("nut"), baselineState.nut);

    // ========================================================================
    // SNR probe: h=1e-3 along D1
    // ========================================================================
    Info<< "\n--- Stage B1 SNR probe h=1e-3 along D1 ---" << endl;
    const scalar hSNR = 1e-3;
    DynamicList<scalar> Jplus, Jminus, DPplus, DPminus;
    for (label rep = 0; rep < stageBSNRRepeats; ++rep)
    {
        restoreStateInto(baselineState);
        forAll(D1, celli)
        {
            x.primitiveFieldRef()[celli] += hSNR*D1[celli];
        }
        recomputeDesignChain();
        solveFrozenPrimalConverged(1000);
        Jplus.append(evaluateObjective());
        DPplus.append(evaluatePressureDrop());

        restoreStateInto(baselineState);
        forAll(D1, celli)
        {
            x.primitiveFieldRef()[celli] -= hSNR*D1[celli];
        }
        recomputeDesignChain();
        solveFrozenPrimalConverged(1000);
        Jminus.append(evaluateObjective());
        DPminus.append(evaluatePressureDrop());
        Info<< "  SNR rep " << rep + 1 << ": Jplus=" << Jplus[rep]
            << " Jminus=" << Jminus[rep]
            << " DPplus=" << DPplus[rep]
            << " DPminus=" << DPminus[rep] << endl;
    }

    auto rangeOf = [&](const DynamicList<scalar>& v)
    {
        scalar mn = GREAT, mx = -GREAT;
        forAll(v, i) { mn = Foam::min(mn, v[i]); mx = Foam::max(mx, v[i]); }
        return mx - mn;
    };
    auto meanOf = [&](const DynamicList<scalar>& v)
    {
        scalar s = 0.0;
        forAll(v, i) s += v[i];
        return s/Foam::max(v.size(), 1);
    };
    scalar baseRangeJ = 0.0;
    scalar baseRangeDP = 0.0;
    if (trialsA.size() > 1)
    {
        scalar mnJ = GREAT, mxJ = -GREAT;
        scalar mnD = GREAT, mxD = -GREAT;
        forAll(trialsA, i)
        {
            mnJ = Foam::min(mnJ, trialsA[i].J);
            mxJ = Foam::max(mxJ, trialsA[i].J);
            mnD = Foam::min(mnD, trialsA[i].DP);
            mxD = Foam::max(mxD, trialsA[i].DP);
        }
        baseRangeJ = mxJ - mnJ;
        baseRangeDP = mxD - mnD;
    }
    const scalar signalJ = mag(meanOf(Jplus) - meanOf(Jminus));
    const scalar signalDP = mag(meanOf(DPplus) - meanOf(DPminus));
    const scalar noiseJ = Foam::max
    (
        Foam::max(rangeOf(Jplus), rangeOf(Jminus)),
        baseRangeJ
    );
    const scalar noiseDP = Foam::max
    (
        Foam::max(rangeOf(DPplus), rangeOf(DPminus)),
        baseRangeDP
    );
    const scalar SNR_J = signalJ/Foam::max(noiseJ, 1e-30);
    const scalar SNR_DP = signalDP/Foam::max(noiseDP, 1e-30);
    Info<< "Stage B1 SNR: signalJ=" << signalJ << " noiseJ=" << noiseJ
        << " SNR_J=" << SNR_J << endl;
    Info<< "Stage B1 SNR: signalDP=" << signalDP << " noiseDP=" << noiseDP
        << " SNR_DP=" << SNR_DP << endl;

    {
        OFstream tsnn(b1OutDir/"stage_b1_noise_floor.tsv");
        tsnn.precision(16);
        tsnn << "Metric\tValue\n";
        tsnn << "baselineJ\t" << Jbase << "\n";
        tsnn << "baselinePressureDropPa\t" << DPbase << "\n";
        tsnn << "baselinegDP\t" << gDPbase << "\n";
        tsnn << "baselinegV\t" << gVbase << "\n";
        tsnn << "baselineRangeJ\t" << baseRangeJ << "\n";
        tsnn << "baselineRangeDP\t" << baseRangeDP << "\n";
        tsnn << "SNRprobeH\t" << hSNR << "\n";
        forAll(Jplus, i) tsnn << "Jplus" << (i + 1) << "\t" << Jplus[i] << "\n";
        forAll(Jminus, i) tsnn << "Jminus" << (i + 1) << "\t" << Jminus[i] << "\n";
        forAll(DPplus, i) tsnn << "DPplus" << (i + 1) << "\t" << DPplus[i] << "\n";
        forAll(DPminus, i) tsnn << "DPminus" << (i + 1) << "\t" << DPminus[i] << "\n";
        tsnn << "signalJ\t" << signalJ << "\n";
        tsnn << "noiseJ\t" << noiseJ << "\n";
        tsnn << "SNR_J\t" << SNR_J << "\n";
        tsnn << "signalDP\t" << signalDP << "\n";
        tsnn << "noiseDP\t" << noiseDP << "\n";
        tsnn << "SNR_DP\t" << SNR_DP << "\n";
    }

    Info<< "Stage B1 outputs written to " << b1OutDir << endl;

    // ---- restore entry state (downstream modules unaffected) ---------------
    restoreStateInto(entryState);
    Info<< "--- Stage B1 complete (state restored) ---\n" << endl;
}
