const label flowCorrectorLimit =
    ransFlowModel ? turbulenceMaxCorrectors : flowInnerIter;
label turbulenceConvergedCount = 0;
turbulenceConverged = !ransFlowModel;

vectorField previousU;
scalarField previousP;
scalarField previousK;
scalarField previousOmega;
scalarField previousNut;

if (ransFlowModel)
{
    previousU = U.primitiveField();
    previousP = p.primitiveField();
    previousK =
        mesh.lookupObject<volScalarField>("k").primitiveField();
    previousOmega =
        mesh.lookupObject<volScalarField>("omega").primitiveField();
    previousNut = nutFrozen.primitiveField();
}

for (i = 0; i < flowCorrectorLimit; ++i)
{
    forwardPrimalOuterIterations = i + 1;
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
    Info<< "NS.H solving, corrector " << i + 1 << endl;
    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        SolverPerformance<vector> perfU =
            solve(UEqn == -fvc::grad(p));
        forwardInitialResidualU = mag(perfU.initialResidual());
        forwardFinalResidualU = mag(perfU.finalResidual());
        fvOptions.correct(U);
    }

    {
        const scalarField& aDiag = UEqn.diag();
        scalar aMin = GREAT, aMax = 0.0, aSum = 0.0;
        forAll(aDiag, celli)
        {
            aSum += aDiag[celli];
            aMin = Foam::min(aMin, aDiag[celli]);
            aMax = Foam::max(aMax, aDiag[celli]);
        }
        reduce(aSum, sumOp<scalar>());
        reduce(aMin, minOp<scalar>());
        reduce(aMax, maxOp<scalar>());
        Info<< "NS UEqn diag: avg=" << aSum/Foam::max(scalar(mesh.nCells()), 1)
            << " min=" << aMin << " max=" << aMax << endl;
        const auto& icList = UEqn.internalCoeffs();
        forAll(icList, patchi)
        {
            if (icList[patchi].size() > 0)
            {
                const scalar icMax = Foam::max(Foam::mag(icList[patchi]));
                Info<< "NS UEqn IC patch " << patchi
                    << " (" << mesh.boundary()[patchi].name() << ")"
                    << " max=" << icMax << endl;
            }
        }
    }
    {
        const scalarField aA = UEqn.A();
        const scalarField dD = UEqn.diag();
        scalar aSum = 0.0, aMin = GREAT, aMax = 0.0;
        scalar dSum = 0.0, dMin = GREAT, dMax = 0.0;
        scalar fvSum = 0.0, fvMin = GREAT, fvMax = 0.0;
        forAll(aA, celli)
        {
            aSum += aA[celli];
            aMin = Foam::min(aMin, aA[celli]);
            aMax = Foam::max(aMax, aA[celli]);
            dSum += dD[celli];
            dMin = Foam::min(dMin, dD[celli]);
            dMax = Foam::max(dMax, dD[celli]);
            (void)celli;
        }
        reduce(aSum, sumOp<scalar>());
        reduce(aMin, minOp<scalar>());
        reduce(aMax, maxOp<scalar>());
        reduce(dSum, sumOp<scalar>());
        reduce(dMin, minOp<scalar>());
        reduce(dMax, maxOp<scalar>());
        (void)fvSum; (void)fvMin; (void)fvMax;
        Info<< "NS UEqn A(): avg=" << aSum/Foam::max(scalar(mesh.nCells()), 1)
            << " min=" << aMin << " max=" << aMax
            << " ; diag() avg=" << dSum/Foam::max(scalar(mesh.nCells()), 1)
            << " min=" << dMin << " max=" << dMax
            << " (fvOptions diag contribution skipped)" << endl;
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

    primalPressureMobility = rAtU().primitiveField();
    // =====================================================================
    // Rx probe (stage B4 diagnostic): design-induced flux derivative
    //   Rx-A: fixed U,p,phi, perturb alpha -> Brinkman residual derivative
    //   Rx-B: fixed U,p, rebuild phi(alpha) -> continuity dRp/dalpha
    // (runs while UEqn is still valid, i.e. before tUEqn.clear())
    // =====================================================================
    static bool rxProbeDone = false;
    if (optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false)
     && !rxProbeDone)
    {
        rxProbeDone = true;
        Info<< "RxProbe: running (U,p,phi fixed at current state)" << endl;

        volScalarField dAlphaDir
        (
            IOobject
            (
                "dAlphaDir",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("da", dimless, 0.0)
        );
        forAll(dAlphaDir.primitiveFieldRef(), celli)
        {
            dAlphaDir.primitiveFieldRef()[celli] =
                Foam::sin(0.371*(celli + 1));
        }
        dAlphaDir.primitiveFieldRef() *= designMask.primitiveField();

        const scalar alphaMaxVal = gMax(mag(alpha.primitiveField()));
        const scalar hAlpha = 1e-6*Foam::max(alphaMaxVal, 1e-12);
        const scalarField alphaP = alpha.primitiveField() + hAlpha*dAlphaDir.primitiveField();
        const scalarField alphaM = alpha.primitiveField() - hAlpha*dAlphaDir.primitiveField();
        const scalarField dA = dAlphaDir.primitiveField();

        // ---------- Rx-A: fixed phi, Brinkman residual derivative ----------
        scalar nRA = 0.0, nRD = 0.0, dotRA = 0.0;
        forAll(U, celli)
        {
            const scalar Vc = mesh.V()[celli];
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                const scalar u = U[celli].component(cmpt);
                const scalar fd = (alphaP[celli] - alphaM[celli])/(2.0*hAlpha)*u*Vc;
                const scalar an = dA[celli]*u*Vc;
                nRA += sqr(an);
                nRD += sqr(fd - an);
                dotRA += fd*an;
            }
        }
        reduce(nRA, sumOp<scalar>());
        reduce(nRD, sumOp<scalar>());
        reduce(dotRA, sumOp<scalar>());
        Info<< "RxProbe Rx-A (fixed phi, Brinkman dR_U/dalpha):"
            << " relL2=" << Foam::sqrt(nRD/Foam::max(nRA, SMALL))
            << " cos=" << dotRA/Foam::max(Foam::sqrt(nRA*nRA), SMALL)
            << " (expect relL2<=1e-6, cos=1)" << endl;

        // ---------- Rx-B: fixed U,p, rebuild phi(alpha) ----------
        volScalarField A0
        (
            IOobject
            (
                "A0",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("A0", UEqn.A().ref().dimensions(), 0.0)
        );
        // --- Sp-diagonal semantics diagnostic ---
        {
            scalar sA = 0.0, mnA = GREAT, mxA = -GREAT;
            scalar sAl = 0.0, mnAl = GREAT, mxAl = -GREAT;
            scalar sAlV = 0.0, mnAlV = GREAT, mxAlV = -GREAT;
            scalar sA0 = 0.0, mnA0 = GREAT, mxA0 = -GREAT;
            scalar sRau = 0.0, mnRau = GREAT, mxRau = -GREAT;
            forAll(UEqn.A().ref(), celli)
            {
                const scalar a = UEqn.A().ref()[celli];
                const scalar al = alpha[celli];
                const scalar alV = al*mesh.V()[celli];
                const scalar a0 = a - al;
                sA += a; mnA = Foam::min(mnA, a); mxA = Foam::max(mxA, a);
                sAl += al; mnAl = Foam::min(mnAl, al); mxAl = Foam::max(mxAl, al);
                sAlV += alV; mnAlV = Foam::min(mnAlV, alV); mxAlV = Foam::max(mxAlV, alV);
                sA0 += a0; mnA0 = Foam::min(mnA0, a0); mxA0 = Foam::max(mxA0, a0);
                const scalar rau = 1.0/Foam::max(a, SMALL);
                sRau += rau; mnRau = Foam::min(mnRau, rau); mxRau = Foam::max(mxRau, rau);
            }
            reduce(sA, sumOp<scalar>()); reduce(mnA, minOp<scalar>()); reduce(mxA, maxOp<scalar>());
            reduce(sAl, sumOp<scalar>()); reduce(mnAl, minOp<scalar>()); reduce(mxAl, maxOp<scalar>());
            reduce(sAlV, sumOp<scalar>()); reduce(mnAlV, minOp<scalar>()); reduce(mxAlV, maxOp<scalar>());
            reduce(sA0, sumOp<scalar>()); reduce(mnA0, minOp<scalar>()); reduce(mxA0, maxOp<scalar>());
            reduce(sRau, sumOp<scalar>()); reduce(mnRau, minOp<scalar>()); reduce(mxRau, maxOp<scalar>());
            const scalar nC = Foam::max(scalar(mesh.nCells()), 1.0);
            Info<< "RxProbe SpDiag: UEqn.A() avg=" << sA/nC
                << " [min=" << mnA << ",max=" << mxA << "]"
                << " ; alpha avg=" << sAl/nC << " [min=" << mnAl << ",max=" << mxAl << "]"
                << " ; alpha*V avg=" << sAlV/nC << " [min=" << mnAlV << ",max=" << mxAlV << "]"
                << " ; A-alpha avg=" << sA0/nC << " [min=" << mnA0 << ",max=" << mxA0 << "]"
                << " ; rAU=1/A avg=" << sRau/nC << " [min=" << mnRau << ",max=" << mxRau << "]"
                << endl;
        }
        // fvm::Sp(alpha, U) assembles alpha on the diagonal (no V factor
        // on the fvMatrix diagonal; V appears only in the residual/source
        // semantics).  So A(alpha) = A0 + alpha, and rAU = 1/(A0 + alpha).
        A0.primitiveFieldRef() =
            UEqn.A().ref().primitiveField() - alpha.primitiveField();
        volVectorField Hbase(UEqn.H().ref());
        volVectorField Ubase(U);
        volScalarField pbase(p);
        volScalarField meshV
        (
            IOobject
            (
                "meshV",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("V", dimVolume, 0.0)
        );
        meshV.primitiveFieldRef() = mesh.V();

        scalarField rpPlus(mesh.nCells(), 0.0);
        scalarField rpMinus(mesh.nCells(), 0.0);
        scalarField phiPlusFace(mesh.nInternalFaces(), 0.0);
        scalarField phiMinusFace(mesh.nInternalFaces(), 0.0);
        for (label sgn = 0; sgn < 2; ++sgn)
        {
            const scalarField alphaModF =
                (sgn == 0) ? alphaP : alphaM;
            volScalarField alphaMod
            (
                IOobject
                (
                    "alphaMod",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("am", alpha.dimensions(), 0.0)
            );
            alphaMod.primitiveFieldRef() = alphaModF;
            alphaMod.correctBoundaryConditions();

            volScalarField rAUmod
            (
                "rAUmod",
                1.0/(A0 + alphaMod)
            );
            rAUmod.correctBoundaryConditions();
            volVectorField HbyAmod
            (
                IOobject
                (
                    "HbyAmod",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rAUmod*Hbase
            );
            constrainHbyA(HbyAmod, Ubase, pbase);
            surfaceScalarField phiHbyAmod(fvc::flux(HbyAmod));
            adjustPhi(phiHbyAmod, Ubase, pbase);
            if (simple.consistent())
            {
                volScalarField rAtUmod
                (
                    "rAtUmod",
                    1.0/(1.0/rAUmod - UEqn.H1())
                );
                phiHbyAmod +=
                    fvc::interpolate(rAtUmod - rAUmod)
                   *fvc::snGrad(pbase)*mesh.magSf();
            }
            // complete SIMPLE phi reconstruction: phi = phiHbyA - pEqn.flux()
            fvScalarMatrix pEqnMod(fvm::laplacian(rAUmod, pbase));
            surfaceScalarField phiMod
            (
                IOobject
                (
                    "phiMod",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phiHbyAmod - pEqnMod.flux()
            );
            // manual continuity residual: internal faces from the rebuilt phi,
            // boundary faces from the converged base phi (U,p fixed semantics)
            scalarField rpMod(mesh.nCells(), 0.0);
            forAll(mesh.owner(), fi)
            {
                rpMod[mesh.owner()[fi]] += phiMod[fi];
                rpMod[mesh.neighbour()[fi]] -= phiMod[fi];
            }
            forAll(phi.boundaryField(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const scalarField& phiB = phi.boundaryField()[patchi];
                forAll(fc, f)
                {
                    rpMod[fc[f]] += phiB[f];
                }
            }
            if (sgn == 0) { rpPlus = rpMod; }
            else          { rpMinus = rpMod; }

            if (sgn == 0) { phiPlusFace = phiMod.primitiveField(); }
            else          { phiMinusFace = phiMod.primitiveField(); }
        }
        scalar nRP = 0.0, maxRP = 0.0;
        label nanP = 0, nanM = 0;
        forAll(rpPlus, celli)
        {
            if (!std::isfinite(rpPlus[celli])) { ++nanP; }
            if (!std::isfinite(rpMinus[celli])) { ++nanM; }
            const scalar fd = (rpPlus[celli] - rpMinus[celli])/(2.0*hAlpha);
            if (std::isfinite(fd))
            {
                nRP += sqr(fd);
                maxRP = Foam::max(maxRP, mag(fd));
            }
        }
        reduce(nRP, sumOp<scalar>());
        reduce(maxRP, maxOp<scalar>());
        reduce(nanP, sumOp<label>());
        reduce(nanM, sumOp<label>());

        // face-level dPhi/dalpha (internal faces)
        scalar nPhiF = 0.0, maxPhiF = 0.0;
        forAll(phiPlusFace, fi)
        {
            const scalar df = (phiPlusFace[fi] - phiMinusFace[fi])/(2.0*hAlpha);
            if (std::isfinite(df))
            {
                nPhiF += sqr(df);
                maxPhiF = Foam::max(maxPhiF, mag(df));
            }
        }
        reduce(nPhiF, sumOp<scalar>());
        reduce(maxPhiF, maxOp<scalar>());

        // boundary-face dPhi/dalpha: with U,p fixed and base phi on the
        // boundary (manual div above), the boundary contribution is zero by
        // construction.  A separate run with reconstructed boundary flux is
        // deferred to the analytic-vs-FD face oracle after implementation.
        Info<< "RxProbe Rx-B raw: rp nan=" << nanP + nanM
            << " |dRp/dalpha|L2=" << Foam::sqrt(nRP)
            << " max=" << maxRP << " hAlpha=" << hAlpha << endl;
        Info<< "RxProbe Rx-B faces: internal |dPhi/dalpha|L2="
            << Foam::sqrt(nPhiF) << " max=" << maxPhiF
            << " (" << mesh.nInternalFaces() << " faces, base-fixed boundary)"
            << endl;
        reduce(nRP, sumOp<scalar>());
        reduce(maxRP, maxOp<scalar>());
        Info<< "RxProbe Rx-B (fixed U,p, phi(alpha) rebuilt, dRp/dalpha):"
            << " |dRp/dalpha|L2=" << Foam::sqrt(nRP)
            << " max|dRp/dalpha|=" << maxRP
            << " (current analytic Rx_P = 0)" << endl;
        if (nRP > SMALL)
        {
            Info<< "RxProbe Rx-B: continuity dRp/dalpha is NONZERO (phi-alpha"
                << " contribution missing in current Rx)" << endl;
        }

        // =================================================================
        // Baseline + frozen-phi column audit using the FULL primal residual
        // (fvMatrix::residual() with RHS = -grad(p) + fvOptions).
        // No SIMPLE solve, no manual diag*U-H reconstruction.
        // =================================================================
        {
            // ---- baseline residual at the converged state ----
            {
                tmp<fvVectorMatrix> tEqBase
                (
                    fvm::div(phi, U)
                  - fvm::laplacian(nuEffFrozen, U)
                  + fvm::Sp(alpha, U)
                );
                fvVectorMatrix& EqBase = tEqBase.ref();
                fvOptions.constrain(EqBase);
                const tmp<Field<vector>> tResBase = EqBase.residual();
                volVectorField resBase
                (
                    IOobject
                    (
                        "resBase",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedVector("rB", dimAcceleration, vector::zero)
                );
                forAll(resBase, celli)
                {
                    resBase[celli] = -tResBase()[celli];
                }
                scalar nRU = 0.0, maxRU = 0.0, nRP = 0.0, maxRP = 0.0;
                forAll(U, celli)
                {
                    const scalar ru = mag(resBase[celli]);
                    nRU += sqr(ru);
                    maxRU = Foam::max(maxRU, ru);
                    nRP += sqr(phi[celli]);
                    maxRP = Foam::max(maxRP, mag(phi[celli]));
                }
                reduce(nRU, sumOp<scalar>());
                reduce(nRP, sumOp<scalar>());
                reduce(maxRU, maxOp<scalar>());
                reduce(maxRP, maxOp<scalar>());
                Info<< "RxProbe baseline: ||RU||L2=" << Foam::sqrt(nRU)
                    << " max=" << maxRU
                    << " ; ||phi||L2=" << Foam::sqrt(nRP)
                    << " max=" << maxRP
                    << " (RU should be ~0 at convergence)" << endl;
            }

            // ---- U-column FD (frozen phi, full residual) ----
            const scalar hU = 1e-3*Foam::max(gMax(mag(U.primitiveField())), SMALL);
            volVectorField dUdir
            (
                IOobject
                (
                    "dUdir",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("du", dimVelocity, vector::zero)
            );
            forAll(dUdir.primitiveFieldRef(), celli)
            {
                const scalar s = Foam::sin(0.271*(celli + 1));
                dUdir.primitiveFieldRef()[celli] = vector(s, 0.5*s, -0.3*s);
            }
            dUdir.correctBoundaryConditions();

            scalarField ruUP(3*mesh.nCells(), 0.0);
            scalarField ruUM(3*mesh.nCells(), 0.0);
            scalarField rpUP(mesh.nCells(), 0.0);
            scalarField rpUM(mesh.nCells(), 0.0);
            for (label sgn = 0; sgn < 2; ++sgn)
            {
                const scalar sg = (sgn == 0) ? 1.0 : -1.0;
                volVectorField Umod(U);
                forAll(Umod.primitiveFieldRef(), celli)
                {
                    Umod.primitiveFieldRef()[celli] += sg*hU*dUdir[celli];
                }
                Umod.correctBoundaryConditions();

                tmp<fvVectorMatrix> tEq
                (
                    fvm::div(phi, Umod)
                  - fvm::laplacian(nuEffFrozen, Umod)
                  + fvm::Sp(alpha, Umod)
                );
                fvVectorMatrix& Eq = tEq.ref();
                fvOptions.constrain(Eq);
                const tmp<Field<vector>> tRes = Eq.residual();
                volVectorField resU
                (
                    IOobject
                    (
                        "resU",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedVector("rU", dimAcceleration, vector::zero)
                );
                forAll(resU, celli)
                {
                    resU[celli] = -tRes()[celli];
                }

                scalarField rpU(mesh.nCells(), 0.0);
                forAll(mesh.owner(), fi)
                {
                    rpU[mesh.owner()[fi]] += phi[fi];
                    rpU[mesh.neighbour()[fi]] -= phi[fi];
                }
                forAll(phi.boundaryField(), patchi)
                {
                    const labelUList& fc = mesh.boundary()[patchi].faceCells();
                    const scalarField& phiB = phi.boundaryField()[patchi];
                    forAll(fc, f) { rpU[fc[f]] += phiB[f]; }
                }
                if (sgn == 0)
                {
                    forAll(U, celli)
                    {
                        rpUP[celli] = rpU[celli];
                        for (direction c = 0; c < 3; ++c)
                        {
                            ruUP[3*celli + c] = resU[celli].component(c);
                        }
                    }
                }
                else
                {
                    forAll(U, celli)
                    {
                        rpUM[celli] = rpU[celli];
                        for (direction c = 0; c < 3; ++c)
                        {
                            ruUM[3*celli + c] = resU[celli].component(c);
                        }
                    }
                }
            }
            OFstream duOut("stageB2_dUdirPhys.mtx");
            OFstream ruOut("stageB2_ruDudU.mtx");
            OFstream rpOut("stageB2_rpDudU.mtx");
            duOut.precision(16); ruOut.precision(16); rpOut.precision(16);
            forAll(U, celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    duOut<< dUdir[celli].component(c) << nl;
                    ruOut<< (ruUP[3*celli + c] - ruUM[3*celli + c])/(2.0*hU) << nl;
                }
                rpOut<< (rpUP[celli] - rpUM[celli])/(2.0*hU) << nl;
            }
            Info<< "RxProbe J-U exported (frozen-phi full residual, hU="
                << hU << ")" << endl;

            // ---- rebuilt-phi P-row FD (continuity): perturb U, rebuild phi
            // with the same path as NS.H (no pEqn.solve), then div(phi)
            scalarField rpRbP(mesh.nCells(), 0.0);
            scalarField rpRbM(mesh.nCells(), 0.0);
            scalarField phiFp(mesh.nInternalFaces(), 0.0);
            scalarField phiFm(mesh.nInternalFaces(), 0.0);
            scalarField APp(mesh.nCells(), 0.0);
            scalarField APm(mesh.nCells(), 0.0);
            vectorField HPp(mesh.nCells(), vector::zero);
            vectorField HPm(mesh.nCells(), vector::zero);
            // baseline: rebuilt phi at the converged state
            {
                volScalarField rAUrb0("rAUrb0", 1.0/UEqn.A().ref());
                rAUrb0.correctBoundaryConditions();
                volVectorField HbyArb0
                (
                    IOobject
                    (
                        "HbyArb0",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    rAUrb0*UEqn.H().ref()
                );
                constrainHbyA(HbyArb0, U, p);
                surfaceScalarField phiHbyArb0(fvc::flux(HbyArb0));
                adjustPhi(phiHbyArb0, U, p);
                fvScalarMatrix pEqnRb0(fvm::laplacian(rAUrb0, p));
                surfaceScalarField phiRb0(phiHbyArb0 - pEqnRb0.flux());
                scalar nR = 0.0, mxR = 0.0;
                forAll(mesh.owner(), fi)
                {
                    const scalar d = phiRb0[fi] - phi[fi];
                    nR += sqr(d);
                    mxR = Foam::max(mxR, mag(d));
                }
                reduce(nR, sumOp<scalar>());
                reduce(mxR, maxOp<scalar>());
                Info<< "RxProbe rebuilt baseline: |phiRb0-phi|L2="
                    << Foam::sqrt(nR) << " max=" << mxR
                    << " (should be small)" << endl;
            }
            for (label sgn = 0; sgn < 2; ++sgn)
            {
                const scalar sg = (sgn == 0) ? 1.0 : -1.0;
                volVectorField Umod(U);
                forAll(Umod.primitiveFieldRef(), celli)
                {
                    Umod.primitiveFieldRef()[celli] += sg*hU*dUdir[celli];
                }
                Umod.correctBoundaryConditions();
                // REASSEMBLE the momentum matrix with Umod so H depends on
                // the perturbed U (frozen H gives dphi/dU = 0).
                tmp<fvVectorMatrix> tUEqnRb
                (
                    fvm::div(phi, Umod)
                  - fvm::laplacian(nuEffFrozen, Umod)
                  + fvm::Sp(alpha, Umod)
                );
                volScalarField rAUrb("rAUrb", 1.0/tUEqnRb.ref().A());
                rAUrb.correctBoundaryConditions();
                volVectorField HbyArb
                (
                    IOobject
                    (
                        "HbyArb",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    rAUrb*tUEqnRb.ref().H()
                );
                constrainHbyA(HbyArb, Umod, p);
                surfaceScalarField phiHbyArb(fvc::flux(HbyArb));
                adjustPhi(phiHbyArb, Umod, p);
                fvScalarMatrix pEqnRb(fvm::laplacian(rAUrb, p));
                surfaceScalarField phiRb
                (
                    IOobject
                    (
                        "phiRb",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    phiHbyArb - pEqnRb.flux()
                );
                scalarField rpRb(mesh.nCells(), 0.0);
                forAll(mesh.owner(), fi)
                {
                    rpRb[mesh.owner()[fi]] += phiRb[fi];
                    rpRb[mesh.neighbour()[fi]] -= phiRb[fi];
                }
                forAll(phi.boundaryField(), patchi)
                {
                    const labelUList& fc = mesh.boundary()[patchi].faceCells();
                    const scalarField& phiB = phi.boundaryField()[patchi];
                    forAll(fc, f) { rpRb[fc[f]] += phiB[f]; }
                }
                // A/H FD snapshots (for deltaA/deltaH comparison)
                if (sgn == 0)
                {
                    rpRbP = rpRb;
                    forAll(phiRb, fi) { phiFp[fi] = phiRb[fi]; }
                    APp = tUEqnRb.ref().A().ref().primitiveField();
                    HPp = tUEqnRb.ref().H();
                }
                else
                {
                    rpRbM = rpRb;
                    forAll(phiRb, fi) { phiFm[fi] = phiRb[fi]; }
                    APm = tUEqnRb.ref().A().ref().primitiveField();
                    HPm = tUEqnRb.ref().H();
                }
            }
            auto hDx2 = [&](const vectorField& P, const vectorField& M, label ci, direction c, scalar h)
            {
                return (P[ci].component(c) - M[ci].component(c))/(2.0*h);
            };
            OFstream rpRbOut("stageB2_rpDudU_rebuilt.mtx");
            OFstream dPhiOut("stageB2_dphiDudU_face.mtx");
            OFstream dAOut("stageB2_dAdU_rebuilt.mtx");
            OFstream dHOut("stageB2_dHdU_rebuilt.mtx");
            rpRbOut.precision(16); dPhiOut.precision(16);
            dAOut.precision(16); dHOut.precision(16);
            forAll(U, celli)
            {
                rpRbOut<< (rpRbP[celli] - rpRbM[celli])/(2.0*hU) << nl;
                dAOut<< (APp[celli] - APm[celli])/(2.0*hU) << nl;
            }
            forAll(phiFp, fi)
            {
                dPhiOut<< (phiFp[fi] - phiFm[fi])/(2.0*hU) << nl;
            }
            // deltaH FD (per component) vs frozen-coefficient analytic
            scalar nHF = 0.0, nHA = 0.0, nHD = 0.0, dotH = 0.0;
            forAll(U, celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    const scalar fd = (HPp[celli].component(c) - HPm[celli].component(c))/(2.0*hU);
                    dHOut<< fd << nl;
                    nHF += sqr(fd);
                    nHD += sqr(fd);
                }
            }
            // frozen analytic: deltaH = -offdiag*dU - alpha*V*dU (diag source)
            scalarField dHx0(mesh.nCells(), 0.0);
            scalarField dHy0(mesh.nCells(), 0.0);
            scalarField dHz0(mesh.nCells(), 0.0);
            forAll(mesh.owner(), fi)
            {
                const label own = mesh.owner()[fi];
                const label nei = mesh.neighbour()[fi];
                const scalar up = UEqn.upper()[fi];
                const scalar lo = UEqn.lower()[fi];
                for (direction c = 0; c < 3; ++c)
                {
                    const scalar duN = hU*dUdir[nei].component(c);
                    const scalar duO = hU*dUdir[own].component(c);
                    if (c == 0) { dHx0[own] -= up*duN; dHx0[nei] -= lo*duO; }
                    if (c == 1) { dHy0[own] -= up*duN; dHy0[nei] -= lo*duO; }
                    if (c == 2) { dHz0[own] -= up*duN; dHz0[nei] -= lo*duO; }
                }
            }
            forAll(U, celli)
            {
                const scalar spT = -alpha[celli]*mesh.V()[celli];
                dHx0[celli] += spT*(hU*dUdir[celli].x());
                dHy0[celli] += spT*(hU*dUdir[celli].y());
                dHz0[celli] += spT*(hU*dUdir[celli].z());
            }
            scalar nA2 = 0.0, nD2 = 0.0, dot2 = 0.0;
            forAll(U, celli)
            {
                const vector fd(hDx2(HPp, HPm, celli, 0, hU), hDx2(HPp, HPm, celli, 1, hU), hDx2(HPp, HPm, celli, 2, hU));
                const vector an(dHx0[celli], dHy0[celli], dHz0[celli]);
                for (direction c = 0; c < 3; ++c)
                {
                    nA2 += sqr(an.component(c));
                    nD2 += sqr(fd.component(c) - an.component(c));
                    dot2 += fd.component(c)*an.component(c);
                }
            }
            reduce(nHF, sumOp<scalar>());
            reduce(nHD, sumOp<scalar>());
            reduce(nA2, sumOp<scalar>());
            reduce(nD2, sumOp<scalar>());
            reduce(dot2, sumOp<scalar>());
            Info<< "RxProbe deltaH: |dH_FD|L2=" << Foam::sqrt(nHF)
                << " ; frozen-analytic relL2="
                << Foam::sqrt(nD2/Foam::max(nA2, SMALL))
                << " cos=" << dot2/Foam::max(Foam::sqrt(nA2*nA2), SMALL)
                << " ; deltaA |dA_FD|L2="
                << Foam::sqrt(Foam::sum(sqr(APp - APm)))/(2.0*hU)
                << endl;
            // ---- decompose dH: div / laplacian / Sp contributions ----
            {
                vectorField HcP(mesh.nCells(), vector::zero);
                vectorField HcM(mesh.nCells(), vector::zero);
                vectorField HlP(mesh.nCells(), vector::zero);
                vectorField HlM(mesh.nCells(), vector::zero);
                vectorField HsP(mesh.nCells(), vector::zero);
                vectorField HsM(mesh.nCells(), vector::zero);
                for (label sgn2 = 0; sgn2 < 2; ++sgn2)
                {
                    const scalar sg2 = (sgn2 == 0) ? 1.0 : -1.0;
                    volVectorField Umod2(U);
                    forAll(Umod2.primitiveFieldRef(), celli)
                    {
                        Umod2.primitiveFieldRef()[celli] += sg2*hU*dUdir[celli];
                    }
                    Umod2.correctBoundaryConditions();
                    tmp<fvVectorMatrix> tC(fvm::div(phi, Umod2));
                    tmp<fvVectorMatrix> tL(fvm::laplacian(nuEffFrozen, Umod2));
                    tmp<fvVectorMatrix> tS(fvm::Sp(alpha, Umod2));
                    const vectorField Hc(tC.ref().H());
                    const vectorField Hl(tL.ref().H());
                    const vectorField Hs(tS.ref().H());
                    if (sgn2 == 0) { HcP = Hc; HlP = Hl; HsP = Hs; }
                    else          { HcM = Hc; HlM = Hl; HsM = Hs; }
                }
                auto l2of = [&](const vectorField& a, const vectorField& b)
                {
                    scalar s = 0.0;
                    forAll(a, i) { s += magSqr(a[i] - b[i]); }
                    return Foam::sqrt(s)/(2.0*hU);
                };
                Info<< "RxProbe dH decomposition: |dH_conv|="
                    << l2of(HcP, HcM)
                    << " |dH_laplacian|=" << l2of(HlP, HlM)
                    << " |dH_Sp|=" << l2of(HsP, HsM) << endl;

                // ---- H() 4-term split on tC (convection) + multi-h scan ----
                {
                    for (label hi = 0; hi < 4; ++hi)
                    {
                        const scalar hh = hU*Foam::pow(0.1, scalar(hi));
                        vectorField Hp(mesh.nCells(), vector::zero);
                        vectorField Hm(mesh.nCells(), vector::zero);
                        vectorField SpS(mesh.nCells(), vector::zero);
                        vectorField SmS(mesh.nCells(), vector::zero);
                        vectorField Ho(mesh.nCells(), vector::zero);
                        vectorField Hmo(mesh.nCells(), vector::zero);
                        vectorField Hbd(mesh.nCells(), vector::zero);
                        vectorField Hbm(mesh.nCells(), vector::zero);
                        for (label sgn3 = 0; sgn3 < 2; ++sgn3)
                        {
                            const scalar sg3 = (sgn3 == 0) ? 1.0 : -1.0;
                            volVectorField Umod3(U);
                            forAll(Umod3.primitiveFieldRef(), celli)
                            {
                                Umod3.primitiveFieldRef()[celli] += sg3*hh*dUdir[celli];
                            }
                            Umod3.correctBoundaryConditions();
                            tmp<fvVectorMatrix> tC3(fvm::div(phi, Umod3));
                            fvVectorMatrix& C3 = tC3.ref();
                            const vectorField src3(C3.source());
                            const vectorField H3(C3.H());
                            // offdiag-only contribution (H without source/boundary)
                            vectorField Hoo(mesh.nCells(), vector::zero);
                            forAll(C3.upper(), fi)
                            {
                                const label own = mesh.owner()[fi];
                                const label nei = mesh.neighbour()[fi];
                                for (direction c = 0; c < 3; ++c)
                                {
                                    Hoo[own].component(c) -= C3.upper()[fi]*Umod3[nei].component(c);
                                    Hoo[nei].component(c) -= C3.lower()[fi]*Umod3[own].component(c);
                                }
                            }
                            // boundary diagonal contribution
                            vectorField Hbdv(mesh.nCells(), vector::zero);
                            forAll(mesh.boundary(), patchi)
                            {
                                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                                if (U.boundaryField()[patchi].fixesValue()) continue;
                                const vectorField icV(C3.internalCoeffs()[patchi]);
                                forAll(fc, f)
                                {
                                    Hbdv[fc[f]].component(0) -= icV[f].component(0)*Umod3[fc[f]].component(0);
                                    Hbdv[fc[f]].component(1) -= icV[f].component(1)*Umod3[fc[f]].component(1);
                                    Hbdv[fc[f]].component(2) -= icV[f].component(2)*Umod3[fc[f]].component(2);
                                }
                            }
                            if (sgn3 == 0) { Hp = H3; SpS = src3; Ho = Hoo; Hbd = Hbdv; }
                            else          { Hm = H3; SmS = src3; Hmo = Hoo; Hbm = Hbdv; }
                        }
                        auto l2v = [&](const vectorField& a, const vectorField& b)
                        {
                            scalar s = 0.0;
                            forAll(a, i) { s += magSqr(a[i] - b[i]); }
                            return Foam::sqrt(s)/(2.0*hh);
                        };
                        auto l2s = [&](const vectorField& a, const vectorField& b)
                        {
                            scalar s = 0.0;
                            forAll(a, i) { s += magSqr(a[i] - b[i]); }
                            return Foam::sqrt(s)/(2.0*hh);
                        };
                        const scalar dHfd = l2v(Hp, Hm);
                        const scalar dHoff = l2v(Ho, Hmo);
                        const scalar dHbnd = l2v(Hbd, Hbm);
                        const scalar dHsrc = l2s(SpS, SmS);
                        {
                    // relax-source prediction: (Drelax-D0)*dU/V vs source FD
                    const scalarField& Drel = UEqn.diag();
                    tmp<fvVectorMatrix> tD0q
                    (
                        fvm::div(phi, U)
                      - fvm::laplacian(nuEffFrozen, U)
                      + fvm::Sp(alpha, U)
                    );
                    const scalarField& D0q = tD0q.ref().diag();
                    scalar nRel = 0.0, mxRel = 0.0;
                    forAll(U, celli)
                    {
                        const vector du(dUdir[celli]);
                        const scalar v = mesh.V()[celli];
                        const vector pred
                        (
                            (Drel[celli] - D0q[celli])*du.x()/v,
                            (Drel[celli] - D0q[celli])*du.y()/v,
                            (Drel[celli] - D0q[celli])*du.z()/v
                        );
                        nRel += magSqr(pred);
                        mxRel = Foam::max(mxRel, mag(pred));
                    }
                    Info<< "RxProbe relaxPred: |(Drel-D0)*dU/V|L2="
                        << Foam::sqrt(nRel) << " max=" << mxRel
                        << " (vs source-FD ~2.4e7 on adjoint basis)" << endl;
                }
                Info<< "RxProbe H4 h=" << hh
                            << " |dH_FD|=" << dHfd
                            << " |dH_offdiag|=" << dHoff
                            << " |dH_boundary|=" << dHbnd
                            << " |dH_source|=" << dHsrc
                            << " ratio(offdiag/FD)=" << dHoff/Foam::max(dHfd, SMALL)
                            << endl;
                    }
                }

                // ---- per-volume statistics + /V-correct analytic dH ----
                {
                    tmp<fvVectorMatrix> tConv0(fvm::div(phi, U));
                    const scalarField& up0 = tConv0.ref().upper();
                    const scalarField& lo0 = tConv0.ref().lower();
                    scalar mnV = GREAT, mxV = -GREAT, sV = 0.0;
                    scalar mnU = GREAT, mxU = -GREAT, sU = 0.0;
                    scalar mnL = GREAT, mxL = -GREAT, sL = 0.0;
                    scalar mnUV = GREAT, mxUV = -GREAT, sUV = 0.0;
                    scalar mnLV = GREAT, mxLV = -GREAT, sLV = 0.0;
                    forAll(up0, fi)
                    {
                        const label own = mesh.owner()[fi];
                        const label nei = mesh.neighbour()[fi];
                        const scalar Vu = mesh.V()[own];
                        const scalar Vn = mesh.V()[nei];
                        sU += up0[fi]; mnU = Foam::min(mnU, up0[fi]); mxU = Foam::max(mxU, up0[fi]);
                        sL += lo0[fi]; mnL = Foam::min(mnL, lo0[fi]); mxL = Foam::max(mxL, lo0[fi]);
                        sUV += up0[fi]/Vu; mnUV = Foam::min(mnUV, up0[fi]/Vu); mxUV = Foam::max(mxUV, up0[fi]/Vu);
                        sLV += lo0[fi]/Vn; mnLV = Foam::min(mnLV, lo0[fi]/Vn); mxLV = Foam::max(mxLV, lo0[fi]/Vn);
                    }
                    forAll(U, celli)
                    {
                        const scalar v = mesh.V()[celli];
                        sV += v; mnV = Foam::min(mnV, v); mxV = Foam::max(mxV, v);
                    }
                    const scalar nF = Foam::max(scalar(mesh.nInternalFaces()), 1.0);
                    const scalar nC = Foam::max(scalar(mesh.nCells()), 1.0);
                    Info<< "RxProbe stats: V[avg=" << sV/nC << " min=" << mnV
                        << " max=" << mxV << "] upper[avg=" << sU/nF
                        << " min=" << mnU << " max=" << mxU
                        << "] lower[avg=" << sL/nF << " min=" << mnL
                        << " max=" << mxL << "] upper/V[avg=" << sUV/nF
                        << " min=" << mnUV << " max=" << mxUV
                        << "] lower/V[avg=" << sLV/nF << " min=" << mnLV
                        << " max=" << mxLV << "]" << endl;

                    // analytic dH_conv with the fvMatrix::H() per-volume semantics
                    vectorField dHcA(mesh.nCells(), vector::zero);
                    forAll(up0, fi)
                    {
                        const label own = mesh.owner()[fi];
                        const label nei = mesh.neighbour()[fi];
                        for (direction c = 0; c < 3; ++c)
                        {
                            const scalar duO = hU*dUdir[own].component(c);
                            const scalar duN = hU*dUdir[nei].component(c);
                            dHcA[own].component(c) -= up0[fi]*duN/mesh.V()[own];
                            dHcA[nei].component(c) -= lo0[fi]*duO/mesh.V()[nei];
                        }
                    }
                    scalar nA3 = 0.0, nD3 = 0.0, dot3 = 0.0;
                    forAll(U, celli)
                    {
                        for (direction c = 0; c < 3; ++c)
                        {
                            const scalar fd = (HcP[celli].component(c) - HcM[celli].component(c))/(2.0*hU);
                            const scalar an = dHcA[celli].component(c);
                            nA3 += sqr(an);
                            nD3 += sqr(fd - an);
                            dot3 += fd*an;
                        }
                    }
                    Info<< "RxProbe dH_conv /V-analytic vs FD: relL2="
                        << Foam::sqrt(nD3/Foam::max(nA3, SMALL))
                        << " cos=" << dot3/Foam::max(Foam::sqrt(nA3*nA3), SMALL)
                        << " |analytic|=" << Foam::sqrt(nA3)
                        << " (H() divides by cell volume)" << endl;
                }
            }

            Info<< "RxProbe J-U rebuilt-phi P-row + face dphi/dU + dA/dH exported"
                << " (hU=" << hU << ")" << endl;
        }
    }
    {
        scalar rSum = 0.0, rMin = GREAT, rMax = 0.0;
        forAll(rAtU(), celli)
        {
            rSum += rAtU()[celli];
            rMin = Foam::min(rMin, rAtU()[celli]);
            rMax = Foam::max(rMax, rAtU()[celli]);
        }
        reduce(rSum, sumOp<scalar>());
        reduce(rMin, minOp<scalar>());
        reduce(rMax, maxOp<scalar>());
        Info<< "NS rAtU: avg=" << rSum/Foam::max(scalar(mesh.nCells()), 1)
            << " min=" << rMin << " max=" << rMax << endl;
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
        SolverPerformance<scalar> perfP = pEqn.solve();
        forwardInitialResidualP = perfP.initialResidual();
        forwardFinalResidualP = perfP.finalResidual();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    p.relax();
    U = HbyA - rAtU()*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    if (ransFlowModel)
    {
        if (!freezeTurbulenceForValidation)
        {
            turbulence->correct();
        }
        #include "updateFrozenTurbulenceFields.H"

        const volScalarField& kField =
            mesh.lookupObject<volScalarField>("k");
        const volScalarField& omegaField =
            mesh.lookupObject<volScalarField>("omega");

        auto coldRelativeScalarChange =
        [&](const scalarField& current, const scalarField& previous)
        {
            scalar numerator = 0.0;
            scalar denominator = 0.0;
            forAll(current, celli)
            {
                if (coldRegionMask[celli] > 0.5)
                {
                    numerator += mesh.V()[celli]
                       *sqr(current[celli] - previous[celli]);
                    denominator += mesh.V()[celli]*sqr(current[celli]);
                }
            }
            reduce(numerator, sumOp<scalar>());
            reduce(denominator, sumOp<scalar>());
            return Foam::sqrt
            (
                numerator/Foam::max(denominator, scalar(VSMALL))
            );
        };
        auto activeTurbulenceRelativeScalarChange =
        [&](const scalarField& current, const scalarField& previous)
        {
            scalar numerator = 0.0;
            scalar denominator = 0.0;
            forAll(current, celli)
            {
                const scalar weight =
                    mesh.V()[celli]*chiTurbulence[celli];
                numerator += weight*sqr(current[celli] - previous[celli]);
                denominator += weight*sqr(current[celli]);
            }
            reduce(numerator, sumOp<scalar>());
            reduce(denominator, sumOp<scalar>());
            return Foam::sqrt
            (
                numerator/Foam::max(denominator, scalar(VSMALL))
            );
        };

        scalar uNumerator = 0.0;
        scalar uDenominator = 0.0;
        forAll(U.primitiveField(), celli)
        {
            if (coldRegionMask[celli] > 0.5)
            {
                uNumerator += mesh.V()[celli]
                   *magSqr(U[celli] - previousU[celli]);
                uDenominator += mesh.V()[celli]*magSqr(U[celli]);
            }
        }
        reduce(uNumerator, sumOp<scalar>());
        reduce(uDenominator, sumOp<scalar>());

        const scalar uChange = Foam::sqrt
        (
            uNumerator/Foam::max(uDenominator, scalar(VSMALL))
        );
        const scalar pChange = coldRelativeScalarChange
        (
            p.primitiveField(),
            previousP
        );
        const scalar kChange = activeTurbulenceRelativeScalarChange
        (
            kField.primitiveField(),
            previousK
        );
        const scalar omegaChange = activeTurbulenceRelativeScalarChange
        (
            omegaField.primitiveField(),
            previousOmega
        );
        const scalar nutChange = activeTurbulenceRelativeScalarChange
        (
            nutFrozen.primitiveField(),
            previousNut
        );
        // The frozen primal/adjoint operators depend on nutFrozen, not on
        // omega directly. Keep omega as a diagnostic while gating the state
        // and the effective frozen transport coefficient that are consumed.
        const scalar maximumRansChange = Foam::max
        (
            Foam::max(uChange, pChange),
            Foam::max(kChange, nutChange)
        );

        // Phase K: expose for CSV logging
        lastRelativeChangeK = kChange;
        lastRelativeChangeOmega = omegaChange;

        if
        (
            i + 1 >= turbulenceMinCorrectors
         && maximumRansChange < turbulenceConvergenceTolerance
        )
        {
            ++turbulenceConvergedCount;
        }
        else
        {
            turbulenceConvergedCount = 0;
        }

        if ((i + 1) % 10 == 0 || turbulenceConvergedCount > 0)
        {
            Info<< "Frozen-RANS relative changes: U=" << uChange
                << ", p=" << pChange
                << ", k=" << kChange
                << ", omega=" << omegaChange
                << ", nut=" << nutChange
                << ", consecutive=" << turbulenceConvergedCount << endl;
        }

        previousU = U.primitiveField();
        previousP = p.primitiveField();
        previousK = kField.primitiveField();
        previousOmega = omegaField.primitiveField();
        previousNut = nutFrozen.primitiveField();

        if (turbulenceConvergedCount >= 5)
        {
            turbulenceConverged = true;
            Info<< "Frozen-RANS flow converged after " << i + 1
                << " correctors." << endl;
            U.storePrevIter();
            p.storePrevIter();
            phi.storePrevIter();
            break;
        }
    }

    U.storePrevIter();
    p.storePrevIter();
    phi.storePrevIter();
}

if (ransFlowModel && !turbulenceConverged)
{
    if (solveFlowAdjoints || mmaUpdateEnabled)
    {
        FatalErrorInFunction
            << "Frozen-RANS flow did not satisfy five consecutive relative-"
            << "change checks below " << turbulenceConvergenceTolerance
            << " after " << turbulenceMaxCorrectors << " correctors. "
            << "Adjoint/design update is blocked." << abort(FatalError);
    }
    else
    {
        WarningInFunction
            << "Forward-only RANS warm-up reached "
            << turbulenceMaxCorrectors << " correctors without satisfying "
            << "the final convergence gate. No adjoint or design update will "
            << "be evaluated; continue the same design at the next outer "
            << "iteration." << endl;
    }
}
