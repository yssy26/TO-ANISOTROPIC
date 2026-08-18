// =============================================================================
// stageB6RxDesignOracle.H  (BFINAL-004, DIAGNOSTIC_ONLY)
// -----------------------------------------------------------------------------
#include <sstream>
// -----------------------------------------------------------------------------
// Guarded diagnostic probe (optProperties switch "stageB6RxDesignOracle",
// default false).  It must be #included ONLY inside src/solveDiscreteFlowAdjoint.H
// AFTER the stageB5 block, and runs only in the pressureDrop adjoint block
// (discreteAdjointLabel == "pressureDrop"), i.e. in the AdjNS_PD.H inclusion
// where Uc/pc are the SOLVED exact pressure-drop adjoint fields.
//
// FIXED-STATE R_x CLOSURE ORACLE (task BFINAL-004 §5-§10).
//
// Residual semantics (LOCKED by BFINAL-003, identical to
// stageB5BoundaryRelaxOracle.H rebuildProduction):
//   R_U = -UEqn.residual()          (relaxation-invariant; full production UEqn
//                                    incl. dev2 when ransFlowModel)
//   R_P = div(phi)                  (raw face-flux sum, NO 1/V; boundary flux
//                                    REBUILT) with phi rebuilt through the exact
//                                    NS.H SIMPLE path:
//     UEqn(phi,U) - laplacian(nuEffFrozen,U) + Sp(alpha,U) == -grad(p)
//        -> relax -> fvOptions.constrain -> residual
//     rAU = 1/A ; HbyA = constrainHbyA(rAU*H, U, p)
//     phiHbyA = fvc::flux(HbyA) ; adjustPhi ; rAtU (consistent=false -> rAU)
//     constrainPressure ; pEqn = laplacian(rAtU,p) == div(phiHbyA)
//     setReference(pRefCell) ; phi = phiHbyA - pEqn.flux()
//
// Gates:
//   RX-A (task §5/§6): fixed U/p/phi/k/omega/nutFrozen/nuEffFrozen, deterministic
//     deltaAlpha (active + unclipped cells only), multi-eps central FD
//     [R(w,alpha+)-R(w,alpha-)]/(2eps) vs a SEPARATELY implemented analytic
//     linearization (no shared formula with the FD side):
//       dR_U/dalpha  = +V*U                       (per cell; -residual() sign)
//       dR_P/dalpha  = div(dphiHbyA_w - dflux_w), the DIRECTIONAL derivative
//         J_P*w with the per-cell weight w (= deltaAlpha for RX-A,
//         dAlphaDxh*deltaXh for RX-B, dAlphaDxh*z for RX-C) applied INSIDE the
//         face-flux assembly, BEFORE the divergence:
//         drAUw   = drAU*w ; dHbyAw = dHbyA*w, with the (unweighted) per-cell
//         drAU    = -rAU^2/alphaRel               (rAU = V/D_rel already carries V)
//         dHbyA   = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)   (V cancels)
//         dphiHbyA_w_f = Sf&(wf*dHbyAw_o+(1-wf)*dHbyAw_n) ; boundary: 0 on
//                      non-assignable U patches (constrainHbyA), extrapolated
//                      on assignable (zeroGradient outlet)
//         dflux_w = flux(fvm::laplacian(drAUw, p)) (drAtU=drAU, consistent=false)
//       R_P = div(phi) is NON-diagonal in alpha (face fluxes interpolate HbyA
//       from own+nei and the laplacian(rAU,p) flux couples neighbours), so the
//       weight must enter BEFORE the divergence: J_P*w, NEVER (J_P*1)*w
//       pointwise (the latter equals the response to a spatially UNIFORM +1
//       perturbation, which differs from the directional derivative when w is
//       space-varying -- cycle-1 RX-A R_P,alpha relL2=1.0115 gap).  R_U is
//       DIAGONAL in alpha (fvm::Sp touches only the own diagonal), so its
//       per-cell product anRUa*w is exact.
//     Per-term norms (div(dphiHbyA_w), -div(dflux_w)) are logged and exported
//     (weighted columns; col1+col2 == col3 == FD oracle).
//   RX-B (task §7): deterministic deltaXh (active, unclipped, away from xh
//     bounds), alpha recomputed via the EXACT production interpolation
//     alpha = max(alphamin, alphaMax*qu*(1-xh)/(qu+xh+SMALL)); compare
//     [R(xh+)-R(xh-)]/(2eps) vs R_alpha*(dAlphaDxh*deltaXh) (RX-A operator).
//   RX-C (task §8): D1/D2/D3 raw-design directions
//     (exact validateStageB2GradientAmplitude.H construction); x +/- eps*d
//     -> Helmholtz filter (exact production solve) -> adaptive-Heaviside
//     projection (bisected eta) -> alpha -> residual; FD vs analytic chain
//     R_alpha*(dAlphaDxh*z) with z = dxh/dx*d via the exact filter tangent
//     solve (M-B)y = -b d and the exact projection Jacobian
//     z = drho*y + dProjectionDeta*[sum(designMask*(1-drho)*V*y)
//         /sum(designMask*dProjectionDeta*V)].
//   §9/§10: ||R_P,x d||/||R_U,x d||; D_momentum = -Uc^T R_U,x d;
//     D_pressure = -pc^T R_P,x d; production projected gsensPressureDrop*d
//     (full production replica: dAlphaDxh -> designMask -> drho + eta
//     correction -> adjoint filter -> designMask).
//   Regression anchor: volume-gradient chain projV_D1 = 0.1555186511144686
//     (production replica of gsenshVol -> drho + eta corr -> adjoint filter).
//
// NEVER modifies U/p/phi/alpha/x/xp/xh/k/omega/nutFrozen/nuEffFrozen; NEVER
// solves UEqn/pEqn (only assembles and reads residual()/A()/H()/flux());
// all perturbed fields are local copies; filter solves act on diagnostic
// local fields only.  No-op when the switch is absent/false.
// =============================================================================
{
    Info<< "--- StageB6: fixed-state R_x closure oracle (BFINAL-004, "
        << "DIAGNOSTIC_ONLY, label=" << discreteAdjointLabel << ") ---" << endl;

    const scalar stageB6Eps[5] = {1e-3, 3e-4, 1e-4, 3e-5, 1e-5};
    const label nIntF = mesh.nInternalFaces();
    const label nBndF = mesh.nFaces() - nIntF;
    const label nFaceTot = nIntF + nBndF;
    const label nC = discreteNCells;

    // ------------------------------------------------------------------
    // 0. State fingerprint (prove no production state is modified).
    // ------------------------------------------------------------------
    auto stageB6Fingerprint = [&]()
    {
        scalar sU = 0.0, sP = 0.0, sPhi = 0.0, sA = 0.0;
        scalar sK = 0.0, sO = 0.0, sNut = 0.0, sNuE = 0.0, sX = 0.0;
        scalar sXh = 0.0, sXp = 0.0;
        forAll(U, celli)
        {
            sU += magSqr(U[celli]);
            sP += sqr(p[celli]);
            sA += sqr(alpha[celli]);
            sX += sqr(x[celli]);
            sXh += sqr(xh[celli]);
            sXp += sqr(xp[celli]);
        }
        forAll(phi, fi) { sPhi += sqr(phi[fi]); }
        const volScalarField& kF = mesh.lookupObject<volScalarField>("k");
        const volScalarField& omegaF = mesh.lookupObject<volScalarField>("omega");
        forAll(kF, celli)
        {
            sK += sqr(kF[celli]);
            sO += sqr(omegaF[celli]);
            sNut += sqr(nutFrozen[celli]);
            sNuE += sqr(nuEffFrozen[celli]);
        }
        reduce(sU, sumOp<scalar>()); reduce(sP, sumOp<scalar>());
        reduce(sPhi, sumOp<scalar>()); reduce(sA, sumOp<scalar>());
        reduce(sK, sumOp<scalar>()); reduce(sO, sumOp<scalar>());
        reduce(sNut, sumOp<scalar>()); reduce(sNuE, sumOp<scalar>());
        reduce(sX, sumOp<scalar>()); reduce(sXh, sumOp<scalar>()); reduce(sXp, sumOp<scalar>());
        Info<< "StageB6 fingerprint: |U|2=" << Foam::sqrt(sU)
            << " |p|2=" << Foam::sqrt(sP)
            << " |phi|2=" << Foam::sqrt(sPhi)
            << " |alpha|2=" << Foam::sqrt(sA)
            << " |k|2=" << Foam::sqrt(sK)
            << " |omega|2=" << Foam::sqrt(sO)
            << " |nutFrozen|2=" << Foam::sqrt(sNut)
            << " |nuEffFrozen|2=" << Foam::sqrt(sNuE)
            << " |x|2=" << Foam::sqrt(sX)
            << " |xp|2=" << Foam::sqrt(sXp)
            << " |xh|2=" << Foam::sqrt(sXh)
            << " eta5=" << eta5 << " del=" << del << endl;
    };
    stageB6Fingerprint();

    // ------------------------------------------------------------------
    // 1. Support classification + production dAlphaDxh (exact sensitivity.H
    //    formula incl. mask/alphamin-clip zeroing).
    // ------------------------------------------------------------------
    boolList activeCell(nC, false);
    boolList unclippedCell(nC, false);
    label nActive = 0, nUnclipped = 0;
    forAll(designMask, celli)
    {
        activeCell[celli] = designMask[celli] > 0.5;
        unclippedCell[celli] =
            activeCell[celli]
         && alpha[celli] > alphamin.value()*(1.0 + 1e-10);
        if (activeCell[celli]) ++nActive;
        if (unclippedCell[celli]) ++nUnclipped;
    }

    scalarField dAlphaDxh(nC, 0.0);
    forAll(xh, celli)
    {
        if (unclippedCell[celli])
        {
            dAlphaDxh[celli] =
               -alphaMax.value()*(1.0 + qu)*qu
               /Foam::sqr(qu + xh[celli] + SMALL);
        }
    }

    // ------------------------------------------------------------------
    // 2. Deterministic perturbation patterns.
    // ------------------------------------------------------------------
    scalarField deltaAlpha(nC, 0.0);
    scalarField deltaXh(nC, 0.0);
    label nDa = 0, nDx = 0;
    forAll(xh, celli)
    {
        const scalar s = Foam::sin(0.271*(celli + 1));
        if (unclippedCell[celli])
        {
            deltaAlpha[celli] = s;
            ++nDa;
        }
        if
        (
            unclippedCell[celli]
         && xh[celli] > 1e-3
         && xh[celli] < 1.0 - 1e-3
        )
        {
            deltaXh[celli] = s;
            ++nDx;
        }
    }
    Info<< "StageB6 support: active=" << nActive << " unclipped=" << nUnclipped
        << " deltaAlpha-support=" << nDa << " deltaXh-support=" << nDx << endl;

    // ------------------------------------------------------------------
    // 3. D1/D2/D3 (exact validateStageB2GradientAmplitude.H construction).
    // ------------------------------------------------------------------
    scalarField D1(nC, 0.0), D2(nC, 0.0), D3(nC, 0.0);
    {
        const vectorField& C = mesh.C().primitiveField();
        const scalar Lx = gMax(C.component(0)) - gMin(C.component(0));
        const scalar Ly = gMax(C.component(1)) - gMin(C.component(1));
        const scalar Lz = gMax(C.component(2)) - gMin(C.component(2));
        forAll(C, celli)
        {
            if (!activeCell[celli]) continue;
            const scalar xx = C[celli].x()/Foam::max(Lx, SMALL);
            const scalar yy = C[celli].y()/Foam::max(Ly, SMALL);
            const scalar zz = C[celli].z()/Foam::max(Lz, SMALL);
            D1[celli] =
                Foam::sin(2.0*Foam::constant::mathematical::pi*xx)
               *Foam::cos(Foam::constant::mathematical::pi*yy)
               + 0.5*Foam::sin(3.0*Foam::constant::mathematical::pi*zz);
            D2[celli] =
                Foam::sin(4.0*Foam::constant::mathematical::pi*xx)
               + 0.3*Foam::cos(2.0*Foam::constant::mathematical::pi*xx);
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
            forAll(d, celli) { mx = Foam::max(mx, mag(d[celli])); }
            reduce(mx, maxOp<scalar>());
            if (mx > 1e-30) d /= mx;
        };
        normalizeDirection(D1);
        normalizeDirection(D2);
        normalizeDirection(D3);
    }
    Info<< "StageB6 directions: |D1|inf=" << gMax(mag(D1))
        << " |D2|inf=" << gMax(mag(D2))
        << " |D3|inf=" << gMax(mag(D3)) << endl;

    // ------------------------------------------------------------------
    // 4. rebuildResidual(alphaMod): THE LOCKED PRODUCTION REDUCED RESIDUAL
    //    (FD side).  Writes rebuilt phi, R_U = -residual, R_P = div(phi),
    //    relaxed rAU, relaxed H(), HbyA cell values.
    // ------------------------------------------------------------------
    auto rebuildResidual =
    [&](const volScalarField& alphaMod,
        surfaceScalarField& phiOut,
        vectorField& resUOut,
        scalarField& rpOut,
        scalarField& rAUrelOut,
        vectorField& HrelOut,
        vectorField& HbyAOut)
    {
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, U)
          - fvm::laplacian(nuEffFrozen, U)
          + fvm::Sp(alphaMod, U)
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

        const tmp<Field<vector>> tRes = UEqn.residual();
        forAll(resUOut, celli)
        {
            resUOut[celli] = -tRes()[celli];
        }

        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
        surfaceScalarField phiHbyA("stageB6_phiHbyA", fvc::flux(HbyA));
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

        constrainPressure(p, U, phiHbyA, rAtU());

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );
        pEqn.setReference(pRefCell, pRefValue);

        phiOut = phiHbyA - pEqn.flux();

        forAll(rAUrelOut, celli) { rAUrelOut[celli] = rAU[celli]; }
        const vectorField Hrel(UEqn.H());
        forAll(HrelOut, celli) { HrelOut[celli] = Hrel[celli]; }
        forAll(HbyAOut, celli) { HbyAOut[celli] = HbyA[celli]; }

        forAll(rpOut, celli) { rpOut[celli] = 0.0; }
        forAll(mesh.owner(), fi)
        {
            rpOut[mesh.owner()[fi]] += phiOut[fi];
            rpOut[mesh.neighbour()[fi]] -= phiOut[fi];
        }
        forAll(phiOut.boundaryField(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const scalarField& pB = phiOut.boundaryField()[patchi];
            forAll(fc, f) { rpOut[fc[f]] += pB[f]; }
        }
    };

    // ------------------------------------------------------------------
    // 5. Metric helper (per-block L2 / relL2 / cosine / max-abs-diff).
    //    relL2 = |a-b|/|b| (b = FD reference), identical to stageB5.
    // ------------------------------------------------------------------
    auto stageB6Metrics =
    [&](const scalarField& a, const scalarField& b, const label nUse,
        const word& labelName)
    {
        scalar nA = 0.0, nB = 0.0, nD = 0.0, dot = 0.0;
        scalar mx = 0.0;
        label arg = -1;
        for (label k = 0; k < nUse; ++k)
        {
            const scalar av = a[k];
            const scalar bv = b[k];
            const scalar d = av - bv;
            nA += av*av; nB += bv*bv; nD += d*d; dot += av*bv;
            if (mag(d) > mx) { mx = mag(d); arg = k; }
        }
        reduce(nA, sumOp<scalar>());
        reduce(nB, sumOp<scalar>());
        reduce(nD, sumOp<scalar>());
        reduce(dot, sumOp<scalar>());
        reduce(mx, maxOp<scalar>());
        reduce(arg, maxOp<label>());
        const scalar nA2 = Foam::sqrt(nA);
        const scalar nB2 = Foam::sqrt(nB);
        const scalar nD2 = Foam::sqrt(nD);
        // Display guard: degenerate only for EXACT zero.  NEVER clamp nB or
        // nA*nB by SMALL: when the reference norm is below SMALL (as in the
        // cycle-1 RX-A R_P,alpha comparison, nB=3.87e-19 < SMALL), a
        // SMALL-clamped denominator printed a geometrically impossible
        // relL2/cos (0.0422 / -1.75e-5 while |a-b|/|b| = 1.0115).  The mtx
        // values are authoritative; the printed relL2/cos must be consistent
        // with the printed |a|,|b|,|err|.
        if (nB2 <= 0.0 || nA2*nB2 <= 0.0)
        {
            Info<< "StageB6 " << labelName << ": |a|L2=" << nA2
                << " |b(FD)|L2=" << nB2
                << " |err|L2=" << nD2
                << " relL2=degenerate(nB=0) cos=degenerate"
                << " maxDiff=" << mx << " @idx=" << arg << endl;
            return 0.0;
        }
        const scalar rel = nD2/nB2;
        const scalar cosv = dot/(nA2*nB2);
        Info<< "StageB6 " << labelName << ": |a|L2=" << nA2
            << " |b(FD)|L2=" << nB2
            << " |err|L2=" << nD2
            << " relL2=" << rel
            << " cos=" << cosv
            << " maxDiff=" << mx << " @idx=" << arg << endl;
        return rel;
    };

    // ------------------------------------------------------------------
    // 6. RX-A: fixed-state alpha perturbation (task §5/§6).
    // ------------------------------------------------------------------
    Info<< "StageB6 RX-A: fixed U/p/phi/turbulence, deterministic deltaAlpha"
        << ", multi-eps central FD vs analytic linearization" << endl;

    // base rebuild at the current alpha (rAU/H/HbyA for the analytic chain)
    surfaceScalarField phiBase("stageB6_phiBase", phi);
    vectorField resUBase(nC, vector::zero);
    scalarField rpBase(nC, 0.0);
    scalarField rAUrelBase(nC, 0.0);
    vectorField HrelBase(nC, vector::zero);
    vectorField HbyABase(nC, vector::zero);
    rebuildResidual(alpha, phiBase, resUBase, rpBase, rAUrelBase, HrelBase, HbyABase);

    // analytic R_U,alpha = +V*U  (R_U = -residual() convention)
    scalarField anRUa(3*nC, 0.0);
    forAll(U, celli)
    {
        const scalar Vc = mesh.V()[celli];
        anRUa[3*celli + 0] = Vc*U[celli].x();
        anRUa[3*celli + 1] = Vc*U[celli].y();
        anRUa[3*celli + 2] = Vc*U[celli].z();
    }

    // analytic R_P,alpha: div(dphiHbyA - dflux) with
    //   drAU = -rAU^2/alphaRel ; dHbyA = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)
    scalarField drAU(nC, 0.0);
    vectorField dHbyA(nC, vector::zero);
    forAll(U, celli)
    {
        const scalar rAUc = rAUrelBase[celli];
        // fvm::Sp(alpha,U) adds alpha*V to the INTEGRATED diag D0 (fvmSup.C:118);
        // relax() -> D_rel = D0/alphaRel (fvMatrix.C:521); A() = D_rel/V
        // (fvMatrix.C:751); so rAU = 1/A() = V/D_rel ALREADY carries V.  Hence
        // drAU/dalpha = -V*(V/alphaRel)/D_rel^2 = -rAU^2/alphaRel (NO extra V).
        // HbyA = H0/D_rel + (1-alphaRel)*U (V cancels through rAU*H() with
        // H() = H0/V) => dHbyA/dalpha = (rAU/alphaRel)*((1-alphaRel)*U - HbyA).
        drAU[celli] = -rAUc*rAUc/alphaRel;
        dHbyA[celli] =
            (rAUc/alphaRel)
           *(
                (1.0 - alphaRel)*U[celli]
              - HbyABase[celli]
            );
    }
    // dphiHbyA per face (internal then boundary)
    scalarField dphiHbyA(nFaceTot, 0.0);
    forAll(mesh.owner(), fi)
    {
        const label own = mesh.owner()[fi];
        const label nei = mesh.neighbour()[fi];
        const scalar wf = mesh.weights()[fi];
        dphiHbyA[fi] =
            mesh.Sf()[fi]
          & (wf*dHbyA[own] + (1.0 - wf)*dHbyA[nei]);
    }
    {
        label bf = nIntF;
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
            const bool uAssignable = U.boundaryField()[patchi].assignable();
            forAll(fc, f)
            {
                const label celli = fc[f];
                const vector dHb = uAssignable ? dHbyA[celli] : vector::zero;
                dphiHbyA[bf] = pSf[f] & dHb;
                ++bf;
            }
        }
    }
    // dflux = flux(fvm::laplacian(drAU, p))  (drAtU = drAU, consistent=false)
    volScalarField drAUfield
    (
        IOobject
        (
            "stageB6_drAU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("d", dimTime, 0.0)
    );
    forAll(drAUfield.primitiveFieldRef(), celli)
    {
        drAUfield.primitiveFieldRef()[celli] = drAU[celli];
    }
    drAUfield.correctBoundaryConditions();
    fvScalarMatrix dpEqn(fvm::laplacian(drAUfield, p));
    surfaceScalarField dflux(dpEqn.flux());
    // Materialize dflux over ALL faces (surfaceScalarField operator[] covers
    // only internal faces; boundary values live in the boundary field).
    scalarField dfluxFaces(nFaceTot, 0.0);
    for (label fi = 0; fi < nIntF; ++fi) { dfluxFaces[fi] = dflux[fi]; }
    {
        label bf = nIntF;
        forAll(dflux.boundaryField(), patchi)
        {
            const scalarField& dfB = dflux.boundaryField()[patchi];
            forAll(dfB, f) { dfluxFaces[bf] = dfB[f]; ++bf; }
        }
    }
    // dphi per face and R_P,alpha per cell (with per-term decomposition)
    scalarField dphiFace(nFaceTot, 0.0);
    scalarField anRPa(nC, 0.0);
    scalarField anRPaHbyA(nC, 0.0);   // div(dphiHbyA)
    scalarField anRPaFlux(nC, 0.0);   // -div(dflux)
    forAll(dphiFace, fi)
    {
        dphiFace[fi] = dphiHbyA[fi] - dfluxFaces[fi];
    }
    forAll(anRPa, celli) { anRPa[celli] = 0.0; anRPaHbyA[celli] = 0.0; anRPaFlux[celli] = 0.0; }
    forAll(mesh.owner(), fi)
    {
        const label own = mesh.owner()[fi];
        const label nei = mesh.neighbour()[fi];
        anRPa[own] += dphiFace[fi];
        anRPa[nei] -= dphiFace[fi];
        anRPaHbyA[own] += dphiHbyA[fi];
        anRPaHbyA[nei] -= dphiHbyA[fi];
        anRPaFlux[own] -= dfluxFaces[fi];
        anRPaFlux[nei] += dfluxFaces[fi];
    }
    {
        label bf = nIntF;
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(fc, f)
            {
                const label celli = fc[f];
                anRPa[celli] += dphiFace[bf];
                anRPaHbyA[celli] += dphiHbyA[bf];
                anRPaFlux[celli] -= dfluxFaces[bf];
                ++bf;
            }
        }
    }
    // per-term L2 norms for the log (UNWEIGHTED row-sum J_P*1 decomposition;
    // the weighted J_P*deltaAlpha decomposition is logged in the RX-A block)
    {
        scalar nH = 0.0, nF = 0.0, nT = 0.0;
        forAll(anRPa, celli)
        {
            nH += sqr(anRPaHbyA[celli]);
            nF += sqr(anRPaFlux[celli]);
            nT += sqr(anRPa[celli]);
        }
        reduce(nH, sumOp<scalar>());
        reduce(nF, sumOp<scalar>());
        reduce(nT, sumOp<scalar>());
        Info<< "StageB6 RX-A analytic R_P,alpha terms (UNWEIGHTED row-sum J_P*1):"
            << " |div(dphiHbyA)|L2=" << Foam::sqrt(nH)
            << " |-div(dflux)|L2=" << Foam::sqrt(nF)
            << " |total|L2=" << Foam::sqrt(nT)
            << " ; alphaRel=" << alphaRel << endl;
    }

    // ------------------------------------------------------------------
    // 6b. assembleWeightedRPa(w): J_P*w -- the DIRECTIONAL P-row derivative.
    //   R_P = div(phi) is NON-diagonal in alpha: phiHbyA faces interpolate
    //   HbyA from own+nei, and the laplacian(rAU,p) flux couples neighbours,
    //   so the directional weight w must be multiplied into the per-cell
    //   derivatives BEFORE the face-flux assembly and the divergence:
    //     dphiHbyA_w[face] = Sf&(wf*(dHbyA*w)[own]+(1-wf)*(dHbyA*w)[nei])
    //     dflux_w          = flux(fvm::laplacian(drAU*w, p))
    //   (both linear in the per-cell derivatives; drAU/dHbyA are the verified
    //   unweighted per-cell values from above).  Returns
    //     anRPwOut      = div(dphiHbyA_w - dflux_w) = J_P*w
    //     anRPaHbyAwOut = div(dphiHbyA_w)           (per-term)
    //     anRPaFluxwOut = -div(dflux_w)             (per-term)
    //   with anRPaHbyAwOut + anRPaFluxwOut == anRPwOut per cell.
    // ------------------------------------------------------------------
    auto assembleWeightedRPa =
    [&](const scalarField& w,
        scalarField& anRPwOut,
        scalarField& anRPaHbyAwOut,
        scalarField& anRPaFluxwOut)
    {
        scalarField drAUw(nC, 0.0);
        vectorField dHbyAw(nC, vector::zero);
        forAll(U, celli)
        {
            drAUw[celli] = drAU[celli]*w[celli];
            dHbyAw[celli] = dHbyA[celli]*w[celli];
        }
        // dphiHbyA_w per face (internal then boundary)
        scalarField dphiHbyAw(nFaceTot, 0.0);
        forAll(mesh.owner(), fi)
        {
            const label own = mesh.owner()[fi];
            const label nei = mesh.neighbour()[fi];
            const scalar wf = mesh.weights()[fi];
            dphiHbyAw[fi] =
                mesh.Sf()[fi]
              & (wf*dHbyAw[own] + (1.0 - wf)*dHbyAw[nei]);
        }
        {
            label bf = nIntF;
            forAll(mesh.boundary(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
                const bool uAssignable = U.boundaryField()[patchi].assignable();
                forAll(fc, f)
                {
                    const label celli = fc[f];
                    const vector dHb = uAssignable ? dHbyAw[celli] : vector::zero;
                    dphiHbyAw[bf] = pSf[f] & dHb;
                    ++bf;
                }
            }
        }
        // dflux_w = flux(fvm::laplacian(drAU*w, p)) (drAtU = drAU, consistent=false)
        volScalarField drAUwField
        (
            IOobject
            (
                "stageB6_drAUw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("d", dimTime, 0.0)
        );
        forAll(drAUwField.primitiveFieldRef(), celli)
        {
            drAUwField.primitiveFieldRef()[celli] = drAUw[celli];
        }
        drAUwField.correctBoundaryConditions();
        fvScalarMatrix dpEqnW(fvm::laplacian(drAUwField, p));
        surfaceScalarField dfluxW(dpEqnW.flux());
        // Materialize dfluxW over ALL faces (surfaceScalarField operator[]
        // covers only internal faces; boundary values live in the boundary field).
        scalarField dfluxWFaces(nFaceTot, 0.0);
        for (label fi = 0; fi < nIntF; ++fi) { dfluxWFaces[fi] = dfluxW[fi]; }
        {
            label bf = nIntF;
            forAll(dfluxW.boundaryField(), patchi)
            {
                const scalarField& dfB = dfluxW.boundaryField()[patchi];
                forAll(dfB, f) { dfluxWFaces[bf] = dfB[f]; ++bf; }
            }
        }
        // divergence with per-term decomposition
        forAll(anRPwOut, celli)
        {
            anRPwOut[celli] = 0.0;
            anRPaHbyAwOut[celli] = 0.0;
            anRPaFluxwOut[celli] = 0.0;
        }
        forAll(mesh.owner(), fi)
        {
            const label own = mesh.owner()[fi];
            const label nei = mesh.neighbour()[fi];
            anRPwOut[own] += dphiHbyAw[fi] - dfluxWFaces[fi];
            anRPwOut[nei] -= dphiHbyAw[fi] - dfluxWFaces[fi];
            anRPaHbyAwOut[own] += dphiHbyAw[fi];
            anRPaHbyAwOut[nei] -= dphiHbyAw[fi];
            anRPaFluxwOut[own] -= dfluxWFaces[fi];
            anRPaFluxwOut[nei] += dfluxWFaces[fi];
        }
        {
            label bf = nIntF;
            forAll(mesh.boundary(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                forAll(fc, f)
                {
                    const label celli = fc[f];
                    anRPwOut[celli] += dphiHbyAw[bf] - dfluxWFaces[bf];
                    anRPaHbyAwOut[celli] += dphiHbyAw[bf];
                    anRPaFluxwOut[celli] -= dfluxWFaces[bf];
                    ++bf;
                }
            }
        }
    };

    // B1: the RX-A FD is [R(alpha+eps*deltaAlpha)-R(alpha-eps*deltaAlpha)]/(2eps)
    // = dR/dalpha * deltaAlpha (directional), so the analytic side must be the
    // DIRECTIONAL derivative J*deltaAlpha.  R_U is DIAGONAL in alpha (fvm::Sp
    // touches only the cell's own diagonal) -> the per-cell product anRUa*da is
    // exact.  R_P = div(phi) is NON-diagonal (face fluxes interpolate HbyA from
    // own+nei and the laplacian(rAU,p) flux couples neighbours) -> the weight
    // must be multiplied into dHbyA/drAU BEFORE the face-flux assembly and the
    // divergence (assembleWeightedRPa): (J_P*1)*da is WRONG -- it equals J_P*1
    // only when da is spatially uniform -- while J_P*da is the true directional
    // derivative.  This fixes the cycle-1 RX-A R_P,alpha relL2=1.0115 gap.
    scalarField anRUaW(3*nC, 0.0);
    forAll(U, celli)
    {
        const scalar da = deltaAlpha[celli];
        anRUaW[3*celli + 0] = anRUa[3*celli + 0]*da;
        anRUaW[3*celli + 1] = anRUa[3*celli + 1]*da;
        anRUaW[3*celli + 2] = anRUa[3*celli + 2]*da;
    }
    scalarField anRPaW(nC, 0.0);
    scalarField anRPaHbyAW(nC, 0.0);
    scalarField anRPaFluxW(nC, 0.0);
    assembleWeightedRPa(deltaAlpha, anRPaW, anRPaHbyAW, anRPaFluxW);
    // weighted per-term norms (sum of the two terms == total == FD oracle)
    {
        scalar nH = 0.0, nF = 0.0, nT = 0.0;
        forAll(anRPaW, celli)
        {
            nH += sqr(anRPaHbyAW[celli]);
            nF += sqr(anRPaFluxW[celli]);
            nT += sqr(anRPaW[celli]);
        }
        reduce(nH, sumOp<scalar>());
        reduce(nF, sumOp<scalar>());
        reduce(nT, sumOp<scalar>());
        Info<< "StageB6 RX-A weighted R_P,alpha terms (J_P*deltaAlpha):"
            << " |div(dphiHbyA_w)|L2=" << Foam::sqrt(nH)
            << " |-div(dflux_w)|L2=" << Foam::sqrt(nF)
            << " |total_w|L2=" << Foam::sqrt(nT)
            << " (two weighted terms sum to total; total == FD oracle if closed)"
            << endl;
    }

    // RX-A FD sweep (scoped OFstreams: flush after each write batch, close on
    // scope exit -> no tail truncation even if a later gate aborts)
    {
        OFstream rxaFD("stageB6_rxa_FD_all_eps.mtx");
        OFstream rxaAn("stageB6_rxa_analytic.mtx");
        OFstream rpaTerms("stageB6_rpa_terms.mtx");
        rxaFD.precision(16); rxaAn.precision(16); rpaTerms.precision(16);
        // weighted analytic fields (compared directly against rxa_FD_all_eps,
        // no further deltaAlpha multiplication needed):
        //   [3N R_U,alpha*deltaAlpha ; N R_P,alpha*deltaAlpha]
        forAll(anRUaW, k) { rxaAn<< anRUaW[k] << nl; }
        forAll(anRPaW, k) { rxaAn<< anRPaW[k] << nl; }
        // weighted per-term decomposition (per cell: col1+col2 == col3 == FD):
        //   [N div(dphiHbyA_w) ; N -div(dflux_w) ; N total_w]
        forAll(anRPaHbyAW, k) { rpaTerms<< anRPaHbyAW[k] << nl; }
        forAll(anRPaFluxW, k) { rpaTerms<< anRPaFluxW[k] << nl; }
        forAll(anRPaW, k) { rpaTerms<< anRPaW[k] << nl; }
        rxaAn.flush(); rpaTerms.flush();

        for (label ei = 0; ei < 5; ++ei)
        {
            const scalar eps = stageB6Eps[ei];
            volScalarField alphaP(alpha);
            volScalarField alphaM(alpha);
            forAll(alphaP.primitiveFieldRef(), celli)
            {
                alphaP.primitiveFieldRef()[celli] += eps*deltaAlpha[celli];
                alphaM.primitiveFieldRef()[celli] -= eps*deltaAlpha[celli];
            }
            alphaP.correctBoundaryConditions();
            alphaM.correctBoundaryConditions();

            surfaceScalarField phiP("stageB6_phiP", phi);
            surfaceScalarField phiM("stageB6_phiM", phi);
            vectorField resUP(nC, vector::zero), resUM(nC, vector::zero);
            scalarField rpP(nC, 0.0), rpM(nC, 0.0);
            scalarField rAUdummy(nC, 0.0);
            vectorField Hdummy(nC, vector::zero), HbyAdummy(nC, vector::zero);
            rebuildResidual(alphaP, phiP, resUP, rpP, rAUdummy, Hdummy, HbyAdummy);
            rebuildResidual(alphaM, phiM, resUM, rpM, rAUdummy, Hdummy, HbyAdummy);

            scalarField fdRUa(3*nC, 0.0);
            scalarField fdRPa(nC, 0.0);
            forAll(U, celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    fdRUa[3*celli + c] =
                        (resUP[celli].component(c) - resUM[celli].component(c))
                       /(2.0*eps);
                }
                fdRPa[celli] = (rpP[celli] - rpM[celli])/(2.0*eps);
            }

            Info<< "StageB6 RX-A eps=" << eps << ":";
            stageB6Metrics(anRUaW, fdRUa, 3*nC, " R_U,alpha");
            stageB6Metrics(anRPaW, fdRPa, nC, " R_P,alpha");

            forAll(fdRUa, k) { rxaFD<< fdRUa[k] << nl; }
            forAll(fdRPa, k) { rxaFD<< fdRPa[k] << nl; }
            rxaFD.flush();
        }
    }

    // ------------------------------------------------------------------
    // 7. RX-B: xh perturbation -> production alpha interpolation (task §7).
    //    analytic = R_alpha*(dAlphaDxh*deltaXh) with the RX-A-verified
    //    operator (anRUa per-cell for U; assembleWeightedRPa for P).
    // ------------------------------------------------------------------
    Info<< "StageB6 RX-B: fixed U/p/phi/turbulence, deterministic deltaXh,"
        << " alpha via production interpolation" << endl;

    scalarField anRUxh(3*nC, 0.0);
    scalarField anRPxh(nC, 0.0);
    scalarField wXh(nC, 0.0);
    forAll(U, celli)
    {
        const scalar f = dAlphaDxh[celli]*deltaXh[celli];
        wXh[celli] = f;
        anRUxh[3*celli + 0] = anRUa[3*celli + 0]*f;
        anRUxh[3*celli + 1] = anRUa[3*celli + 1]*f;
        anRUxh[3*celli + 2] = anRUa[3*celli + 2]*f;
    }
    // P-row: alpha is DIAGONAL in xh (alpha[celli] depends only on xh[celli]),
    // but R_P is NON-diagonal in alpha -> assemble J_P*(dAlphaDxh*deltaXh)
    // with the weight inside the face-flux assembly (not (J_P*1)*f pointwise).
    scalarField anRPaHbyAxh(nC, 0.0);
    scalarField anRPaFluxxh(nC, 0.0);
    assembleWeightedRPa(wXh, anRPxh, anRPaHbyAxh, anRPaFluxxh);

    {
        OFstream rxbFD("stageB6_rxb_FD_all_eps.mtx");
        OFstream rxbAn("stageB6_rxb_analytic.mtx");
        rxbFD.precision(16); rxbAn.precision(16);
        forAll(anRUxh, k) { rxbAn<< anRUxh[k] << nl; }
        forAll(anRPxh, k) { rxbAn<< anRPxh[k] << nl; }
        rxbAn.flush();

        for (label ei = 0; ei < 5; ++ei)
        {
            const scalar eps = stageB6Eps[ei];
            volScalarField xhP(xh);
            volScalarField xhM(xh);
            forAll(xhP.primitiveFieldRef(), celli)
            {
                xhP.primitiveFieldRef()[celli] += eps*deltaXh[celli];
                xhM.primitiveFieldRef()[celli] -= eps*deltaXh[celli];
            }
            xhP.correctBoundaryConditions();
            xhM.correctBoundaryConditions();

            volScalarField alphaP
            (
                IOobject("stageB6_alphaP", runTime.timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("a", alpha.dimensions(), 0.0)
            );
            volScalarField alphaM
            (
                IOobject("stageB6_alphaM", runTime.timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("a", alpha.dimensions(), 0.0)
            );
            forAll(alphaP.primitiveFieldRef(), celli)
            {
                const scalar vP = alphaMax.value()*qu*(1.0 - xhP[celli])
                                / (qu + xhP[celli] + SMALL);
                const scalar vM = alphaMax.value()*qu*(1.0 - xhM[celli])
                                / (qu + xhM[celli] + SMALL);
                alphaP.primitiveFieldRef()[celli] =
                    Foam::max(alphamin.value(), vP);
                alphaM.primitiveFieldRef()[celli] =
                    Foam::max(alphamin.value(), vM);
            }
            alphaP.correctBoundaryConditions();
            alphaM.correctBoundaryConditions();

            surfaceScalarField phiP("stageB6_phiP2", phi);
            surfaceScalarField phiM("stageB6_phiM2", phi);
            vectorField resUP(nC, vector::zero), resUM(nC, vector::zero);
            scalarField rpP(nC, 0.0), rpM(nC, 0.0);
            scalarField rAUdummy(nC, 0.0);
            vectorField Hdummy(nC, vector::zero), HbyAdummy(nC, vector::zero);
            rebuildResidual(alphaP, phiP, resUP, rpP, rAUdummy, Hdummy, HbyAdummy);
            rebuildResidual(alphaM, phiM, resUM, rpM, rAUdummy, Hdummy, HbyAdummy);

            scalarField fdRUxh(3*nC, 0.0);
            scalarField fdRPxh(nC, 0.0);
            forAll(U, celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    fdRUxh[3*celli + c] =
                        (resUP[celli].component(c) - resUM[celli].component(c))
                       /(2.0*eps);
                }
                fdRPxh[celli] = (rpP[celli] - rpM[celli])/(2.0*eps);
            }

            Info<< "StageB6 RX-B eps=" << eps << ":";
            stageB6Metrics(anRUxh, fdRUxh, 3*nC, " R_U,xh");
            stageB6Metrics(anRPxh, fdRPxh, nC, " R_P,xh");

            forAll(fdRUxh, k) { rxbFD<< fdRUxh[k] << nl; }
            forAll(fdRPxh, k) { rxbFD<< fdRPxh[k] << nl; }
            rxbFD.flush();
        }
    }

    // ------------------------------------------------------------------
    // 8. RX-C: raw-design directions via filter+projection+alpha (task §8).
    //    FD side: exact production chain (filter solve, bisected eta,
    //    projection, alpha interpolation), U/p fixed.
    //    Analytic side: R_alpha*(dAlphaDxh*z), z = dxh/dx*d via the exact
    //    filter tangent solve and the exact projection Jacobian (incl. the
    //    implicit adaptive-eta derivative).
    // ------------------------------------------------------------------
    Info<< "StageB6 RX-C: D1/D2/D3 raw-design directions, exact "
        << "filter->projection->alpha chain, U/p fixed" << endl;

    // drho / dProjectionDeta at the CURRENT xp/eta5 (exact production
    // formulas, filter_chainrule.H L35-71).
    scalarField drhoF(nC, 0.0);
    scalarField dProjEta(nC, 0.0);
    scalar projEtaDenom = 0.0;
    if (del > SMALL)
    {
        forAll(xp, celli)
        {
            if (!activeCell[celli]) continue;
            if (xp[celli] <= eta5)
            {
                const scalar pe = Foam::exp(-del*(1.0 - xp[celli]/eta5));
                drhoF[celli] = del*pe + Foam::exp(-del);
                dProjEta[celli] =
                    pe*(1.0 - del*xp[celli]/eta5) - Foam::exp(-del);
            }
            else
            {
                const scalar pe = Foam::exp
                (
                    -del*(xp[celli] - eta5)/(1.0 - eta5)
                );
                drhoF[celli] = del*pe + Foam::exp(-del);
                dProjEta[celli] =
                    pe*(1.0 - del*(1.0 - xp[celli])/(1.0 - eta5))
                  - Foam::exp(-del);
            }
            projEtaDenom += designMask[celli]*dProjEta[celli]*mesh.V()[celli];
        }
        reduce(projEtaDenom, sumOp<scalar>());
    }
    Info<< "StageB6 RX-C projection: eta5=" << eta5 << " del=" << del
        << " projectionEtaDenominator=" << projEtaDenom << endl;

    // D4 fix: production divides by the RAW (possibly negative) eta denominator
    // with a mag()<=SMALL singularity guard (filter_chainrule.H L81/L108-111).
    // The old SMALL-clamping max() guard was WRONG: with a negative denominator
    // it returned +SMALL=1e-15 (since 1e-15 > -3.14e-6), a ~3.14e9 sign-flipping
    // amplification that poisoned RX-C/§9/§10 (cycle-1).  Mirror the production
    // guard: divide by the raw value; singular -> 0 + Info (diagnostic
    // fallback, production would FatalError).
    auto safeEtaDivide =
    [&](const scalar num) -> scalar
    {
        if (mag(projEtaDenom) <= SMALL)
        {
            Info<< "StageB6 warning: projectionEtaDenominator singular ("
                << projEtaDenom << "), eta response treated as 0" << endl;
            return 0.0;
        }
        return num/projEtaDenom;
    };

    // Filter tangent solve: (M-B) y = -b d.
    auto filterTangent =
    [&](const scalarField& d, scalarField& y)
    {
        volScalarField dfield
        (
            IOobject("stageB6_d", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("d", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        volScalarField yfield
        (
            IOobject("stageB6_y", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("y", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(dfield.primitiveFieldRef(), celli)
        {
            dfield.primitiveFieldRef()[celli] = d[celli];
        }
        dfield.correctBoundaryConditions();
        yfield.correctBoundaryConditions();
        dictionary filterSolver;
        filterSolver.add("solver", word("PCG"));
        filterSolver.add("preconditioner", word("DIC"));
        filterSolver.add("tolerance", scalar(1e-9));
        filterSolver.add("relTol", scalar(0.0));
        fvScalarMatrix filterEqn
        (
            fvm::laplacian(designFilterFaceMask, yfield)
          - fvm::Sp(b, yfield)
          + dfield*b
        );
        filterEqn.solve(filterSolver);
        forAll(y, celli) { y[celli] = yfield[celli]; }
    };

    // Projection tangent: z = drho*y + dProjEta*[sum(designMask*(1-drho)*V*y)
    //   /sum(designMask*dProjEta*V)]   (true tangent incl. implicit eta)
    auto projectionTangent =
    [&](const scalarField& y, scalarField& z)
    {
        scalar num = 0.0;
        forAll(y, celli)
        {
            num += designMask[celli]*(1.0 - drhoF[celli])*mesh.V()[celli]*y[celli];
        }
        reduce(num, sumOp<scalar>());
        const scalar etaResp = safeEtaDivide(num);
        forAll(z, celli)
        {
            z[celli] = drhoF[celli]*y[celli] + dProjEta[celli]*etaResp;
        }
    };

    // Adaptive-eta bisection + projection on a LOCAL xpMod (exact production
    // path, filter_x.H), writing xhMod.
    auto projectLocal =
    [&](const volScalarField& xpMod, volScalarField& xhMod)
    {
        scalar e0 = 0.0001;
        scalar e1 = 0.9999;
        scalar etaB = 0.5;
        scalar y0 = diff(xpMod, mesh.V(), designMask, del, e0, n);
        scalar y1 = diff(xpMod, mesh.V(), designMask, del, e1, n);
        reduce(y0, sumOp<scalar>());
        reduce(y1, sumOp<scalar>());
        if (y0*y1 > 0.0)
        {
            etaB = 0.5;
        }
        else
        {
            do
            {
                etaB = 0.5*(e0 + e1);
                scalar y5 = diff(xpMod, mesh.V(), designMask, del, etaB, n);
                reduce(y5, sumOp<scalar>());
                if (y0*y5 < 0.0) { e1 = etaB; y1 = y5; }
                else             { e0 = etaB; y0 = y5; }
            }
            while ((e1 - e0) > projectionEtaTolerance);
            etaB = 0.5*(e0 + e1);
        }
        forAll(xhMod.primitiveFieldRef(), celli)
        {
            if (designMask[celli] <= 0.5) continue;
            if (xpMod[celli] <= etaB)
            {
                xhMod.primitiveFieldRef()[celli] = etaB*
                (
                    Foam::exp(-del*(1.0 - xpMod[celli]/etaB))
                  - (1.0 - xpMod[celli]/etaB)*Foam::exp(-del)
                );
            }
            else
            {
                xhMod.primitiveFieldRef()[celli] = etaB + (1.0 - etaB)*
                (
                    1.0
                  - Foam::exp(-del*(xpMod[celli] - etaB)/(1.0 - etaB))
                  + (xpMod[celli] - etaB)*Foam::exp(-del)/(1.0 - etaB)
                );
            }
        }
        xhMod.correctBoundaryConditions();
    };

    // Analytic chain per direction: y = filterTangent(d); z = projectionTangent(y);
    // dR/dx*d = R_alpha*(dAlphaDxh*z) (U per-cell; P via assembleWeightedRPa).
    auto analyticRxDirection =
    [&](const scalarField& d, scalarField& anRUd, scalarField& anRPd,
        scalarField& zOut)
    {
        scalarField y(nC, 0.0);
        filterTangent(d, y);
        projectionTangent(y, zOut);
        scalarField w(nC, 0.0);
        forAll(U, celli)
        {
            const scalar f = dAlphaDxh[celli]*zOut[celli];
            w[celli] = f;
            anRUd[3*celli + 0] = anRUa[3*celli + 0]*f;
            anRUd[3*celli + 1] = anRUa[3*celli + 1]*f;
            anRUd[3*celli + 2] = anRUa[3*celli + 2]*f;
        }
        // P-row: assemble J_P*(dAlphaDxh*z) with the weight inside the
        // face-flux assembly (RX-A-verified operator; not (J_P*1)*f pointwise).
        scalarField anRPaHbyAd(nC, 0.0);
        scalarField anRPaFluxd(nC, 0.0);
        assembleWeightedRPa(w, anRPd, anRPaHbyAd, anRPaFluxd);
    };

    const scalarField* dirs[3] = {&D1, &D2, &D3};
    const char* dirNames[3] = {"D1", "D2", "D3"};

    {
        OFstream rxcFD("stageB6_rxc_FD_all_eps.mtx");
        OFstream rxcAn("stageB6_rxc_analytic.mtx");
        OFstream rxcZ("stageB6_rxc_z_analytic.mtx");
        OFstream rxcXhFD("stageB6_rxc_xh_FD_all_eps.mtx");
        OFstream rxcXhAn("stageB6_rxc_xh_analytic.mtx");
        OFstream rxcAlphFD("stageB6_rxc_alpha_FD_all_eps.mtx");
        rxcFD.precision(16); rxcAn.precision(16); rxcZ.precision(16);
        rxcXhFD.precision(16); rxcXhAn.precision(16); rxcAlphFD.precision(16);

        for (label di = 0; di < 3; ++di)
        {
            const scalarField& d = *dirs[di];
            scalarField anRUd(3*nC, 0.0), anRPd(nC, 0.0), zAn(nC, 0.0);
            analyticRxDirection(d, anRUd, anRPd, zAn);
            forAll(anRUd, k) { rxcAn<< anRUd[k] << nl; }
            forAll(anRPd, k) { rxcAn<< anRPd[k] << nl; }
            forAll(zAn, k) { rxcZ<< zAn[k] << nl; }
            forAll(zAn, k) { rxcXhAn<< zAn[k] << nl; }
            rxcAn.flush(); rxcZ.flush(); rxcXhAn.flush();

            for (label ei = 0; ei < 5; ++ei)
            {
                const scalar eps = stageB6Eps[ei];

                volScalarField xP(x);
                volScalarField xM(x);
                forAll(xP.primitiveFieldRef(), celli)
                {
                    xP.primitiveFieldRef()[celli] += eps*d[celli];
                    xM.primitiveFieldRef()[celli] -= eps*d[celli];
                }
                xP.correctBoundaryConditions();
                xM.correctBoundaryConditions();

                // Helmholtz filter (exact production equation)
                dictionary filterSolver;
                filterSolver.add("solver", word("PCG"));
                filterSolver.add("preconditioner", word("DIC"));
                filterSolver.add("tolerance", scalar(1e-9));
                filterSolver.add("relTol", scalar(0.0));
                volScalarField xpP(xp);
                volScalarField xpM(xp);
                {
                    fvScalarMatrix filterEqnP
                    (
                        fvm::laplacian(designFilterFaceMask, xpP)
                      - fvm::Sp(b, xpP)
                      + xP*b
                    );
                    filterEqnP.solve(filterSolver);
                    fvScalarMatrix filterEqnM
                    (
                        fvm::laplacian(designFilterFaceMask, xpM)
                      - fvm::Sp(b, xpM)
                      + xM*b
                    );
                    filterEqnM.solve(filterSolver);
                }

                // adaptive-Heaviside projection
                volScalarField xhP(xh);
                volScalarField xhM(xh);
                projectLocal(xpP, xhP);
                projectLocal(xpM, xhM);

                // alpha via the exact production interpolation
                volScalarField alphaP
                (
                    IOobject("stageB6_alphaP3", runTime.timeName(), mesh,
                             IOobject::NO_READ, IOobject::NO_WRITE),
                    mesh,
                    dimensionedScalar("a", alpha.dimensions(), 0.0)
                );
                volScalarField alphaM
                (
                    IOobject("stageB6_alphaM3", runTime.timeName(), mesh,
                             IOobject::NO_READ, IOobject::NO_WRITE),
                    mesh,
                    dimensionedScalar("a", alpha.dimensions(), 0.0)
                );
                forAll(alphaP.primitiveFieldRef(), celli)
                {
                    const scalar vP = alphaMax.value()*qu*(1.0 - xhP[celli])
                                    / (qu + xhP[celli] + SMALL);
                    const scalar vM = alphaMax.value()*qu*(1.0 - xhM[celli])
                                    / (qu + xhM[celli] + SMALL);
                    alphaP.primitiveFieldRef()[celli] =
                        Foam::max(alphamin.value(), vP);
                    alphaM.primitiveFieldRef()[celli] =
                        Foam::max(alphamin.value(), vM);
                }
                alphaP.correctBoundaryConditions();
                alphaM.correctBoundaryConditions();

                // residual rebuild (U/p fixed)
                surfaceScalarField phiP("stageB6_phiP4", phi);
                surfaceScalarField phiM("stageB6_phiM4", phi);
                vectorField resUP(nC, vector::zero), resUM(nC, vector::zero);
                scalarField rpP(nC, 0.0), rpM(nC, 0.0);
                scalarField rAUdummy(nC, 0.0);
                vectorField Hdummy(nC, vector::zero), HbyAdummy(nC, vector::zero);
                rebuildResidual(alphaP, phiP, resUP, rpP, rAUdummy, Hdummy, HbyAdummy);
                rebuildResidual(alphaM, phiM, resUM, rpM, rAUdummy, Hdummy, HbyAdummy);

                scalarField fdRUd(3*nC, 0.0);
                scalarField fdRPd(nC, 0.0);
                scalarField zFD(nC, 0.0);
                scalarField aFD(nC, 0.0);
                forAll(U, celli)
                {
                    for (direction c = 0; c < 3; ++c)
                    {
                        fdRUd[3*celli + c] =
                            (resUP[celli].component(c) - resUM[celli].component(c))
                           /(2.0*eps);
                    }
                    fdRPd[celli] = (rpP[celli] - rpM[celli])/(2.0*eps);
                    zFD[celli] = (xhP[celli] - xhM[celli])/(2.0*eps);
                    aFD[celli] = (alphaP[celli] - alphaM[celli])/(2.0*eps);
                }

                // internal gates: xh tangent and alpha tangent closure
                {
                    scalarField zAnAl(nC, 0.0);
                    forAll(zAnAl, celli) { zAnAl[celli] = zAn[celli]; }
                    Info<< "StageB6 RX-C " << dirNames[di] << " eps=" << eps << ":";
                    stageB6Metrics(zAnAl, zFD, nC, " xh-tangent");
                    scalarField aAn(nC, 0.0);
                    forAll(aAn, celli) { aAn[celli] = dAlphaDxh[celli]*zAn[celli]; }
                    stageB6Metrics(aAn, aFD, nC, " alpha-tangent");
                    stageB6Metrics(anRUd, fdRUd, 3*nC, " R_U,xd");
                    stageB6Metrics(anRPd, fdRPd, nC, " R_P,xd");
                }

                forAll(fdRUd, k) { rxcFD<< fdRUd[k] << nl; }
                forAll(fdRPd, k) { rxcFD<< fdRPd[k] << nl; }
                forAll(zFD, k) { rxcXhFD<< zFD[k] << nl; }
                forAll(aFD, k) { rxcAlphFD<< aFD[k] << nl; }
                rxcFD.flush(); rxcXhFD.flush(); rxcAlphFD.flush();
            }
        }
    }

    // ------------------------------------------------------------------
    // 9. Volume-gradient chain regression anchor (production replica).
    //    gsenshVol = designMask*V/designVolume -> drho + eta correction ->
    //    adjoint filter -> gsensVol ; projV_Dk = sum_active gsensVol*Dk.
    // ------------------------------------------------------------------
    {
        volScalarField gsenshVolR
        (
            IOobject("stageB6_gsenshVol", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("g", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(gsenshVolR.primitiveFieldRef(), celli)
        {
            gsenshVolR.primitiveFieldRef()[celli] =
                designMask[celli]*mesh.V()[celli]/Foam::max(designVolume, SMALL);
        }
        // eta correction (production filter_chainrule.H order)
        scalar volEtaCorr = 0.0;
        forAll(gsenshVolR, celli)
        {
            volEtaCorr +=
                designMask[celli]*dProjEta[celli]*gsenshVolR[celli];
        }
        reduce(volEtaCorr, sumOp<scalar>());
        volEtaCorr = safeEtaDivide(volEtaCorr);
        forAll(gsenshVolR.primitiveFieldRef(), celli)
        {
            gsenshVolR.primitiveFieldRef()[celli] *= drhoF[celli];
            gsenshVolR.primitiveFieldRef()[celli] +=
                designMask[celli]*mesh.V()[celli]*(1.0 - drhoF[celli])
               *volEtaCorr;
        }
        gsenshVolR.correctBoundaryConditions();

        volScalarField gsensVolR
        (
            IOobject("stageB6_gsensVol", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("g", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        gsensVolR.correctBoundaryConditions();
        dictionary filterSolver;
        filterSolver.add("solver", word("PCG"));
        filterSolver.add("preconditioner", word("DIC"));
        filterSolver.add("tolerance", scalar(1e-9));
        filterSolver.add("relTol", scalar(0.0));
        fvScalarMatrix volFilterEqn
        (
            fvm::laplacian(designFilterFaceMask, gsensVolR)
          - fvm::Sp(b, gsensVolR)
          + gsenshVolR*b
        );
        volFilterEqn.solve(filterSolver);
        gsensVolR.primitiveFieldRef() *= designMask.primitiveField();

        const scalarField* dirsV[3] = {&D1, &D2, &D3};
        const char* dirNamesV[3] = {"D1", "D2", "D3"};
        for (label di = 0; di < 3; ++di)
        {
            scalar projV = 0.0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                projV += gsensVolR[celli]*(*dirsV[di])[celli];
            }
            reduce(projV, sumOp<scalar>());
            std::ostringstream projVStr;
            projVStr.precision(16);
            projVStr << projV;
            Info<< "StageB6 volume-anchor " << dirNamesV[di]
                << ": projV=" << projVStr.str() << endl;
        }
        {
            OFstream vOut("stageB6_prod_gsensVol.mtx");
            vOut.precision(16);
            forAll(gsensVolR, celli) { vOut<< gsensVolR[celli] << nl; }
            vOut.flush();
        }
    }

    // ------------------------------------------------------------------
    // 10. Production gsenshPressureDrop chain replica (momentum-only) and
    //     the §9/§10 weighted contributions.
    // ------------------------------------------------------------------
    volScalarField gsenshR
    (
        IOobject("stageB6_gsenshPressureDrop", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("g", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(gsenshR.primitiveFieldRef(), celli)
    {
        gsenshR.primitiveFieldRef()[celli] =
           -dAlphaDxh[celli]*(U[celli] & Uc[celli])*mesh.V()[celli];
    }
    // designMask (filter_chainrule.H L12)
    gsenshR.primitiveFieldRef() *= designMask.primitiveField();
    // eta correction (production order: numerator uses the masked source)
    {
        scalar pressEtaCorr = 0.0;
        forAll(gsenshR, celli)
        {
            pressEtaCorr +=
                designMask[celli]*dProjEta[celli]*gsenshR[celli];
        }
        reduce(pressEtaCorr, sumOp<scalar>());
        pressEtaCorr = safeEtaDivide(pressEtaCorr);
        forAll(gsenshR.primitiveFieldRef(), celli)
        {
            gsenshR.primitiveFieldRef()[celli] *= drhoF[celli];
            gsenshR.primitiveFieldRef()[celli] +=
                designMask[celli]*mesh.V()[celli]*(1.0 - drhoF[celli])
               *pressEtaCorr;
        }
    }
    gsenshR.correctBoundaryConditions();

    volScalarField gsensR
    (
        IOobject("stageB6_gsensPressureDrop", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("g", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
    gsensR.correctBoundaryConditions();
    {
        dictionary filterSolver;
        filterSolver.add("solver", word("PCG"));
        filterSolver.add("preconditioner", word("DIC"));
        filterSolver.add("tolerance", scalar(1e-9));
        filterSolver.add("relTol", scalar(0.0));
        fvScalarMatrix pressFilterEqn
        (
            fvm::laplacian(designFilterFaceMask, gsensR)
          - fvm::Sp(b, gsensR)
          + gsenshR*b
        );
        pressFilterEqn.solve(filterSolver);
    }
    gsensR.primitiveFieldRef() *= designMask.primitiveField();
    {
        OFstream ghOut("stageB6_prod_gsensh.mtx");
        OFstream gOut("stageB6_prod_gsens.mtx");
        ghOut.precision(16); gOut.precision(16);
        forAll(gsenshR, celli) { ghOut<< gsenshR[celli] << nl; }
        forAll(gsensR, celli) { gOut<< gsensR[celli] << nl; }
        ghOut.flush(); gOut.flush();
    }

    // ---- §9/§10 weighted contributions (analytic R_x, eps-independent) ----
    Info<< "StageB6 §9/§10: weighted design-derivative contributions "
        << "(analytic chain; pressureGradientScale=" << pressureGradientScale
        << ")" << endl;
    {
        // lambda_U = Uc (3N), lambda_P = pc (N)
        OFstream lamOut("stageB6_lambda.mtx");
        lamOut.precision(16);
        forAll(Uc, celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                lamOut<< Uc[celli].component(c) << nl;
            }
        }
        forAll(pc, celli) { lamOut<< pc[celli] << nl; }
        lamOut.flush();

        OFstream daOut("stageB6_deltaAlpha.mtx");
        OFstream dxOut("stageB6_deltaXh.mtx");
        OFstream dadxOut("stageB6_dalphadxh.mtx");
        OFstream dirOut("stageB6_dirs.mtx");
        OFstream ctOut("stageB6_celltype.mtx");
        OFstream prOut("stageB6_proj.mtx");
        daOut.precision(16); dxOut.precision(16); dadxOut.precision(16);
        dirOut.precision(16); prOut.precision(16);
        forAll(deltaAlpha, celli)
        {
            daOut<< deltaAlpha[celli] << nl;
            dxOut<< deltaXh[celli] << nl;
            dadxOut<< dAlphaDxh[celli] << nl;
            ctOut<< (activeCell[celli] ? (unclippedCell[celli] ? 2 : 1) : 0) << nl;
            prOut<< drhoF[celli] << nl;
        }
        forAll(D1, celli)
        {
            dirOut<< D1[celli] << nl;
        }
        forAll(D2, celli)
        {
            dirOut<< D2[celli] << nl;
        }
        forAll(D3, celli)
        {
            dirOut<< D3[celli] << nl;
        }
        daOut.flush(); dxOut.flush(); dadxOut.flush();
        dirOut.flush(); ctOut.flush(); prOut.flush();

        for (label di = 0; di < 3; ++di)
        {
            const scalarField& d = *dirs[di];
            scalarField anRUd(3*nC, 0.0), anRPd(nC, 0.0), zAn(nC, 0.0);
            analyticRxDirection(d, anRUd, anRPd, zAn);

            // norms
            scalar nRU = 0.0, nRP = 0.0;
            forAll(U, celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    nRU += sqr(anRUd[3*celli + c]);
                }
                nRP += sqr(anRPd[celli]);
            }
            reduce(nRU, sumOp<scalar>());
            reduce(nRP, sumOp<scalar>());
            const scalar ratio = Foam::sqrt(nRP/Foam::max(nRU, SMALL));

            // lambda-weighted
            scalar lU = 0.0, lP = 0.0;
            forAll(U, celli)
            {
                for (direction c = 0; c < 3; ++c)
                {
                    lU += Uc[celli].component(c)*anRUd[3*celli + c];
                }
                lP += pc[celli]*anRPd[celli];
            }
            reduce(lU, sumOp<scalar>());
            reduce(lP, sumOp<scalar>());
            const scalar Dmom = -lU;
            const scalar Dpres = -lP;
            const scalar Dtot = Dmom + Dpres;

            // production projected contribution
            scalar prodContrib = 0.0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                prodContrib += gsensR[celli]*d[celli];
            }
            reduce(prodContrib, sumOp<scalar>());
            const scalar prodScaled = pressureGradientScale*prodContrib;

            Info<< "StageB6 weighted " << dirNames[di] << ":"
                << " ||R_U,x d||L2=" << Foam::sqrt(nRU)
                << " ||R_P,x d||L2=" << Foam::sqrt(nRP)
                << " ratio||R_P||/||R_U||=" << ratio
                << " lambda_U^T R_U,x d=" << lU
                << " lambda_P^T R_P,x d=" << lP
                << " D_momentum=" << Dmom
                << " D_pressure=" << Dpres
                << " D_total=" << Dtot
                << " prodGsenDPressDrop*d=" << prodContrib
                << " prodDgdx1*d(scaled)=" << prodScaled
                << " |D_mom-prod|/|prod|="
                << (mag(prodContrib) > SMALL
                      ? mag(Dmom - prodContrib)/mag(prodContrib) : 0.0)
                << endl;
        }
    }

    stageB6Fingerprint();
    Info<< "--- StageB6 oracle complete (artifacts: stageB6_*.mtx) ---" << endl;
}
