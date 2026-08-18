// =============================================================================
// stageB8JPPActualResidualFD.H  (BFINAL-008, PATCH-validation probe)
// -----------------------------------------------------------------------------
// Switch-guarded (optProperties "stageB8JPPActualResidualFD", default false).
// Must be #included ONLY at the end of src/solveDiscreteFlowAdjoint.H, at a
// point where the reduced-J state (applyDiscreteFlowJ, rAUAdj,
// discretePrimalMomentum, discreteOwner/Neighbour/Weights/Sf/Phi,
// discreteUIndex/PIndex, primalPressureMobility, U, p, phi, alpha,
// nuEffFrozen, mesh, runTime, optProperties, simple, pRefCell, pRefValue,
// fvOptions, ransFlowModel) is in scope and the flow state is the CONVERGED
// frozen production state.  No-op when the switch is false (the caller guards
// this block with the dictionary lookup).
//
// The probe rebuilds the ACTUAL frozen-primal residual pair
//   R(w) = [R_U(w); R_P(w)],  w = [U; p]
// from production semantics and central-FD's it along deterministic
// directions:
//   U rows:  UEqn = fvm::div(phi,U) - fvm::laplacian(nuEffFrozen,U)
//            + fvm::Sp(alpha,U) == -fvc::grad(p) + fvOptions(U)
//            (+ dev2 when ransFlowModel), relax, constrain;
//            resU = -UEqn.residual()  (NS.H sign, stageB5-validated pairing).
//   P rows:  rAU = 1/A(); HbyA = constrainHbyA(rAU*H(), U, p);
//            phiHbyA = fvc::flux(HbyA); adjustPhi; rAtU (consistent branch);
//            constrainPressure(p, U, phiHbyA, rAtU());
//            pEqn = fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA);
//            pEqn.setReference(pRefCell, pRefValue);
//            phi = phiHbyA - pEqn.flux()   (production flux() semantics incl.
//            the fixedValue outlet internalCoeffs boundary contribution);
//            resP = div(phi)  (raw face-flux sum, owner +, neighbour -,
//            boundary +, NO 1/V).
//   FD:      [R(w+eps*v) - R(w-eps*v)]/(2 eps), eps sweep.
//   P1:      v = [0; n_P] (dp-only), n_P = normalized P block of the
//            BFINAL-006 null vector wPrime_TAN_D1.mtx: compare J_corrected*v
//            (applyDiscreteFlowJ) vs FD; report full-P / outlet-adjacent /
//            off-outlet relL2 + sign/structure.
//   P3:      deterministic mixed directions (n_P; dp-heavy outlet-supported;
//            seeded smooth dU+dp): report U rows / P-internal /
//            P-outlet-adjacent / P-other-boundary / P-total with relL2 /
//            cosine / normRatio / max / eps-plateau.
// Artifacts (Jv, FD per eps, directions, n_P, cell classification) are written
// under evidence/agent-group/BFINAL-008/cycle-1/artifacts/ only.
//
// NO production field/matrix is modified: every perturbation acts on LOCAL
// copies (Umod/pmod and their rebuilt UEqn/pEqn); the probe never solves
// UEqn/pEqn and never writes into the case.
// =============================================================================
#include <fstream>
#include <sstream>

{
    Info<< "--- stageB8JPPActualResidualFD: start ---" << endl;

    const label s8N = mesh.nCells();
    const label s8NIntF = mesh.nInternalFaces();
    const label s8NUnk = 4*s8N;
    const label s8NVel = 3*s8N;
    const fileName s8ArtDir
    (
        "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
        "BFINAL-008/cycle-1/artifacts"
    );
    Foam::mkDir(s8ArtDir);

    // ---- (0) fixesValue-p patches -> outlet-adjacent cell set ----
    boolList s8OutletCell(s8N, false);
    label s8OutletPatchi = -1;
    forAll(mesh.boundary(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            s8OutletPatchi = patchi;
            Info<< "  stageB8: fixesValue(p) patch = "
                << mesh.boundary()[patchi].name()
                << "  nFaces=" << mesh.boundary()[patchi].size() << endl;
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(fc, f) { s8OutletCell[fc[f]] = true; }
        }
    }
    if (s8OutletPatchi < 0)
    {
        FatalErrorInFunction
            << "stageB8: no fixesValue pressure patch found"
            << abort(FatalError);
    }

    // ---- (1) cell classification: internal / other-boundary / outlet ----
    boolList s8HasBnd(s8N, false);
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        forAll(fc, f) { s8HasBnd[fc[f]] = true; }
    }
    labelList s8IdxPint, s8IdxPout, s8IdxPother, s8IdxPtot(s8N);
    for (label c = 0; c < s8N; ++c)
    {
        s8IdxPtot[c] = c;
        if (s8OutletCell[c])          s8IdxPout.append(c);
        else if (s8HasBnd[c])         s8IdxPother.append(c);
        else                          s8IdxPint.append(c);
    }
    labelList s8IdxMom(s8NVel);
    for (label i = 0; i < s8NVel; ++i) s8IdxMom[i] = i;
    Info<< "  stageB8: cells N=" << s8N
        << " internal=" << s8IdxPint.size()
        << " otherBoundary=" << s8IdxPother.size()
        << " outletAdjacent=" << s8IdxPout.size() << endl;

    // ---- (2) read the BFINAL-006 null vector, extract normalized P block ----
    scalarField s8W(s8NUnk, 0.0);
    {
        std::ifstream s8In
        (
            "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
            "BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx"
        );
        if (!s8In.good())
        {
            FatalErrorInFunction
                << "stageB8: cannot open wPrime_TAN_D1.mtx"
                << abort(FatalError);
        }
        for (label i = 0; i < s8NUnk; ++i)
        {
            double s8V = 0.0;
            s8In >> s8V;
            if (!s8In.good())
            {
                FatalErrorInFunction
                    << "stageB8: wPrime_TAN_D1.mtx truncated at " << i
                    << " of " << s8NUnk << abort(FatalError);
            }
            s8W[i] = scalar(s8V);
        }
    }
    const scalar s8Wnorm = Foam::sqrt(gSumSqr(s8W));
    scalarField s8nP(s8N, 0.0);
    for (label c = 0; c < s8N; ++c)
    {
        s8nP[c] = s8W[3*s8N + c]/s8Wnorm;
    }
    Info<< "  stageB8: null vector ||w||=" << s8Wnorm
        << "  ||n_P||=" << Foam::sqrt(gSumSqr(s8nP))
        << "  n_P[PREF=cell0]=" << s8nP[0] << endl;

    // ---- (3) ACTUAL frozen-primal residual rebuild (production semantics) ----
    // Takes local perturbed copies (Umod, pmod); returns resU (3N, NS.H sign)
    // and resP (N, raw div of the rebuilt phi incl. fixedValue outlet
    // internalCoeffs through pEqn.flux()).
    //
    // SEMANTICS (BFINAL-007 Q1-verified, and what the corrected J implements):
    //   * the U-row residual R_U responds to the pressure through the UEqn
    //     source -fvc::grad(pmod)  (this is the J's U-row dp contribution,
    //     the "pressureContribution" term);
    //   * the P-row phiHbyA is built from the PRODUCTION p (p-independent in
    //     value: constrainHbyA/adjustPhi only use the p BC *type*), so the
    //     only p-dependence of R_P = div(phiHbyA) - div(pEqn.flux()) enters
    //     through pEqn built on psi = pmod.  Rebuilding phiHbyA from pmod
    //     would inject a d(phiHbyA)/dp term that the J does NOT contain and
    //     that BFINAL-007's verified FD does not have (measured: such a
    //     rebuild mismatches the verified FD at relL2=0.37).
    auto s8Rebuild = [&](volVectorField& Umod, volScalarField& pmod,
                         vectorField& resU, scalarField& resP)
    {
        // U-row residual: UEqn with the perturbed pressure source
        tmp<fvVectorMatrix> tUEqn8
        (
            fvm::div(phi, Umod)
          - fvm::laplacian(nuEffFrozen, Umod)
          + fvm::Sp(alpha, Umod)
         ==
           -fvc::grad(pmod)
          + fvOptions(Umod)
        );
        fvVectorMatrix& UEqn8 = tUEqn8.ref();

        if (ransFlowModel)
        {
            UEqn8 -= fvc::div
            (
                nuEffFrozen*dev2(Foam::T(fvc::grad(Umod)))
            );
        }
        UEqn8.relax();
        fvOptions.constrain(UEqn8);
        const tmp<Field<vector>> tRes8 = UEqn8.residual();
        forAll(resU, celli)
        {
            resU[celli] = -tRes8()[celli];
        }

        // P-row phi path: rebuilt from the PRODUCTION p (p-independent
        // phiHbyA), exactly the BFINAL-007 Q1 semantics.
        tmp<fvVectorMatrix> tUEqnP8
        (
            fvm::div(phi, Umod)
          - fvm::laplacian(nuEffFrozen, Umod)
          + fvm::Sp(alpha, Umod)
         ==
           -fvc::grad(p)
          + fvOptions(Umod)
        );
        fvVectorMatrix& UEqnP8 = tUEqnP8.ref();
        if (ransFlowModel)
        {
            UEqnP8 -= fvc::div
            (
                nuEffFrozen*dev2(Foam::T(fvc::grad(Umod)))
            );
        }
        UEqnP8.relax();
        fvOptions.constrain(UEqnP8);

        volScalarField rAU8(1.0/UEqnP8.A());
        volVectorField HbyA8(constrainHbyA(rAU8*UEqnP8.H(), Umod, p));
        surfaceScalarField phiHbyA8("stageB8_phiHbyA", fvc::flux(HbyA8));
        adjustPhi(phiHbyA8, Umod, p);

        tmp<volScalarField> rAtU8(rAU8);
        if (simple.consistent())
        {
            rAtU8 = 1.0/(1.0/rAU8 - UEqnP8.H1());
            phiHbyA8 +=
                fvc::interpolate(rAtU8() - rAU8)
               *fvc::snGrad(p)*mesh.magSf();
            HbyA8 -= (rAU8 - rAtU8())*fvc::grad(p);
        }

        constrainPressure(p, Umod, phiHbyA8, rAtU8());

        // P-row pEqn: psi = pmod so flux() responds to dp (the P-P block,
        // incl. the fixedValue outlet internalCoeffs boundary contribution).
        fvScalarMatrix pEqn8
        (
            fvm::laplacian(rAtU8(), pmod) == fvc::div(phiHbyA8)
        );
        pEqn8.setReference(pRefCell, pRefValue);

        const surfaceScalarField phi8(phiHbyA8 - pEqn8.flux());

        resP = scalar(0.0);
        for (label fi = 0; fi < s8NIntF; ++fi)
        {
            resP[mesh.owner()[fi]] += phi8[fi];
            resP[mesh.neighbour()[fi]] -= phi8[fi];
        }
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const scalarField& pb = phi8.boundaryField()[patchi];
            forAll(fc, f) { resP[fc[f]] += pb[f]; }
        }
    };

    // ---- (4) central FD of the actual residual pair along v ----
    auto s8CentralFD = [&](const scalarField& v, const scalar eps,
                           scalarField& fd)
    {
        volVectorField UmodP(U);
        volScalarField pmodP(p);
        for (label c = 0; c < s8N; ++c)
        {
            for (direction d = 0; d < 3; ++d)
            {
                UmodP.primitiveFieldRef()[c].component(d) +=
                    eps*v[discreteUIndex(c, d)];
            }
            pmodP.primitiveFieldRef()[c] += eps*v[discretePIndex(c)];
        }
        UmodP.correctBoundaryConditions();
        pmodP.correctBoundaryConditions();
        vectorField resUP(s8N, vector::zero);
        scalarField resPP(s8N, 0.0);
        s8Rebuild(UmodP, pmodP, resUP, resPP);

        volVectorField UmodM(U);
        volScalarField pmodM(p);
        for (label c = 0; c < s8N; ++c)
        {
            for (direction d = 0; d < 3; ++d)
            {
                UmodM.primitiveFieldRef()[c].component(d) -=
                    eps*v[discreteUIndex(c, d)];
            }
            pmodM.primitiveFieldRef()[c] -= eps*v[discretePIndex(c)];
        }
        UmodM.correctBoundaryConditions();
        pmodM.correctBoundaryConditions();
        vectorField resUM(s8N, vector::zero);
        scalarField resPM(s8N, 0.0);
        s8Rebuild(UmodM, pmodM, resUM, resPM);

        for (label c = 0; c < s8N; ++c)
        {
            for (direction d = 0; d < 3; ++d)
            {
                fd[discreteUIndex(c, d)] =
                    (resUP[c].component(d) - resUM[c].component(d))
                   /(2.0*eps);
            }
            fd[discretePIndex(c)] = (resPP[c] - resPM[c])/(2.0*eps);
        }
    };

    // ---- (5) block metrics helper (P blocks are cell-indexed here) ----
    auto s8BlockMetrics = [&](const scalarField& a, const scalarField& b,
                              const labelList& idx, const word& name,
                              const scalar eps, const bool isU)
    {
        scalar nA = 0.0, nB = 0.0, nD = 0.0, dot = 0.0;
        scalar mx = 0.0;
        label arg = -1;
        forAll(idx, k)
        {
            const label i = idx[k];
            const label slot = isU ? i : discretePIndex(i);
            const scalar av = a[slot];
            const scalar bv = b[slot];
            const scalar d = av - bv;
            nA += av*av; nB += bv*bv; nD += d*d; dot += av*bv;
            if (mag(d) > mx) { mx = mag(d); arg = i; }
        }
        reduce(nA, sumOp<scalar>());
        reduce(nB, sumOp<scalar>());
        reduce(nD, sumOp<scalar>());
        reduce(dot, sumOp<scalar>());
        reduce(mx, maxOp<scalar>());
        reduce(arg, maxOp<label>());
        // NOTE: do NOT clamp the denominators with Foam::SMALL (1e-15): the
        // P-block squared norms here are O(1e-18) (FD values O(1e-9)), far
        // below SMALL, so max(nB,SMALL) would corrupt relL2/cos/normRatio.
        // Use exact denominators (VSMALL only as a divide-by-zero guard).
        const scalar rel = Foam::sqrt(nD/Foam::max(nB, VSMALL));
        const scalar cosv =
            dot/Foam::max(Foam::sqrt(nA*nB), VSMALL);
        const scalar nr = Foam::sqrt(nA/Foam::max(nB, VSMALL));
        Info<< "  stageB8 " << name << " eps=" << eps
            << ": relL2=" << rel
            << " cos=" << cosv
            << " normRatio=" << nr
            << " |Jv|=" << Foam::sqrt(nA)
            << " |FD|=" << Foam::sqrt(nB)
            << " maxDiff=" << mx << " @"
            << (isU ? word("Ucell") : word("Pcell")) << "=" << arg << endl;
        return rel;
    };

    // ---- (6) directions ----
    // dir1 (P1): dp-only null-vector P block
    scalarField v1(s8NUnk, 0.0);
    for (label c = 0; c < s8N; ++c) { v1[discretePIndex(c)] = s8nP[c]; }
    // dir2 (P3): dp-heavy, outlet-supported (deterministic sin on outlet
    // cells; small seeded smooth dU component)
    scalarField v2(s8NUnk, 0.0);
    {
        label k = 0;
        forAll(mesh.boundary(), patchi)
        {
            if (!p.boundaryField()[patchi].fixesValue()) continue;
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(fc, f)
            {
                const label c = fc[f];
                if (s8OutletCell[c])
                {
                    v2[discretePIndex(c)] = Foam::sin(0.7*(k + 1));
                    ++k;
                }
            }
        }
        for (label c = 0; c < s8N; ++c)
        {
            const scalar s = Foam::sin(0.271*(c + 1));
            v2[discreteUIndex(c, 0)] = 0.05*s;
            v2[discreteUIndex(c, 1)] = 0.025*s;
            v2[discreteUIndex(c, 2)] = -0.015*s;
        }
    }
    // dir3 (P3): seeded smooth dU + dp
    scalarField v3(s8NUnk, 0.0);
    for (label c = 0; c < s8N; ++c)
    {
        const scalar s = Foam::sin(0.271*(c + 1));
        v3[discreteUIndex(c, 0)] = s;
        v3[discreteUIndex(c, 1)] = 0.5*s;
        v3[discreteUIndex(c, 2)] = -0.3*s;
        v3[discretePIndex(c)] = Foam::sin(0.173*(c + 1));
    }

    // ---- (7) Jv (corrected matrix-free forward J) ----
    scalarField Jv1(s8NUnk, 0.0), Jv2(s8NUnk, 0.0), Jv3(s8NUnk, 0.0);
    applyDiscreteFlowJ(v1, Jv1);
    applyDiscreteFlowJ(v2, Jv2);
    applyDiscreteFlowJ(v3, Jv3);

    // ---- (8) eps sweep + P1/P3 reporting ----
    const scalar epsList[4] = {1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5};
    OFstream mOut(s8ArtDir + "/stageB8_metrics.txt");
    mOut.precision(16);

    scalarField fd1(s8NUnk, 0.0), fd2(s8NUnk, 0.0), fd3(s8NUnk, 0.0);

    for (label ei = 0; ei < 4; ++ei)
    {
        const scalar eps = epsList[ei];
        Info<< "  stageB8 eps=" << eps << ":" << endl;

        // ---- P1 (dir1): dp-only null vector ----
        s8CentralFD(v1, eps, fd1);
        {
            // full-P / outlet-adjacent / off-outlet
            labelList idxPoff;
            for (label c = 0; c < s8N; ++c)
            {
                if (!s8OutletCell[c]) idxPoff.append(c);
            }
            const scalar rFull =
                s8BlockMetrics(Jv1, fd1, s8IdxPtot, "P1 full-P", eps, false);
            const scalar rOut =
                s8BlockMetrics(Jv1, fd1, s8IdxPout, "P1 outlet-adjacent", eps, false);
            const scalar rOff =
                s8BlockMetrics(Jv1, fd1, idxPoff, "P1 off-outlet", eps, false);
            // sign/structure: dominant-sign agreement + outlet mass fraction
            scalar sPos = 0.0, sNeg = 0.0, fdOutMass = 0.0, fdTotMass = 0.0;
            scalar jvOutMass = 0.0, jvTotMass = 0.0;
            for (label c = 0; c < s8N; ++c)
            {
                const scalar fdv = fd1[discretePIndex(c)];
                const scalar jvv = Jv1[discretePIndex(c)];
                if (fdv*jvv > 0) sPos += 1.0; else sNeg += 1.0;
                fdTotMass += mag(fdv);
                jvTotMass += mag(jvv);
                if (s8OutletCell[c])
                {
                    fdOutMass += mag(fdv);
                    jvOutMass += mag(jvv);
                }
            }
            reduce(sPos, sumOp<scalar>());
            reduce(sNeg, sumOp<scalar>());
            reduce(fdOutMass, sumOp<scalar>());
            reduce(fdTotMass, sumOp<scalar>());
            reduce(jvOutMass, sumOp<scalar>());
            reduce(jvTotMass, sumOp<scalar>());
            Info<< "  stageB8 P1 sign: sameSignFraction="
                << sPos/Foam::max(sPos + sNeg, SMALL)
                << " (pos=" << sPos << " neg=" << sNeg << ")"
                << "  |FD| mass outlet=" << fdOutMass/Foam::max(fdTotMass, SMALL)
                << "  |Jv| mass outlet=" << jvOutMass/Foam::max(jvTotMass, SMALL)
                << endl;
            mOut<< "P1 " << eps << ' ' << rFull << ' ' << rOut << ' ' << rOff
                << ' ' << sPos/Foam::max(sPos + sNeg, SMALL)
                << ' ' << fdOutMass/Foam::max(fdTotMass, SMALL)
                << ' ' << jvOutMass/Foam::max(jvTotMass, SMALL) << '\n';
        }

        // ---- P3 dir1 (n_P): full block decomposition ----
        {
            const scalar rU = s8BlockMetrics(Jv1, fd1, s8IdxMom, "P3d1 U rows", eps, true);
            const scalar rPi = s8BlockMetrics(Jv1, fd1, s8IdxPint, "P3d1 P-internal", eps, false);
            const scalar rPo = s8BlockMetrics(Jv1, fd1, s8IdxPout, "P3d1 P-outlet-adjacent", eps, false);
            const scalar rPb = s8BlockMetrics(Jv1, fd1, s8IdxPother, "P3d1 P-other-boundary", eps, false);
            const scalar rPt = s8BlockMetrics(Jv1, fd1, s8IdxPtot, "P3d1 P-total", eps, false);
            mOut<< "P3d1 " << eps << ' ' << rU << ' ' << rPi << ' ' << rPo
                << ' ' << rPb << ' ' << rPt << '\n';
        }

        // ---- P3 dir2 (dp-heavy outlet-supported) ----
        s8CentralFD(v2, eps, fd2);
        {
            const scalar rU = s8BlockMetrics(Jv2, fd2, s8IdxMom, "P3d2 U rows", eps, true);
            const scalar rPi = s8BlockMetrics(Jv2, fd2, s8IdxPint, "P3d2 P-internal", eps, false);
            const scalar rPo = s8BlockMetrics(Jv2, fd2, s8IdxPout, "P3d2 P-outlet-adjacent", eps, false);
            const scalar rPb = s8BlockMetrics(Jv2, fd2, s8IdxPother, "P3d2 P-other-boundary", eps, false);
            const scalar rPt = s8BlockMetrics(Jv2, fd2, s8IdxPtot, "P3d2 P-total", eps, false);
            mOut<< "P3d2 " << eps << ' ' << rU << ' ' << rPi << ' ' << rPo
                << ' ' << rPb << ' ' << rPt << '\n';
        }

        // ---- P3 dir3 (seeded smooth dU+dp) ----
        s8CentralFD(v3, eps, fd3);
        {
            const scalar rU = s8BlockMetrics(Jv3, fd3, s8IdxMom, "P3d3 U rows", eps, true);
            const scalar rPi = s8BlockMetrics(Jv3, fd3, s8IdxPint, "P3d3 P-internal", eps, false);
            const scalar rPo = s8BlockMetrics(Jv3, fd3, s8IdxPout, "P3d3 P-outlet-adjacent", eps, false);
            const scalar rPb = s8BlockMetrics(Jv3, fd3, s8IdxPother, "P3d3 P-other-boundary", eps, false);
            const scalar rPt = s8BlockMetrics(Jv3, fd3, s8IdxPtot, "P3d3 P-total", eps, false);
            mOut<< "P3d3 " << eps << ' ' << rU << ' ' << rPi << ' ' << rPo
                << ' ' << rPb << ' ' << rPt << '\n';
        }

        // ---- per-eps FD artifacts ----
        const word epsName =
            ei == 0 ? "1e-2" : (ei == 1 ? "1e-3" : (ei == 2 ? "1e-4" : "1e-5"));
        OFstream fdO1(s8ArtDir + "/stageB8_FD1_eps_" + epsName + ".mtx");
        OFstream fdO2(s8ArtDir + "/stageB8_FD2_eps_" + epsName + ".mtx");
        OFstream fdO3(s8ArtDir + "/stageB8_FD3_eps_" + epsName + ".mtx");
        fdO1.precision(16); fdO2.precision(16); fdO3.precision(16);
        for (label i = 0; i < s8NUnk; ++i)
        {
            fdO1<< fd1[i] << '\n';
            fdO2<< fd2[i] << '\n';
            fdO3<< fd3[i] << '\n';
        }
    }

    // ---- (9) static exports ----
    {
        OFstream d1(s8ArtDir + "/stageB8_dir1.mtx");
        OFstream d2(s8ArtDir + "/stageB8_dir2.mtx");
        OFstream d3(s8ArtDir + "/stageB8_dir3.mtx");
        OFstream j1(s8ArtDir + "/stageB8_Jv1.mtx");
        OFstream j2(s8ArtDir + "/stageB8_Jv2.mtx");
        OFstream j3(s8ArtDir + "/stageB8_Jv3.mtx");
        OFstream nP(s8ArtDir + "/stageB8_nP.mtx");
        OFstream ct(s8ArtDir + "/stageB8_celltype.mtx");
        OFstream oc(s8ArtDir + "/stageB8_outlet_cells.mtx");
        d1.precision(16); d2.precision(16); d3.precision(16);
        j1.precision(16); j2.precision(16); j3.precision(16);
        nP.precision(16);
        for (label i = 0; i < s8NUnk; ++i)
        {
            d1<< v1[i] << '\n';
            d2<< v2[i] << '\n';
            d3<< v3[i] << '\n';
            j1<< Jv1[i] << '\n';
            j2<< Jv2[i] << '\n';
            j3<< Jv3[i] << '\n';
        }
        for (label c = 0; c < s8N; ++c)
        {
            nP<< s8nP[c] << '\n';
            ct<< (s8OutletCell[c] ? 2 : (s8HasBnd[c] ? 1 : 0)) << '\n';
            if (s8OutletCell[c]) oc<< c << '\n';
        }
    }

    Info<< "--- stageB8JPPActualResidualFD: complete (artifacts in "
        << s8ArtDir << ") ---" << endl;
}
