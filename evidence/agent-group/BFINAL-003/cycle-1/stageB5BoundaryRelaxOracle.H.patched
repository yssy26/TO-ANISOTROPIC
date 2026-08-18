// =============================================================================
// stageB5BoundaryRelaxOracle.H  (BFINAL-002, DIAGNOSTIC_ONLY)
// -----------------------------------------------------------------------------
// Guarded diagnostic probe (optProperties switch "stageB5BoundaryRelaxOracle",
// default false).  It must be #included ONLY inside src/solveDiscreteFlowAdjoint.H
// (which is itself only compiled when stageB4JacobianProbe=true), at a point
// where the reduced-J state (applyDiscreteFlowJ, rAUAdj, discretePrimalMomentum,
// discreteOwner/Neighbour/Weights/Sf/Phi, discreteUIndex/PIndex,
// primalPressureMobility, U, p, phi, alpha, nuEffFrozen, mesh, runTime,
// optProperties, simple, pRefCell, pRefValue, fvOptions, ransFlowModel) is in
// scope and the flow state is the CONVERGED production state.
//
// No-op when disabled: the caller guards this block with
//   optProperties.lookupOrDefault<Switch>("stageB5BoundaryRelaxOracle", false)
// so when the switch is absent/false nothing here executes, no field is
// created/modified and no file is written.  This probe never modifies U, p,
// phi, alpha, k, omega, nutFrozen, nuEffFrozen, primalPressureMobility,
// rAUAdj or any production value/residual.
//
// Gate A (Q-A):  J*v (applyDiscreteFlowJ, dp=0) vs central FD
//   [R(U+eps*dU)-R(U-eps*dU)]/(2eps) of the FULL production reduced residual,
//   with phi rebuilt through the EXACT NS.H SIMPLE sequence
//   (fresh UEqn -> relax -> rAU=1/A -> constrainHbyA -> fvc::flux -> adjustPhi
//    -> rAtU(consistent=false) -> constrainPressure -> pEqn=laplacian(rAtU,p)
//    -> setReference(pRefCell) -> phi = phiHbyA - pEqn.flux(), boundary flux
//    REBUILT, dp = 0).  Reported per block: momentum / P-internal /
//   P-boundary / P-total with |Jv|, |FD|, |err|, relL2, cosine,
//   max-abs-diff+location, epsilon table.  Also reproduces the BFINAL-001
//   frozen-phi momentum anchor (relL2 ~ 1.65e-4).
//
// Gate B (Q-B): isolated dphi/dU oracle.  p/design/alpha/k/omega/nutFrozen/
//   nuEffFrozen fixed, deterministic dU, phi(U+/-eps dU) via the production
//   path.  Candidate A = current J dphi tangent (unrelaxed rAUAdj x unrelaxed
//   frozen-coefficient dH, internal faces from collectedDPhiJ, boundary faces
//   from the J boundary term).  Candidate B = production relaxed mobility x
//   production relaxed-H tangent (rAU_rel x dH_rel with the relax source and
//   dev2, constrainHbyA boundary pinning).  True FD = rebuild central diff.
//   Per candidate per eps: |analytic|, |FD|, norm ratio, relL2, cosine,
//   internal-face rel err, boundary-face rel err.  Plus relax factor and
//   rAUAdj vs relaxed-rAU avg/min/max/norm-ratio.
// =============================================================================
{
    Info<< "--- StageB5: boundary-consistent P-row + dphi/dU relaxation "
        << "oracle (BFINAL-002, DIAGNOSTIC_ONLY) ---" << endl;

    // ------------------------------------------------------------------
    // 0. Deterministic direction dU (identical construction to the
    //    BFINAL-001 stageB2_dUdirPhys.mtx, sin(0.271*(celli+1)) pattern).
    // ------------------------------------------------------------------
    volVectorField dUdir
    (
        IOobject
        (
            "stageB5_dUdir",
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

    // state vector v = [dU; dp=0] (4N unknowns)
    scalarField stateJ(discreteNUnknowns, 0.0);
    forAll(U, celli)
    {
        stateJ[discreteUIndex(celli, 0)] = dUdir[celli].x();
        stateJ[discreteUIndex(celli, 1)] = dUdir[celli].y();
        stateJ[discreteUIndex(celli, 2)] = dUdir[celli].z();
    }

    // ------------------------------------------------------------------
    // 1. J*v via the existing matrix-free reduced Jacobian (dp = 0).
    // ------------------------------------------------------------------
    scalarField outJ(discreteNUnknowns, 0.0);
    collectDPhiJ = true;
    applyDiscreteFlowJ(stateJ, outJ);
    collectDPhiJ = false;
    // collectedDPhiJ (internal faces) is now J's dphi tangent (Candidate A
    // internal faces).  outJ[P] is J's P-row = div(dphi_J) + boundary dphi.
    // P-only view of the J output (P components live at 3N + cell).
    scalarField outJP(discreteNCells, 0.0);
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        outJP[celli] = outJ[discretePIndex(celli)];
    }

    // Recompute J's dH (copy of applyDiscreteFlowJ L502-582, frozen
    // unrelaxed coefficients) so the per-boundary-face Candidate-A dphi can
    // be extracted; verified below against outJ[P] by self-consistency.
    scalarField dHxJ(discreteNCells, 0.0);
    scalarField dHyJ(discreteNCells, 0.0);
    scalarField dHzJ(discreteNCells, 0.0);
    forAll(discreteUpper, facei)
    {
        const label own = discreteOwner[facei];
        const label nei = discreteNeighbour[facei];
        dHxJ[own] -= discreteUpper[facei]*stateJ[discreteUIndex(nei, 0)];
        dHyJ[own] -= discreteUpper[facei]*stateJ[discreteUIndex(nei, 1)];
        dHzJ[own] -= discreteUpper[facei]*stateJ[discreteUIndex(nei, 2)];
        dHxJ[nei] -= discreteLower[facei]*stateJ[discreteUIndex(own, 0)];
        dHyJ[nei] -= discreteLower[facei]*stateJ[discreteUIndex(own, 1)];
        dHzJ[nei] -= discreteLower[facei]*stateJ[discreteUIndex(own, 2)];
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const vectorField icV(discretePrimalMomentum.internalCoeffs()[patchi]);
        forAll(fc, f)
        {
            const label celli = fc[f];
            const scalar cav = (icV[f].x() + icV[f].y() + icV[f].z())/3.0;
            dHxJ[celli] +=
                (-icV[f].x() + cav)*stateJ[discreteUIndex(celli, 0)];
            dHyJ[celli] +=
                (-icV[f].y() + cav)*stateJ[discreteUIndex(celli, 1)];
            dHzJ[celli] +=
                (-icV[f].z() + cav)*stateJ[discreteUIndex(celli, 2)];
        }
    }
    for (label celli = 0; celli < discreteNCells; ++celli)
    {
        const scalar Vinv = 1.0/Foam::max(mesh.V()[celli], SMALL);
        dHxJ[celli] *= Vinv;
        dHyJ[celli] *= Vinv;
        dHzJ[celli] *= Vinv;
    }

    // Candidate-A boundary dphi per boundary face: J's corrected boundary
    // term (BFINAL-003 P<-U relaxation patch, applyDiscreteFlowJ boundary
    // loop): exactly 0 on non-assignable U patches (production constrainHbyA
    // pins HbyA_b = U_b) and Sf_b & (alphaRel*rAUAdj*dH + (1-alphaRel)*dU_cell)
    // on assignable (zeroGradient) U patches.
    // Per-boundary-face arrays are indexed by (internal faces, then boundary
    // faces in patch order).  Also accumulate per-cell boundary dphi for the
    // self-consistency check against outJ[P].
    const label nIntF = mesh.nInternalFaces();
    const label nBndF = mesh.nFaces() - nIntF;
    const label nFaceTot = nIntF + nBndF;
    scalarField dphiA_faces(nFaceTot, 0.0);
    scalarField dphiA_cell(discreteNCells, 0.0);
    // Candidate-A internal-face dphi = J's collectedDPhiJ (from the
    // applyDiscreteFlowJ call above; now the corrected relaxed P-U tangent).
    forAll(collectedDPhiJ, fi)
    {
        dphiA_faces[fi] = collectedDPhiJ[fi];
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
                const vector dHb(dHxJ[celli], dHyJ[celli], dHzJ[celli]);
                const scalar db =
                    uAssignable
                  ? (pSf[f]
                     & (alphaRel*rAUAdj[celli]*dHb
                        + (1.0 - alphaRel)*dUdir[celli]))
                  : scalar(0.0);
                dphiA_faces[bf] = db;
                dphiA_cell[celli] += db;
                ++bf;
            }
        }
    }
    // self-consistency: outJ[P(c)] == internal-face div(collectedDPhiJ) +
    // dphiA_cell(c)
    {
        scalarField Pcheck(discreteNCells, 0.0);
        forAll(discreteOwner, fi)
        {
            Pcheck[discreteOwner[fi]] += collectedDPhiJ[fi];
            Pcheck[discreteNeighbour[fi]] -= collectedDPhiJ[fi];
        }
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            Pcheck[celli] += dphiA_cell[celli];
        }
        scalar nSC = 0.0, nSP = 0.0;
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            const scalar d = Pcheck[celli] - outJ[discretePIndex(celli)];
            nSC += sqr(d);
            nSP += sqr(outJ[discretePIndex(celli)]);
        }
        reduce(nSC, sumOp<scalar>());
        reduce(nSP, sumOp<scalar>());
        Info<< "StageB5 selfcheck (recomputed dH boundary + internal div vs "
            << "outJ[P]): relL2="
            << Foam::sqrt(nSC/Foam::max(nSP, SMALL)) << endl;
    }

    // ------------------------------------------------------------------
    // 2. Boundary-cell classification + patch of first boundary face.
    // ------------------------------------------------------------------
    boolList cellHasBnd(discreteNCells, false);
    labelList cellBndPatch(discreteNCells, -1);
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        forAll(fc, f)
        {
            cellHasBnd[fc[f]] = true;
            if (cellBndPatch[fc[f]] < 0) cellBndPatch[fc[f]] = patchi;
        }
    }

    // ------------------------------------------------------------------
    // 3. Production phi-rebuild lambda (EXACT NS.H SIMPLE sequence).
    //    Writes: phiOut (rebuilt phi), resUOut (production momentum
    //    residual, NS.H sign: resU = -Eq.residual()), rAUrelOut (relaxed
    //    primal rAU), HrelOut (relaxed production H incl. relax source).
    // ------------------------------------------------------------------
    auto rebuildProduction =
    [&](volVectorField& Umod,
        surfaceScalarField& phiOut,
        vectorField& resUOut,
        scalarField& rAUrelOut,
        vectorField& HrelOut)
    {
        tmp<fvVectorMatrix> tUEqn
        (
            fvm::div(phi, Umod)
          - fvm::laplacian(nuEffFrozen, Umod)
          + fvm::Sp(alpha, Umod)
         ==
           -fvc::grad(p)
          + fvOptions(Umod)
        );
        fvVectorMatrix& UEqn = tUEqn.ref();

        if (ransFlowModel)
        {
            UEqn -= fvc::div
            (
                nuEffFrozen*dev2(Foam::T(fvc::grad(Umod)))
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
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), Umod, p));
        surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
        adjustPhi(phiHbyA, Umod, p);

        tmp<volScalarField> rAtU(rAU);
        if (simple.consistent())
        {
            rAtU = 1.0/(1.0/rAU - UEqn.H1());
            phiHbyA +=
                fvc::interpolate(rAtU() - rAU)
               *fvc::snGrad(p)*mesh.magSf();
            HbyA -= (rAU - rAtU())*fvc::grad(p);
        }

        constrainPressure(p, Umod, phiHbyA, rAtU());

        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );
        pEqn.setReference(pRefCell, pRefValue);

        phiOut = phiHbyA - pEqn.flux();

        forAll(rAUrelOut, celli) { rAUrelOut[celli] = rAU[celli]; }
        const vectorField Hrel(UEqn.H());
        forAll(HrelOut, celli) { HrelOut[celli] = Hrel[celli]; }
    };

    // ------------------------------------------------------------------
    // 4. Metric helper (per-block L2 / relL2 / cosine / max-abs-diff).
    //    relL2 convention identical to BFINAL-001 crosscheck_bf2.py:
    //    relL2 = |a-b|_L2 / |b|_L2 (b = reference/FD).
    // ------------------------------------------------------------------
    auto blockMetrics =
    [&](const scalarField& a, const scalarField& b,
        const labelList& idx, const word& labelName,
        const bool wantLoc, const labelList& locCell,
        const wordList& patchNames)
    {
        scalar nA = 0.0, nB = 0.0, nD = 0.0, dot = 0.0;
        scalar mx = 0.0;
        label arg = -1;
        forAll(idx, k)
        {
            const label i = idx[k];
            const scalar av = a[i];
            const scalar bv = b[i];
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
        const scalar rel = Foam::sqrt(nD/Foam::max(nB, SMALL));
        const scalar cosv = dot/Foam::max(Foam::sqrt(nA*nB), SMALL);
        Info<< "StageB5 " << labelName << ": |a|L2=" << Foam::sqrt(nA)
            << " |b(FD)|L2=" << Foam::sqrt(nB)
            << " |err|L2=" << Foam::sqrt(nD)
            << " relL2=" << rel
            << " cos=" << cosv
            << " maxDiff=" << mx << " @idx=" << arg;
        if (wantLoc && arg >= 0)
        {
            // For P blocks `a` is cell-indexed (size discreteNCells); for the
            // momentum block `a` is U-slot-indexed (size 4N or 3N).
            if (a.size() == discreteNCells)
            {
                const label c = arg;
                if (locCell.size() > 0 && locCell[c] >= 0)
                {
                    Info<< " (P cell " << c << " patch "
                        << patchNames[locCell[c]] << ")";
                }
                else
                {
                    Info<< " (P cell " << c << ")";
                }
            }
            else
            {
                Info<< " (U cell " << arg/3 << " cmpt " << arg%3 << ")";
            }
        }
        Info<< endl;
        return rel;
    };

    // indices for the four Gate-A blocks
    labelList idxMom(discreteNVelocity);
    for (label i = 0; i < discreteNVelocity; ++i) idxMom[i] = i;
    labelList idxPint, idxPbnd, idxPtot(discreteNCells);
    for (label c = 0; c < discreteNCells; ++c)
    {
        idxPtot[c] = c;
        if (cellHasBnd[c]) idxPbnd.append(c);
        else               idxPint.append(c);
    }

    wordList patchNames(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        patchNames[patchi] = mesh.boundary()[patchi].name();
    }

    // ------------------------------------------------------------------
    // 5. Anchor: frozen-phi momentum FD (exact BFINAL-001 NS.H recipe,
    //    hU = 1e-3*max|U|), J*v momentum rows vs FD.
    // ------------------------------------------------------------------
    {
        const scalar hU = 1e-3*Foam::max(gMax(mag(U.primitiveField())), SMALL);
        vectorField ruP(discreteNCells, vector::zero);
        vectorField ruM(discreteNCells, vector::zero);
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
            if (sgn == 0)
            {
                forAll(Umod, celli) { ruP[celli] = -tRes()[celli]; }
            }
            else
            {
                forAll(Umod, celli) { ruM[celli] = -tRes()[celli]; }
            }
        }
        scalarField fdU(discreteNVelocity, 0.0);
        scalarField jvU(discreteNVelocity, 0.0);
        forAll(U, celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                fdU[discreteUIndex(celli, c)] =
                    (ruP[celli].component(c) - ruM[celli].component(c))/(2.0*hU);
                jvU[discreteUIndex(celli, c)] =
                    outJ[discreteUIndex(celli, c)];
            }
        }
        Info<< "StageB5 anchor (frozen-phi momentum FD, hU=" << hU << "):";
        blockMetrics(jvU, fdU, idxMom, "momentum", false, labelList(), wordList());
    }

    // ------------------------------------------------------------------
    // 6. Gate A + Gate B epsilon sweeps (shared rebuilds).
    // ------------------------------------------------------------------
    const scalar gateAEps[5] = {1e-3, 3e-4, 1e-4, 3e-5, 1e-5};
    // per-eps FD state vectors (exported for the Post-Reviewer)
    // index: eps table row; U part 3N, P part N
    // dphi FD per face (internal+boundary); dphi candidate B per face
    // (candidate A is eps-independent)
    OFstream fdAll("stageB5_gateA_FD_all_eps.mtx");
    OFstream dphiFDA("stageB5_gateB_dphi_FD_all_eps.mtx");
    OFstream dphiBA("stageB5_gateB_dphi_B_all_eps.mtx");
    fdAll.precision(16); dphiFDA.precision(16); dphiBA.precision(16);

    for (label ei = 0; ei < 5; ++ei)
    {
        const scalar eps = gateAEps[ei];

        surfaceScalarField phiP("stageB5_phiP", phi);
        surfaceScalarField phiM("stageB5_phiM", phi);
        vectorField resUP(discreteNCells, vector::zero);
        vectorField resUM(discreteNCells, vector::zero);
        scalarField rAUP(discreteNCells, 0.0);
        scalarField rAUM(discreteNCells, 0.0);
        vectorField HP(discreteNCells, vector::zero);
        vectorField HM(discreteNCells, vector::zero);

        volVectorField UmodP(U);
        forAll(UmodP.primitiveFieldRef(), celli)
        {
            UmodP.primitiveFieldRef()[celli] += eps*dUdir[celli];
        }
        UmodP.correctBoundaryConditions();
        rebuildProduction(UmodP, phiP, resUP, rAUP, HP);

        volVectorField UmodM(U);
        forAll(UmodM.primitiveFieldRef(), celli)
        {
            UmodM.primitiveFieldRef()[celli] -= eps*dUdir[celli];
        }
        UmodM.correctBoundaryConditions();
        rebuildProduction(UmodM, phiM, resUM, rAUM, HM);

        // ---- Gate A FD vectors (this eps) ----
        scalarField fdU(discreteNVelocity, 0.0);
        scalarField fdP(discreteNCells, 0.0);
        forAll(U, celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                fdU[discreteUIndex(celli, c)] =
                    (resUP[celli].component(c) - resUM[celli].component(c))
                   /(2.0*eps);
            }
        }
        // P residual = div(rebuilt phi), boundary flux REBUILT
        {
            scalarField rpP(discreteNCells, 0.0);
            scalarField rpM(discreteNCells, 0.0);
            forAll(mesh.owner(), fi)
            {
                rpP[mesh.owner()[fi]] += phiP[fi];
                rpP[mesh.neighbour()[fi]] -= phiP[fi];
                rpM[mesh.owner()[fi]] += phiM[fi];
                rpM[mesh.neighbour()[fi]] -= phiM[fi];
            }
            forAll(phiP.boundaryField(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const scalarField& pB = phiP.boundaryField()[patchi];
                const scalarField& mB = phiM.boundaryField()[patchi];
                forAll(fc, f)
                {
                    rpP[fc[f]] += pB[f];
                    rpM[fc[f]] += mB[f];
                }
            }
            forAll(fdP, celli)
            {
                fdP[celli] = (rpP[celli] - rpM[celli])/(2.0*eps);
            }
        }

        // ---- Gate A block metrics ----
        Info<< "StageB5 GateA eps=" << eps << ":";
        blockMetrics
        (
            outJ, fdU, idxMom, " momentum", false, labelList(), wordList()
        );
        blockMetrics
        (
            outJP, fdP, idxPint, " P-internal", false, labelList(), wordList()
        );
        blockMetrics
        (
            outJP, fdP, idxPbnd, " P-boundary", true, cellBndPatch, patchNames
        );
        blockMetrics
        (
            outJP, fdP, idxPtot, " P-total", true, cellBndPatch, patchNames
        );

        // ---- Gate B dphi vectors (this eps) ----
        // FD dphi per face (internal then boundary)
        scalarField dphiFD(nFaceTot, 0.0);
        forAll(phiP, fi)  // internal faces
        {
            dphiFD[fi] = (phiP[fi] - phiM[fi])/(2.0*eps);
        }
        {
            label bf = nIntF;
            forAll(phiP.boundaryField(), patchi)
            {
                const scalarField& pB = phiP.boundaryField()[patchi];
                const scalarField& mB = phiM.boundaryField()[patchi];
                forAll(pB, f)
                {
                    dphiFD[bf] = (pB[f] - mB[f])/(2.0*eps);
                    ++bf;
                }
            }
        }

        // Candidate B: relaxed mobility x relaxed-H tangent (production)
        //   dH_rel = (HP - HM)/(2eps)  (affine in U => exact at any eps)
        //   dphiB_f = Sf & (wf*rAUrel_o*dHrel_o + (1-wf)*rAUrel_n*dHrel_n)
        //   boundary: 0 at !assignable (constrainHbyA pinning), extrapolated
        //   at assignable U patches (outlet).
        scalarField dHxR(discreteNCells, 0.0);
        scalarField dHyR(discreteNCells, 0.0);
        scalarField dHzR(discreteNCells, 0.0);
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            dHxR[celli] = (HP[celli].x() - HM[celli].x())/(2.0*eps);
            dHyR[celli] = (HP[celli].y() - HM[celli].y())/(2.0*eps);
            dHzR[celli] = (HP[celli].z() - HM[celli].z())/(2.0*eps);
        }
        // rAU_rel = relaxed production rAU (use the plus-side; A is
        // U-independent so rAUP == rAUM to roundoff)
        const scalarField& rAUrel = rAUP;

        scalarField dphiB(nFaceTot, 0.0);
        forAll(discreteOwner, fi)
        {
            const label own = discreteOwner[fi];
            const label nei = discreteNeighbour[fi];
            const scalar wf = discreteWeights[fi];
            const vector dHo(dHxR[own], dHyR[own], dHzR[own]);
            const vector dHn(dHxR[nei], dHyR[nei], dHzR[nei]);
            dphiB[fi] =
                discreteSf[fi]
              & (wf*rAUrel[own]*dHo + (1.0 - wf)*rAUrel[nei]*dHn);
        }
        {
            label bf = nIntF;
            forAll(mesh.boundary(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
                const bool uAssignable =
                    U.boundaryField()[patchi].assignable();
                forAll(fc, f)
                {
                    const label celli = fc[f];
                    if (!uAssignable)
                    {
                        dphiB[bf] = 0.0;   // constrainHbyA pinning
                    }
                    else
                    {
                        const vector dHc(dHxR[celli], dHyR[celli], dHzR[celli]);
                        dphiB[bf] = (pSf[f] & (rAUrel[celli]*dHc));
                    }
                    ++bf;
                }
            }
        }

        // ---- Gate B metrics ----
        // internal / boundary / all face index subsets
        labelList idxIntF(nIntF), idxBndF(nBndF), idxAllF(nFaceTot);
        for (label fi = 0; fi < nIntF; ++fi) idxIntF[fi] = fi;
        for (label b = 0; b < nBndF; ++b) idxBndF[b] = nIntF + b;
        for (label fi = 0; fi < nFaceTot; ++fi) idxAllF[fi] = fi;

        Info<< "StageB5 GateB eps=" << eps
            << " Candidate A (J unrelaxed):";
        blockMetrics(dphiA_faces, dphiFD, idxIntF, " A-internal-faces",
                     false, labelList(), wordList());
        blockMetrics(dphiA_faces, dphiFD, idxBndF, " A-boundary-faces",
                     false, labelList(), wordList());
        blockMetrics(dphiA_faces, dphiFD, idxAllF, " A-all-faces",
                     false, labelList(), wordList());
        // per-candidate summary numbers (for the S4 tables)
        {
            scalar nA = 0.0, nB = 0.0, nD = 0.0, dot = 0.0;
            scalar nAI = 0.0, nBI = 0.0, nDI = 0.0, dotI = 0.0;
            scalar nAB = 0.0, nBB = 0.0, nDB = 0.0, dotB = 0.0;
            for (label fi = 0; fi < nIntF; ++fi)
            {
                const scalar av = dphiA_faces[fi];
                const scalar bv = dphiFD[fi];
                nA += av*av; nB += bv*bv; nD += (av-bv)*(av-bv); dot += av*bv;
                nAI += av*av; nBI += bv*bv; nDI += (av-bv)*(av-bv); dotI += av*bv;
            }
            for (label b = 0; b < nBndF; ++b)
            {
                const scalar av = dphiA_faces[nIntF + b];
                const scalar bv = dphiFD[nIntF + b];
                nA += av*av; nB += bv*bv; nD += (av-bv)*(av-bv); dot += av*bv;
                nAB += av*av; nBB += bv*bv; nDB += (av-bv)*(av-bv); dotB += av*bv;
            }
            Info<< "StageB5 GateB eps=" << eps
                << " candA: |A|L2=" << Foam::sqrt(nA)
                << " |FD|L2=" << Foam::sqrt(nB)
                << " normRatio=" << Foam::sqrt(nA/Foam::max(nB, SMALL))
                << " relL2=" << Foam::sqrt(nD/Foam::max(nB, SMALL))
                << " cos=" << dot/Foam::max(Foam::sqrt(nA*nB), SMALL)
                << " intRel=" << Foam::sqrt(nDI/Foam::max(nBI, SMALL))
                << " bndRel=" << Foam::sqrt(nDB/Foam::max(nBB, SMALL))
                << endl;
        }
        Info<< "StageB5 GateB eps=" << eps
            << " Candidate B (production relaxed):";
        blockMetrics(dphiB, dphiFD, idxIntF, " B-internal-faces",
                     false, labelList(), wordList());
        blockMetrics(dphiB, dphiFD, idxBndF, " B-boundary-faces",
                     false, labelList(), wordList());
        blockMetrics(dphiB, dphiFD, idxAllF, " B-all-faces",
                     false, labelList(), wordList());
        {
            scalar nA = 0.0, nB = 0.0, nD = 0.0, dot = 0.0;
            scalar nAI = 0.0, nBI = 0.0, nDI = 0.0, dotI = 0.0;
            scalar nAB = 0.0, nBB = 0.0, nDB = 0.0, dotB = 0.0;
            for (label fi = 0; fi < nIntF; ++fi)
            {
                const scalar av = dphiB[fi];
                const scalar bv = dphiFD[fi];
                nA += av*av; nB += bv*bv; nD += (av-bv)*(av-bv); dot += av*bv;
                nAI += av*av; nBI += bv*bv; nDI += (av-bv)*(av-bv); dotI += av*bv;
            }
            for (label b = 0; b < nBndF; ++b)
            {
                const scalar av = dphiB[nIntF + b];
                const scalar bv = dphiFD[nIntF + b];
                nA += av*av; nB += bv*bv; nD += (av-bv)*(av-bv); dot += av*bv;
                nAB += av*av; nBB += bv*bv; nDB += (av-bv)*(av-bv); dotB += av*bv;
            }
            Info<< "StageB5 GateB eps=" << eps
                << " candB: |B|L2=" << Foam::sqrt(nA)
                << " |FD|L2=" << Foam::sqrt(nB)
                << " normRatio=" << Foam::sqrt(nA/Foam::max(nB, SMALL))
                << " relL2=" << Foam::sqrt(nD/Foam::max(nB, SMALL))
                << " cos=" << dot/Foam::max(Foam::sqrt(nA*nB), SMALL)
                << " intRel=" << Foam::sqrt(nDI/Foam::max(nBI, SMALL))
                << " bndRel=" << Foam::sqrt(nDB/Foam::max(nBB, SMALL))
                << endl;
        }

        // ---- exports (per eps) ----
        forAll(U, celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                fdAll<< fdU[discreteUIndex(celli, c)] << nl;
            }
            fdAll<< fdP[celli] << nl;
        }
        for (label fi = 0; fi < nFaceTot; ++fi)
        {
            dphiFDA<< dphiFD[fi] << nl;
            dphiBA<< dphiB[fi] << nl;
        }
    }
    // OFstream objects are auto-closed on destruction (no close() in OF7).

    // ------------------------------------------------------------------
    // 7. Exports: J*v, dUdir, candidate A dphi, cell/face classification,
    //    rAU stats.
    // ------------------------------------------------------------------
    {
        OFstream jvOut("stageB5_gateA_Jv.mtx");
        jvOut.precision(16);
        for (label i = 0; i < discreteNUnknowns; ++i)
        {
            jvOut<< outJ[i] << nl;
        }
    }
    {
        OFstream duOut("stageB5_dUdir.mtx");
        duOut.precision(16);
        forAll(U, celli)
        {
            for (direction c = 0; c < 3; ++c)
            {
                duOut<< dUdir[celli].component(c) << nl;
            }
        }
    }
    {
        OFstream daOut("stageB5_gateB_dphi_A.mtx");
        daOut.precision(16);
        for (label fi = 0; fi < nFaceTot; ++fi)
        {
            daOut<< dphiA_faces[fi] << nl;
        }
    }
    {
        OFstream ftOut("stageB5_face_type.mtx");
        for (label fi = 0; fi < nIntF; ++fi)
        {
            ftOut<< 0 << nl;
        }
        forAll(mesh.boundary(), patchi)
        {
            for (label f = 0; f < mesh.boundary()[patchi].size(); ++f)
            {
                ftOut<< patchi + 1 << nl;
            }
        }
    }
    {
        OFstream ctOut("stageB5_cell_type.mtx");
        for (label c = 0; c < discreteNCells; ++c)
        {
            ctOut<< (cellHasBnd[c] ? 1 : 0) << nl;
        }
    }

    // ---- rAU statistics (relax factor, rAUAdj vs relaxed rAU) ----
    {
        scalar sA = 0.0, mnA = GREAT, mxA = 0.0;
        scalar sR = 0.0, mnR = GREAT, mxR = 0.0;
        scalar nA2 = 0.0, nR2 = 0.0;
        for (label celli = 0; celli < discreteNCells; ++celli)
        {
            const scalar a = rAUAdj[celli];
            const scalar r = primalPressureMobility[celli];
            sA += a; mnA = Foam::min(mnA, a); mxA = Foam::max(mxA, a);
            sR += r; mnR = Foam::min(mnR, r); mxR = Foam::max(mxR, r);
            nA2 += sqr(a); nR2 += sqr(r);
        }
        reduce(sA, sumOp<scalar>()); reduce(mnA, minOp<scalar>()); reduce(mxA, maxOp<scalar>());
        reduce(sR, sumOp<scalar>()); reduce(mnR, minOp<scalar>()); reduce(mxR, maxOp<scalar>());
        reduce(nA2, sumOp<scalar>()); reduce(nR2, sumOp<scalar>());
        const scalar nC = Foam::max(scalar(discreteNCells), 1.0);
        scalar uRelaxVal = -1.0;
        if (mesh.solutionDict().found("relaxationFactors"))
        {
            const dictionary& rf = mesh.solutionDict().subDict("relaxationFactors");
            if (rf.found("equations"))
            {
                uRelaxVal =
                    rf.subDict("equations").lookupOrDefault<scalar>("U", -1.0);
            }
        }
        Info<< "StageB5 rAU: relaxFactorU(fvSolution)=" << uRelaxVal
            << " ; rAUAdj(unrelaxed)[avg=" << sA/nC
            << " min=" << mnA << " max=" << mxA << "]"
            << " ; rAUrel(primalPressureMobility)[avg=" << sR/nC
            << " min=" << mnR << " max=" << mxR << "]"
            << " ; avgRatio(rAUAdj/rAUrel)="
            << (sA/Foam::max(sR, SMALL))
            << " ; normRatio(|rAUAdj|/|rAUrel|)="
            << Foam::sqrt(nA2/Foam::max(nR2, SMALL))
            << endl;
    }

    Info<< "--- StageB5 oracle complete (artifacts: stageB5_*.mtx) ---" << endl;
}
