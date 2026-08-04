//-----------------------------------------------------------------------------
// validateCommon.H
// Shared helper lambdas for the Gate 4/6 validation modules.
//
// Extracted from validateFrozenGradient.H (Stage F) and validateSSTDirection.H
// (Stage G) so that the new Gate 6 module (validateGate6SST.H) can reuse them.
// [v2 R4] Extraction only performed after byte-identical diff of the two
// copies of recomputeDesignChain and evaluateObjective (both verified).
//
// All lambdas use [&] capture; every referenced field/parameter is declared
// in main scope BEFORE this file is included (MTO_HF.C, inside the while
// loop, before the Stage E/F/G/6 includes).
//-----------------------------------------------------------------------------

// -- Helper: evaluate objective J = -MeanT / T_ref ---------------------------
auto evaluateObjective = [&]() -> scalar
{
    scalar mFlow = sum
    (
        phi.boundaryField()[coldOutletPatchID]
    );
    scalar mT = sum
    (
        phi.boundaryField()[coldOutletPatchID]
       *T.boundaryField()[coldOutletPatchID]
    );
    reduce(mFlow, sumOp<scalar>());
    reduce(mT, sumOp<scalar>());

    if (mag(mFlow) < SMALL) return GREAT;

    scalar meanT = mT/mFlow;
    return -meanT/Foam::max(thermalReferenceTemperature, SMALL);
};

// -- Helper: evaluate pressure drop (Pa) -------------------------------------
auto evaluatePressureDrop = [&]() -> scalar
{
    if (nObjPatch < 2) return 0.0;

    const label inletPatchI = conPatchList[0];
    const label outletPatchI = conPatchList[1];

    scalar inletArea = sum(mesh.magSf().boundaryField()[inletPatchI]);
    scalar outletArea = sum(mesh.magSf().boundaryField()[outletPatchI]);
    scalar inletPIntegral = sum
    (
        p.boundaryField()[inletPatchI]
       *mesh.magSf().boundaryField()[inletPatchI]
    );
    scalar outletPIntegral = sum
    (
        p.boundaryField()[outletPatchI]
       *mesh.magSf().boundaryField()[outletPatchI]
    );

    reduce(inletArea, sumOp<scalar>());
    reduce(outletArea, sumOp<scalar>());
    reduce(inletPIntegral, sumOp<scalar>());
    reduce(outletPIntegral, sumOp<scalar>());

    scalar pIn = inletPIntegral/Foam::max(inletArea, SMALL);
    scalar pOut = outletPIntegral/Foam::max(outletArea, SMALL);
    return rhoFluid*(pIn - pOut);
};

// -- Helper: evaluate volume constraint V ------------------------------------
auto evaluateVolumeConstraint = [&]() -> scalar
{
    scalar designFluidFraction =
        fvc::domainIntegrate(xh*designMask).value()
       /Foam::max(designVolume, SMALL);
    return solidFractionMin - (1.0 - designFluidFraction);
};

// -- Helper: recompute design chain from x -----------------------------------
auto recomputeDesignChain = [&]()
{
    #include "enforceFixedDesignRegions.H"
    fvScalarMatrix filterEqn
    (
        fvm::laplacian(designFilterFaceMask, xp)
      - fvm::Sp(b, xp)
      + x*b
    );
    filterEqn.solve();
    #include "enforceFixedDesignRegions.H"

    {
        scalar e0 = 0.0001;
        scalar e1 = 0.9999;
        scalar y0 = diff(xp, mesh.V(), designMask, del, e0, n);
        scalar y1 = diff(xp, mesh.V(), designMask, del, e1, n);
        reduce(y0, sumOp<scalar>());
        reduce(y1, sumOp<scalar>());
        if (y0*y1 > 0.0)
        {
            eta5 = 0.5;
        }
        else
        {
            do
            {
                eta5 = 0.5*(e0 + e1);
                scalar y5 = diff(xp, mesh.V(), designMask, del, eta5, n);
                reduce(y5, sumOp<scalar>());
                if (y0*y5 < 0.0) { e1 = eta5; y1 = y5; }
                else { e0 = eta5; y0 = y5; }
            }
            while ((e1 - e0) > projectionEtaTolerance);
            eta5 = 0.5*(e0 + e1);
        }
    }

    for (i = 0; i < n; ++i)
    {
        if (designMask[i] <= 0.5) continue;
        if (xp[i] <= eta5)
        {
            xh[i] = eta5*
            (
                Foam::exp(-del*(1.0 - xp[i]/eta5))
              - (1.0 - xp[i]/eta5)*Foam::exp(-del)
            );
        }
        else
        {
            xh[i] = eta5 + (1.0 - eta5)*
            (
                1.0
              - Foam::exp(-del*(xp[i] - eta5)/(1.0 - eta5))
              + (xp[i] - eta5)*Foam::exp(-del)/(1.0 - eta5)
            );
        }
    }
    #include "enforceFixedDesignRegions.H"

    alpha = Foam::max
    (
        alphamin,
        alphaMax*qu*(1 - xh)/(qu + xh + SMALL)
    );
    forAll(DTMolecular.primitiveFieldRef(), celli)
    {
        const scalar g = xh[celli]*(1.0 + qu)
                       / (qu + xh[celli] + SMALL);
        DTMolecular.primitiveFieldRef()[celli] =
            (KsTensor + (kfIdentity - KsTensor)*g) / rhoc.value();
    }
    DTMolecular.correctBoundaryConditions();
    Q = dimensionedScalar("zeroQ", Q.dimensions(), 0.0);
    if (cells_solid.size() > 0)
    {
        const scalar Qsolid_solver = 5e7/rhoc.value();
        forAll(cells_solid, ci)
        {
            Q[cells_solid[ci]] = Qsolid_solver;
        }
    }
    #include "applyFrozenHotRegionProperties.H"

    DTEffective = DTMolecular + alphaTurbulent * I_dim;
    DTEffective.correctBoundaryConditions();
};

// -- Helper: solve full SST primal with convergence checking -----------------
// Returns number of iterations, or -1 if not converged.
auto solveFullSST = [&](label maxIter, scalar tol) -> label
{
    freezeTurbulenceForValidation = false;

    label consecutiveConverged = 0;
    vectorField prevU = U.primitiveField();
    scalarField prevNut = nutFrozen.primitiveField();
    scalar uRel = GREAT;
    scalar nutRel = GREAT;

    for (label iter = 0; iter < maxIter; ++iter)
    {
        // --- Momentum equation ---
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

        // --- Pressure equation ---
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

        // --- Full SST turbulence correction ---
        turbulence->correct();
        #include "updateFrozenTurbulenceFields.H"

        // --- Thermal equation ---
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

        // --- Convergence check ---
        scalar uNum = 0.0;
        scalar uDen = 0.0;
        forAll(U.primitiveField(), celli)
        {
            if (coldRegionMask[celli] > 0.5)
            {
                uNum += mesh.V()[celli]
                       *magSqr(U[celli] - prevU[celli]);
                uDen += mesh.V()[celli]*magSqr(U[celli]);
            }
        }
        reduce(uNum, sumOp<scalar>());
        reduce(uDen, sumOp<scalar>());
        uRel = Foam::sqrt(uNum/Foam::max(uDen, scalar(VSMALL)));

        scalar nutNum = 0.0;
        scalar nutDen = 0.0;
        forAll(nutFrozen.primitiveField(), celli)
        {
            const scalar w = mesh.V()[celli]*chiTurbulence[celli];
            nutNum += w*sqr(nutFrozen[celli] - prevNut[celli]);
            nutDen += w*sqr(nutFrozen[celli]);
        }
        reduce(nutNum, sumOp<scalar>());
        reduce(nutDen, sumOp<scalar>());
        nutRel = Foam::sqrt
        (
            nutNum/Foam::max(nutDen, scalar(VSMALL))
        );

        prevU = U.primitiveField();
        prevNut = nutFrozen.primitiveField();

        if (uRel < tol && nutRel < tol)
        {
            ++consecutiveConverged;
        }
        else
        {
            consecutiveConverged = 0;
        }

        if ((iter + 1) % 10 == 0 || consecutiveConverged > 0)
        {
            Info<< "    Full-SST iter " << iter + 1
                << ": Urel=" << uRel
                << ", nutRel=" << nutRel
                << ", consecutive=" << consecutiveConverged
                << endl;
        }

        if (consecutiveConverged >= 3)
        {
            Info<< "    Full-SST converged after " << iter + 1
                << " iterations." << endl;
            U.storePrevIter();
            p.storePrevIter();
            phi.storePrevIter();
            return iter + 1;
        }
    }

    Info<< "    Full-SST did NOT converge after " << maxIter
        << " iterations (Urel=" << uRel
        << ", nutRel=" << nutRel << ")" << endl;
    return -1;
};

// -- Helper: evaluate cold-side mass imbalance (relative) --------------------
// Same formula as costfunction.H (coldMassImbalance).
auto evaluateMassBalance = [&]() -> scalar
{
    scalar coldInletFlow = -sum(phi.boundaryField()[conPatchList[0]]);
    scalar coldOutletFlow = sum(phi.boundaryField()[conPatchList[1]]);
    reduce(coldInletFlow, sumOp<scalar>());
    reduce(coldOutletFlow, sumOp<scalar>());
    return mag(coldOutletFlow - coldInletFlow)
      / Foam::max(Foam::max(mag(coldInletFlow), mag(coldOutletFlow)), SMALL);
};

// -- Helper: evaluate thermal balance error (relative) -----------------------
// Convective enthalpy balance from the thermal flux and temperature at the
// cold/hot inlet-outlet patches.  Equivalent to the costfunction.H flux-based
// diagnostic for the validation cases (T is fixedValue/zeroGradient at the
// inlet/outlet patches, so the diffusive flux contribution there is
// negligible).  [v2 Gate 6] The fvMatrix::flux()-based variant crashed when
// applied to an UNSOLVED matrix in parallel, so the direct flux form is used.
auto evaluateHeatBalance = [&]() -> scalar
{
    if (!frozenHotRegionEnabled)
    {
        return 0.0;   // only defined for frozen-hot mode
    }

    #include "updateThermalFlux.H"
    const surfaceScalarField& activeFlux =
        frozenHotRegionEnabled ? phiThermal : phi;

    scalar coldEnthalpyIn = -sum
    (
        activeFlux.boundaryField()[coldInletPatchID]
       *T.boundaryField()[coldInletPatchID]
    );
    scalar coldEnthalpyOut = sum
    (
        activeFlux.boundaryField()[coldOutletPatchID]
       *T.boundaryField()[coldOutletPatchID]
    );
    scalar hotEnthalpyIn = -sum
    (
        activeFlux.boundaryField()[hotInletPatchID]
       *T.boundaryField()[hotInletPatchID]
    );
    scalar hotEnthalpyOut = sum
    (
        activeFlux.boundaryField()[hotOutletPatchID]
       *T.boundaryField()[hotOutletPatchID]
    );

    reduce(coldEnthalpyIn, sumOp<scalar>());
    reduce(coldEnthalpyOut, sumOp<scalar>());
    reduce(hotEnthalpyIn, sumOp<scalar>());
    reduce(hotEnthalpyOut, sumOp<scalar>());

    const scalar Qcold = rhoc.value()*(coldEnthalpyOut - coldEnthalpyIn);
    const scalar Qhot = rhoc.value()*(hotEnthalpyIn - hotEnthalpyOut);
    return mag(Qhot - Qcold)
      / Foam::max(Foam::max(mag(Qhot), mag(Qcold)), SMALL);
};
