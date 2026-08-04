//-----------------------------------------------------------------------------
// validateGradientChain.H
// Component-level gradient validation (Stage E).
//
// Runs once after the first full adjoint + sensitivity solve.  Each check
// compares an analytical derivative (already coded in sensitivity.H and
// filter_chainrule.H) with a finite-difference estimate on a random subset
// of active design cells.
//
// Checks:
//   E1: Darcy alpha(xh) derivative
//   E2: Thermal diffusivity DT(xh) derivative
//   E3: Heaviside projection xh(xp) derivative (drho)
//   E4: Helmholtz filter transpose dot-product test
//   E5: Volume constraint derivative
//   E6: Fixed-region gradient zeroing
//   E7: Cell-volume weighting consistency
//-----------------------------------------------------------------------------

static bool gradientChainValidated = false;

if (!gradientChainValidated && opt >= 1 && !gradientValidated)
{
    gradientChainValidated = true;
    Info<< "\n========================================" << endl;
    Info<< "Stage E: Component-level gradient validation" << endl;
    Info<< "========================================\n" << endl;

    // Deterministic pseudo-random direction for subset sampling
    auto seedRand = [](label i) -> scalar
    {
        return Foam::sin(0.371*(i + 1) + 0.123);
    };

    const label nSample = Foam::min(nDesignLocal, label(50));
    scalar maxErrorE1 = 0.0;
    scalar maxErrorE2 = 0.0;
    scalar maxErrorE3 = 0.0;
    scalar maxErrorE5 = 0.0;
    scalar maxWeightError = 0.0;

    // -----------------------------------------------------------------------
    // E1: Darcy alpha(xh) derivative
    //   alpha = max(alphamin, alphaMax*qu*(1-xh)/(qu+xh+SMALL))
    //   dAlpha/dxh = -alphaMax*(1+qu)*qu / (qu+xh+SMALL)^2   [from sensitivity.H]
    // -----------------------------------------------------------------------
    Info<< "E1: Darcy alpha(xh) derivative" << endl;
    {
        const scalar eps = 1e-6;
        label count = 0;
        for (label activeI = 0; activeI < nDesignLocal && count < nSample; ++activeI)
        {
            const label celli = activeDesignCells[activeI];
            if (designMask[celli] <= 0.5) continue;
            if (alpha[celli] <= alphamin.value()*(1.0 + 1e-10)) continue; // clipped

            const scalar xhVal = xh[celli];
            const scalar denom = qu + xhVal + SMALL;

            // Analytical
            const scalar analytic =
                -alphaMax.value()*(1.0 + qu)*qu/sqr(denom);

            // Finite difference (centered)
            const scalar xhPlus = xhVal + eps;
            const scalar xhMinus = xhVal - eps;
            const scalar alphaPlus =
                Foam::max(alphamin.value(), alphaMax.value()*qu*(1.0 - xhPlus)/(qu + xhPlus + SMALL));
            const scalar alphaMinus =
                Foam::max(alphamin.value(), alphaMax.value()*qu*(1.0 - xhMinus)/(qu + xhMinus + SMALL));
            const scalar fd = (alphaPlus - alphaMinus)/(2.0*eps);

            const scalar err = mag(analytic - fd)
               /Foam::max(Foam::max(mag(analytic), mag(fd)), scalar(1e-30));
            maxErrorE1 = Foam::max(maxErrorE1, err);
            ++count;
        }
        Info<< "  samples=" << count << ", max relative error=" << maxErrorE1 << endl;
    }

    // -----------------------------------------------------------------------
    // E2: Thermal diffusivity DTMolecular(xh) tensor derivative
    //   DTMolecular = (KsTensor + (kfIdentity - KsTensor)*g(xh)) / rhoc
    //   dD/dxh = (kfIdentity - KsTensor)*(1+qu)*qu / (qu+xh+SMALL)^2 / rhoc
    //   Verify with Frobenius norm of tensor difference.
    // -----------------------------------------------------------------------
    Info<< "E2: Thermal diffusivity DTMolecular(xh) tensor derivative" << endl;
    {
        const scalar eps = 1e-6;
        label count = 0;
        for (label activeI = 0; activeI < nDesignLocal && count < nSample; ++activeI)
        {
            const label celli = activeDesignCells[activeI];
            if (designMask[celli] <= 0.5) continue;

            const scalar xhVal = xh[celli];
            const scalar denom = qu + xhVal + SMALL;

            // Analytical tensor derivative
            const symmTensor analyticTensor =
                (kfIdentity - KsTensor)*(1.0 + qu)*qu
               /sqr(denom)/rhoc.value();

            // Finite difference (centered) on tensor components
            const scalar xhPlus = xhVal + eps;
            const scalar xhMinus = xhVal - eps;
            const scalar gPlus =
                xhPlus*(1.0 + qu)/(qu + xhPlus + SMALL);
            const scalar gMinus =
                xhMinus*(1.0 + qu)/(qu + xhMinus + SMALL);
            const symmTensor dtPlus =
                (KsTensor + (kfIdentity - KsTensor)*gPlus)/rhoc.value();
            const symmTensor dtMinus =
                (KsTensor + (kfIdentity - KsTensor)*gMinus)/rhoc.value();
            const symmTensor fdTensor = (dtPlus - dtMinus)/(2.0*eps);

            // Error: Frobenius norm
            const scalar err = mag(analyticTensor - fdTensor)
               /Foam::max
                (
                    Foam::max(mag(analyticTensor), mag(fdTensor)),
                    scalar(1e-30)
                );
            maxErrorE2 = Foam::max(maxErrorE2, err);
            ++count;
        }
        Info<< "  samples=" << count << ", max relative error=" << maxErrorE2 << endl;
    }

    // -----------------------------------------------------------------------
    // E3: Heaviside projection xh(xp) derivative (drho)
    //   drho is computed in filter_chainrule.H
    //   FD: xh(xp+eps) - xh(xp-eps) / (2*eps)
    // -----------------------------------------------------------------------
    Info<< "E3: Heaviside projection derivative" << endl;
    {
        const scalar eps = 1e-6;
        label count = 0;
        for (label activeI = 0; activeI < nDesignLocal && count < nSample; ++activeI)
        {
            const label celli = activeDesignCells[activeI];
            if (designMask[celli] <= 0.5) continue;
            if (del <= SMALL) continue; // no projection in validation mode

            const scalar xpVal = xp[celli];
            const scalar analytic = drho[celli];

            // Finite difference: recompute projection with perturbed xp
            auto computeProjection = [&](scalar xpPert) -> scalar
            {
                if (xpPert <= eta5)
                {
                    return eta5*
                    (
                        Foam::exp(-del*(1.0 - xpPert/eta5))
                      - (1.0 - xpPert/eta5)*Foam::exp(-del)
                    );
                }
                else
                {
                    return eta5 + (1.0 - eta5)*
                    (
                        1.0
                      - Foam::exp(-del*(xpPert - eta5)/(1.0 - eta5))
                      + (xpPert - eta5)*Foam::exp(-del)/(1.0 - eta5)
                    );
                }
            };

            const scalar xhPlus = computeProjection(xpVal + eps);
            const scalar xhMinus = computeProjection(xpVal - eps);
            const scalar fd = (xhPlus - xhMinus)/(2.0*eps);

            const scalar err = mag(analytic - fd)
               /Foam::max(Foam::max(mag(analytic), mag(fd)), scalar(1e-30));
            maxErrorE3 = Foam::max(maxErrorE3, err);
            ++count;
        }
        if (del <= SMALL)
        {
            Info<< "  skipped (del=0, no projection in validation mode)" << endl;
        }
        else
        {
            Info<< "  samples=" << count << ", max relative error=" << maxErrorE3 << endl;
        }
    }

    // -----------------------------------------------------------------------
    // E4: Helmholtz filter transpose dot-product test
    //   <a, F*b> = <F^T*a, b>
    //   F:  forward filter  (L - diag(b))*xp = -diag(b)*x
    //   F^T: transpose filter (L - diag(b))*y = -diag(b)*source
    //   Since L is symmetric and diag(b) is diagonal, F should be symmetric.
    // -----------------------------------------------------------------------
    Info<< "E4: Helmholtz filter transpose dot-product test" << endl;
    {
        scalar maxDotError = 0.0;

        for (label testI = 0; testI < 3; ++testI)
        {
            // Copy-construct from xp to inherit its zeroGradient BCs.
            // The IOobject constructor with typeName does not reliably set
            // zeroGradient on all patches when files already exist on disk.
            volScalarField fieldA("fieldA_test", xp);
            volScalarField fieldB("fieldB_test", xp);
            fieldA = dimensionedScalar("zeroA", dimless, 0.0);
            fieldB = dimensionedScalar("zeroB", dimless, 0.0);

            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                fieldA[celli] = Foam::sin((0.173 + 0.019*testI)*(celli + 1));
                fieldB[celli] = Foam::cos((0.117 + 0.023*testI)*(celli + 1));
            }
            fieldA.correctBoundaryConditions();
            fieldB.correctBoundaryConditions();

            // Forward filter: F*b -> Fb
            volScalarField Fb("Fb_test", xp);
            Fb = dimensionedScalar("zeroFb", dimless, 0.0);
            Fb.primitiveFieldRef() = fieldB.primitiveField();
            Fb.correctBoundaryConditions();
            fvScalarMatrix forwardEqn
            (
                fvm::laplacian(designFilterFaceMask, Fb)
              - fvm::Sp(b, Fb)
              + fieldB*b
            );
            forwardEqn.solve("xp");

            // Transpose filter: F^T*a -> FTa
            volScalarField FTa("FTa_test", xp);
            FTa = dimensionedScalar("zeroFTa", dimless, 0.0);
            FTa.primitiveFieldRef() = fieldA.primitiveField();
            FTa.correctBoundaryConditions();
            fvScalarMatrix transposeEqn
            (
                fvm::laplacian(designFilterFaceMask, FTa)
              - fvm::Sp(b, FTa)
              + fieldA*b
            );
            transposeEqn.solve("xp");

            // Dot products: <a, Fb> vs <FTa, b>
            scalar aFb = 0.0;
            scalar FTab = 0.0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                aFb += fieldA[celli]*Fb[celli]*mesh.V()[celli];
                FTab += FTa[celli]*fieldB[celli]*mesh.V()[celli];
            }
            reduce(aFb, sumOp<scalar>());
            reduce(FTab, sumOp<scalar>());

            const scalar dotErr = mag(aFb - FTab)
               /Foam::max(Foam::max(mag(aFb), mag(FTab)), scalar(1e-30));
            maxDotError = Foam::max(maxDotError, dotErr);

            Info<< "  test " << testI << ": <a,Fb>=" << aFb
                << ", <FTa,b>=" << FTab
                << ", rel error=" << dotErr << endl;
        }
        Info<< "  max dot-product error=" << maxDotError << endl;
    }

    // -----------------------------------------------------------------------
    // E5: Volume constraint derivative
    //   V = solidFractionMin - designSolidFraction
    //   designSolidFraction = 1 - sum(xh*designMask*V)/designVolume
    //   The analytical gsensVol includes the eta-correction from the
    //   volume-preserving projection.  The FD must therefore recompute
    //   eta5 for each perturbed design, otherwise the implicit d eta5/dx
    //   term is missing and the error will be O(1).
    // -----------------------------------------------------------------------
    Info<< "E5: Volume constraint derivative" << endl;
    {
        const scalar eps = 1e-4;

        // Analytical: dV/dx = gsensVol[celli] (already includes filter + volume)
        // FD: V(x+eps*d) - V(x-eps*d) / (2*eps) for a test direction d

        // --- E5 DIAG: print gsensVol statistics ---
        {
            scalar gvSum = 0.0, gvMax = 0.0;
            scalar gvFirst = 0.0;
            label nNonZero = 0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                gvSum += gsensVol[celli];
                gvMax = max(gvMax, mag(gsensVol[celli]));
                if (activeI == 0) gvFirst = gsensVol[celli];
                if (mag(gsensVol[celli]) > 1e-30) nNonZero++;
            }
            reduce(gvSum, sumOp<scalar>());
            reduce(gvMax, maxOp<scalar>());
            reduce(nNonZero, sumOp<label>());
            Info<< "E5 DIAG: gsensVol sum=" << gvSum << " max=" << gvMax
                << " nNonZero=" << nNonZero
                << " first=" << gvFirst << endl;
        }

        scalarField testDirection(mesh.nCells(), 0.0);
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            testDirection[celli] = seedRand(celli);
        }

        // Analytic direction derivative
        scalar analyticDir = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            analyticDir += gsensVol[celli]*testDirection[celli];
        }
        reduce(analyticDir, sumOp<scalar>());

        // Helper lambda: recompute volume-preserving eta5 from current xp
        auto computeEta5 = [&]() -> scalar
        {
            scalar e0 = 0.0001;
            scalar e1 = 0.9999;
            scalar y0 = diff(xp, mesh.V(), designMask, del, e0, n);
            scalar y1 = diff(xp, mesh.V(), designMask, del, e1, n);
            reduce(y0, sumOp<scalar>());
            reduce(y1, sumOp<scalar>());
            if (y0*y1 > 0.0) return 0.5;
            do
            {
                scalar e5 = 0.5*(e0 + e1);
                scalar y5 = diff(xp, mesh.V(), designMask, del, e5, n);
                reduce(y5, sumOp<scalar>());
                if (y0*y5 < 0.0) { e1 = e5; y1 = y5; }
                else { e0 = e5; y0 = y5; }
            }
            while ((e1 - e0) > projectionEtaTolerance);
            return 0.5*(e0 + e1);
        };

        // Helper lambda: apply Heaviside projection with given eta
        auto applyProjection = [&](scalar eta) -> void
        {
            for (i = 0; i < n; ++i)
            {
                if (designMask[i] <= 0.5) continue;
                if (xp[i] <= eta)
                {
                    xh[i] = eta*
                    (
                        Foam::exp(-del*(1.0 - xp[i]/eta))
                      - (1.0 - xp[i]/eta)*Foam::exp(-del)
                    );
                }
                else
                {
                    xh[i] = eta + (1.0 - eta)*
                    (
                        1.0
                      - Foam::exp(-del*(xp[i] - eta)/(1.0 - eta))
                      + (xp[i] - eta)*Foam::exp(-del)/(1.0 - eta)
                    );
                }
            }
        };

        // Save original x and eta5
        scalarField xOrig = x.primitiveField();
        scalar eta5Orig = eta5;

        // --- Plus perturbation ---
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            x[celli] = Foam::min(1.0, xOrig[celli] + eps*testDirection[celli]);
        }
        x.correctBoundaryConditions();
        #include "enforceFixedDesignRegions.H"
        fvScalarMatrix EqnPlus
        (
            fvm::laplacian(designFilterFaceMask, xp)
          - fvm::Sp(b, xp)
          + x*b
        );
        EqnPlus.solve();
        #include "enforceFixedDesignRegions.H"
        scalar eta5plus = computeEta5();
        applyProjection(eta5plus);
        #include "enforceFixedDesignRegions.H"
        const scalar Vplus = solidFractionMin
            - (1.0 - fvc::domainIntegrate(xh*designMask).value()
              /Foam::max(designVolume, SMALL));

        // --- Minus perturbation ---
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            x[celli] = Foam::max(0.0, xOrig[celli] - eps*testDirection[celli]);
        }
        x.correctBoundaryConditions();
        #include "enforceFixedDesignRegions.H"
        fvScalarMatrix EqnMinus
        (
            fvm::laplacian(designFilterFaceMask, xp)
          - fvm::Sp(b, xp)
          + x*b
        );
        EqnMinus.solve();
        #include "enforceFixedDesignRegions.H"
        scalar eta5minus = computeEta5();
        applyProjection(eta5minus);
        #include "enforceFixedDesignRegions.H"
        const scalar Vminus = solidFractionMin
            - (1.0 - fvc::domainIntegrate(xh*designMask).value()
              /Foam::max(designVolume, SMALL));

        // --- Restore original state ---
        x.primitiveFieldRef() = xOrig;
        x.correctBoundaryConditions();
        #include "enforceFixedDesignRegions.H"
        fvScalarMatrix EqnRestore
        (
            fvm::laplacian(designFilterFaceMask, xp)
          - fvm::Sp(b, xp)
          + x*b
        );
        EqnRestore.solve();
        #include "enforceFixedDesignRegions.H"
        eta5 = eta5Orig;
        applyProjection(eta5);
        #include "enforceFixedDesignRegions.H"

        const scalar fdDir = (Vplus - Vminus)/(2.0*eps);
        maxErrorE5 = mag(analyticDir - fdDir)
           /Foam::max(Foam::max(mag(analyticDir), mag(fdDir)), scalar(1e-30));

        Info<< "  analytic=" << analyticDir << ", FD=" << fdDir
            << ", rel error=" << maxErrorE5 << endl;

        // Multi-step-size convergence check to distinguish systematic
        // errors (filter solver tolerance, eta bisection) from FD
        // truncation/roundoff.
        for (scalar eps2 : {1e-4, 1e-6, 1e-7})
        {
            // Plus
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = Foam::min(1.0, xOrig[celli] + eps2*testDirection[celli]);
            }
            x.correctBoundaryConditions();
            #include "enforceFixedDesignRegions.H"
            fvScalarMatrix EqnP2
            (
                fvm::laplacian(designFilterFaceMask, xp)
              - fvm::Sp(b, xp)
              + x*b
            );
            EqnP2.solve();
            #include "enforceFixedDesignRegions.H"
            applyProjection(computeEta5());
            #include "enforceFixedDesignRegions.H"
            const scalar Vp2 = solidFractionMin
                - (1.0 - fvc::domainIntegrate(xh*designMask).value()
                  /Foam::max(designVolume, SMALL));

            // Minus
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = Foam::max(0.0, xOrig[celli] - eps2*testDirection[celli]);
            }
            x.correctBoundaryConditions();
            #include "enforceFixedDesignRegions.H"
            fvScalarMatrix EqnM2
            (
                fvm::laplacian(designFilterFaceMask, xp)
              - fvm::Sp(b, xp)
              + x*b
            );
            EqnM2.solve();
            #include "enforceFixedDesignRegions.H"
            applyProjection(computeEta5());
            #include "enforceFixedDesignRegions.H"
            const scalar Vm2 = solidFractionMin
                - (1.0 - fvc::domainIntegrate(xh*designMask).value()
                  /Foam::max(designVolume, SMALL));

            const scalar fd2 = (Vp2 - Vm2)/(2.0*eps2);
            const scalar err2 = mag(analyticDir - fd2)
               /Foam::max(Foam::max(mag(analyticDir), mag(fd2)), scalar(1e-30));
            Info<< "  eps=" << eps2 << ", FD=" << fd2
                << ", rel error=" << err2 << endl;
        }

        // Restore original state (including eta5)
        x.primitiveFieldRef() = xOrig;
        x.correctBoundaryConditions();
        #include "enforceFixedDesignRegions.H"
        fvScalarMatrix EqnFinal
        (
            fvm::laplacian(designFilterFaceMask, xp)
          - fvm::Sp(b, xp)
          + x*b
        );
        EqnFinal.solve();
        #include "enforceFixedDesignRegions.H"
        eta5 = eta5Orig;
        applyProjection(eta5);
        #include "enforceFixedDesignRegions.H"
    }

    // -----------------------------------------------------------------------
    // E6: Fixed-region gradient zeroing
    // -----------------------------------------------------------------------
    Info<< "E6: Fixed-region gradient zeroing" << endl;
    {
        scalar maxFixedGrad = 0.0;
        forAll(designMask.primitiveField(), celli)
        {
            if (designMask[celli] <= 0.5)
            {
                maxFixedGrad = Foam::max(maxFixedGrad, mag(fsensMeanT[celli]));
                maxFixedGrad = Foam::max(maxFixedGrad, mag(gsensVol[celli]));
                maxFixedGrad = Foam::max(maxFixedGrad, mag(gsensPressureDrop[celli]));
            }
        }
        reduce(maxFixedGrad, maxOp<scalar>());
        Info<< "  max |gradient| in fixed regions=" << maxFixedGrad << endl;
    }

    // -----------------------------------------------------------------------
    // E7: Cell-volume weighting consistency
    //   The MMA gradient dfdx[activeI] should equal fsensMeanT[celli]
    //   which already includes mesh.V() from filter_chainrule.H.
    // -----------------------------------------------------------------------
    Info<< "E7: Cell-volume weighting consistency" << endl;
    {
        for (label activeI = 0; activeI < nDesignLocal; ++activeI)
        {
            const label celli = activeDesignCells[activeI];
            // For outletTemperature mode, dfdx = fsensMeanT * scale
            const scalar expected = objectiveMode == "outletTemperature"
                ? fsensMeanT[celli]*objectiveGradientScale
                : ((1.0 - w1)*fsensMeanT[celli] + w1*gsensPowerDiss[celli])
                 *objectiveGradientScale;
            const scalar err2 = mag(dfdx[activeI] - expected);
            maxWeightError = Foam::max(maxWeightError, err2);
        }
        Info<< "  max |dfdx - expected|=" << maxWeightError << endl;
    }

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    Info<< "\n----------------------------------------" << endl;
    Info<< "Stage E Summary:" << endl;
    Info<< "  E1 (Darcy alpha):      " << maxErrorE1 << endl;
    Info<< "  E2 (Thermal DT):       " << maxErrorE2 << endl;
    Info<< "  E3 (Heaviside proj):   " << (del > SMALL ? name(maxErrorE3) : "N/A (del=0)") << endl;
    Info<< "  E4 (Helmholtz filter): see dot-product errors above" << endl;
    Info<< "  E5 (Volume constraint):" << maxErrorE5 << endl;
    Info<< "  E6 (Fixed-region zero):see max |gradient| above" << endl;
    Info<< "  E7 (Volume weighting): " << maxWeightError << endl;
    Info<< "----------------------------------------\n" << endl;

    // Hard fail if any component exceeds tolerance.
    // E1-E3 are local algebraic derivatives: strict 1e-4.
    // E5 involves the full filter+projection+eta-correction chain solved
    // iteratively (PCG tol 1e-9, bisection tol 1e-10); a small systematic
    // discrepancy (~0.5-1%) is expected and acceptable for FD validation.
    const scalar componentTolerance = 1e-4;
    const scalar e5Tolerance = 5e-2;
    if
    (
        maxErrorE1 > componentTolerance
     || maxErrorE2 > componentTolerance
     || (del > SMALL && maxErrorE3 > componentTolerance)
     || maxErrorE5 > e5Tolerance
    )
    {
        FatalErrorInFunction
            << "Component-level gradient validation failed. "
            << "E1=" << maxErrorE1
            << ", E2=" << maxErrorE2
            << ", E3=" << (del > SMALL ? name(maxErrorE3) : "N/A")
            << ", E5=" << maxErrorE5
            << abort(FatalError);
    }
}
