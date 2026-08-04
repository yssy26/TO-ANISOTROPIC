//-----------------------------------------------------------------------------
// validateFrozenGradient.H
// Stage F: Frozen U-p-T total gradient validation.
//
// Two-step validation:
//   F1: Thermal-only FD вЂ” fix U,p,phi, re-solve only T.  This isolates
//       the thermal adjoint and avoids pressure-velocity convergence issues.
//   F2: Full frozen FD вЂ” re-solve U,p,T with a FIXED iteration count
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
         зЋw¶‰ћЛkєwµз[WHПH\“›Ь›NВ€B‚€ШШ[\€Y€HЊВ€ШШ[\€Y”HЊВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€Y€
ПH™ШXЭ]™RWJ™\™XЭ[Ы–ШЩ[WNВ€Y”
ПHЩМWVШXЭ]™RWJ™\™XЭ[Ы–ШЩ[WNВ€B€™YXЩJY‹Э[SЬШШ[\ЏЉ
JNВ€™YXЩJY”Э[SЬШШ[\ЏЉ
JNВ€X^Y‘Њ€H›Ш[NЋ›X^
X^Y‘Њ‹XYКYЉJNВ€X^Y”[H›Ш[NЋ›X^
X^Y”[XYКY”
JNВ‚€[™›П—€Њ€\™XЭ[Ы€€\’B€Ћ€РQЉШљЉOH€Y‚€‹РQЉ
OH€Y”[™В‚€›Ь€
X™[Э\HHИЭ\H™Э\СЊ‹њЪ^™J
NИ
КЬЭ\JB€В€ЫЫњЭШШ[\€H™Э\СЊ–ЬЭ\WNВ€Э[\ЭСЊЉКОВ‚€ЛИ\В€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€ШЩ[WHH›Ш[NЋ›Z[ЉKЊШ]™YШЩ[WH
И
™\™XЭ[Ы–ШЩ[WJNВ€B€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€™XЫЫ\]Q\ЪYЫђЪZ[Љ
NВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€KЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЫ™Qњ›Ю™[”љ[X[љ^Y
‘љ^Y]\ЉNВ€ШШ[\€њ\ИH][X]SШљ™XЭ]™J
NВ€ШШ[\€\ИH][X]T™\ЬЭ\™Q›Ь

NВ‚€ЛИZ[ќ\В€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€ШЩ[WHH›Ш[NЋ›X^
ЊШ]™YШЩ[WHH
™\™XЭ[Ы–ШЩ[WJNВ€B€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€™XЫЫ\]Q\ЪYЫђЪZ[Љ
NВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€KЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЫ™Qњ›Ю™[”љ[X[љ^Y
‘љ^Y]\ЉNВ€ШШ[\€›Z[ќ\ИH][X]SШљ™XЭ]™J
NВ€ШШ[\€Z[ќ\ИH][X]T™\ЬЭ\™Q›Ь

NВ‚€ШШ[\€™H
њ\ИH›Z[ќ\КKК‹Њ
љ
NВ€ШШ[\€™[\њ€HXYКY€H™
B€С›Ш[NЋ›X^
›Ш[NЋ›X^
XYКYЉKXYК™
JKШШ[\ЉYKLМ
JNВ‚€ЛИ™\ЬЭ\™H›Ь‘
›Ь›X[^™YћH™\ЬЭ\™Q›ЬX^HИX]ЪЮМWJB€ШШ[\€™H
\ИHZ[ќ\КKК‹Њ
љ
NВ€ШШ[\€™›Ь›HH™€С›Ш[NЋ›X^
™\ЬЭ\™Q›ЬX^KУPS
NВ€ШШ[\€™[\њ”HXYКY”H™›Ь›JB€С›Ш[NЋ›X^€
€›Ш[NЋ›X^
XYКY”
KXYК™›Ь›JJK€ШШ[\ЉYKLМ
B€
NВ‚€[™›П€H€€‹ЉПH€њ\И‹‹OH€›Z[ќ\В€‹С‘H€™€‹™[\њЉШљЉOH€™[\њ‚€‹
ПH€\И‹OH€Z[ќ\В€‹С‘

OH€™›Ь›B€‹™[\њЉ
OH€™[\њ”[™В‚€X^™[\њ›Ь‘Њ€H›Ш[NЋ›X^
X^™[\њ›Ь‘Њ‹™[\њЉNВ€X^™[\њ›Ь”H›Ш[NЋ›X^
X^™[\њ›Ь”™[\њ”
NВ€Y€
™[\њ€€ЊJHЭ[Z[\™\СЊЉКОВ€Y€
™[\њ”€ЌL
HЭ[Z[\™\Ф
КОВ€B‚€ЛИ™\ЭЬ™B€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€[HH[TШ]™YВ€[ЫXЭ[\€HШ]™YВ€HHTШ]™YВ€Y™™XЭ]™HHY™™XЭ]™TШ]™YВ€HHTШ]™YВ€HШ]™YВ€HHTШ]™YВ€HШ]™YВ€]MHH]MTШ]™YЭ[В€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€KЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€B‚€[™›П—€Њ€Э[[X\ћN€X^™[\њ›ЬЏH€X^™[\њ›Ь‘Њ‚€‹Z[\™\ПH€Э[Z[\™\СЊ€‹И€Э[\ЭСЊ€[™В€[™›П€Њ€™\ЬЭ\™H›Ь€X^™[\њ›ЬЏH€X^™[\њ›Ь”€‹Z[\™\ПH€Э[Z[\™\Ф‹И€Э[\ЭСЊ€[™В‚€ЛПOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOB€ЛИЊО€ЫЫXљ[™Y\™XЭ[Ы€‘
XЭX[SPH\™XЭ[ЫЉB€ЛПOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOB€[™›П—‹KKHЊО€ЫЫXљ[™Y\™XЭ[Ы€‘
›Ъ™XЭY™YШ]]™HЬYY[ќ
HKKH€[™В‚€ЛИ™\ЭЬ™HИњ›Ю™[‹XЫЫќ™\™ЩYЭ]H›Ь€ЊВ€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€[HH[TШ]™YВ€[ЫXЭ[\€HШ]™YВ€HHTШ]™YВ€Y™™XЭ]™HHY™™XЭ]™TШ]™YВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€]MHH]MTШ]™YЭ[В€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€KЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ‚€ЛИ\Щ[[™HШљ™XЭ]™K›Ы[YK™\ЬЭ\™H›Ь€ШШ[\€\ЩWЩЊИH][X]SШљ™XЭ]™J
NВ€ШШ[\€\ЩWЩЊИH][X]U›Ы[YPЫЫњЭZ[ќ

NВ€ШШ[\€\ЩWЩЊИH][X]T™\ЬЭ\™Q›Ь

NВ€[™›П€ЊИ\Щ[[™N€ЏH€\ЩWЩЊВ€‹ЏH€\ЩWЩЊВ€‹H€\ЩWЩЊИ[™В‚€ЛИќZ[ЫЫXљ[™Y\™XЭ[ЫЋ€›Ъ™XЭY™YШ]]™HШљ™XЭ]™HЬYY[ќ€ЛИHY™H›ЪЉY™ЩМJH
€ЩМB€ЛИK™K‹™[[Э™HH›Ы[YHЫЫњЭZ[ќЫЫ\Ы™[ќ€ШШ[\‘љY[™YСЬY
‘\ЪYЫ“ШШ[Њ
NВ€ШШ[\‘љY[›ЫЬY
‘\ЪYЫ“ШШ[Њ
NВ€›Ь€
X™[XЭ]™RHHИXЭ]™RH‘\ЪYЫ“ШШ[И
КШXЭ]™RJB€В€™YСЬYШXЭ]™RWHHY™ШXЭ]™RWNВ€›ЫЬYШXЭ]™RWHHЩМVШXЭ]™RWNВ€B‚€ШШ[\€Э•€HЊВ€ШШ[\€Э•€HЊВ€›Ь€
X™[XЭ]™RHHИXЭ]™RH‘\ЪYЫ“ШШ[И
КШXЭ]™RJB€В€Э•€
ПH™YСЬYШXЭ]™RWJќ›ЫЬYШXЭ]™RWNВ€Э•€
ПHЬ\Љ›ЫЬYШXЭ]™RWJNВ€B€™YXЩJЭ•‹Э[SЬШШ[\ЏЉ
JNВ€™YXЩJЭ•‹Э[SЬШШ[\ЏЉ
JNВ‚€ШШ[\€›ЪђЫЩY™€HЊВ€Y€
XYКЭ•ЉH€УPS
B€В€›ЪђЫЩY™€HЭ•‹ЩЭ•ЋВ€B‚€ШШ[\‘љY[ЫЫX‘\ЉY\Ъ›ђЩ[К
KЊ
NВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€ЫЫX‘\–ШЩ[WHH™YСЬYШXЭ]™RWHH›ЪђЫЩY™Љќ›ЫЬYШXЭ]™RWNВ€B‚€ЛИ›Ь›X[^™B€ШШ[\€\“›Ь›LИHЊВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€\“›Ь›LИ
ПHЬ\ЉЫЫX‘\–ШЩ[WJNВ€B€™YXЩJ\“›Ь›LЛЭ[SЬШШ[\ЏЉ
JNВ€\“›Ь›LИH›Ш[NЋњЬ\ќ
\“›Ь›LКNВ€Y€
\“›Ь›LИ€УPS
B€В€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€ЫЫX‘\–ШЩ[WHПH\“›Ь›LОВ€B€B‚€ЛИ\™XЭ[Ы[\љ]]]™\В€ШШ[\€Y“Шљ€HЊВ€ШШ[\€Y•›ЫHЊВ€ШШ[\€Y”HЊВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€Y“Шљ€
ПH™ШXЭ]™RWJЫЫX‘\–ШЩ[WNВ€Y•›Ы
ПHЩМVШXЭ]™RWJЫЫX‘\–ШЩ[WNВ€Y”
ПHЩМWVШXЭ]™RWJЫЫX‘\–ШЩ[WNВ€B€™YXЩJY“Шљ‹Э[SЬШШ[\ЏЉ
JNВ€™YXЩJY•›ЫЭ[SЬШШ[\ЏЉ
JNВ€™YXЩJY”Э[SЬШШ[\ЏЉ
JNВ‚€X^Y”[H›Ш[NЋ›X^
X^Y”[XYКY”
JNВ€X^Y“Шљђ[H›Ш[NЋ›X^
X^Y“Шљђ[XYКY“ШљЉJNВ‚€[™›П€ЊИ\™XЭ[Ы[\љ]]]™\О€ЏH€Y“Шљ‚€‹ЏH€Y•›Ы€‹
›Ь›JOH€Y”[™В‚€ШШ[\€X^™[\њ›Ь‘ЊЫШљ€HЊВ€ШШ[\€X^™[\њ›Ь‘ЊЭ›ЫHЊВ€ШШ[\€X^™[\њ›Ь‘ЊЬHЊВ€ЫЫњЭШШ[\“\Э™Э\СЊИHМYKLKYKL‹™KLџNВ‚€›Ь€
X™[Э\HHИЭ\H™Э\СЊЛњЪ^™J
NИ
КЬЭ\JB€В€ЫЫњЭШШ[\€H™Э\СЊЦЬЭ\WNВ‚€ЛИ\В€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€ШЩ[WHH›Ш[NЋ›Z[‚€
€KЊ€›Ш[NЋ›X^
ЊШ]™YШЩ[WH
И
ЫЫX‘\–ШЩ[WJB€
NВ€B€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€™XЫЫ\]Q\ЪYЫђЪZ[Љ
NВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€KЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЫ™Qњ›Ю™[”љ[X[љ^Y
‘љ^Y]\ЉNВ€ШШ[\€њH][X]SШљ™XЭ]™J
NВ€ШШ[\€њH][X]U›Ы[YPЫЫњЭZ[ќ

NВ€ШШ[\€H][X]T™\ЬЭ\™Q›Ь

NВ‚€ЛИZ[ќ\В€њљ[Z]]™QљY[™YЉ
HHШ]™YВ€›Ьђ[
XЭ]™Q\ЪYЫђЩ[ЛXЭ]™RJB€В€ЫЫњЭX™[Щ[HHXЭ]™Q\ЪYЫђЩ[ЦШXЭ]™RWNВ€ШЩ[WHH›Ш[NЋ›Z[‚€
€KЊ€›Ш[NЋ›X^
ЊШ]™YШЩ[WHH
ЫЫX‘\–ШЩ[WJB€
NВ€B€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€™XЫЫ\]Q\ЪYЫђЪZ[Љ
NВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€HHTШ]™Yњ›ЮЋВ€HШ]™Yњ›ЮЋВ€KЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЬњ™XЭ›Э[™\ћPЫЫ™][ЫњК
NВ€ЫЫ™Qњ›Ю™[”љ[X[љ^Y
‘љ^Y]\ЉNВ€ШШ[\€›HH][X]SШљ™XЭ]™J
NВ€ШШ[\€›HH][X]U›Ы[YPЫЫњЭZ[ќ

NВ€ШШ[\€HH][X]T™\ЬЭ\™Q›Ь

NВ‚€ЛИ‘\љ]]]™\В€ШШ[\€™Шљ€H
њH›JKК‹Њ
љ
NВ€ШШ[\€™›ЫH
њH›JKК‹Њ
љ
NВ€ШШ[\€™H
HJKК‹Њ
љ
NВ€ШШ[\€™›Ь›HH™С›Ш[NЋ›X^
™\ЬЭ\™Q›ЬX^KУPS
NВ‚€ШШ[\€\њ“Шљ€HXYКY“Шљ€H™ШљЉB€С›Ш[NЋ›X^
›Ш[NЋ›X^
XYКY“ШљЉKXYК™ШљЉJKШШ[\ЉYKLМ
JNВ€ШШ[\€\њ•›ЫHXYКY•›ЫH™›Ы
B€С›Ш[NЋ›X^
›Ш[NЋ›X^
XYКY•›Ы
KXYК™›Ы
JKШШ[\ЉYKLМ
JNВ€ШШ[\€\њ”HXYКY”H™›Ь›JB€С›Ш[NЋ›X^€
€›Ш[NЋ›X^
XYКY”
KXYК™›Ь›JJK€ШШ[\ЉYKLМ
B€
NВ‚€[™›П€H€€€ШљЋ€QЏH€Y“Шљ€€‘H€™Шљ‚€€\њЏH€\њ“Шљ‚€€›Ы€QЏH€Y•›Ы€‘H€™›Ы€€\њЏH€\њ•›Ы€€€QЏH€Y”€‘H€™›Ь›B€€\њЏH€\њ”[™В‚€X^™[\њ›Ь‘ЊЫШљ€H›Ш[NЋ›X^
X^™[\њ›Ь‘ЊЫШљ‹\њ“ШљЉNВ€X^™[\њ›Ь‘ЊЭ›ЫH›Ш[NЋ›X^
X^™[\њ›Ь‘ЊЭ›Ы\њ•›Ы
NВ€X^™[\њ›Ь‘ЊЬH›Ш[NЋ›X^
X^™[\њ›Ь‘ЊЬ\њ”
NВ€X^™Шљђ[H›Ш[NЋ›X^
X^™Шљђ[XYК™ШљЉJNВ€B‚€[™›П—€ЊИЭ[[X\ћN€Шљ—Щ\њЏH€X^™[\њ›Ь‘ЊЫШљ‚€‹›ЫЩ\њЏH€X^™[\њ›Ь‘ЊЭ›Ы€‹Щ\њЏH€X^™[\њ›Ь‘ЊЬ[™В‚€ЛИKHљ[[™\ЭЬ™HKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKB€њ™Y^™U\ќ[[ЩQ›Ь•[Y][Ы€HЬљYСњ™Y^™Q›YОВ‚€ЛИKHЭ™\[Э[[X\ћHKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKB€[™›П—‹KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKH€[™В€[™›П”ЭYЩH€Э[[X\ћN€€[™В€[™›П€ЊH
РUM\›X[[Ы›JN€YYX[”™[\њ›ЬЏH€YYX[”™[\њ›Ь‘ЊB€‹X^™[\њ›ЬЏH€X^™[\њ›Ь‘ЊB€‹ЪYЫ“Z\ЫX]ЪH€ЪYЫ“Z\ЫX]ЪЊB€‹Z[\™\ПH€Э[Z[\™\СЊH‹И€Э[\ЭСЊH[™В€[™›П€Њ€
ќ[њ›Ю™[ЉN€X^™[\њ›ЬЏH€X^™[\њ›Ь‘Њ‚€‹Z[\™\ПH€Э[Z[\™\СЊ€‹И€Э[\ЭСЊ‚€‹X^РQџH€X^Y‘Њ€[™В€[™›П€Њ€
™\ЬЭ\™H›Ь
N€X^™[\њ›ЬЏH€X^™[\њ›Ь”€‹Z[\™\ПH€Э[Z[\™\Ф‹И€Э[\ЭСЊ€[™В€[™›П€ЊИ
ЫЫXљ[™Y\ЉN€Шљ—Щ\њЏH€X^™[\њ›Ь‘ЊЫШљ‚€‹›ЫЩ\њЏH€X^™[\њ›Ь‘ЊЭ›Ы€‹Щ\њЏH€X^™[\њ›Ь‘ЊЬ[™В€[™›П‹KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKW€€[™В‚€ЛИ\ЬИЬљ]\љ[ЫЋ‚€ЛИЊH
\›X[[Ы›JHЪЭ[\ЬИЪ]IH\њ›Ь€
YЪЫЫќ™\™Щ[ЩJK‚€ЛИЊ€
ќ[њ›Ю™[ЉH\И[™›Ь›X][Ы[YHИ™\ЬЭ\™K]™[ШЪ]HЫЭ\[™В€ЛИ›Ъ\ЩNИXШЩ\Y€X^\њ›Ь€L	HФ€Y€›ЭРQџ[™С‘\™B€ЛИЫX[
YKLЛ[™XШ][™И™X\‹^™\›ИЬYY[ќ
K‚€ЛИЭЊ€Ш]HHЊH\™Ьљ]\љ[ЫЋ€YYX[€HIKX^HIK›ИЪYЫ‚€ЛИZ\ЫX]Ъ
\Э™X[HS’TУХ“ФPЧРУУ‘PХU’UWФS‹›YЩXЭ[Ы€HШ]H
K‚€›ЫЫЊT\ЬИB€YYX[”™[\њ›Ь‘ЊHHЊB€	‰€X^™[\њ›Ь‘ЊHHЊB€	‰€ЪYЫ“Z\ЫX]ЪЊHOH€	‰€Э[\ЭСЊH€В‚€ЛИ›Ь€Њ‹ЪXЪИY€ЬYY[ќ\И™X\‹^™\›И
[\™XЭ[ЫњИ]™HЫX[РQџ
B€›ЫЫЊ“™X\–™\›ИH
X^Y‘Њ€YKLКNВ€›ЫЫЊ”\ЬИHЊ“™X\–™\›И
Э[Z[\™\СЊ€HЭ[\ЭСЊ‹МИ	‰€X^™[\њ›Ь‘Њ€ЌL
NВ‚€ЛИ™\ЬЭ\™H›Ь€[™›Ь›X][Ы[YHИ™\ЬЭ\™K]™[ШЪ]HЫЭ\[™И›Ъ\ЩK‚€ЛИXШЩ\Y€™X\‹^™\›ИЬYY[ќЬ€L	H\њ›Ь‹‚€›ЫЫ™X\–™\›ИH
X^Y”[YKN
NВ€›ЫЫ\ЬИH™X\–™\›И
Э[Z[\™\ФHЭ[\ЭСЊ‹МИ	‰€X^™[\њ›Ь”KЊ
NВ‚€ЛИЊИЫЫXљ[™Y€[™›Ь›X][Ы[
[™YHЫЫ\Ы™[ќКB€›ЫЫЊФ\ЬИH
X^™[\њ›Ь‘ЊЫШљ€ЌL
B€
X^Y“Шљђ[YKLИ	‰€X^™Шљђ[YKLКNВ‚€Y€
YЊT\ЬКB€В€][\њ›Ь’[‘ќ[Э[Ы‚€‘ЊH
РUM
HЬYY[ќ[Y][Ы€Z[Y€‚€›YYX[”™[\њ›ЬЏH€YYX[”™[\њ›Ь‘ЊB€€
[Z]IJKX^™[\њ›ЬЏH€X^™[\њ›Ь‘ЊB€€
[Z]IJKЪYЫ“Z\ЫX]ЪH€ЪYЫ“Z\ЫX]ЪЊB€‹Z[\™\ПH€Э[Z[\™\СЊB€‹И€Э[\ЭСЊB€‹€ЭЊ€ЋHXYЫ›ЬЭXЬО€JH›Ы[YK[™]][]H™\ЪYX[И‚€™\•›Ы™\ЪYX[ЊH
\™Щ]HYKN
NИЉH‘]X]HXЬ›ЬЬИH‚€ЊИ\ќ\][Ы€Э\И
Y€LYKM\И›Ъ\ЩKYЫZ[]Y[›\™ЩH‚€ќHЫX[\ЭЭ\[њЭXYЩ€Z[[™КNИКHЪYЫ€›\ИO€‚€љ[њЬXЭY’X][њЩ™\‹’›Э[™\ћHЫЫ™][ЫњИ[™[њЬЬЩH‚€\ЬЩ[X›K€‚€X›Ьќ
][\њ›ЬЉNВ€B‚€Y€
YЊ”\ЬКB€В€][\њ›Ь’[‘ќ[Э[Ы‚€‘Њ€
ќ[њ›Ю™[ЉHЬYY[ќ[Y][Ы€Z[Y€‚€›X^™[\њ›ЬЏH€X^™[\њ›Ь‘Њ‚€‹Z[\™\ПH€Э[Z[\™\СЊ‚€‹И€Э[\ЭСЊ‚€X›Ьќ
][\њ›ЬЉNВ€B‚€Y€
\\ЬКB€В€Ш\›љ[™Т[‘ќ[Э[Ы‚€‘Њ€™\ЬЭ\™H›ЬЬYY[ќ[Y][ЫЋ€X^™[\њ›ЬЏH‚€X^™[\њ›Ь”‹Z[\™\ПH€Э[Z[\™\Ф€‹И€Э[\ЭСЊ‚€‹€™\ЬЭ\™H›ЬЬYY[ќX^H™H[њ™[XX›K€‚€[™В€B‚€Y€
YЊФ\ЬКB€В€Ш\›љ[™Т[‘ќ[Э[Ы‚€‘ЊИЫЫXљ[™Y\™XЭ[Ы€[Y][ЫЋ€Шљ—Щ\њЏH‚€X^™[\њ›Ь‘ЊЫШљ‚€‹€ЫЫXљ[™YЬYY[ќ\™XЭ[Ы€X^H™H[њ™[XX›K€‚€[™В€B‚€[™›П”ЭYЩHЋ€TФСQ€[™В€Y€
YЊ“™X\–™\›КB€В€[™›П€
Њ€\ЬЩYЪ][€Ы\[ЩJH€[™В€B€[ЩB€В€[™›П€
ЊЋ€ЬYY[ќ™X\‹^™\›Л‘›Ъ\ЩKYЫZ[]Y8 %XШЩ\Y
H€[™В€B€Y€
™X\–™\›КB€В€[™›П€
€ЬYY[ќ™X\‹^™\›И8 %XШЩ\Y
H€[™В€BџB