//-----------------------------------------------------------------------------
// validateGate6SSTStrict.H
// Strict full-SST directional-derivative validation gate.
//
// Corrections relative to the legacy Gate 6 implementation:
//   1. dgdx[1] is already the derivative of the normalized pressure-drop
//      constraint, so AD_DP is NOT divided by pressureDropMaxPa again.
//   2. Both sign and magnitude agreement are required.
//   3. At least four non-collapsed directions must be tested.
//   4. A repeated full-SST baseline solve checks objective/constraint noise.
//   5. dDTDxh and all turbulence-derived fields are restored transactionally.
//
// This gate remains a validation tool. It does not add k/omega/nut adjoints.
//-----------------------------------------------------------------------------

static bool strictGate6Executed = false;
if
(
   !strictGate6Executed
 && opt >= 1
 && ransFlowModel
 && !ransDirectionValidated
)
{
    strictGate6Executed = true;

    Info<< "\n========================================" << endl;
    Info<< "GATE6-STRICT: full-SST derivative validation" << endl;
    Info<< "========================================\n" << endl;

    const Switch originalFreezeFlag = freezeTurbulenceForValidation;

    const scalar plateauTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "gate6PlateauTolerance",
            0.25
        );
    const scalar magnitudeTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "gate6RelativeGradientTolerance",
            0.35
        );
    const scalar repeatabilityTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "gate6RepeatabilityTolerance",
            1e-5
        );
    const scalar fullSSTTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "gate6FullSSTTolerance",
            1e-4
        );
    const scalar massBalanceTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "gate6MassBalanceTolerance",
            1e-3
        );
    const scalar thermalBalanceTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "gate6ThermalBalanceTolerance",
            1e-2
        );
    const label minimumValidDirections =
        optProperties.lookupOrDefault<label>
        (
            "gate6MinimumValidDirections",
            4
        );
    const label baselineMaxIterations =
        optProperties.lookupOrDefault<label>
        (
            "gate6BaselineMaxIterations",
            300
        );
    const label perturbedMaxIterations =
        optProperties.lookupOrDefault<label>
        (
            "gate6PerturbedMaxIterations",
            200
        );

    if
    (
        plateauTolerance <= 0.0
     || magnitudeTolerance <= 0.0
     || repeatabilityTolerance <= 0.0
     || fullSSTTolerance <= 0.0
     || minimumValidDirections < 4
     || baselineMaxIterations < 1
     || perturbedMaxIterations < 1
    )
    {
        FatalErrorInFunction
            << "Invalid GATE6-STRICT controls." << abort(FatalError);
    }

    const scalarList fdSteps({1e-3, 3e-4, 1e-4});
    const label directionCandidates =
        Foam::max(minimumValidDirections + 4, label(8));
    const scalar nearZero = 1e-12;

    // Save the complete state that may be changed by the validation solves.
    scalarField xSaved = x.primitiveField();
    scalarField xpSaved = xp.primitiveField();
    scalarField xhSaved = xh.primitiveField();
    volScalarField alphaSaved("alpha_gate6_saved", alpha);
    volSymmTensorField DTMolecularSaved("DT_gate6_saved", DTMolecular);
    volSymmTensorField DTEffectiveSaved("DTEff_gate6_saved", DTEffective);
    volSymmTensorField dDTDxhSaved("dDT_gate6_saved", dDTDxh);
    volScalarField QSaved("Q_gate6_saved", Q);
    volVectorField USaved("U_gate6_saved", U);
    volScalarField pSaved("p_gate6_saved", p);
    surfaceScalarField phiSaved("phi_gate6_saved", phi);
    volScalarField TSaved("T_gate6_saved", T);
    volScalarField nutFrozenSaved("nutF_gate6_saved", nutFrozen);
    volScalarField nuEffFrozenSaved("nuEff_gate6_saved", nuEffFrozen);
    volScalarField alphaTurbulentSaved("alphaT_gate6_saved", alphaTurbulent);
    volScalarField fluidIndicatorSaved("fluidI_gate6_saved", fluidIndicator);
    volScalarField chiTurbulenceSaved("chiT_gate6_saved", chiTurbulence);
    volScalarField kSaved
    (
        "k_gate6_saved",
        mesh.lookupObject<volScalarField>("k")
    );
    volScalarField omegaSaved
    (
        "omega_gate6_saved",
        mesh.lookupObject<volScalarField>("omega")
    );
    volScalarField nutSaved
    (
        "nut_gate6_saved",
        mesh.lookupObject<volScalarField>("nut")
    );
    const scalar eta5Saved = eta5;

    auto relativeDifference = [&](const scalar a, const scalar b)
    {
        return mag(a - b)
          / Foam::max
            (
                Foam::max(mag(a), mag(b)),
                scalar(1e-30)
            );
    };

    auto medianOfThree = [&](const scalarList& values)
    {
        const scalar minimumValue = Foam::min
        (
            Foam::min(values[0], values[1]),
            values[2]
        );
        const scalar maximumValue = Foam::max
        (
            Foam::max(values[0], values[1]),
            values[2]
        );
        return values[0] + values[1] + values[2]
             - minimumValue - maximumValue;
    };

    auto sameNonZeroSign = [&](const scalarList& values)
    {
        if
        (
            mag(values[0]) <= nearZero
         || mag(values[1]) <= nearZero
         || mag(values[2]) <= nearZero
        )
        {
            return false;
        }
        return
        (
            (values[0] > 0.0 && values[1] > 0.0 && values[2] > 0.0)
         || (values[0] < 0.0 && values[1] < 0.0 && values[2] < 0.0)
        );
    };

    auto sampleSpread = [&](const scalarList& values)
    {
        scalar minimumMagnitude = GREAT;
        scalar maximumMagnitude = 0.0;
        scalar meanMagnitude = 0.0;
        forAll(values, sampleI)
        {
            const scalar sampleMagnitude = mag(values[sampleI]);
            minimumMagnitude = Foam::min(minimumMagnitude, sampleMagnitude);
            maximumMagnitude = Foam::max(maximumMagnitude, sampleMagnitude);
            meanMagnitude += sampleMagnitude;
        }
        meanMagnitude /= Foam::max(values.size(), label(1));
        return (maximumMagnitude - minimumMagnitude)
             / Foam::max(meanMagnitude, scalar(1e-30));
    };

    // Full-SST baseline plus a second solve for repeatability/noise checking.
    const label baselineIterations1 =
        solveFullSST(baselineMaxIterations, fullSSTTolerance);
    const scalar baselineJ1 = evaluateObjective();
    const scalar baselineDP1 =
        evaluatePressureDrop()/Foam::max(pressureDropMaxPa, SMALL) - 1.0;

    const label baselineIterations2 =
        solveFullSST(baselineMaxIterations, fullSSTTolerance);
    const scalar baselineJ2 = evaluateObjective();
    const scalar baselineDP2 =
        evaluatePressureDrop()/Foam::max(pressureDropMaxPa, SMALL) - 1.0;

    const scalar baselineJNoise = relativeDifference(baselineJ1, baselineJ2);
    const scalar baselineDPNoise = relativeDifference(baselineDP1, baselineDP2);
    const bool repeatabilityPass =
           baselineIterations1 > 0
        && baselineIterations2 > 0
        && baselineJNoise <= repeatabilityTolerance
        && baselineDPNoise <= repeatabilityTolerance;

    Info<< "GATE6-STRICT baseline repeatability: Jnoise="
        << baselineJNoise << ", DPnoise=" << baselineDPNoise
        << ", tolerance=" << repeatabilityTolerance
        << (repeatabilityPass ? " PASS" : " FAIL") << endl;

    // Save the repeated, fully converged baseline for identical warm starts.
    scalarField xBaseline = x.primitiveField();
    scalarField xpBaseline = xp.primitiveField();
    scalarField xhBaseline = xh.primitiveField();
    volScalarField alphaBaseline("alpha_gate6_base", alpha);
    volSymmTensorField DTMolecularBaseline("DT_gate6_base", DTMolecular);
    volSymmTensorField DTEffectiveBaseline("DTEff_gate6_base", DTEffective);
    volSymmTensorField dDTDxhBaseline("dDT_gate6_base", dDTDxh);
    volScalarField QBaseline("Q_gate6_base", Q);
    volVectorField UBaseline("U_gate6_base", U);
    volScalarField pBaseline("p_gate6_base", p);
    surfaceScalarField phiBaseline("phi_gate6_base", phi);
    volScalarField TBaseline("T_gate6_base", T);
    volScalarField nutFrozenBaseline("nutF_gate6_base", nutFrozen);
    volScalarField nuEffFrozenBaseline("nuEff_gate6_base", nuEffFrozen);
    volScalarField alphaTurbulentBaseline("alphaT_gate6_base", alphaTurbulent);
    volScalarField fluidIndicatorBaseline("fluidI_gate6_base", fluidIndicator);
    volScalarField chiTurbulenceBaseline("chiT_gate6_base", chiTurbulence);
    volScalarField kBaseline
    (
        "k_gate6_base",
        mesh.lookupObject<volScalarField>("k")
    );
    volScalarField omegaBaseline
    (
        "omega_gate6_base",
        mesh.lookupObject<volScalarField>("omega")
    );
    volScalarField nutBaseline
    (
        "nut_gate6_base",
        mesh.lookupObject<volScalarField>("nut")
    );
    const scalar eta5Baseline = eta5;

    auto restoreBaseline = [&]()
    {
        x.primitiveFieldRef() = xBaseline;
        xp.primitiveFieldRef() = xpBaseline;
        xh.primitiveFieldRef() = xhBaseline;
        alpha = alphaBaseline;
        DTMolecular = DTMolecularBaseline;
        DTEffective = DTEffectiveBaseline;
        dDTDxh = dDTDxhBaseline;
        Q = QBaseline;
        U = UBaseline;
        p = pBaseline;
        phi = phiBaseline;
        T = TBaseline;
        nutFrozen = nutFrozenBaseline;
        nuEffFrozen = nuEffFrozenBaseline;
        alphaTurbulent = alphaTurbulentBaseline;
        fluidIndicator = fluidIndicatorBaseline;
        chiTurbulence = chiTurbulenceBaseline;
        mesh.lookupObjectRef<volScalarField>("k") = kBaseline;
        mesh.lookupObjectRef<volScalarField>("omega") = omegaBaseline;
        mesh.lookupObjectRef<volScalarField>("nut") = nutBaseline;
        eta5 = eta5Baseline;
        x.correctBoundaryConditions();
        U.correctBoundaryConditions();
        T.correctBoundaryConditions();
    };

    bool allConverged = repeatabilityPass;
    bool allBalancesPass = true;
    bool allPlateausPass = true;
    bool allSignsPass = true;
    bool allMagnitudesPass = true;
    label validDirections = 0;
    label collapsedDirections = 0;
    label signFailures = 0;
    label magnitudeFailures = 0;
    label plateauFailures = 0;
    label convergenceFailures = 0;
    label balanceFailures = 0;

    for (label candidateI = 0; candidateI < directionCandidates; ++candidateI)
    {
        scalarField direction(mesh.nCells(), 0.0);
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            if (candidateI == 0)
            {
                direction[celli] = -dfdx[activeI];
            }
            else
            {
                direction[celli] =
                    Foam::sin
                    (
                        (0.173 + 0.031*candidateI)*(celli + 1)
                      + 0.37*candidateI
                    )
                  + 0.25*Foam::cos
                    (
                        (0.119 + 0.017*candidateI)*(celli + 1)
                    );
            }
        }

        // Project out the volume-gradient component twice, while excluding
        // cells too close to the [0,1] design bounds for the largest step.
        for (label projectionI = 0; projectionI < 2; ++projectionI)
        {
            scalar dotDirectionVolume = 0.0;
            scalar dotVolumeVolume = 0.0;
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                dotDirectionVolume += direction[celli]*gsensVol[celli];
                dotVolumeVolume += sqr(gsensVol[celli]);
            }
            reduce(dotDirectionVolume, sumOp<scalar>());
            reduce(dotVolumeVolume, sumOp<scalar>());

            if (mag(dotVolumeVolume) > SMALL)
            {
                const scalar coefficient =
                    dotDirectionVolume/dotVolumeVolume;
                forAll(activeDesignCells, activeI)
                {
                    const label celli = activeDesignCells[activeI];
                    direction[celli] -= coefficient*gsensVol[celli];
                }
            }

            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                const scalar margin = 2.0*fdSteps[0]*mag(direction[celli]);
                if
                (
                    xBaseline[celli] <= margin
                 || xBaseline[celli] >= 1.0 - margin
                )
                {
                    direction[celli] = 0.0;
                }
            }
        }

        scalar directionNorm2 = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            directionNorm2 += sqr(direction[celli]);
        }
        reduce(directionNorm2, sumOp<scalar>());
        const scalar directionNorm = Foam::sqrt(directionNorm2);

        if (directionNorm <= SMALL)
        {
            ++collapsedDirections;
            continue;
        }

        forAll(activeDesignCells, activeI)
        {
            direction[activeDesignCells[activeI]] /= directionNorm;
        }
        ++validDirections;

        scalar adjointJ = 0.0;
        scalar adjointDP = 0.0;
        forAll(activeDesignCells, activeI)
        {
            const label celli = activeDesignCells[activeI];
            adjointJ += dfdx[activeI]*direction[celli];
            // dgdx[1] already differentiates DP/DPmax - 1.
            adjointDP += dgdx[1][activeI]*direction[celli];
        }
        reduce(adjointJ, sumOp<scalar>());
        reduce(adjointDP, sumOp<scalar>());

        scalarList fdJ(fdSteps.size(), 0.0);
        scalarList fdDP(fdSteps.size(), 0.0);
        bool directionConverged = true;
        bool directionBalancePass = true;

        forAll(fdSteps, stepI)
        {
            const scalar h = fdSteps[stepI];

            restoreBaseline();
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = xBaseline[celli] + h*direction[celli];
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            const label plusIterations =
                solveFullSST(perturbedMaxIterations, fullSSTTolerance);
            const scalar Jplus = evaluateObjective();
            const scalar DPplus =
                evaluatePressureDrop()/Foam::max(pressureDropMaxPa, SMALL)
              - 1.0;
            const scalar massPlus = evaluateMassBalance();
            const scalar heatPlus = evaluateHeatBalance();

            restoreBaseline();
            forAll(activeDesignCells, activeI)
            {
                const label celli = activeDesignCells[activeI];
                x[celli] = xBaseline[celli] - h*direction[celli];
            }
            x.correctBoundaryConditions();
            recomputeDesignChain();
            const label minusIterations =
                solveFullSST(perturbedMaxIterations, fullSSTTolerance);
            const scalar Jminus = evaluateObjective();
            const scalar DPminus =
                evaluatePressureDrop()/Foam::max(pressureDropMaxPa, SMALL)
              - 1.0;
            const scalar massMinus = evaluateMassBalance();
            const scalar heatMinus = evaluateHeatBalance();

            fdJ[stepI] = (Jplus - Jminus)/(2.0*h);
            fdDP[stepI] = (DPplus - DPminus)/(2.0*h);

            const bool convergencePass =
                plusIterations > 0 && minusIterations > 0;
            const bool balancePass =
                   massPlus <= massBalanceTolerance
                && massMinus <= massBalanceTolerance
                && heatPlus <= thermalBalanceTolerance
                && heatMinus <= thermalBalanceTolerance;

            directionConverged = directionConverged && convergencePass;
            directionBalancePass = directionBalancePass && balancePass;

            Info<< "  GATE6-STRICT direction=" << validDirections - 1
                << ", h=" << h
                << ", FD_J=" << fdJ[stepI]
                << ", AD_J=" << adjointJ
                << ", FD_DP=" << fdDP[stepI]
                << ", AD_DP=" << adjointDP
                << ", iterations(+/-)=" << plusIterations
                << '/' << minusIterations
                << ", mass(+/-)=" << massPlus << '/' << massMinus
                << ", heat(+/-)=" << heatPlus << '/' << heatMinus
                << endl;
        }

        const scalar medianJ = medianOfThree(fdJ);
        const scalar medianDP = medianOfThree(fdDP);
        const scalar spreadJ = sampleSpread(fdJ);
        const scalar spreadDP = sampleSpread(fdDP);

        const bool plateauPass =
               sameNonZeroSign(fdJ)
            && sameNonZeroSign(fdDP)
            && spreadJ <= plateauTolerance
            && spreadDP <= plateauTolerance;

        const bool signJPass =
               (mag(adjointJ) <= nearZero && mag(medianJ) <= nearZero)
            || adjointJ*medianJ > 0.0;
        const bool signDPPass =
               (mag(adjointDP) <= nearZero && mag(medianDP) <= nearZero)
            || adjointDP*medianDP > 0.0;
        const bool signPass = signJPass && signDPPass;

        const scalar magnitudeErrorJ =
            mag(adjointJ - medianJ)
          / Foam::max
            (
                Foam::max(mag(adjointJ), mag(medianJ)),
                scalar(1e-30)
            );
        const scalar magnitudeErrorDP =
            mag(adjointDP - medianDP)
          / Foam::max
            (
                Foam::max(mag(adjointDP), mag(medianDP)),
                scalar(1e-30)
            );
        const bool magnitudePass =
               magnitudeErrorJ <= magnitudeTolerance
            && magnitudeErrorDP <= magnitudeTolerance;

        if (!directionConverged) ++convergenceFailures;
        if (!directionBalancePass) ++balanceFailures;
        if (!plateauPass) ++plateauFailures;
        if (!signPass) ++signFailures;
        if (!magnitudePass) ++magnitudeFailures;

        allConverged = allConverged && directionConverged;
        allBalancesPass = allBalancesPass && directionBalancePass;
        allPlateausPass = allPlateausPass && plateauPass;
        allSignsPass = allSignsPass && signPass;
        allMagnitudesPass = allMagnitudesPass && magnitudePass;

        Info<< "  Direction summary: spread(J/DP)="
            << spreadJ << '/' << spreadDP
            << ", magnitudeError(J/DP)="
            << magnitudeErrorJ << '/' << magnitudeErrorDP
            << ", plateau=" << (plateauPass ? "PASS" : "FAIL")
            << ", sign=" << (signPass ? "PASS" : "FAIL")
            << ", magnitude=" << (magnitudePass ? "PASS" : "FAIL")
            << endl;
    }

    const bool directionCountPass =
        validDirections >= minimumValidDirections;
    const bool strictGatePass =
           repeatabilityPass
        && directionCountPass
        && allConverged
        && allBalancesPass
        && allPlateausPass
        && allSignsPass
        && allMagnitudesPass;

    Info<< "\n----------------------------------------" << endl;
    Info<< "GATE6-STRICT Summary" << endl;
    Info<< "  repeatability:       "
        << (repeatabilityPass ? "PASS" : "FAIL") << endl;
    Info<< "  valid directions:    " << validDirections
        << " (required " << minimumValidDirections
        << ", collapsed " << collapsedDirections << ") "
        << (directionCountPass ? "PASS" : "FAIL") << endl;
    Info<< "  convergence failures:" << convergenceFailures << endl;
    Info<< "  balance failures:    " << balanceFailures << endl;
    Info<< "  plateau failures:    " << plateauFailures << endl;
    Info<< "  sign failures:       " << signFailures << endl;
    Info<< "  magnitude failures:  " << magnitudeFailures << endl;
    Info<< "----------------------------------------" << endl;

    if (strictGatePass)
    {
        ransDirectionValidated = true;
        gradientValidated = true;
        Info<< "GATE6-STRICT: PASSED" << endl;
    }
    else
    {
        ransDirectionValidated = false;
        gradientValidated = false;
        WarningInFunction
            << "GATE6-STRICT: FAILED. Full-SST topology updates remain "
            << "blocked unless the explicit direction-only safeguard and "
            << "full-SST candidate acceptance are enabled." << endl;
    }

    // Restore the state that existed before the validation gate.
    freezeTurbulenceForValidation = originalFreezeFlag;
    x.primitiveFieldRef() = xSaved;
    xp.primitiveFieldRef() = xpSaved;
    xh.primitiveFieldRef() = xhSaved;
    alpha = alphaSaved;
    DTMolecular = DTMolecularSaved;
    DTEffective = DTEffectiveSaved;
    dDTDxh = dDTDxhSaved;
    Q = QSaved;
    U = USaved;
    p = pSaved;
    phi = phiSaved;
    T = TSaved;
    nutFrozen = nutFrozenSaved;
    nuEffFrozen = nuEffFrozenSaved;
    alphaTurbulent = alphaTurbulentSaved;
    fluidIndicator = fluidIndicatorSaved;
    chiTurbulence = chiTurbulenceSaved;
    mesh.lookupObjectRef<volScalarField>("k") = kSaved;
    mesh.lookupObjectRef<volScalarField>("omega") = omegaSaved;
    mesh.lookupObjectRef<volScalarField>("nut") = nutSaved;
    eta5 = eta5Saved;
    x.correctBoundaryConditions();
    U.correctBoundaryConditions();
    T.correctBoundaryConditions();

    Info<< "\n========================================\n" << endl;
}
