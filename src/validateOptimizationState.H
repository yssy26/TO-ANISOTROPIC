const scalar minTemperature = gMin(T.primitiveField());
const scalar maxTemperature = gMax(T.primitiveField());
const scalar maxColdVelocity = gMax(mag(U.primitiveField()));

if
(
   !std::isfinite(MeanT)
 || !std::isfinite(Objective)
 || !std::isfinite(PowerDiss)
 || !std::isfinite(PressureDropPa)
 || !std::isfinite(V)
 || !std::isfinite(minTemperature)
 || !std::isfinite(maxTemperature)
 || !std::isfinite(maxColdVelocity)
)
{
    FatalErrorInFunction
        << "Non-finite primal/objective state detected: " << nl
        << "  thermalMetric=" << MeanT << nl
        << "  objective=" << Objective << nl
        << "  powerDiss=" << PowerDiss << nl
        << "  pressureDropPa=" << PressureDropPa << nl
        << "  volumeConstraint=" << V << nl
        << "  Tmin/Tmax=" << minTemperature << '/' << maxTemperature << nl
        << "  max|U|=" << maxColdVelocity
        << abort(FatalError);
}

scalar fixedDesignError = 0.0;
if (solid_area)
{
    forAll(cells_solid, i)
    {
        const label celli = cells_solid[i];
        fixedDesignError = Foam::max
        (
            fixedDesignError,
            Foam::max(mag(x[celli]), Foam::max(mag(xp[celli]), mag(xh[celli])))
        );
    }
}
if (fluid_area)
{
    forAll(cells_fluid, i)
    {
        const label celli = cells_fluid[i];
        fixedDesignError = Foam::max
        (
            fixedDesignError,
            Foam::max
            (
                mag(x[celli] - 0.999),
                Foam::max(mag(xp[celli] - 0.999), mag(xh[celli] - 0.999))
            )
        );
    }
}
if (frozenHotRegionEnabled)
{
    forAll(hotRegionCells, i)
    {
        const label celli = hotRegionCells[i];
        fixedDesignError = Foam::max
        (
            fixedDesignError,
            Foam::max
            (
                mag(x[celli] - 0.999),
                Foam::max(mag(xp[celli] - 0.999), mag(xh[celli] - 0.999))
            )
        );
    }
    forAll(wallRegionCells, i)
    {
        const label celli = wallRegionCells[i];
        fixedDesignError = Foam::max
        (
            fixedDesignError,
            Foam::max(mag(x[celli]), Foam::max(mag(xp[celli]), mag(xh[celli])))
        );
    }
}
reduce(fixedDesignError, maxOp<scalar>());
Info<< "Fixed-region design error=" << fixedDesignError << endl << flush;
if (fixedDesignError > 1e-14)
{
    FatalErrorInFunction
        << "A prescribed design region changed: max error="
        << fixedDesignError << abort(FatalError);
}
if (ransFlowModel)
{
    const volScalarField& kRANS =
        mesh.lookupObject<volScalarField>("k");
    const volScalarField& omegaRANS =
        mesh.lookupObject<volScalarField>("omega");
    const volScalarField& nutRANS =
        mesh.lookupObject<volScalarField>("nut");

    const scalar minK = gMin(kRANS.primitiveField());
    const scalar minOmega = gMin(omegaRANS.primitiveField());
    const scalar minNut = gMin(nutRANS.primitiveField());
    const scalar maxNutToNu = gMax
    (
        nutFrozen.primitiveField()/Foam::max(nu.value(), SMALL)
    );
    scalar maxFixedNutFrozen = 0.0;
    forAll(nutFrozen.primitiveField(), celli)
    {
        if (coldRegionMask[celli] <= 0.5)
        {
            maxFixedNutFrozen = Foam::max
            (
                maxFixedNutFrozen,
                mag(nutFrozen[celli])
            );
        }
    }
    reduce(maxFixedNutFrozen, maxOp<scalar>());

    Info<< "Frozen-RANS diagnostics: min(k)=" << minK
        << ", min(omega)=" << minOmega
        << ", min(nut)=" << minNut
        << ", max(nutFrozen/nu)=" << maxNutToNu
        << ", max fixed-region nutFrozen=" << maxFixedNutFrozen
        << ", cold mass imbalance=" << coldMassImbalance
        << ", converged=" << turbulenceConverged << endl << flush;

    if
    (
       !std::isfinite(minK)
     || !std::isfinite(minOmega)
     || !std::isfinite(minNut)
     || !std::isfinite(maxNutToNu)
     || minK <= 0.0
     || minOmega <= 0.0
     || minNut < -SMALL
    )
    {
        FatalErrorInFunction
            << "Invalid SST state: min(k)=" << minK
            << ", min(omega)=" << minOmega
            << ", min(nut)=" << minNut
            << ", max(nutFrozen/nu)=" << maxNutToNu
            << abort(FatalError);
    }
    if (maxFixedNutFrozen > 1e-14)
    {
        FatalErrorInFunction
            << "nutFrozen is non-zero in a prescribed non-cold region: "
            << maxFixedNutFrozen << abort(FatalError);
    }
    if (coldMassImbalance > 1e-8)
    {
        FatalErrorInFunction
            << "Cold-side mass imbalance " << coldMassImbalance
            << " exceeds the RANS validation limit 1e-8."
            << abort(FatalError);
    }
    const scalar thermalBalanceToleranceVal =
        optProperties.lookupOrDefault<scalar>("thermalBalanceTolerance", 1e-2);
    if (thermalBalanceError > thermalBalanceToleranceVal)
    {
        FatalErrorInFunction
            << "Thermal balance error " << thermalBalanceError
            << " exceeds the RANS validation limit "
            << thermalBalanceToleranceVal << "."
            << abort(FatalError);
    }
    if
    (
        runTime.value() >= runTime.endTime().value() - SMALL
     && !freezeColdFlowForValidation
     && !turbulenceConverged
    )
    {
        FatalErrorInFunction
            << "The final forward-only RANS iteration is not converged."
            << abort(FatalError);
    }
}
