// Do not normalize individual gradient fields by their current maxima. The
// objective values already use fixed reference scales, so dynamic gradient
// normalization would make dfdx inconsistent with the reported objective.
// Volume constraint: V = solidFractionMin - (1 - domainIntegrate(xh*designMask)/designVolume)
// dV/dxh = designMask * mesh.V() / designVolume (per-cell, includes V from domainIntegrate)
gsenshVol.primitiveFieldRef() =
    designMask.primitiveField()*mesh.V()/Foam::max(designVolume, SMALL);

// Remove all physical-density derivatives outside the active design region.
fsenshMeanT.primitiveFieldRef() *= designMask.primitiveField();
gsenshPowerDiss.primitiveFieldRef() *= designMask.primitiveField();
gsenshPressureDrop.primitiveFieldRef() *= designMask.primitiveField();
gsenshVol.primitiveFieldRef() *= designMask.primitiveField();

// Heaviside projection chain rule. eta5 is selected from a global
// volume-preserving equation, so its implicit derivative is included below.
drho = dimensionedScalar("zeroDrho", dimless, 0.0);
volScalarField dProjectionDeta
(
    IOobject
    (
        "dProjectionDeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroDProjectionDeta", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

if (del > SMALL)
{
    for (i = 0; i < n; ++i)
    {
        if (designMask[i] <= 0.5)
        {
            drho[i] = 0.0;
            dProjectionDeta[i] = 0.0;
            continue;
        }

        if (xp[i] <= eta5)
        {
            const scalar projectionExponential =
                Foam::exp(-del*(1.0 - xp[i]/eta5));
            drho[i] =
                del*projectionExponential
              + Foam::exp(-del);
            dProjectionDeta[i] =
                projectionExponential*(1.0 - del*xp[i]/eta5)
              - Foam::exp(-del);
        }
        else
        {
            const scalar projectionExponential = Foam::exp
            (
                -del*(xp[i] - eta5)/(1.0 - eta5)
            );
            drho[i] =
                del*projectionExponential
              + Foam::exp(-del);
            dProjectionDeta[i] =
                projectionExponential
               *(
                    1.0
                  - del*(1.0 - xp[i])/(1.0 - eta5)
                )
              - Foam::exp(-del);
        }
    }

    if (projectionThresholdAdaptive)
    {
        const scalar projectionEtaDenominator = fvc::domainIntegrate
        (
            designMask*dProjectionDeta
        ).value();

        if (mag(projectionEtaDenominator) <= SMALL)
        {
            FatalErrorInFunction
                << "The volume-preserving projection has a singular eta "
                << "derivative." << abort(FatalError);
        }

        // All sensitivity sources are per-cell quantities (not densities).
        // Use plain sum (not domainIntegrate) for the numerator, since
        // domainIntegrate would add an extra mesh.V() factor.
        scalar meanTEtaCorrection = 0.0;
        scalar powerEtaCorrection = 0.0;
        scalar pressureEtaCorrection = 0.0;
        scalar volumeEtaCorrection = 0.0;
        forAll(designMask, celli)
        {
            const scalar maskEta =
                designMask[celli]*dProjectionDeta[celli];
            meanTEtaCorrection    += maskEta*fsenshMeanT[celli];
            powerEtaCorrection    += maskEta*gsenshPowerDiss[celli];
            pressureEtaCorrection += maskEta*gsenshPressureDrop[celli];
            volumeEtaCorrection   += maskEta*gsenshVol[celli];
        }
        reduce(meanTEtaCorrection, sumOp<scalar>());
        reduce(powerEtaCorrection, sumOp<scalar>());
        reduce(pressureEtaCorrection, sumOp<scalar>());
        reduce(volumeEtaCorrection, sumOp<scalar>());
        meanTEtaCorrection    /= projectionEtaDenominator;
        powerEtaCorrection    /= projectionEtaDenominator;
        pressureEtaCorrection /= projectionEtaDenominator;
        volumeEtaCorrection   /= projectionEtaDenominator;

        // --- E5 DIAG: print intermediate values ---
        {
            scalar gsenshVolSumPre = 0.0, gsenshVolMaxPre = 0.0;
            scalar drhoSum = 0.0, drhoMax = 0.0;
            scalar dPDESum = 0.0, dPDEMax = 0.0;
            forAll(designMask, celli)
            {
                if (designMask[celli] > 0.5)
                {
                    gsenshVolSumPre += gsenshVol[celli];
                    gsenshVolMaxPre = max(gsenshVolMaxPre, mag(gsenshVol[celli]));
                    drhoSum += drho[celli];
                    drhoMax = max(drhoMax, mag(drho[celli]));
                    dPDESum += dProjectionDeta[celli];
                    dPDEMax = max(dPDEMax, mag(dProjectionDeta[celli]));
                }
            }
            reduce(gsenshVolSumPre, sumOp<scalar>());
            reduce(gsenshVolMaxPre, maxOp<scalar>());
            reduce(drhoSum, sumOp<scalar>());
            reduce(drhoMax, maxOp<scalar>());
            reduce(dPDESum, sumOp<scalar>());
            reduce(dPDEMax, maxOp<scalar>());
            Info<< "E5 DIAG: BEFORE drho:" << endl
                << "  gsenshVol sum=" << gsenshVolSumPre << " max=" << gsenshVolMaxPre << endl
                << "  drho sum=" << drhoSum << " max=" << drhoMax << endl
                << "  dProjectionDeta sum=" << dPDESum << " max=" << dPDEMax << endl
                << "  projectionEtaDenominator=" << projectionEtaDenominator << endl
                << "  volumeEtaCorrection=" << volumeEtaCorrection << endl;
        }

        // Correct eta correction from volume-preserving projection:
        //   dJ/dxp[j] = dJ/dxh[j]*drho[j] + C_J * V[j] * (1-drho[j])
        // where C_J = sum(dJ/dxh * dProjectionDeta) / sum(dProjectionDeta * V)
        // (dProjectionDeta is zero outside designMask, so designMask is implicit)
        fsenshMeanT *= drho;
        gsenshPowerDiss *= drho;
        gsenshPressureDrop *= drho;
        gsenshVol *= drho;

        forAll(designMask, celli)
        {
            const scalar etaCorrFactor =
                mesh.V()[celli]*(1.0 - drho[celli]);
            fsenshMeanT[celli]        += designMask[celli]*etaCorrFactor*meanTEtaCorrection;
            gsenshPowerDiss[celli]    += designMask[celli]*etaCorrFactor*powerEtaCorrection;
            gsenshPressureDrop[celli] += designMask[celli]*etaCorrFactor*pressureEtaCorrection;
            gsenshVol[celli]          += designMask[celli]*etaCorrFactor*volumeEtaCorrection;
        }

        // --- E5 DIAG: print after eta correction ---
        {
            scalar gsenshVolSumPost = 0.0, gsenshVolMaxPost = 0.0;
            forAll(designMask, celli)
            {
                if (designMask[celli] > 0.5)
                {
                    gsenshVolSumPost += gsenshVol[celli];
                    gsenshVolMaxPost = max(gsenshVolMaxPost, mag(gsenshVol[celli]));
                }
            }
            reduce(gsenshVolSumPost, sumOp<scalar>());
            reduce(gsenshVolMaxPost, maxOp<scalar>());
            Info<< "E5 DIAG: AFTER eta correction:" << endl
                << "  gsenshVol sum=" << gsenshVolSumPost << " max=" << gsenshVolMaxPost << endl;
        }
    }
    else
    {
        fsenshMeanT *= drho;
        gsenshPowerDiss *= drho;
        gsenshPressureDrop *= drho;
        gsenshVol *= drho;
    }
}

// Adjoint PDE filter.
solve
(
    fvm::laplacian(designFilterFaceMask, fsensMeanT)
  - fvm::Sp(b, fsensMeanT)
  + fsenshMeanT*b
);
solve
(
    fvm::laplacian(designFilterFaceMask, gsensPowerDiss)
  - fvm::Sp(b, gsensPowerDiss)
  + gsenshPowerDiss*b
);
solve
(
    fvm::laplacian(designFilterFaceMask, gsensPressureDrop)
  - fvm::Sp(b, gsensPressureDrop)
  + gsenshPressureDrop*b
);
solve
(
    fvm::laplacian(designFilterFaceMask, gsensVol)
  - fvm::Sp(b, gsensVol)
  + gsenshVol*b
);

// Fixed cells solve a local zero-source equation and remain exactly zero.
fsensMeanT.primitiveFieldRef() *= designMask.primitiveField();
gsensPowerDiss.primitiveFieldRef() *= designMask.primitiveField();
gsensPressureDrop.primitiveFieldRef() *= designMask.primitiveField();
gsensVol.primitiveFieldRef() *= designMask.primitiveField();

// All sensitivity sources are per-cell quantities:
//   - fsenshMeanT: Brinkman (-dAlphaDxh*(U&Ub)*V) + thermal (face-sum)
//   - gsenshPowerDiss: dAlphaDxh*(U&Ua)*V/PowerDiss0
//   - gsenshPressureDrop: -dAlphaDxh*(U&Uc)*V
//   - gsenshVol: designMask*V/designVolume
// The Helmholtz filter adjoint produces per-cell output directly — no
// mesh.V() multiplication needed on the filter output.

fsensMeanT.primitiveFieldRef() *= designMask.primitiveField();
gsensPowerDiss.primitiveFieldRef() *= designMask.primitiveField();
gsensPressureDrop.primitiveFieldRef() *= designMask.primitiveField();
gsensVol.primitiveFieldRef() *= designMask.primitiveField();
