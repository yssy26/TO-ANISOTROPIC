//-----------------------------------------------------------------------------
// stageB13ContractionProbe.H
// BFINAL-013 P2: per-term contraction decomposition, included INSIDE
// validateStageB2GradientAmplitude.H right after the baseline adjoint
// projections so it runs at exactly the B2-baseline linearization point
// that produced the BFINAL-012 ADJ numbers.  Switch "stageB13GradientProbe"
// (same switch as the P1 source probe), read-only.
//
// Decomposition along the in-scope D1/D2/D3:
//   gDP: momentum row (recomputed -dAlphaDxh*(U&Uc)*V), pressure row
//        (xh-total field minus momentum), xh-total (gsenshPressureDrop),
//        x-total (gsensPressureDrop, post filter_chainrule).
//   J  : Brinkman (-dAlphaDxh*(U&Ub)*V), thermal-diffusion
//        (thermalDiffusionDerivativeDTCell), xh-total (fsenshMeanT),
//        x-total (fsensMeanT).
// Assembly self-consistency: component sums vs the assembled totals.
//-----------------------------------------------------------------------------

{
    volScalarField b13dAlphaDxh
    (
        IOobject
        (
            "b13dAlphaDxh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroB13", dimless, 0.0)
    );
    b13dAlphaDxh.primitiveFieldRef() =
        -alphaMax*(1.0 + qu)*qu/sqr(qu + xh + SMALL);
    forAll(b13dAlphaDxh.primitiveField(), celli)
    {
        if
        (
            designMask[celli] <= 0.5
         || alpha[celli] <= alphamin.value()*(1.0 + 1e-10)
        )
        {
            b13dAlphaDxh[celli] = 0.0;
        }
    }

    Info<< "B13CONTRACT P2 decomposition (B2-baseline state):" << endl;
    for (label d = 0; d < 3; ++d)
    {
        const scalarField& D = d == 0 ? D1 : (d == 1 ? D2 : D3);

        // gDP terms
        scalar gMomentum = 0.0;
        forAll(D, celli)
        {
            gMomentum +=
               -b13dAlphaDxh[celli]*(U[celli] & Uc[celli])
               *mesh.V()[celli]*D[celli];
        }
        reduce(gMomentum, sumOp<scalar>());
        scalar gTotalXh = 0.0;
        forAll(D, celli)
        {
            gTotalXh += gsenshPressureDrop[celli]*D[celli];
        }
        reduce(gTotalXh, sumOp<scalar>());
        const scalar gPressureRow = gTotalXh - gMomentum;
        scalar gTotalX = 0.0;
        forAll(activeDesignCells, activeI)
        {
            gTotalX += gsensPressureDrop[activeDesignCells[activeI]]
                      *D[activeDesignCells[activeI]];
        }
        reduce(gTotalX, sumOp<scalar>());

        // J terms
        scalar jBrinkman = 0.0;
        scalar jThermal = 0.0;
        forAll(D, celli)
        {
            jBrinkman +=
               -b13dAlphaDxh[celli]*(U[celli] & Ub[celli])
               *mesh.V()[celli]*D[celli];
            jThermal +=
                thermalDiffusionDerivativeDTCell[celli]*D[celli];
        }
        reduce(jBrinkman, sumOp<scalar>());
        reduce(jThermal, sumOp<scalar>());
        scalar jTotalXh = 0.0;
        forAll(D, celli)
        {
            jTotalXh += fsenshMeanT[celli]*D[celli];
        }
        reduce(jTotalXh, sumOp<scalar>());
        scalar jTotalX = 0.0;
        forAll(activeDesignCells, activeI)
        {
            jTotalX += fsensMeanT[activeDesignCells[activeI]]
                      *D[activeDesignCells[activeI]];
        }
        reduce(jTotalX, sumOp<scalar>());

        Info<< "  B13CONTRACT D" << d + 1
            << " gDP: momentum=" << gMomentum
            << " pressureRow=" << gPressureRow
            << " xhTotal=" << gTotalXh
            << " xTotal=" << gTotalX
            << " | J: brinkman=" << jBrinkman
            << " thermal=" << jThermal
            << " xhTotal=" << jTotalXh
            << " xTotal=" << jTotalX
            << " | scales: obj=" << objectiveGradientScale
            << " pd=" << pressureGradientScale << endl;
    }
}
