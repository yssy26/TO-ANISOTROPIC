if (frozenHotRegionEnabled)
{
    // U and p represent the cold stream only. Block that system from the
    // prescribed hot-fluid cells and from the separator wall.
    const symmTensor hotKIdentity
    (
        hotK.value(), 0, 0,
        hotK.value(), 0,
        hotK.value()
    );

    forAll(hotRegionCells, i)
    {
        const label celli = hotRegionCells[i];
        alpha[celli] = coldBlockAlpha.value();
        DTMolecular[celli] = hotKIdentity
                           / Foam::max(rhoc.value(), SMALL);
        Q[celli] = 0.0;
    }

    // Separator wall: x_h=0 makes the interpolation formula give K_s/(rho*c)
    // automatically, so no explicit DTMolecular override is needed.
    forAll(wallRegionCells, i)
    {
        const label celli = wallRegionCells[i];
        alpha[celli] = coldBlockAlpha.value();
        Q[celli] = 0.0;
    }

    // Heat enters through the prescribed hot-stream inlet enthalpy. The
    // legacy volumetric source must not remain active in this mode.
    Q = dimensionedScalar("zeroQ", Q.dimensions(), 0.0);

    #include "enforceFixedDesignRegions.H"

    alpha.correctBoundaryConditions();
    DTMolecular.correctBoundaryConditions();
    Q.correctBoundaryConditions();

    designVolume = fvc::domainIntegrate(designMask).value();
    if (designVolume <= SMALL)
    {
        FatalErrorInFunction
            << "Design volume is zero after excluding hot/wall regions."
            << abort(FatalError);
    }
}
