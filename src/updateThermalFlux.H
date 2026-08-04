// Assemble the thermal convection flux only after NS.H has updated phi for
// the current design. The hot-side contribution remains prescribed.
if (frozenHotRegionEnabled)
{
    phiThermal =
        coldFaceMask*phi
      + hotRhoCpRatio*hotFaceMask*phiHotFrozen;
}
else
{
    phiThermal = phi;
}
