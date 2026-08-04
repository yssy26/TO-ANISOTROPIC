// Write fields only after the primal and adjoint equations have been solved for
// the current design and before MMA advances to the next design.
const bool finalSolvedIteration =
    runTime.value() >= runTime.endTime().value() - SMALL;

if (runTime.writeTime() || finalSolvedIteration)
{
    #include "enforceFixedDesignRegions.H"
    runTime.write();

    designMask.write();
    hotRegionMask.write();
    wallRegionMask.write();
    coldRegionMask.write();
    thermalObjectiveMask.write();

    if (frozenHotRegionEnabled)
    {
        hotFaceMask.write();
        coldFaceMask.write();
        UHotFrozen.write();
        phiHotFrozen.write();
        phiThermal.write();
    }
}
