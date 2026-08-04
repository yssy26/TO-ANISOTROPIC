// Keep all prescribed regions outside the optimizer and projection chain.
// TO-HF convention: x=1 is fluid and x=0 is solid.

if (solid_area)
{
    forAll(cells_solid, i)
    {
        const label celli = cells_solid[i];
        x[celli] = 0.0;
        xp[celli] = 0.0;
        xh[celli] = 0.0;
    }
}

if (fluid_area)
{
    forAll(cells_fluid, i)
    {
        const label celli = cells_fluid[i];
        x[celli] = 0.999;
        xp[celli] = 0.999;
        xh[celli] = 0.999;
    }
}

if (frozenHotRegionEnabled)
{
    forAll(hotRegionCells, i)
    {
        const label celli = hotRegionCells[i];
        x[celli] = 0.999;
        xp[celli] = 0.999;
        xh[celli] = 0.999;
    }

    forAll(wallRegionCells, i)
    {
        const label celli = wallRegionCells[i];
        x[celli] = 0.0;
        xp[celli] = 0.0;
        xh[celli] = 0.0;
    }
}

x.correctBoundaryConditions();
xp.correctBoundaryConditions();
xh.correctBoundaryConditions();
