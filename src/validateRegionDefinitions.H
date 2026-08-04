if (frozenHotRegionEnabled)
{
    label hotFixedSolidOverlap = 0;
    label wallFixedFluidOverlap = 0;

    if (solid_area)
    {
        forAll(cells_solid, i)
        {
            if (hotRegionMask[cells_solid[i]] > 0.5)
            {
                ++hotFixedSolidOverlap;
            }
        }
    }

    if (fluid_area)
    {
        forAll(cells_fluid, i)
        {
            if (wallRegionMask[cells_fluid[i]] > 0.5)
            {
                ++wallFixedFluidOverlap;
            }
        }
    }

    reduce(hotFixedSolidOverlap, sumOp<label>());
    reduce(wallFixedFluidOverlap, sumOp<label>());

    if (hotFixedSolidOverlap > 0)
    {
        FatalErrorInFunction
            << hotFixedSolidOverlap
            << " cells are defined as both hot fluid and legacy fixed solid."
            << abort(FatalError);
    }

    if (wallFixedFluidOverlap > 0)
    {
        FatalErrorInFunction
            << wallFixedFluidOverlap
            << " cells are defined as both separator wall and legacy fixed "
            << "fluid." << abort(FatalError);
    }
}
