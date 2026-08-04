opt++;

// The optimizer, PDE filter and projection all operate on full local arrays.
// Restore prescribed regions before any material interpolation is evaluated.
#include "enforceFixedDesignRegions.H"

// Continuation is an optimizer operation. A validation run must evaluate one
// fixed material model while the primal and adjoint systems converge.
if (mmaUpdateEnabled && continuationEnabled)
{
    rex = 10;
    if (opt <= rex)
    {
        alphaMax = alphamax*(opt/7.0 + 1.0);
    }
    else
    {
        alphaMax = alphaMax*1.15;
    }

    alphaMax.value() = Foam::min(alphaMax.value(), alphaMAX.value());

    if (opt >= rex)
    {
        qu = Foam::max(0.1 - (opt - rex)*9e-3, 0.01);
    }
}

// Cold design-region interpolation. Fixed hot/wall properties are applied
// afterwards by updateMaterialProperties.H (which includes
// applyFrozenHotRegionProperties.H).
#include "updateMaterialProperties.H"
