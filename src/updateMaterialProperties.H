// updateMaterialProperties.H
// Unified material property update from design variable xh.
// Called from: update.H, evaluateCandidate.H, and validation modules.
//
// Computes:
//   1. DTMolecular  — molecular thermal diffusivity (volSymmTensorField)
//   2. alpha        — Brinkman penalty (volScalarField, unchanged)
//   3. Q            — volumetric heat source (volScalarField, unchanged)
//   4. Frozen hot region overrides (via applyFrozenHotRegionProperties.H)
//   5. dDTDxh       — derivative of DTMolecular w.r.t. xh (volSymmTensorField)
//
// Requires in scope: solidK/kfIdentity/KsTensor (from readThermalProperties.H),
//   xh, alpha, alphaMax, alphamin, qu, rhoc, Q, cells_solid, designMask,
//   DTMolecular, DTEffective, dDTDxh, frozenHotRegionEnabled, etc.

// --- Molecular thermal diffusivity (tensor) ---
forAll(DTMolecular.primitiveFieldRef(), celli)
{
    const scalar g = xh[celli]*(1.0 + qu)
                   / (qu + xh[celli] + SMALL);
    DTMolecular.primitiveFieldRef()[celli] =
        (KsTensor + (kfIdentity - KsTensor)*g) / rhoc.value();
}

// --- Brinkman penalty (unchanged, scalar) ---
alpha = Foam::max
(
    alphamin,
    alphaMax*qu*(1.0 - xh)/(qu + xh + SMALL)
);

// --- Volumetric heat source ---
Q = dimensionedScalar("zeroQ", Q.dimensions(), 0.0);
if (cells_solid.size() > 0)
{
    const scalar Qsolid_solver = 5e7/rhoc.value();
    forAll(cells_solid, i)
    {
        Q[cells_solid[i]] = Qsolid_solver;
    }
}

// --- Frozen hot region overrides ---
#include "applyFrozenHotRegionProperties.H"

// --- Boundary conditions ---
DTMolecular.correctBoundaryConditions();
alpha.correctBoundaryConditions();
Q.correctBoundaryConditions();

// --- Thermal diffusivity derivative dD/dxh (tensor) ---
// Computed here because AdjHeatTransfer.H (face loop) needs dDTDxh
// BEFORE sensitivity.H runs. The include order in MTO_HF.C is:
//   update.H -> updateFrozenTurbulenceFields.H -> ... -> AdjHeatTransfer.H
//   -> ... -> sensitivity.H
//
// Cell-loop approach: OpenFOAM 7 lacks operator/ for
// (dimensioned)SymmTensor divided by volScalarField, so we compute
// per-cell like the DTMolecular loop above.
const symmTensor dDTDxhTensor =
    (kfIdentity - KsTensor)*(1.0 + qu)*qu/rhoc.value();

forAll(dDTDxh.primitiveFieldRef(), celli)
{
    const scalar denom = sqr(qu + xh[celli] + SMALL);
    dDTDxh.primitiveFieldRef()[celli] = dDTDxhTensor/denom;

    if (designMask[celli] <= 0.5)
    {
        dDTDxh.primitiveFieldRef()[celli] = symmTensor::zero;
    }
}
dDTDxh.correctBoundaryConditions();
