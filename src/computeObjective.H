//-----------------------------------------------------------------------------
// computeObjective.H
// Recompute objective and constraint derivatives based on the current forward
// state (U, p, phi, T).  Must be included AFTER the forward solve (NS.H +
// HeatTransfer.H) and BEFORE the adjoint (AdjHeatTransfer.H).
//
// The thermal objective for frozen-hot mode is the mass-flow-weighted mixed
// outlet temperature:
//
//   Tmix = sum(phi_out * T_out) / sum(phi_out)
//   J    = -Tmix / Tref                         (minimization form)
//
// Derivative with respect to the cell temperature:
//
//   dJ/dT[celli] = -phi[facei] / (sum(phi_out) * Tref)
//
// for cells adjacent to the cold outlet, zero elsewhere.
//
// Because phi changes every outer iteration, this derivative MUST be
// refreshed before each adjoint solve.  The legacy volume-integral objective
// (legacyMeanTemperature) does not depend on phi and is left unchanged.
//-----------------------------------------------------------------------------

// --- Thermal objective derivative dJ/dT -------------------------------------

thermalObjectiveDerivative = 0.0;

if (thermalObjectiveType == "maximizeColdOutletTemperature")
{
    const scalarField& phiOutlet =
        phi.boundaryField()[coldOutletPatchID];
    const labelUList& outletCells =
        mesh.boundary()[coldOutletPatchID].faceCells();

    scalar coldOutletMassFlow = sum(phiOutlet);
    reduce(coldOutletMassFlow, sumOp<scalar>());

    if (mag(coldOutletMassFlow) < SMALL)
    {
        FatalErrorInFunction
            << "Cold outlet mass flow is near zero ("
            << coldOutletMassFlow << "). Cannot compute "
            << "mass-flow-weighted objective derivative."
            << abort(FatalError);
    }

    forAll(outletCells, facei)
    {
        thermalObjectiveDerivative[outletCells[facei]] -=
            phiOutlet[facei]
           /Foam::max
            (
                coldOutletMassFlow*thermalReferenceTemperature,
                SMALL
            );
    }
}
// legacyMeanTemperature: J = domainIntegrate(T*mask) / (V_obj * T_ref)
// dJ/dT[celli] = mask[celli] * V[celli] / (V_obj * T_ref)
else if (thermalObjectiveType == "legacyMeanTemperature")
{
    forAll(thermalObjectiveDerivative, celli)
    {
        thermalObjectiveDerivative[celli] =
            thermalObjectiveMask[celli]*mesh.V()[celli]
           /Foam::max
            (
                thermalObjectiveVolume*thermalReferenceTemperature,
                SMALL
            );
    }
}
