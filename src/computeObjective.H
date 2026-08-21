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

// --- BFINAL-025: M_frozen snapshot (maximizeTotalHeatTransfer only) ----------
// Q = N - T_in*M, J_Q = -Q/(Tref*M_frozen).  M_frozen = sum(phi_out) is
// snapshotted ONCE per outer state HERE and kept CONSTANT for every J
// evaluation in that state: probeFD (validateStageB2GradientAmplitude.H) and
// the F1/F2/F3 probes (validateFrozenGradient.H) evaluate J via
// evaluateObjective() and never re-include this file, so all FD +/-h points
// and all adjoint projections share the same snapshot (reviewer point 1:
// no FD point silently re-normalizes by the current M).
if (thermalObjectiveType == "maximizeTotalHeatTransfer")
{
    const scalarField& phiOutlet =
        phi.boundaryField()[coldOutletPatchID];
    thermalObjectiveMassFlowFrozen = sum(phiOutlet);
    reduce(thermalObjectiveMassFlowFrozen, sumOp<scalar>());

    if (mag(thermalObjectiveMassFlowFrozen) < SMALL)
    {
        FatalErrorInFunction
            << "Cold outlet mass flow is near zero ("
            << thermalObjectiveMassFlowFrozen
            << "). Cannot compute total-heat-transfer objective."
            << abort(FatalError);
    }
}

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
// BFINAL-025: Q = N - T_in*M, J_Q = -Q/(Tref*M_frozen).
//   dJ_Q/dT[celli] = -phi_f/(Tref*M_frozen)   (cells at cold outlet)
// Same form as the old type with M -> M_frozen (snapshot, constant in state).
else if (thermalObjectiveType == "maximizeTotalHeatTransfer")
{
    const scalarField& phiOutlet =
        phi.boundaryField()[coldOutletPatchID];
    const labelUList& outletCells =
        mesh.boundary()[coldOutletPatchID].faceCells();

    forAll(outletCells, facei)
    {
        thermalObjectiveDerivative[outletCells[facei]] -=
            phiOutlet[facei]
           /Foam::max
            (
                thermalObjectiveMassFlowFrozen*thermalReferenceTemperature,
                SMALL
            );
    }
}
// -------------------------------------------------------------------------
// Thermal objective flux derivative dJ/dphi on the cold outlet boundary.
// J = -N/(M*Tref) with M=sum(phi_out), N=sum(phi_out*T_out), Tmix=N/M.
//   dJ/dphi_f = -(T_f - Tmix)/(M*Tref)
// This is the objective's EXPLICIT dependence on the outlet face flux phi.
// It feeds the discrete flow adjoint through discreteExternalFaceFluxAdjoint
// (AdjNS_HT.H); the thermal-residual boundary flux term is zero for outflow
// zeroGradient T under the fvm discretization (verified by B0.1 FD tests),
// so the combined flow-adjoint boundary flux RHS is exactly dJ/dphi.
// -------------------------------------------------------------------------
thermalObjectiveFluxDerivative.primitiveFieldRef() = 0.0;
forAll(mesh.boundary(), patchi)
{
    thermalObjectiveFluxDerivative.boundaryFieldRef()[patchi] = 0.0;
}
if (thermalObjectiveType == "maximizeColdOutletTemperature")
{
    const scalarField& phiOutlet =
        phi.boundaryField()[coldOutletPatchID];
    const scalarField TOutlet =
        T.boundaryField()[coldOutletPatchID].patchInternalField();
    scalar M = sum(phiOutlet);
    scalar N = sum(phiOutlet*TOutlet);
    reduce(M, sumOp<scalar>());
    reduce(N, sumOp<scalar>());
    if (mag(M) < SMALL)
    {
        FatalErrorInFunction
            << "Cold outlet mass flow is near zero (" << M
            << "). Cannot compute mass-flow-weighted objective flux derivative."
            << abort(FatalError);
    }
    const scalar Tmix = N/M;
    const labelUList& outletCells =
        mesh.boundary()[coldOutletPatchID].faceCells();
    forAll(outletCells, facei)
    {
        thermalObjectiveFluxDerivative.boundaryFieldRef()[coldOutletPatchID][facei] =
           -(TOutlet[facei] - Tmix)
           /Foam::max(M*thermalReferenceTemperature, SMALL);
    }
}
// BFINAL-025: Q = N - T_in*M, J_Q = -Q/(Tref*M_frozen).
//   dJ_Q/dphi_f = -(T_f - T_in)/(Tref*M_frozen)
// The old type's near-cancellation (T_f - Tmix) ~ 0 is replaced by the
// full temperature rise (T_f - T_in) ~ O(heating), which is non-small.
// T_in is read from the cold-inlet case BC (gAverage of the fixedValue
// boundary field), never hardcoded.
else if (thermalObjectiveType == "maximizeTotalHeatTransfer")
{
    const scalarField& phiOutlet =
        phi.boundaryField()[coldOutletPatchID];
    const scalarField TOutlet =
        T.boundaryField()[coldOutletPatchID].patchInternalField();
    const scalar TIn =
        gAverage(T.boundaryField()[coldInletPatchID]);
    const labelUList& outletCells =
        mesh.boundary()[coldOutletPatchID].faceCells();
    forAll(outletCells, facei)
    {
        thermalObjectiveFluxDerivative.boundaryFieldRef()[coldOutletPatchID][facei] =
           -(TOutlet[facei] - TIn)
           /Foam::max
            (
                thermalObjectiveMassFlowFrozen*thermalReferenceTemperature,
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
