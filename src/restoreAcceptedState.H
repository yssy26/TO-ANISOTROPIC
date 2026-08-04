//-----------------------------------------------------------------------------
// restoreAcceptedState.H
// Phase I: Restore the accepted design + primal state from backup storage.
//
// Used when a candidate design is rejected by the SST acceptance check.
// After restoration the fields are identical to the last accepted state,
// allowing the optimizer to retry with a reduced move limit or simply
// proceed to the next outer iteration with the unchanged accepted design.
//-----------------------------------------------------------------------------

// --- Design variables -------------------------------------------------------
x.primitiveFieldRef() = xTxBackup;
x.correctBoundaryConditions();
xp.primitiveFieldRef() = xpTxBackup;
xp.correctBoundaryConditions();
xh.primitiveFieldRef() = xhTxBackup;
xh.correctBoundaryConditions();
eta5 = eta5TxBackup;

// --- Material properties ----------------------------------------------------
alpha.primitiveFieldRef() = alphaTxBackup;
alpha.correctBoundaryConditions();
DTMolecular.primitiveFieldRef() = DTTxBackup;
DTMolecular.correctBoundaryConditions();
DTEffective.primitiveFieldRef() = DTEffectiveTxBackup;
DTEffective.correctBoundaryConditions();

// --- Flow fields ------------------------------------------------------------
U.primitiveFieldRef() = UTxBackup;
U.correctBoundaryConditions();
p.primitiveFieldRef() = pTxBackup;
p.correctBoundaryConditions();
phi = phiTxBackup;

// --- Temperature ------------------------------------------------------------
T.primitiveFieldRef() = TTxBackup;
T.correctBoundaryConditions();

// --- Turbulence (managed by the turbulence model) --------------------------
// Only restore internal fields.  correctBoundaryConditions() is intentionally
// skipped for k, omega, nut because omegaWallFunction::updateCoeffs() looks up
// kOmegaSST:G from the objectRegistry, which is only guaranteed to exist
// during turbulence->correct().  Boundary conditions will be refreshed by the
// next turbulence->correct() call in evaluateCandidate.H or the next outer
// iteration's NS.H.
if (ransFlowModel)
{
    mesh.lookupObjectRef<volScalarField>("k").primitiveFieldRef() =
        kTxBackup;
    mesh.lookupObjectRef<volScalarField>("omega").primitiveFieldRef() =
        omegaTxBackup;
    mesh.lookupObjectRef<volScalarField>("nut").primitiveFieldRef() =
        nutTxBackup;
}

// --- Frozen turbulence derived fields --------------------------------------
nutFrozen.primitiveFieldRef() = nutFrozenTxBackup;
nutFrozen.correctBoundaryConditions();
nuEffFrozen.primitiveFieldRef() = nuEffFrozenTxBackup;
nuEffFrozen.correctBoundaryConditions();
alphaTurbulent.primitiveFieldRef() = alphaTurbulentTxBackup;
alphaTurbulent.correctBoundaryConditions();

// --- Objective and constraint values ---------------------------------------
Objective = objectiveTxBackup;
V = volumeConstraintTxBackup;
PressureDropPa = pressureDropPaTxBackup;
MeanT = meanTTxBackup;
PowerDiss = powerDissTxBackup;

// --- Solver flags -----------------------------------------------------------
freezeTurbulenceForValidation = freezeTurbulenceForValidationTxBackup;
