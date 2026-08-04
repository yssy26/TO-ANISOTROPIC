//-----------------------------------------------------------------------------
// saveAcceptedState.H
// Phase I: Snapshot the current accepted design + primal state into backup
// storage so that a rejected candidate can be fully rolled back.
//
// Must be included AFTER the forward solve (NS.H + HeatTransfer.H) and
// AFTER costfunction.H so that Objective, V, PressureDropPa, MeanT reflect
// the accepted design.
//
// The backup fields are declared in opt_initialization.H and are only
// allocated when enableSSTAcceptanceCheck is true.
//-----------------------------------------------------------------------------

// --- Design variables -------------------------------------------------------
xTxBackup = x.primitiveField();
xpTxBackup = xp.primitiveField();
xhTxBackup = xh.primitiveField();
eta5TxBackup = eta5;

// --- Material properties ----------------------------------------------------
alphaTxBackup = alpha.primitiveField();
DTTxBackup = DTMolecular.primitiveField();
DTEffectiveTxBackup = DTEffective.primitiveField();

// --- Flow fields ------------------------------------------------------------
UTxBackup = U.primitiveField();
pTxBackup = p.primitiveField();
phiTxBackup = phi;

// --- Temperature ------------------------------------------------------------
TTxBackup = T.primitiveField();

// --- Turbulence (managed by the turbulence model) --------------------------
if (ransFlowModel)
{
    kTxBackup =
        mesh.lookupObject<volScalarField>("k").primitiveField();
    omegaTxBackup =
        mesh.lookupObject<volScalarField>("omega").primitiveField();
    nutTxBackup =
        mesh.lookupObject<volScalarField>("nut").primitiveField();
}

// --- Frozen turbulence derived fields --------------------------------------
nutFrozenTxBackup = nutFrozen.primitiveField();
nuEffFrozenTxBackup = nuEffFrozen.primitiveField();
alphaTurbulentTxBackup = alphaTurbulent.primitiveField();

// --- Objective and constraint values ---------------------------------------
objectiveTxBackup = Objective;
volumeConstraintTxBackup = V;
pressureDropTxBackup =
    PressureDropPa/Foam::max(pressureDropMaxPa, SMALL) - 1.0;
pressureDropPaTxBackup = PressureDropPa;
meanTTxBackup = MeanT;
powerDissTxBackup = PowerDiss;

// --- Solver flags -----------------------------------------------------------
freezeTurbulenceForValidationTxBackup = freezeTurbulenceForValidation;
