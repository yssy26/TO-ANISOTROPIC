Info<< "Initialization\n" << endl;
std::remove("Voluse.txt");
std::remove("Time.txt");
std::remove("MeanT.txt");
std::remove("PowerDiss.txt");

int i = 0;
int opt = runTime.timeIndex();
int rex;
double del = 0.0;
scalar eta0 = 0.0;
scalar eta5 = 0.5;
scalar eta1 = 1.0;
scalar y0 = 0.0;
scalar y5 = 0.0;
bool projectionThresholdAdaptive = false;

// n remains the number of local cells for field/filter loops.
const label n = x.size();

DynamicList<label> activeDesignCellBuffer;
forAll(designMask.primitiveField(), celli)
{
    if (designMask[celli] > 0.5)
    {
        activeDesignCellBuffer.append(celli);
    }
}

labelList activeDesignCells(activeDesignCellBuffer.size());
forAll(activeDesignCellBuffer, activeI)
{
    activeDesignCells[activeI] = activeDesignCellBuffer[activeI];
}

const label nDesignLocal = activeDesignCells.size();
label nDesignGlobal = nDesignLocal;
reduce(nDesignGlobal, sumOp<label>());

if (nDesignGlobal <= 0)
{
    FatalErrorInFunction
        << "No active design cells are available for MMA."
        << abort(FatalError);
}

scalar V = 0.0;
scalar MeanT = 0.0;
scalar PowerDiss = 0.0;
scalar Objective = 0.0;
scalar PressureDropPa = 0.0;

// Forward solver residual diagnostics (captured from the last corrector
// iteration of NS.H and the last solve of HeatTransfer.H).
scalar forwardInitialResidualU = 0.0;
scalar forwardFinalResidualU = 0.0;
scalar forwardInitialResidualP = 0.0;
scalar forwardFinalResidualP = 0.0;
scalar forwardInitialResidualT = 0.0;
scalar forwardFinalResidualT = 0.0;
label forwardPrimalOuterIterations = 0;

// ---------------------------------------------------------------------------
// Phase K: CSV logging diagnostics
// ---------------------------------------------------------------------------
// SST relative changes (exposed from NS.H for the last corrector iteration)
scalar lastRelativeChangeK = 0.0;
scalar lastRelativeChangeOmega = 0.0;

// Adjoint solver residuals (captured from AdjHeatTransfer.H and AdjNS_HT.H)
scalar adjointResidualT = 0.0;
scalar adjointResidualU = 0.0;
scalar adjointResidualP = 0.0;

// Transaction tracking for CSV rows
label csvAcceptedIteration = 0;   // last accepted iteration number
label csvCandidateTrial = 0;      // 0=accepted, 1+=candidate attempt
label csvAcceptedFlag = 1;        // 1 if this row is accepted, 0 if rejected
word rejectionReason("");         // empty if accepted

// CSV file path (computed once, reused for every write)
const fileName csvLogPath = runTime.path()/"optimization_history.csv";

const double len = std::pow
(
    designVolume/Foam::max(scalar(nDesignGlobal), SMALL),
    1.0/3.0
);
b.value() =
    1.0
   /Foam::max
    (
        sqr(filterR*len/3.464),
        SMALL
    );

Info<< "Active MMA design cells: local=" << nDesignLocal
    << ", global=" << nDesignGlobal
    << ", designVolume=" << designVolume << endl;

// MMA uses active design variables only. Fixed hot-fluid, separator-wall and
// legacy prescribed zones never enter the MMA subproblem.
const int m = 2;
MMA mma(nDesignLocal, m);
mma.asyminit = mma_init;
mma.asymdec = mma_dec;
mma.asyminc = mma_inc;
mma.raa0 = raa0;
mma.movelimit = movlim;

std::vector<double> xmma(nDesignLocal);
std::vector<double> dfdx(nDesignLocal);
std::vector<double> gx(m);
std::vector<std::vector<double>> dgdx(m);
for (int constraintI = 0; constraintI < m; ++constraintI)
{
    dgdx[constraintI].resize(nDesignLocal);
}

// ---------------------------------------------------------------------------
// Phase I: Transaction mechanism backup storage.
// Only used when enableSSTAcceptanceCheck is true.  Declared unconditionally
// so that saveAcceptedState.H / restoreAcceptedState.H can always compile.
// ---------------------------------------------------------------------------
scalarField xTxBackup(n);
scalarField xpTxBackup(n);
scalarField xhTxBackup(n);
scalarField alphaTxBackup(n);
Field<symmTensor> DTTxBackup(n);
Field<symmTensor> DTEffectiveTxBackup(n);
vectorField UTxBackup(n);
scalarField pTxBackup(n);
surfaceScalarField phiTxBackup
(
    IOobject
    (
        "phiTxBackup",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi
);
scalarField TTxBackup(n);
scalarField kTxBackup;
scalarField omegaTxBackup;
scalarField nutTxBackup;
scalarField nutFrozenTxBackup(n);
scalarField nuEffFrozenTxBackup(n);
scalarField alphaTurbulentTxBackup(n);
scalar eta5TxBackup = 0.0;
scalar objectiveTxBackup = 0.0;
scalar volumeConstraintTxBackup = 0.0;
scalar pressureDropTxBackup = 0.0;
scalar pressureDropPaTxBackup = 0.0;
scalar meanTTxBackup = 0.0;
scalar powerDissTxBackup = 0.0;
Switch freezeTurbulenceForValidationTxBackup = false;

if (!mmaUpdateEnabled && recomputeDesignMappingWhenMmaDisabled)
{
    Info<< "Validation mode: recomputing xp/xh from raw x with fixed "
        << "projection beta=" << validationProjectionBeta << endl;
    #include "filter_x.H"
}

// --- Phase J: Load optimizer state for restart ---
// Must come after MMA object creation and all state variables are declared.
#include "loadOptimizerState.H"
