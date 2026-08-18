//-----------------------------------------------------------------------------
// validateMmaUnlockGate.H
// Hard safety gate for the frozen-turbulence optimization route.
//
// Frozen-gradient amplitude validation (Stage B2/B3) and full-SST directional
// consistency are deliberately separate concepts. Production MMA is unlocked
// by an explicit case acknowledgement that the frozen-model amplitude gate has
// passed. Until full-SST consistency is also validated, every topology update
// must use the full-SST candidate acceptance/rollback layer.
//-----------------------------------------------------------------------------

const Switch frozenGradientValidated =
    optProperties.lookupOrDefault<Switch>
    (
        "frozenGradientValidated",
        false
    );
const Switch stageB4JacobianProbeEnabled =
    optProperties.lookupOrDefault<Switch>("stageB4JacobianProbe", false);
const Switch discreteUseExplicitSolutionEnabled =
    optProperties.lookupOrDefault<Switch>("discreteUseExplicitSolution", false);
const Switch discreteExplicitSolutionIncludesDeviatoric =
    optProperties.lookupOrDefault<Switch>
    (
        "discreteExplicitSolutionIncludesDeviatoric",
        false
    );

if
(
    mmaUpdateEnabled
 && (stageBEnabled || stageB2Enabled || stageB4JacobianProbeEnabled)
)
{
    FatalErrorInFunction
        << "MMA updates are forbidden while Stage-B validation/probe modes are "
        << "enabled. Run validation with mmaUpdateEnabled=false."
        << abort(FatalError);
}

if (stageB2Enabled && stageB4JacobianProbeEnabled)
{
    FatalErrorInFunction
        << "stageB2Enabled and stageB4JacobianProbe must not be enabled in the "
        << "same run. B2/B3 uses the production adjoint solver; B4 is the "
        << "diagnostic/oracle path."
        << abort(FatalError);
}

// In frozen RANS the exported explicitJT.mtx from the legacy B4 diagnostic
// path is a core CSR matrix plus a matrix-free deviatoric-stress transpose.
// A solution obtained from that file alone is not a full J^T solution.
if
(
    ransFlowModel
 && stageB4JacobianProbeEnabled
 && discreteUseExplicitSolutionEnabled
 && !discreteExplicitSolutionIncludesDeviatoric
)
{
    FatalErrorInFunction
        << "discreteUseExplicitSolution=true is unsafe for frozen RANS unless "
        << "the external solution was built with the FULL J^T, including the "
        << "matrix-free deviatoric-stress transpose. Set "
        << "discreteExplicitSolutionIncludesDeviatoric=true only for a "
        << "verified full-operator external solve."
        << abort(FatalError);
}

if (mmaUpdateEnabled && !frozenGradientValidated)
{
    FatalErrorInFunction
        << "MMA is blocked until frozen-turbulence Stage B2/B3 gradient "
        << "directions and amplitudes pass on the current source/configuration. "
        << "Keep mmaUpdateEnabled=false during validation; only set "
        << "frozenGradientValidated=true after reviewing a current-head PASS "
        << "report. Direction-only approval cannot bypass this gate."
        << abort(FatalError);
}

// A validated frozen gradient is sufficient to define a trustworthy search
// direction for the frozen model, but not to assume that the omitted SST
// derivatives are negligible. Until the strict full-SST consistency gate has
// passed, require transactional full-SST candidate validation on every update.
if
(
    mmaUpdateEnabled
 && ransFlowModel
 && frozenGradientValidated
 && !gradientValidated
 && !enableSSTAcceptanceCheck
)
{
    FatalErrorInFunction
        << "Frozen gradients are validated, but full-SST consistency is not. "
        << "Enable enableSSTAcceptanceCheck=true so each MMA candidate is "
        << "re-converged with full SST and can be rejected/rolled back."
        << abort(FatalError);
}

Info<< "Frozen-gradient MMA gate: frozenGradientValidated="
    << frozenGradientValidated
    << ", fullSSTGradientValidated=" << gradientValidated
    << ", SSTcandidateAcceptance=" << enableSSTAcceptanceCheck
    << ", mmaUpdateEnabled=" << mmaUpdateEnabled
    << ", stageB4JacobianProbe=" << stageB4JacobianProbeEnabled
    << endl;
