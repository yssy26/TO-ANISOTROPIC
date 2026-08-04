// The legacy model uses phi. Frozen-hot mode uses the conservative combined
// cold/hot thermal flux. Reuse the existing div(phi,T) scheme entry so old
// cases do not require an additional div(phiThermal,T) dictionary key.
const surfaceScalarField& activeThermalFlux =
    frozenHotRegionEnabled ? phiThermal : phi;

fvScalarMatrix TEqn
(
    fvm::div(activeThermalFlux, T, "div(phi,T)")
  - fvm::laplacian(DTEffective, T)
 ==
    Q
);
SolverPerformance<scalar> perfT = TEqn.solve();

forwardInitialResidualT = perfT.initialResidual();
forwardFinalResidualT = perfT.finalResidual();

Info<< "TEqn: initialResidual=" << perfT.initialResidual()
    << ", finalResidual=" << perfT.finalResidual()
    << ", nIterations=" << perfT.nIterations() << endl;
