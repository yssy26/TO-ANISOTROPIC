if (adjointMode == "discrete")
{
    vectorField discreteFlowVelocityRhs(mesh.nCells(), vector::zero);
    scalarField discreteFlowPressureRhs(pressureConstraintDerivative);
    scalarField discreteExternalFaceFluxAdjoint(mesh.nFaces(), 0.0);
    volVectorField& discreteAdjointVelocity = Uc;
    volScalarField& discreteAdjointPressure = pc;
    const word discreteAdjointLabel("pressureDrop");
    #include "solveDiscreteFlowAdjoint.H"
}
else
{volScalarField pressureAdjointContinuitySource
(
    IOobject
    (
        "pressureAdjointContinuitySource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(phic)*scalar(0.0)
);

forAll(pressureConstraintDerivative, celli)
{
    pressureAdjointContinuitySource[celli] =
        pressureConstraintDerivative[celli]
       /Foam::max(mesh.V()[celli], SMALL);
}
pressureAdjointContinuitySource.correctBoundaryConditions();
for(i=0;i<adjointInnerIter;i++)
{
  volVectorField adjointTransposeConvectionc((Uc & fvc::grad(U)));

  tmp<fvVectorMatrix> tUcEqn
  (
      fvm::div(-phi, Uc)
    + adjointTransposeConvectionc
    - fvm::laplacian(nuEffFrozen,Uc)
    + fvm::Sp(alpha, Uc)
   ==
      fvOptions(Uc)
  );
  fvVectorMatrix& UcEqn = tUcEqn.ref();

  UcEqn.relax();

  fvOptions.constrain(UcEqn);
  Info << "AdjNS_PD.H - Attempting to solve UcEqn for Uc (iter " << i << ")." << endl;
  solve(UcEqn == -fvc::grad(pc));
  Info << "AdjNS_PD.H - UcEqn for Uc solved." << endl;
  fvOptions.correct(Uc);

  volScalarField rAUc(1.0/UcEqn.A());
  Info << "AdjNS_PD.H - Min/Max rAUc: " << gMin(rAUc) << " / " << gMax(rAUc) << endl;

  if (gMin(rAUc) <= 0) {
      FatalErrorInFunction << "rAUc contains non-positive values!" << abort(FatalError);
  }

  volVectorField HbyAc("HbyAc", Uc);
  HbyAc = rAUc*UcEqn.H();

  surfaceScalarField phiHbyAc("phiHbyAc", fvc::flux(HbyAc));
  adjustPhi(phiHbyAc, Uc, pc);

  tmp<volScalarField> rAtUc(rAUc);

  if (simple.consistent())
  {
      rAtUc = 1.0/(1.0/rAUc - UcEqn.H1());
      phiHbyAc +=
          fvc::interpolate(rAtUc() - rAUc)*fvc::snGrad(pc)*mesh.magSf();
      HbyAc -= (rAUc - rAtUc())*fvc::grad(pc);
  }

  tUcEqn.clear();

  constrainPressure(pc, Uc, phiHbyAc, rAtUc());

  while (simple.correctNonOrthogonal())
  {
      Info << "AdjNS_PD.H - Iter: " << i << ", INSIDE while loop for pc solve." << endl;
      fvScalarMatrix pcEqn
      (
          fvm::laplacian(rAtUc(), pc)
       == fvc::div(phiHbyAc) + pressureAdjointContinuitySource
      );

      pcEqn.setReference(pcRefCell, pcRefValue);
      Info << "AdjNS_PD.H: Attempting to solve pcEqn for pc. Outer iter: " << opt << ", Inner iter: " << i << endl;
      pcEqn.solve();
      Info << "AdjNS_PD.H: pcEqn for pc solved. Outer iter: " << opt << ", Inner iter: " << i << endl;

      if (simple.finalNonOrthogonalIter())
      {
          phic = phiHbyAc - pcEqn.flux();
      }
  }

  pc.relax();

  Uc = HbyAc - rAtUc()*fvc::grad(pc);
  Uc.correctBoundaryConditions();
  fvOptions.correct(Uc);

  Uc.storePrevIter();
  pc.storePrevIter();
  phic.storePrevIter();
}

}
