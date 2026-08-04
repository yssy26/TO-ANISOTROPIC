for(i=0;i<adjointInnerIter;i++)
{ 
  volVectorField adjointTransposeConvectiona((Ua & fvc::grad(U)));

  tmp<fvVectorMatrix> tUaEqn
  (
      fvm::div(-phi, Ua)
    + adjointTransposeConvectiona
    - fvm::laplacian(nuEffFrozen,Ua)
    + fvm::Sp(alpha, Ua)
   ==
      fvOptions(Ua)
  );
  fvVectorMatrix& UaEqn = tUaEqn.ref();

  UaEqn.relax();

  fvOptions.constrain(UaEqn);
Info << "AdjNS_FF.H - Attempting to solve UaEqn for Ua (iter " << i << ")." << endl;
  solve(UaEqn == -fvc::grad(pa));
Info << "AdjNS_FF.H - UaEqn for Ua solved. "<< endl;
  fvOptions.correct(Ua);
//****************************************
 
volScalarField rAUa(1.0/UaEqn.A());
Info << "Min/Max rAUa: " << gMin(rAUa) << " / " << gMax(rAUa) << endl;
// 检查 rAUa 是否有非正值
if (gMin(rAUa) <= 0) {
    FatalErrorInFunction << "rAUa contains non-positive values!" << abort(FatalError);
}
  volVectorField HbyAa("HbyAa", Ua);
  HbyAa = rAUa*UaEqn.H();
       
  surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
  adjustPhi(phiHbyAa, Ua, pa);
  
  tmp<volScalarField> rAtUa(rAUa);
  
   if (simple.consistent())
  {
      rAtUa = 1.0/(1.0/rAUa - UaEqn.H1());
      phiHbyAa +=
          fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh.magSf();
      HbyAa -= (rAUa - rAtUa())*fvc::grad(pa);
  }

  tUaEqn.clear();

  // Update the pressure BCs to ensure flux consistency
  constrainPressure(pa, Ua, phiHbyAa, rAtUa());           
  // Non-orthogonal pressure corrector loop
  while (simple.correctNonOrthogonal())
  {
 Info << "AdjNS_FF.H - Iter: " << i << ", INSIDE while loop for pa solve." << endl;
      fvScalarMatrix paEqn
      (
          fvm::laplacian(rAtUa(), pa) == fvc::div(phiHbyAa)
      );

      paEqn.setReference(paRefCell, paRefValue);
Info << "AdjNS_FF.H: Attempting to solve paEqn for pa. Outer iter: " << opt << ", Inner iter: " << i << ", NonOrth iter" << endl;
      paEqn.solve();
Info << "AdjNS_FF.H: paEqn for pa solved. Outer iter: " << opt << ", Inner iter: " << i << ", NonOrth iter" << endl;
      if (simple.finalNonOrthogonalIter())
      {
          phia = phiHbyAa - paEqn.flux();
      }
  }
Info << "AdjNS_FF.H - Iter: " << i << ", After while loop." << endl;
  //#include "adjointContinuityErrs.H"

  // Explicitly relax pressure for adjoint momentum corrector
  pa.relax();

  // Adjoint momentum corrector
  Ua = HbyAa - rAtUa()*fvc::grad(pa);
  Ua.correctBoundaryConditions();
  fvOptions.correct(Ua);
  
  Ua.storePrevIter();
  pa.storePrevIter();
  phia.storePrevIter();
}
