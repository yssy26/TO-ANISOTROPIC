if (adjointMode == "discrete")
{
    vectorField discreteFlowVelocityRhs(mesh.nCells(), vector::zero);
    scalarField discreteFlowPressureRhs(mesh.nCells(), 0.0);
    scalarField discreteExternalFaceFluxAdjoint(mesh.nFaces(), 0.0);
    const labelUList& thermalFlowOwner = mesh.owner();
    const labelUList& thermalFlowNeighbour = mesh.neighbour();
    const scalarField& activeThermalPhi = phiThermal.primitiveField();
    const scalarField& coldThermalFaceMask = coldFaceMask.primitiveField();

    forAll(thermalFlowOwner, facei)
    {
        if (coldThermalFaceMask[facei] <= 0.0)
        {
            continue;
        }

        const label own = thermalFlowOwner[facei];
        const label nei = thermalFlowNeighbour[facei];
        const scalar temperatureJump = T[nei] - T[own];
        const scalar adjointDownwind =
            activeThermalPhi[facei] >= 0.0 ? Tb[nei] : Tb[own];
        discreteExternalFaceFluxAdjoint[facei] =
           -coldThermalFaceMask[facei]
           *adjointDownwind*temperatureJump;
    }
    volVectorField& discreteAdjointVelocity = Ub;
    volScalarField& discreteAdjointPressure = pb;
    const word discreteAdjointLabel("thermalCoupling");
    #include "solveDiscreteFlowAdjoint.H"

    // Phase K: capture discrete adjoint residual for CSV logging
    adjointResidualU = discreteRelativeResidual;
    adjointResidualP = 0.0;  // solved jointly in discrete mode
}
else
{for (i = 0; i < adjointInnerIter; ++i)
{
    volVectorField adjointTransposeConvectionb
    (
        (Ub & fvc::grad(U))
    );

    // Only the solved cold-side velocity depends on topology. The prescribed
    // hot-side velocity has no momentum adjoint. thermalObjectiveMask is one
    // in the legacy model and equals coldRegionMask in frozen-hot mode.
    //
    // Coupling source from the thermal adjoint: the primal thermal convection
    // term div(phi, T) depends on U through phi.  Its adjoint w.r.t. U gives
    // a source -T*grad(Tb) (derived via integration by parts of the weak
    // form <Tb, div(u*T)>).  Placing +T*grad(Tb) on the LHS yields the
    // correct effective RHS source of -T*grad(Tb).
    tmp<fvVectorMatrix> tUbEqn
    (
        fvm::div(-phi, Ub)
      + adjointTransposeConvectionb
      - fvm::laplacian(nuEffFrozen, Ub)
      + fvm::Sp(alpha, Ub)
      + thermalObjectiveMask*T*fvc::grad(Tb)
     ==
        fvOptions(Ub)
    );
    fvVectorMatrix& UbEqn = tUbEqn.ref();

    UbEqn.relax();
    fvOptions.constrain(UbEqn);
    {
        SolverPerformance<vector> perfUb =
            solve(UbEqn == -fvc::grad(pb));
        adjointResidualU = mag(perfUb.finalResidual());
    }
    fvOptions.correct(Ub);

    volScalarField rAUb(1.0/UbEqn.A());
    if (gMin(rAUb) <= 0)
    {
        FatalErrorInFunction
            << "rAUb contains non-positive values."
            << abort(FatalError);
    }

    volVectorField HbyAb("HbyAb", Ub);
    HbyAb = rAUb*UbEqn.H();
    surfaceScalarField phiHbyAb("phiHbyAb", fvc::flux(HbyAb));
    adjustPhi(phiHbyAb, Ub, pb);

    tmp<volScalarField> rAtUb(rAUb);
    if (simple.consistent())
    {
        rAtUb = 1.0/(1.0/rAUb - UbEqn.H1());
        phiHbyAb +=
            fvc::interpolate(rAtUb() - rAUb)
           *fvc::snGrad(pb)*mesh.magSf();
        HbyAb -= (rAUb - rAtUb())*fvc::grad(pb);
    }

    tUbEqn.clear();
    constrainPressure(pb, Ub, phiHbyAb, rAtUb());

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pbEqn
        (
            fvm::laplacian(rAtUb(), pb) == fvc::div(phiHbyAb)
        );

        pbEqn.setReference(pbRefCell, pbRefValue);
        {
            SolverPerformance<scalar> perfPb = pbEqn.solve();
            adjointResidualP = perfPb.finalResidual();
        }

        if (simple.finalNonOrthogonalIter())
        {
            phib = phiHbyAb - pbEqn.flux();
        }
    }

    pb.relax();
    Ub = HbyAb - rAtUb()*fvc::grad(pb);
    Ub.correctBoundaryConditions();
    fvOptions.correct(Ub);

    Ub.storePrevIter();
    pb.storePrevIter();
    phib.storePrevIter();
}

}
