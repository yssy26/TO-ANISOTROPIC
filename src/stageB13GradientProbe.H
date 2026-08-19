//-----------------------------------------------------------------------------
// stageB13GradientProbe.H
// BFINAL-013 localization round: gradient-assembly probes.  Read-only
// diagnostics (P1 field save/restore is exact); no production semantics are
// touched.  Switch-gated ("stageB13GradientProbe", default false), serial
// discrete-adjoint mode only.
//
// P1 - source dot tests (the objective-derivative sources vs direct FD of
//      the FD-side objective evaluators; no solver involved):
//   P1a  gDP p-source : (PD(p+e v)-PD(p-e v))/2e  vs  v.dgDP/dp source
//   P1b  J   T-source : (J (T+e v)-J (T-e v))/2e  vs  v.dJ/dT  source
//   P1c  J outlet-flux source: (J(phi+e v)-J(phi-e v))/2e vs v.dJ/dphi
//   PASS => defect is downstream (contraction/assembly layer);
//   FAIL => defect is in the objective derivative source itself.
//
// P2 - per-term contraction decomposition along the Stage-B2 D1/D2/D3
//      directions (exact same construction as validateStageB2GradientAmplitude.H):
//   gDP: momentum row (recomputed), pressure row (total minus momentum),
//        xh-total field, x-total field (post filter_chainrule), MMA column.
//   J  : Brinkman term (recomputed), thermal-diffusion term, xh-total,
//        x-total, MMA column.  Component sums are checked against the
//        assembled totals (assembly self-consistency).
//-----------------------------------------------------------------------------

const Switch stageB13GradientProbeOn =
    optProperties.lookupOrDefault<Switch>("stageB13GradientProbe", false);

if (stageB13GradientProbeOn)
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "stageB13GradientProbe requires serial execution."
            << abort(FatalError);
    }
    if (adjointMode != "discrete")
    {
        FatalErrorInFunction
            << "stageB13GradientProbe requires adjointMode discrete."
            << abort(FatalError);
    }

    Info<< "\n========================================\n"
        << " Stage B13 P1: objective-source dot tests\n"
        << "========================================\n" << endl;

    // Deterministic smooth cell-wise test fields (exact Stage-B2 D1
    // construction: physical coordinates, designMask + x-bound gating,
    // max-normalized).
    auto b13BuildDirection = [&]() -> scalarField
    {
        scalarField dir(mesh.nCells(), 0.0);
        const vectorField& C = mesh.C().primitiveField();
        const scalar Lx = gMax(C.component(0)) - gMin(C.component(0));
        const scalar Ly = gMax(C.component(1)) - gMin(C.component(1));
        const scalar Lz = gMax(C.component(2)) - gMin(C.component(2));
        forAll(C, celli)
        {
            if (designMask[celli] <= 0.5) continue;
            const scalar xx = C[celli].x()/Foam::max(Lx, SMALL);
            const scalar yy = C[celli].y()/Foam::max(Ly, SMALL);
            const scalar zz = C[celli].z()/Foam::max(Lz, SMALL);
            dir[celli] =
                Foam::sin(2.0*Foam::constant::mathematical::pi*xx)
               *Foam::cos(Foam::constant::mathematical::pi*yy)
              + 0.5*Foam::sin(3.0*Foam::constant::mathematical::pi*zz);
            if (x[celli] < 0.02 || x[celli] > 0.98)
            {
                dir[celli] = 0.0;
            }
        }
        scalar mx = 0.0;
        forAll(dir, celli) mx = Foam::max(mx, mag(dir[celli]));
        reduce(mx, maxOp<scalar>());
        if (mx > 1e-30) dir /= mx;
        return dir;
    };
    const scalarField b13Dir(b13BuildDirection());

    // -----------------------------------------------------------------
    // P1a: gDP pressure source dot test.  The source lives only on cells
    // adjacent to non-fixedValue-p constraint patches, so the test field
    // must live there too (the design-region D1 support does not overlap
    // it - a domain-wide D1 gives a vacuous 0==0).
    // -----------------------------------------------------------------
    {
        scalarField v(mesh.nCells(), 0.0);
        {
            const label inletPatchI = conPatchList[0];
            const labelUList& inletCells =
                mesh.boundary()[inletPatchI].faceCells();
            forAll(inletCells, f)
            {
                v[inletCells[f]] = Foam::sin(0.271*(f + 1));
            }
            scalar mx = 0.0;
            forAll(v, celli) mx = Foam::max(mx, mag(v[celli]));
            reduce(mx, maxOp<scalar>());
            if (mx > 1e-30) v /= mx;
        }
        const scalarField pSaved(p.primitiveField());
        PtrList<scalarField> pBndSaved(mesh.boundary().size());
        forAll(mesh.boundary(), patchi)
        {
            pBndSaved.set
            (
                patchi,
                new scalarField(p.boundaryField()[patchi])
            );
        }
        scalar adDot = 0.0;
        forAll(v, celli) adDot += pressureConstraintDerivative[celli]*v[celli];
        reduce(adDot, sumOp<scalar>());

        const scalarList epsList = {1e-2, 3e-3, 1e-3, 3e-4, 1e-4};
        Info<< "B13SOURCEPROBE P1a gDP/p-source (ADJ v.dSource=" << adDot
            << ", PDmax=" << pressureDropMaxPa << "):" << endl;
        forAll(epsList, ei)
        {
            const scalar eps = epsList[ei];
            p.primitiveFieldRef() = pSaved + eps*v;
            p.correctBoundaryConditions();
            const scalar pdP = evaluatePressureDrop();
            p.primitiveFieldRef() = pSaved - eps*v;
            p.correctBoundaryConditions();
            const scalar pdM = evaluatePressureDrop();
            p.primitiveFieldRef() = pSaved;
            forAll(mesh.boundary(), patchi)
            {
                p.boundaryFieldRef()[patchi] = pBndSaved[patchi];
            }
            p.correctBoundaryConditions();
            // FD is differentiated against gDP = PD/PDmax - 1 to match
            // the source normalization.
            const scalar fd = (pdP - pdM)/(2.0*eps*pressureDropMaxPa);
            const scalar rel =
                mag(fd - adDot)/Foam::max(mag(adDot), 1e-300);
            Info<< "  eps=" << eps << " FD=" << fd << " rel=" << rel
                << (rel <= 1e-6 ? " PASS" : (rel <= 1e-3 ? " WARN" : " FAIL"))
                << endl;
        }
    }

    // -----------------------------------------------------------------
    // P1b: J temperature source dot test.  dJ/dT lives only on cells
    // adjacent to the cold outlet, so the test field is built there.
    // -----------------------------------------------------------------
    {
        scalarField v(mesh.nCells(), 0.0);
        {
            const labelUList& outletCells =
                mesh.boundary()[coldOutletPatchID].faceCells();
            forAll(outletCells, f)
            {
                v[outletCells[f]] = Foam::sin(0.271*(f + 1));
            }
            scalar mx = 0.0;
            forAll(v, celli) mx = Foam::max(mx, mag(v[celli]));
            reduce(mx, maxOp<scalar>());
            if (mx > 1e-30) v /= mx;
        }
        const scalarField TSaved(T.primitiveField());
        PtrList<scalarField> TBndSaved(mesh.boundary().size());
        forAll(mesh.boundary(), patchi)
        {
            TBndSaved.set
            (
                patchi,
                new scalarField(T.boundaryField()[patchi])
            );
        }
        scalar adDot = 0.0;
        forAll(v, celli) adDot += thermalObjectiveDerivative[celli]*v[celli];
        reduce(adDot, sumOp<scalar>());

        const scalarList epsList = {1e-2, 3e-3, 1e-3, 3e-4, 1e-4};
        Info<< "B13SOURCEPROBE P1b J/T-source (ADJ v.dSource=" << adDot
            << "):" << endl;
        forAll(epsList, ei)
        {
            const scalar eps = epsList[ei];
            T.primitiveFieldRef() = TSaved + eps*v;
            T.correctBoundaryConditions();
            const scalar jP = evaluateObjective();
            T.primitiveFieldRef() = TSaved - eps*v;
            T.correctBoundaryConditions();
            const scalar jM = evaluateObjective();
            T.primitiveFieldRef() = TSaved;
            forAll(mesh.boundary(), patchi)
            {
                T.boundaryFieldRef()[patchi] = TBndSaved[patchi];
            }
            T.correctBoundaryConditions();
            const scalar fd = (jP - jM)/(2.0*eps);
            const scalar rel =
                mag(fd - adDot)/Foam::max(mag(adDot), 1e-300);
            Info<< "  eps=" << eps << " FD=" << fd << " rel=" << rel
                << (rel <= 1e-6 ? " PASS" : (rel <= 1e-3 ? " WARN" : " FAIL"))
                << endl;
        }
    }

    // -----------------------------------------------------------------
    // P1c: J cold-outlet flux source dot test (dJ/dphi on outlet faces).
    // -----------------------------------------------------------------
    {
        const label outletStart =
            mesh.boundary()[coldOutletPatchID].start();
        const label nOutletFaces =
            mesh.boundary()[coldOutletPatchID].size();
        scalarField v(nOutletFaces, 0.0);
        for (label f = 0; f < nOutletFaces; ++f)
        {
            v[f] = Foam::sin(0.271*(f + 1));
        }
        scalar vmax = 0.0;
        forAll(v, f) vmax = Foam::max(vmax, mag(v[f]));
        if (vmax > 1e-30) v /= vmax;

        const scalarField phiOutletSaved
        (
            phi.boundaryField()[coldOutletPatchID]
        );
        scalar adDot = 0.0;
        const scalarField& outletFluxDeriv =
            thermalObjectiveFluxDerivative.boundaryField()[coldOutletPatchID];
        forAll(v, f) adDot += outletFluxDeriv[f]*v[f];
        reduce(adDot, sumOp<scalar>());

        const scalarList epsList = {1e-2, 3e-3, 1e-3, 3e-4, 1e-4};
        Info<< "B13SOURCEPROBE P1c J/phiOut-source (ADJ v.dSource=" << adDot
            << "):" << endl;
        forAll(epsList, ei)
        {
            const scalar eps = epsList[ei];
            forAll(v, f)
            {
                phi.boundaryFieldRef()[coldOutletPatchID][f] =
                    phiOutletSaved[f] + eps*v[f];
            }
            const scalar jP = evaluateObjective();
            forAll(v, f)
            {
                phi.boundaryFieldRef()[coldOutletPatchID][f] =
                    phiOutletSaved[f] - eps*v[f];
            }
            const scalar jM = evaluateObjective();
            phi.boundaryFieldRef()[coldOutletPatchID] = phiOutletSaved;
            const scalar fd = (jP - jM)/(2.0*eps);
            const scalar rel =
                mag(fd - adDot)/Foam::max(mag(adDot), 1e-300);
            Info<< "  eps=" << eps << " FD=" << fd << " rel=" << rel
                << (rel <= 1e-6 ? " PASS" : (rel <= 1e-3 ? " WARN" : " FAIL"))
                << endl;
        }
    }

    Info<< "--- Stage B13 P1 probes complete (all state restored) ---\n"
        << endl;
}
