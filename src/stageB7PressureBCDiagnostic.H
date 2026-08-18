// =============================================================================
// stageB7PressureBCDiagnostic.H  (BFINAL-007, DIAGNOSTIC_ONLY, switch-guarded)
// -----------------------------------------------------------------------------
// Read-only diagnostic probe gated by the optProperties switch
//     stageB7PressureBCDiagnostic   (default false)
// When disabled the whole block is a complete no-op (no field creation, no
// Info output, no file writes) => production outputs are byte-identical to
// the pre-probe build.
//
// Answers the BFINAL-007 data needs (Q1/Q2/Q3/Q4b/Q5 raw evidence):
//   (a) export the ACTUAL production pressure-equation matrix
//       pEqnDbg = fvm::laplacian(mobilityField, p)  (identical schemes/BCs to
//       the NS.H pEqn, incl. the "corrected" laplacian scheme and the
//       fluxRequired face-flux correction): diag / upper / source, the
//       outlet fixedValue patch internalCoeffs / boundaryCoeffs /
//       deltaCoeffs / magSf / faceCells, and the addBoundaryDiag-augmented
//       diagonal (what solveSegregated would assemble);
//   (b) central FD of the ACTUAL frozen-primal pressure residual
//       R_P(p) = div(phiHbyA) - div(pEqnDbg.flux())   [production pEqn.flux()
//       semantics, raw face-flux sum, NO 1/V, boundary faces owner +]
//       along p +/- eps*n_P, where n_P is the BFINAL-006 null-vector
//       pressure block (read from
//       evidence/agent-group/BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx,
//       normalized); eps sweep exported (exact linearity => eps-independent);
//   (c) export the candidate outlet boundary diagonal
//       +rAtU_c*deltaCoeffs_b*|Sf|_b per outlet cell (residual-Jacobian
//       convention) and the raw internalCoeffs_b (laplacian convention).
//
// Mobility field: internal values = persisted primalPressureMobility
// (production NS.H rAtU internal field, recorded at NS.H L129); boundary
// values = rebuilt 1/UEqn.A() boundary values from the read-only
// reduced-SIMPLE rebuild (same pattern as rxPressureRowTranspose.H).  The
// rebuilt rAU internal field is exported for a fidelity check.
//
// NO production field/matrix is modified: every perturbation acts on LOCAL
// copies; the only writes are diagnostic artifacts under
// evidence/agent-group/BFINAL-007/artifacts/.
// =============================================================================
#include <fstream>
#include <sstream>

if (optProperties.lookupOrDefault<Switch>("stageB7PressureBCDiagnostic", false))
{
    Info<< "--- stageB7PressureBCDiagnostic: start ---" << endl;

    const label s7N = mesh.nCells();
    const label s7NIntF = mesh.nInternalFaces();
    const label s7NUnk = 4*s7N;
    const fileName s7ArtDir
    (
        "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-007/artifacts"
    );
    Foam::mkDir(s7ArtDir);

    // ---- (0) locate the outlet = the only fixesValue patch of p ----
    label s7OutletPatchi = -1;
    forAll(mesh.boundary(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            s7OutletPatchi = patchi;
            Info<< "  stageB7: fixesValue(p) patch = "
                << mesh.boundary()[patchi].name()
                << "  nFaces=" << mesh.boundary()[patchi].size() << endl;
        }
    }
    if (s7OutletPatchi < 0)
    {
        FatalErrorInFunction
            << "stageB7: no fixesValue pressure patch found"
            << abort(FatalError);
    }
    const label s7OB = s7OutletPatchi;
    const label s7NB = mesh.boundary()[s7OB].size();

    // ---- (1) read-only reduced-SIMPLE rebuild (rxPressureRowTranspose.H
    //          pattern) for phiHbyA and the rAU boundary values ----
    tmp<fvVectorMatrix> s7TUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nuEffFrozen, U)
      + fvm::Sp(alpha, U)
     ==
       -fvc::grad(p)
      + fvOptions(U)
    );
    fvVectorMatrix& s7UEqn = s7TUEqn.ref();
    if (ransFlowModel)
    {
        s7UEqn -= fvc::div(nuEffFrozen*dev2(Foam::T(fvc::grad(U))));
    }
    s7UEqn.relax();
    fvOptions.constrain(s7UEqn);
    volScalarField s7rAU(1.0/s7UEqn.A());
    volVectorField s7HbyA(constrainHbyA(s7rAU*s7UEqn.H(), U, p));
    surfaceScalarField s7phiHbyA("stageB7_phiHbyA", fvc::flux(s7HbyA));
    adjustPhi(s7phiHbyA, U, p);

    // rebuild fidelity vs the persisted production rAtU (internal)
    {
        scalar s7Num = 0.0, s7Den = 0.0, s7MaxRel = 0.0;
        forAll(s7rAU, celli)
        {
            const scalar s7D =
                Foam::mag(s7rAU[celli] - primalPressureMobility[celli]);
            const scalar s7E = Foam::mag(primalPressureMobility[celli]);
            s7Num += s7D*s7D;
            s7Den += s7E*s7E;
            s7MaxRel = Foam::max(s7MaxRel, s7E > SMALL ? s7D/s7E : 0.0);
        }
        reduce(s7Num, sumOp<scalar>());
        reduce(s7Den, sumOp<scalar>());
        reduce(s7MaxRel, maxOp<scalar>());
        Info<< "  stageB7: rebuilt rAU vs primalPressureMobility: relL2="
            << Foam::sqrt(s7Num/Foam::max(s7Den, SMALL))
            << " maxRel=" << s7MaxRel << endl;
    }

    // ---- (2) mobility field: production internal (persisted) + rebuilt
    //          boundary (rAtU boundary semantics from 1/UEqn.A()) ----
    volScalarField s7Mobility
    (
        IOobject("stageB7_mobility", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("stageB7_mobility", dimTime, scalar(0.0))
    );
    forAll(s7Mobility, celli)
    {
        s7Mobility[celli] = primalPressureMobility[celli];
    }
    forAll(mesh.boundary(), patchi)
    {
        s7Mobility.boundaryFieldRef()[patchi] == s7rAU.boundaryField()[patchi];
    }

    // ---- (3) production-equivalent pEqn matrix ----
    fvScalarMatrix s7PEqn(fvm::laplacian(s7Mobility, p));
    const scalarField& s7Diag = s7PEqn.diag();
    const scalarField& s7Upper = s7PEqn.upper();
    const scalarField& s7Source = s7PEqn.source();
    const scalarField& s7IcB = s7PEqn.internalCoeffs()[s7OB];
    const scalarField& s7BcB = s7PEqn.boundaryCoeffs()[s7OB];
    const scalarField& s7DeltaB = mesh.deltaCoeffs().boundaryField()[s7OB];
    const scalarField& s7MagSfB = mesh.magSf().boundaryField()[s7OB];
    const labelUList& s7FcB = mesh.boundary()[s7OB].faceCells();

    // addBoundaryDiag-augmented diagonal (what solveSegregated would use)
    scalarField s7DiagBnd(s7Diag);
    forAll(s7FcB, f)
    {
        s7DiagBnd[s7FcB[f]] += s7IcB[f];
    }

    // candidate boundary diagonal (residual-Jacobian convention):
    // +rAtU_c*deltaCoeffs_b*|Sf|_b on each outlet cell
    scalarField s7Cand(s7N, 0.0);
    scalarField s7CandB(s7NB, 0.0);
    forAll(s7FcB, f)
    {
        const scalar s7T =
            primalPressureMobility[s7FcB[f]]*s7DeltaB[f]*s7MagSfB[f];
        s7Cand[s7FcB[f]] = s7T;
        s7CandB[f] = s7T;
    }

    // ---- (4) null-vector pressure block (BFINAL-006, D1 tangent solution
    //          = normalized numerical null vector) ----
    scalarField s7W(s7NUnk, 0.0);
    {
        std::ifstream s7In
        (
            "/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/"
            "BFINAL-006/cycle-2/artifacts/wPrime_TAN_D1.mtx"
        );
        if (!s7In.good())
        {
            FatalErrorInFunction
                << "stageB7: cannot open wPrime_TAN_D1.mtx" << abort(FatalError);
        }
        for (label i = 0; i < s7NUnk; ++i)
        {
            double s7V = 0.0;
            s7In >> s7V;
            if (!s7In.good())
            {
                FatalErrorInFunction
                    << "stageB7: wPrime_TAN_D1.mtx truncated at " << i
                    << " of " << s7NUnk << abort(FatalError);
            }
            s7W[i] = scalar(s7V);
        }
    }
    const scalar s7Wnorm = Foam::sqrt(gSumSqr(s7W));
    scalarField s7nP(s7N, 0.0);
    for (label c = 0; c < s7N; ++c)
    {
        s7nP[c] = s7W[3*s7N + c]/s7Wnorm;
    }
    Info<< "  stageB7: null vector ||w||=" << s7Wnorm
        << "  ||n_P||=" << Foam::sqrt(gSumSqr(s7nP))
        << "  n_P[PREF=cell0]=" << s7nP[0] << endl;

    // ---- (5) actual-residual central FD: R_P(p) = div(phiHbyA) - div(flux)
    //          with production pEqn.flux() semantics (incl. the corrected
    //          face-flux correction), raw face-flux sum (owner +, nei -).
    //          FD = (R_P(p+eps n) - R_P(p-eps n))/(2 eps)
    //             = -(div(flux+) - div(flux-))/(2 eps)   (div(phiHbyA) cancels)
    //          A direct evaluation with psi = n_P (exact linearity, no
    //          central-difference roundoff) is exported as FD_lin = -div(flux).
    // ----
    scalarField s7divPhiHbyA(s7N, 0.0);
    for (label fi = 0; fi < s7NIntF; ++fi)
    {
        const label s7O = mesh.owner()[fi];
        const label s7M = mesh.neighbour()[fi];
        s7divPhiHbyA[s7O] += s7phiHbyA[fi];
        s7divPhiHbyA[s7M] -= s7phiHbyA[fi];
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& s7FC = mesh.boundary()[patchi].faceCells();
        const scalarField& s7PB = s7phiHbyA.boundaryField()[patchi];
        forAll(s7FC, f) { s7divPhiHbyA[s7FC[f]] += s7PB[f]; }
    }

    auto s7AssembleRPFlux = [&](volScalarField& s7Pp, scalarField& s7Out)
    {
        fvScalarMatrix s7PDbg(fvm::laplacian(s7Mobility, s7Pp));
        const surfaceScalarField s7FluxDbg(s7PDbg.flux());
        s7Out = scalar(0.0);
        for (label fi = 0; fi < s7NIntF; ++fi)
        {
            const label s7O = mesh.owner()[fi];
            const label s7M = mesh.neighbour()[fi];
            s7Out[s7O] += s7FluxDbg[fi];
            s7Out[s7M] -= s7FluxDbg[fi];
        }
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& s7FC = mesh.boundary()[patchi].faceCells();
            const scalarField& s7FB = s7FluxDbg.boundaryField()[patchi];
            forAll(s7FC, f) { s7Out[s7FC[f]] += s7FB[f]; }
        }
    };

    // unperturbed residual
    scalarField s7RP0flux(s7N, 0.0);
    s7AssembleRPFlux(p, s7RP0flux);

    // eps sweep (exact linearity => FD eps-independent up to roundoff)
    const scalar s7EpsA = 1.0e-2;
    const scalar s7EpsB = 1.0e-4;
    const scalar s7EpsC = 1.0e-6;
    scalarField s7FDa(s7N, 0.0), s7FDb(s7N, 0.0), s7FDc(s7N, 0.0);
    scalarField s7RPp(s7N, 0.0), s7RPm(s7N, 0.0);

    auto s7FDPerturb = [&](const scalar s7Eps, scalarField& s7FDout)
    {
        volScalarField s7pP(p);
        volScalarField s7pM(p);
        forAll(s7pP, c)
        {
            s7pP.primitiveFieldRef()[c] = p[c] + s7Eps*s7nP[c];
            s7pM.primitiveFieldRef()[c] = p[c] - s7Eps*s7nP[c];
        }
        s7AssembleRPFlux(s7pP, s7RPp);
        s7AssembleRPFlux(s7pM, s7RPm);
        forAll(s7FDout, c)
        {
            // R_P = div(phiHbyA) - div(flux) -> FD = -(divFlux+ - divFlux-)/2eps
            s7FDout[c] = -(s7RPp[c] - s7RPm[c])/(2.0*s7Eps);
        }
    };
    s7FDPerturb(s7EpsA, s7FDa);
    s7FDPerturb(s7EpsB, s7FDb);
    s7FDPerturb(s7EpsC, s7FDc);

    // full residual at the unperturbed state (informational)
    scalarField s7RP0(s7N, 0.0);
    forAll(s7RP0, c) { s7RP0[c] = s7divPhiHbyA[c] - s7RP0flux[c]; }

    // direct (exact) derivative: FD_lin = -div(flux(psi = n_P))
    scalarField s7divFluxN(s7N, 0.0);
    scalarField s7FDlin(s7N, 0.0);
    {
        volScalarField s7pN(p);
        s7pN.primitiveFieldRef() = scalar(0.0);
        forAll(s7pN, c) { s7pN.primitiveFieldRef()[c] = s7nP[c]; }
        s7AssembleRPFlux(s7pN, s7divFluxN);
        forAll(s7FDlin, c) { s7FDlin[c] = -s7divFluxN[c]; }
    }

    // ---- in-probe self-check: FD_lin == -div(flux(n_P)) == -(L_prod . n_P)
    //      where L_prod = diag (with addBoundaryDiag) + upper.  On this case
    //      zeroGradient patches have gradientInternalCoeffs = 0 and the mesh
    //      is perfectly orthogonal (non-orthogonal correction = 0), so the
    //      laplacian boundary coefficient exists ONLY on the fixedValue
    //      outlet -> diagBnd is the complete production diagonal. ----
    scalarField s7Lnp(s7N, 0.0);   // (L_prod . n_P)
    forAll(s7Lnp, c) { s7Lnp[c] = s7DiagBnd[c]*s7nP[c]; }
    for (label fi = 0; fi < s7NIntF; ++fi)
    {
        const label s7O = mesh.owner()[fi];
        const label s7M = mesh.neighbour()[fi];
        s7Lnp[s7O] += s7Upper[fi]*s7nP[s7M];
        s7Lnp[s7M] += s7Upper[fi]*s7nP[s7O];
    }
    {
        scalar s7Num = 0.0, s7Den = 0.0;
        forAll(s7Lnp, c)
        {
            const scalar s7D = s7FDlin[c] + s7Lnp[c];   // FD_lin - (-L.n)
            s7Num += s7D*s7D;
            s7Den += s7Lnp[c]*s7Lnp[c];
        }
        reduce(s7Num, sumOp<scalar>());
        reduce(s7Den, sumOp<scalar>());
        Info<< "  stageB7: self-check |FD_lin - (-L_prod.n_P)|/|L_prod.n_P| = "
            << Foam::sqrt(s7Num/Foam::max(s7Den, SMALL)) << endl;
    }

    // ---- outlet concentration + closure of FD vs candidate.n_P ----
    {
        // candidate.n_P (residual-convention boundary term action)
        scalarField s7CandNP(s7N, 0.0);
        forAll(s7CandNP, c) { s7CandNP[c] = s7Cand[c]*s7nP[c]; }

        scalar s7FDmax = 0.0, s7FDsum = 0.0, s7FDoutSum = 0.0;
        scalar s7OutNum = 0.0, s7OutDen = 0.0;
        scalar s7Dot = 0.0, s7NormFD = 0.0, s7NormCand = 0.0;
        forAll(s7FDlin, c)
        {
            const scalar s7V = Foam::mag(s7FDlin[c]);
            s7FDmax = Foam::max(s7FDmax, s7V);
            s7FDsum += s7V;
        }
        forAll(s7FcB, f)
        {
            const label s7C = s7FcB[f];
            s7FDoutSum += Foam::mag(s7FDlin[s7C]);
            s7OutNum += Foam::sqr(s7FDlin[s7C] - s7CandNP[s7C]);
            s7OutDen += Foam::sqr(s7CandNP[s7C]);
            s7Dot += s7FDlin[s7C]*s7CandNP[s7C];
            s7NormFD += Foam::sqr(s7FDlin[s7C]);
            s7NormCand += Foam::sqr(s7CandNP[s7C]);
        }
        reduce(s7FDmax, maxOp<scalar>());
        reduce(s7FDsum, sumOp<scalar>());
        reduce(s7FDoutSum, sumOp<scalar>());
        reduce(s7OutNum, sumOp<scalar>());
        reduce(s7OutDen, sumOp<scalar>());
        reduce(s7Dot, sumOp<scalar>());
        reduce(s7NormFD, sumOp<scalar>());
        reduce(s7NormCand, sumOp<scalar>());
        const scalar s7Cos =
            Foam::sqrt(Foam::max(s7NormFD, SMALL))
           *Foam::sqrt(Foam::max(s7NormCand, SMALL)) > SMALL
          ? s7Dot/(Foam::sqrt(s7NormFD)*Foam::sqrt(s7NormCand)) : 0.0;
        Info<< "  stageB7: FD_lin max|.|=" << s7FDmax
            << "  outlet-mass-fraction=" << (s7FDsum > SMALL ? s7FDoutSum/s7FDsum : 0.0)
            << "  |FD-cand.n_P|_outlet/|cand.n_P|_outlet="
            << Foam::sqrt(s7OutNum/Foam::max(s7OutDen, SMALL))
            << "  cos(FD, cand.n_P)_outlet=" << s7Cos << endl;
    }

    // ---- eps robustness (sign-corrected central FD) ----
    {
        scalar s7Num = 0.0, s7Den = 0.0;
        forAll(s7FDa, c)
        {
            const scalar s7V = s7FDb[c];
            s7Num += Foam::sqr(s7FDa[c] - s7V) + Foam::sqr(s7FDc[c] - s7V);
            s7Den += Foam::sqr(s7V);
        }
        reduce(s7Num, sumOp<scalar>());
        reduce(s7Den, sumOp<scalar>());
        Info<< "  stageB7: FD(central-diff) eps-spread relL2(1e-2,1e-6 vs 1e-4) = "
            << Foam::sqrt(s7Num/Foam::max(s7Den, SMALL))
            << "  ; |FD_cd(1e-4) - FD_lin|/|FD_lin| = ";
        scalar s7N2 = 0.0, s7D2 = 0.0;
        forAll(s7FDb, c)
        {
            const scalar s7W = s7FDb[c] - s7FDlin[c];
            s7N2 += s7W*s7W;
            s7D2 += s7FDlin[c]*s7FDlin[c];
        }
        reduce(s7N2, sumOp<scalar>());
        reduce(s7D2, sumOp<scalar>());
        Info<< Foam::sqrt(s7N2/Foam::max(s7D2, SMALL)) << endl;
    }

    // ---- (a/b/c) exports (plain text, loadtxt-friendly) ----
    {
        OFstream s7Out(s7ArtDir + "/stageB7_outlet_patch.mtx");
        s7Out.precision(16);
        forAll(s7FcB, f)
        {
            s7Out << s7FcB[f] << ' ' << s7DeltaB[f] << ' ' << s7MagSfB[f]
                  << ' ' << s7IcB[f] << ' ' << s7BcB[f] << ' '
                  << s7rAU.boundaryField()[s7OB][f] << ' '
                  << primalPressureMobility[s7FcB[f]] << ' ' << s7CandB[f]
                  << '\n';
        }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_laplacian_diag.mtx");
        s7Out.precision(16);
        forAll(s7Diag, c) { s7Out << s7Diag[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_laplacian_upper.mtx");
        s7Out.precision(16);
        forAll(s7Upper, f) { s7Out << s7Upper[f] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_laplacian_diag_bnd.mtx");
        s7Out.precision(16);
        forAll(s7DiagBnd, c) { s7Out << s7DiagBnd[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_laplacian_source.mtx");
        s7Out.precision(16);
        forAll(s7Source, c) { s7Out << s7Source[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_nP.mtx");
        s7Out.precision(16);
        forAll(s7nP, c) { s7Out << s7nP[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_divPhiHbyA.mtx");
        s7Out.precision(16);
        forAll(s7divPhiHbyA, c) { s7Out << s7divPhiHbyA[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_RP_full.mtx");
        s7Out.precision(16);
        forAll(s7RP0, c)
        {
            s7Out << s7RP0[c] << ' ' << s7RPp[c] << ' ' << s7RPm[c] << '\n';
        }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_RP_FD.mtx");
        s7Out.precision(16);
        forAll(s7FDb, c) { s7Out << s7FDb[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_RP_FD_lin.mtx");
        s7Out.precision(16);
        forAll(s7FDlin, c) { s7Out << s7FDlin[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_RP_divFluxN.mtx");
        s7Out.precision(16);
        forAll(s7divFluxN, c) { s7Out << s7divFluxN[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_RP_FD_all_eps.mtx");
        s7Out.precision(16);
        forAll(s7FDa, c)
        {
            s7Out << s7FDa[c] << ' ' << s7FDb[c] << ' ' << s7FDc[c] << '\n';
        }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_candidate_diag.mtx");
        s7Out.precision(16);
        forAll(s7Cand, c) { s7Out << s7Cand[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_rebuilt_rAU.mtx");
        s7Out.precision(16);
        forAll(s7rAU, c) { s7Out << s7rAU[c] << '\n'; }
    }
    {
        OFstream s7Out(s7ArtDir + "/stageB7_laplacian_Lnp.mtx");
        s7Out.precision(16);
        forAll(s7Lnp, c) { s7Out << s7Lnp[c] << '\n'; }
    }
    // full addBoundaryDiag diagonal over ALL patches (complete production
    // diagonal after solveSegregated; on this case only the fixedValue outlet
    // contributes, but export the complete assembly for the audit)
    {
        scalarField s7DiagBndFull(s7Diag);
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& s7FCp = mesh.boundary()[patchi].faceCells();
            const scalarField& s7IcP = s7PEqn.internalCoeffs()[patchi];
            forAll(s7FCp, f) { s7DiagBndFull[s7FCp[f]] += s7IcP[f]; }
        }
        OFstream s7Out(s7ArtDir + "/stageB7_laplacian_diag_bnd_full.mtx");
        s7Out.precision(16);
        forAll(s7DiagBndFull, c) { s7Out << s7DiagBndFull[c] << '\n'; }
    }
    // per-boundary-face boundary coefficients over ALL patches
    // (patchIndex cell deltaCoeffs magSf internalCoeffs boundaryCoeffs)
    {
        OFstream s7Out(s7ArtDir + "/stageB7_boundary_all_patches.mtx");
        s7Out.precision(16);
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& s7FCp = mesh.boundary()[patchi].faceCells();
            const scalarField& s7IcP = s7PEqn.internalCoeffs()[patchi];
            const scalarField& s7BcP = s7PEqn.boundaryCoeffs()[patchi];
            const scalarField& s7DP =
                mesh.deltaCoeffs().boundaryField()[patchi];
            const scalarField& s7MS =
                mesh.magSf().boundaryField()[patchi];
            forAll(s7FCp, f)
            {
                s7Out << patchi << ' ' << s7FCp[f] << ' ' << s7DP[f] << ' '
                      << s7MS[f] << ' ' << s7IcP[f] << ' ' << s7BcP[f]
                      << '\n';
            }
        }
    }

    Info<< "--- stageB7PressureBCDiagnostic: complete (artifacts in "
        << s7ArtDir << ") ---" << endl;
}
