// =============================================================================
// rxPressureRowTranspose.H  (BFINAL-005, PATCH_RX_PRESSURE_ROW;
//                            BFINAL-019 parameterized)
// -----------------------------------------------------------------------------
// Computes the per-cell field RX_OUTPUT_FIELD = J_P^T <adjoint pressure>: the
// EXACT transpose of the BFINAL-004-validated forward pressure/continuity-row
// design operator R_P,x = J_P*(dAlphaDxh*z) with R_P = div(phi) (raw face-flux
// sum, NO 1/V).  Historical/default contraction: pc (pressure-drop adjoint);
// BFINAL-019 adds the pb (thermal-coupling adjoint) contraction, see the
// parameterization block below the header.
//
// Forward operator (stageB6RxDesignOracle.H assembleWeightedRPa, BFINAL-004
// CLOSED; identical reduced-SIMPLE semantics as the task §7):
//   J_P*w = div(dphiHbyA_w - dflux_w),   w = per-cell directional weight
//     dphiHbyA_w[f] = Sf&(wf*(dHbyA*w)_own + (1-wf)*(dHbyA*w)_nei)  (internal)
//     dphiHbyA_w[b] = Sf&(uAssignable ? (dHbyA*w)_cell : 0)         (boundary)
//     dflux_w       = flux(fvm::laplacian(drAU*w, p))
// with the BFINAL-004-LOCKED per-cell derivatives (NO extra V):
//     drAU  = -rAU^2/alphaRel
//     dHbyA = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)
// rebuilt READ-ONLY from the reduced-SIMPLE layer (no solve, no state write):
//     UEqn = div(phi,U) - laplacian(nuEffFrozen,U) + Sp(alpha,U)
//          == -grad(p) + fvOptions(U)
//          [- div(nuEffFrozen*dev2(T(grad U))) if ransFlowModel]
//     UEqn.relax(); fvOptions.constrain(UEqn)
//     alphaRel = relaxEquation("U") ? equationRelaxationFactor("U") : 1.0
//     rAU = 1/UEqn.A(); HbyA = constrainHbyA(rAU*UEqn.H(), U, p)
//
// Transpose (per cell, accumulated over faces; T = T1 + T2):
//   T1[own] += wf*(Sf&dHbyA[own])*(pc[own]-pc[nei])
//   T1[nei] += (1-wf)*(Sf&dHbyA[nei])*(pc[own]-pc[nei])
//   T1[b]   += (Sf&dHbyA[cell])*pc[cell]            (assignable U patches only,
//              mirrors the constrainHbyA pinning in the forward)
//   g0 = flux(fvm::laplacian(one,p))  (one = uniform 1.0; scheme-exact face
//        gradient; fvm::laplacian is LINEAR in the coefficient field and
//        fvMatrix::flux() is linear in the matrix coefficients, so
//        dflux_w[f] = g0[f]*(wf*drAUw[own]+(1-wf)*drAUw[nei]))
//   T2[own] += -wf*g0[f]*(pc[own]-pc[nei])*drAU[own]
//   T2[nei] += -(1-wf)*g0[f]*(pc[own]-pc[nei])*drAU[nei]
//   T2[b]   += 0   (the forward coefficient field drAU*w carries the default
//              calculated BC with boundary value 0, so the interpolated
//              boundary gamma_b = 0 and dflux_b == 0 exactly in this case;
//              the runtime dot test below arbitrates this claim)
//
// PRODUCTION SEMANTICS (BFINAL-005 Fixer round, reviewer-mandated): the
// reduced-SIMPLE rebuild, drAU/dHbyA, g0, T1/T2 face loops, and the
// rxPressureRowT assignment are UNCONDITIONAL — the pressure/continuity-row
// term is part of the production pressure-drop sensitivity with the
// stageB6RxDesignOracle switch OFF as well (task §6/§22).  ONLY the runtime
// dot-test self-check (forward re-assembly + FatalError) and its diagnostic
// Info are gated behind "stageB6RxDesignOracle".
//
// Runtime self-check (dot test, gated): with the deterministic weight
//     w = sin(0.271*(celli+1)) on active unclipped cells
// (identical pattern to the BFINAL-004 deltaAlpha support), re-assemble
// J_P*w locally (exact stageB6 assembleWeightedRPa replica) and verify
//     sum(T*w) == sum(RX_ADJ_PRESSURE*(J_P*w))    relative |diff| <= 1e-8
// FatalError on mismatch -> pins the sign/transpose BEFORE the production
// gates consume T.  Expected relErr ~ 1e-14 (machine precision).
//
// Include-safe: no global state, no project-header includes.  Include inside
// src/sensitivity.H DISCRETE branch, once per contraction, before the
// corresponding sensitivity term is assembled.
//
// BFINAL-019 parameterization: the helper can be contracted with either
// pressure adjoint (production pressure-drop chain uses pc; the thermal-
// coupling J chain uses pb) and can additionally produce the Gx direct-term
// basis.  The caller selects the contraction via the following macros; with
// NONE of them defined the expansion is token-identical to the historical
// (BFINAL-005) text, so the production pc path is numerically unchanged:
//   RX_ADJ_PRESSURE   adjoint pressure field          (default pc)
//   RX_OUTPUT_FIELD   output field variable name      (default rxPressureRowT)
//   RX_OUTPUT_NAME    registered IOobject name string (default "rxPressureRowT")
//   RX_LABEL          log label string                (default "pc")
//   RX_FLUX_DIRECT_SOURCE  (OPTIONAL) name of a const scalarField of size
//                     mesh.nFaces() holding the thermal-coupling face
//                     functional g_f = dJ/dphi_f - (dR_T/dphi_f)^T Tb (the
//                     SAME functional the BFINAL-018 b_TC folding routed into
//                     the (U,p) rows).  When defined, the additional per-cell
//                     basis field rxFluxDirectT is produced such that
//                         sum(rxFluxDirectT*dAlphaDxh*z) == g^T phi_x
//                     i.e. the Gx direct term of DERIVATION_FIX3,
//                         phi_x[f] = Sf&(interp(dHbyA)*wc)
//                                  - interp(drAU)*wc*g0[f]        (internal)
//                         phi_x[b] = Sf_b&dHbyA[cell]*wc[cell]    (U-assignable
//                                    boundary faces only; dflux_b == 0 exactly
//                                    because drAU*w carries a zero boundary
//                                    value -> interpolated gamma_b = 0)
//                     with wc = dAlphaDxh*z and the SAME BFINAL-004-LOCKED
//                     per-cell derivatives used by the pressure-row transpose.
// =============================================================================
#include <sstream>

#ifndef RX_ADJ_PRESSURE
#define RX_ADJ_PRESSURE pc
#endif
#ifndef RX_OUTPUT_FIELD
#define RX_OUTPUT_FIELD rxPressureRowT
#endif
#ifndef RX_OUTPUT_NAME
#define RX_OUTPUT_NAME "rxPressureRowT"
#endif
#ifndef RX_LABEL
#define RX_LABEL "pc"
#endif

volScalarField RX_OUTPUT_FIELD
(
    IOobject
    (
        RX_OUTPUT_NAME,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rxPressureRowT", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

#ifdef RX_FLUX_DIRECT_SOURCE
// BFINAL-019 Gx basis field (only produced for the contraction that defines
// RX_FLUX_DIRECT_SOURCE).  Values are assigned through primitiveFieldRef()
// from the raw scalarField contraction below; the (0, 2, -2) dimension is the
// same per-cell adjoint-pressure scale as RX_OUTPUT_FIELD so that the
// production term gsenshMeanTFluxDirect = +rxFluxDirectT*dAlphaDxh carries
// the (0, 2, -3) dimensions of fsenshMeanT.
volScalarField rxFluxDirectT
(
    IOobject
    (
        "rxFluxDirectT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rxFluxDirectT", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);
#endif

{
    const label rxNCells = mesh.nCells();
    const label rxNIntF = mesh.nInternalFaces();
    const label rxNFaceTot = mesh.nFaces();

    // ---- (a) read-only reduced-SIMPLE layer rebuild (LOCKED form) ----
    tmp<fvVectorMatrix> rxTUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nuEffFrozen, U)
      + fvm::Sp(alpha, U)
     ==
       -fvc::grad(p)
      + fvOptions(U)
    );
    fvVectorMatrix& rxUEqn = rxTUEqn.ref();
    if (ransFlowModel)
    {
        rxUEqn -= fvc::div
        (
            nuEffFrozen*dev2(Foam::T(fvc::grad(U)))
        );
    }
    rxUEqn.relax();
    fvOptions.constrain(rxUEqn);

    const scalar rxAlphaRel =
        mesh.relaxEquation("U")
      ? mesh.equationRelaxationFactor("U")
      : scalar(1.0);
    volScalarField rxrAU(1.0/rxUEqn.A());
    volVectorField rxHbyA(constrainHbyA(rxrAU*rxUEqn.H(), U, p));

    // ---- (b) BFINAL-004-locked per-cell derivatives (NO extra V) ----
    scalarField rxDrAU(rxNCells, 0.0);
    vectorField rxDHbyA(rxNCells, vector::zero);
    forAll(U, celli)
    {
        const scalar rAUc = rxrAU[celli];
        rxDrAU[celli] = -rAUc*rAUc/rxAlphaRel;
        rxDHbyA[celli] =
            (rAUc/rxAlphaRel)
           *((1.0 - rxAlphaRel)*U[celli] - rxHbyA[celli]);
    }

    // ---- (c) g0 = flux(fvm::laplacian(one,p)) over ALL faces ----
    volScalarField rxOne
    (
        IOobject
        (
            "rxPressureRowOne",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );
    fvScalarMatrix rxOneEqn(fvm::laplacian(rxOne, p));
    surfaceScalarField rxG0Field(rxOneEqn.flux());
    scalarField rxG0(rxNFaceTot, 0.0);
    for (label fi = 0; fi < rxNIntF; ++fi) { rxG0[fi] = rxG0Field[fi]; }
    scalar rxG0BndMax = 0.0;
    {
        label bf = rxNIntF;
        forAll(rxG0Field.boundaryField(), patchi)
        {
            const scalarField& gB = rxG0Field.boundaryField()[patchi];
            forAll(gB, f)
            {
                rxG0[bf] = gB[f];
                rxG0BndMax = Foam::max(rxG0BndMax, Foam::mag(gB[f]));
                ++bf;
            }
        }
    }

    // ---- (c) transpose face loops: T1 (HbyA) + T2 (flux) ----
    scalarField rxT1(rxNCells, 0.0);
    scalarField rxT2(rxNCells, 0.0);
    forAll(mesh.owner(), fi)
    {
        const label own = mesh.owner()[fi];
        const label nei = mesh.neighbour()[fi];
        const scalar wf = mesh.weights()[fi];
        const scalar pcDiff = RX_ADJ_PRESSURE[own] - RX_ADJ_PRESSURE[nei];
        rxT1[own] += wf*(mesh.Sf()[fi] & rxDHbyA[own])*pcDiff;
        rxT1[nei] += (1.0 - wf)*(mesh.Sf()[fi] & rxDHbyA[nei])*pcDiff;
        rxT2[own] += -wf*rxG0[fi]*pcDiff*rxDrAU[own];
        rxT2[nei] += -(1.0 - wf)*rxG0[fi]*pcDiff*rxDrAU[nei];
    }
    {
        label bf = rxNIntF;
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
            const bool uAssignable = U.boundaryField()[patchi].assignable();
            forAll(fc, f)
            {
                const label celli = fc[f];
                if (uAssignable)
                {
                    rxT1[celli] +=
                        (pSf[f] & rxDHbyA[celli])*RX_ADJ_PRESSURE[celli];
                }
                ++bf;
            }
        }
    }
    scalarField rxT(rxNCells, 0.0);
    forAll(rxT, celli) { rxT[celli] = rxT1[celli] + rxT2[celli]; }

    // ---- PRODUCTION assignment: RX_OUTPUT_FIELD = T (UNCONDITIONAL) ----
    // This runs with stageB6RxDesignOracle OFF as well, so the production
    // sensitivities contain the pressure/continuity-row term
    // (pressure-drop chain) gsenshPressureDropPressureRow
    //   = -rxPressureRowT*dAlphaDxh        (BFINAL-005, pc)
    // and (thermal-coupling J chain, BFINAL-019) gsenshMeanTPressureRow
    //   = -rxPressureRowTb*dAlphaDxh       (pb).
    forAll(RX_OUTPUT_FIELD.primitiveFieldRef(), celli)
    {
        RX_OUTPUT_FIELD.primitiveFieldRef()[celli] = rxT[celli];
    }

    // BFINAL-032: gated export of the internal-face segment (the raw
    // J_P operator ingredients) so the offline B32 decomposition can
    // compare them against the coordinator reconstruction.  The values
    // are identical between the pc and pb contractions (same primal
    // rebuild), so the writes fire only on the pc contraction (the pc
    // include is the only one WITHOUT RX_FLUX_DIRECT_SOURCE) to avoid
    // duplication.  Names carry the deterministic _r1/_r2 round tag (see
    // sensitivity.H top); rxrAU/rxG0Field are written as named copies
    // because their IOobjects were tmp-constructed.
#ifndef RX_FLUX_DIRECT_SOURCE
    if (b32WriteDecomposition)
    {
        volScalarField rxrAUR
        (
            IOobject
            (
                "rxrAU" + b32RoundTag,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rxrAU
        );
        rxrAUR.write();
        surfaceScalarField rxG0FieldR
        (
            IOobject
            (
                "rxG0Field" + b32RoundTag,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rxG0Field
        );
        rxG0FieldR.write();
        // rxDHbyA is a raw vectorField; wrap it in a volVectorField for
        // a regular ascii field write (b31b-style parser reads ascii).
        volVectorField rxDHbyAField
        (
            IOobject
            (
                "rxDHbyA" + b32RoundTag,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "rxDHbyA",
                dimensionSet(0, 1, -1, 0, 0, 0, 0),
                vector::zero
            ),
            zeroGradientFvPatchVectorField::typeName
        );
        rxDHbyAField.primitiveFieldRef() = rxDHbyA;
        rxDHbyAField.write();
    }
#endif  // !RX_FLUX_DIRECT_SOURCE (pc-only internal-face write)

#ifdef RX_FLUX_DIRECT_SOURCE
    // ---- (e) BFINAL-019 Gx direct-term basis: per-cell field whose
    // contraction with dAlphaDxh*z equals g^T phi_x (see file header).
    // Face channels mirror the (dphiHbyA_w - dflux_w) channels of the
    // forward operator above, with the face functional g in place of the
    // cell-centred adjoint pressure difference:
    //   internal: own += wf*((Sf&dHbyA[own]) - g0[f]*drAU[own])*g[f]
    //             nei += (1-wf)*((Sf&dHbyA[nei]) - g0[f]*drAU[nei])*g[f]
    //   boundary (U-assignable): += (Sf_b&dHbyA[cell])*g[b]
    //     (the g0/drAU channel has exactly zero boundary contribution,
    //      identical to T2[b] above)
    scalarField rxFx(rxNCells, 0.0);
    forAll(mesh.owner(), fi)
    {
        const label own = mesh.owner()[fi];
        const label nei = mesh.neighbour()[fi];
        const scalar wf = mesh.weights()[fi];
        const scalar gFace = RX_FLUX_DIRECT_SOURCE[fi];
        rxFx[own] +=
            wf
           *( (mesh.Sf()[fi] & rxDHbyA[own]) - rxG0[fi]*rxDrAU[own] )
           *gFace;
        rxFx[nei] +=
            (1.0 - wf)
           *( (mesh.Sf()[fi] & rxDHbyA[nei]) - rxG0[fi]*rxDrAU[nei] )
           *gFace;
    }
    {
        label bf = rxNIntF;
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
            const bool uAssignable = U.boundaryField()[patchi].assignable();
            forAll(fc, f)
            {
                const label celli = fc[f];
                if (uAssignable)
                {
                    rxFx[celli] +=
                        (pSf[f] & rxDHbyA[celli])*RX_FLUX_DIRECT_SOURCE[bf];
                }
                ++bf;
            }
        }
    }
    forAll(rxFluxDirectT.primitiveFieldRef(), celli)
    {
        rxFluxDirectT.primitiveFieldRef()[celli] = rxFx[celli];
    }
#endif

    // ---- (d) runtime dot-test: sum(T*w) vs sum(adjoint*(J_P*w)) ----
    // DIAGNOSTIC ONLY (gated): re-assembles the forward operator to pin the
    // sign/transpose at machine precision before the gates consume T.
    if (optProperties.lookupOrDefault<Switch>("stageB6RxDesignOracle", false))
    {
        // deterministic weight (BFINAL-004 deltaAlpha support pattern)
        scalarField rxW(rxNCells, 0.0);
        label rxWSupport = 0;
        forAll(xh, celli)
        {
            if
            (
                designMask[celli] > 0.5
             && alpha[celli] > alphamin.value()*(1.0 + 1e-10)
            )
            {
                rxW[celli] = Foam::sin(0.271*(celli + 1));
                ++rxWSupport;
            }
        }
        // forward J_P*w (exact stageB6 assembleWeightedRPa replica)
        scalarField rxDrAUw(rxNCells, 0.0);
        vectorField rxDHbyAw(rxNCells, vector::zero);
        forAll(U, celli)
        {
            rxDrAUw[celli] = rxDrAU[celli]*rxW[celli];
            rxDHbyAw[celli] = rxDHbyA[celli]*rxW[celli];
        }
        scalarField rxDphiHbyAw(rxNFaceTot, 0.0);
        forAll(mesh.owner(), fi)
        {
            const label own = mesh.owner()[fi];
            const label nei = mesh.neighbour()[fi];
            const scalar wf = mesh.weights()[fi];
            rxDphiHbyAw[fi] =
                mesh.Sf()[fi]
              & (wf*rxDHbyAw[own] + (1.0 - wf)*rxDHbyAw[nei]);
        }
        {
            label bf = rxNIntF;
            forAll(mesh.boundary(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
                const bool uAssignable = U.boundaryField()[patchi].assignable();
                forAll(fc, f)
                {
                    const label celli = fc[f];
                    const vector dHb =
                        uAssignable ? rxDHbyAw[celli] : vector::zero;
                    rxDphiHbyAw[bf] = pSf[f] & dHb;
                    ++bf;
                }
            }
        }
        volScalarField rxDrAUwField
        (
            IOobject
            (
                "rxPressureRowDrAUw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("d", dimTime, 0.0)
        );
        forAll(rxDrAUwField.primitiveFieldRef(), celli)
        {
            rxDrAUwField.primitiveFieldRef()[celli] = rxDrAUw[celli];
        }
        rxDrAUwField.correctBoundaryConditions();
        fvScalarMatrix rxDPEqn(fvm::laplacian(rxDrAUwField, p));
        surfaceScalarField rxDfluxW(rxDPEqn.flux());
        scalarField rxDfluxWFaces(rxNFaceTot, 0.0);
        for (label fi = 0; fi < rxNIntF; ++fi) { rxDfluxWFaces[fi] = rxDfluxW[fi]; }
        {
            label bf = rxNIntF;
            forAll(rxDfluxW.boundaryField(), patchi)
            {
                const scalarField& dfB = rxDfluxW.boundaryField()[patchi];
                forAll(dfB, f) { rxDfluxWFaces[bf] = dfB[f]; ++bf; }
            }
        }
        scalarField rxAnRPw(rxNCells, 0.0);
        forAll(mesh.owner(), fi)
        {
            const label own = mesh.owner()[fi];
            const label nei = mesh.neighbour()[fi];
            rxAnRPw[own] += rxDphiHbyAw[fi] - rxDfluxWFaces[fi];
            rxAnRPw[nei] -= rxDphiHbyAw[fi] - rxDfluxWFaces[fi];
        }
        {
            label bf = rxNIntF;
            forAll(mesh.boundary(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                forAll(fc, f)
                {
                    const label celli = fc[f];
                    rxAnRPw[celli] += rxDphiHbyAw[bf] - rxDfluxWFaces[bf];
                    ++bf;
                }
            }
        }

        scalar rxDotFwd = 0.0, rxDotT = 0.0;
        scalar rxNormT = 0.0, rxNormFwd = 0.0;
        forAll(rxT, celli)
        {
            rxDotFwd += RX_ADJ_PRESSURE[celli]*rxAnRPw[celli];
            rxDotT += rxT[celli]*rxW[celli];
            rxNormT += sqr(rxT[celli]);
            rxNormFwd += sqr(rxAnRPw[celli]);
        }
        reduce(rxDotFwd, sumOp<scalar>());
        reduce(rxDotT, sumOp<scalar>());
        reduce(rxNormT, sumOp<scalar>());
        reduce(rxNormFwd, sumOp<scalar>());
        reduce(rxWSupport, sumOp<label>());
        const scalar rxDotScale =
            Foam::max(Foam::mag(rxDotFwd), Foam::mag(rxDotT));
        const scalar rxRelDot =
            rxDotScale > SMALL ? Foam::mag(rxDotFwd - rxDotT)/rxDotScale : 0.0;
        std::ostringstream rxDotFwdStr;
        rxDotFwdStr.precision(16);
        rxDotFwdStr << rxDotFwd;
        std::ostringstream rxDotTStr;
        rxDotTStr.precision(16);
        rxDotTStr << rxDotT;
        Info<< "RxPressureRowTranspose[" RX_LABEL "]: sum(T*w)=" << rxDotTStr.str()
            << " sum(" RX_LABEL "*(J_P*w))=" << rxDotFwdStr.str()
            << " relErr=" << rxRelDot
            << " alphaRel=" << rxAlphaRel
            << " |T|L2=" << Foam::sqrt(rxNormT)
            << " |J_P*w|L2=" << Foam::sqrt(rxNormFwd)
            << " |g0_bnd|max=" << rxG0BndMax
            << " w-support=" << rxWSupport << endl;
        if (rxRelDot > 1e-8)
        {
            FatalErrorInFunction
                << "RxPressureRowTranspose dot-test FAILED: relErr="
                << rxRelDot << " > 1e-8. rxPressureRowT is NOT the transpose "
                << "of the validated J_P operator." << abort(FatalError);
        }
    }
}

// Undefine the BFINAL-019 parameterization macros so each inclusion of this
// file is self-contained (defaults re-apply unless the caller re-defines).
#undef RX_ADJ_PRESSURE
#undef RX_OUTPUT_FIELD
#undef RX_OUTPUT_NAME
#undef RX_LABEL
#undef RX_FLUX_DIRECT_SOURCE
