scalarField thermalDiffusionDerivativeDTCell(mesh.nCells(), 0.0);

// Adjoint thermal convection flux (frozen-hot mode uses the combined
// cold/hot flux).  Declared here so both the discrete-transpose branch and
// the continuous branch can use it.
const surfaceScalarField& activeAdjointThermalFlux =
    frozenHotRegionEnabled ? phiThermal : phi;

if (adjointMode == "discrete")
{
    scalarField primalDiag(TEqn.diag());
    forAll(TEqn.internalCoeffs(), patchi)
    {
        const labelUList& patchCells = mesh.boundary()[patchi].faceCells();
        const scalarField& patchDiag = TEqn.internalCoeffs()[patchi];
        forAll(patchCells, facei)
        {
            primalDiag[patchCells[facei]] += patchDiag[facei];
        }
    }

    const scalarField primalLower(TEqn.lower());
    const scalarField primalUpper(TEqn.upper());
    const dimensionSet transposeDimensions =
        TEqn.dimensions()/T.dimensions()*Tb.dimensions();
    fvScalarMatrix thermalTransposeEqn(Tb, transposeDimensions);
    thermalTransposeEqn.diag() = primalDiag;
    thermalTransposeEqn.lower() = primalUpper;
    thermalTransposeEqn.upper() = primalLower;
    thermalTransposeEqn.source() = thermalObjectiveDerivative;

    forAll(thermalTransposeEqn.internalCoeffs(), patchi)
    {
        if (mesh.boundary()[patchi].coupled())
        {
            // Parallel discrete adjoint fix: the transposed operator's
            // processor/cyclic coupling is obtained by SWAPPING the roles of
            // the forward operator's internalCoeffs (diagonal contribution)
            // and boundaryCoeffs (off-diagonal coupling to the remote cells).
            // This is exact for the symmetric laplacian and for the upwind
            // divergence (upwind(-phi) is the transpose of upwind(phi)).
            // Previously all patches were zeroed, which decoupled the
            // subdomains and gave Tb=0 outside the objective subdomains in
            // parallel runs.
            const labelUList& patchCells =
                mesh.boundary()[patchi].faceCells();
            const scalarField& fInt = TEqn.internalCoeffs()[patchi];
            const scalarField& fBou = TEqn.boundaryCoeffs()[patchi];
            // addBoundaryDiag() (inside fvMatrix::solve) adds internalCoeffs
            // of ALL patches to the diagonal; primalDiag already contains the
            // forward internalCoeffs (true transpose diagonal), so pre-
            // compensate the swapped-in boundary coefficients.
            forAll(patchCells, facei)
            {
                thermalTransposeEqn.diag()[patchCells[facei]] -=
                    fBou[facei];
            }
            thermalTransposeEqn.internalCoeffs()[patchi] = fBou;
            thermalTransposeEqn.boundaryCoeffs()[patchi] = fInt;
        }
        else
        {
            // Physical boundaries: homogeneous (zero boundary contribution)
            // in the transposed system, as validated in serial by Gate 4.
            thermalTransposeEqn.internalCoeffs()[patchi] = scalar(0.0);
            thermalTransposeEqn.boundaryCoeffs()[patchi] = scalar(0.0);
        }
    }

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    scalar maxTransposeError = 0.0;

    for (label testI = 0; testI < 3; ++testI)
    {
        scalarField a(mesh.nCells(), 0.0);
        scalarField bTest(mesh.nCells(), 0.0);
        scalarField Ab(mesh.nCells(), 0.0);
        scalarField ATa(mesh.nCells(), 0.0);

        forAll(a, celli)
        {
            a[celli] = Foam::sin((0.173 + 0.019*testI)*(celli + 1));
            bTest[celli] = Foam::cos((0.117 + 0.023*testI)*(celli + 1));
            Ab[celli] = primalDiag[celli]*bTest[celli];
            ATa[celli] = primalDiag[celli]*a[celli];
        }

        forAll(primalUpper, facei)
        {
            const label own = owner[facei];
            const label nei = neighbour[facei];
            Ab[own] += primalUpper[facei]*bTest[nei];
            Ab[nei] += primalLower[facei]*bTest[own];
            ATa[own] += primalLower[facei]*a[nei];
            ATa[nei] += primalUpper[facei]*a[own];
        }

        scalar aAb = sum(a*Ab);
        scalar bATa = sum(bTest*ATa);
        reduce(aAb, sumOp<scalar>());
        reduce(bATa, sumOp<scalar>());
        const scalar transposeError =
            mag(aAb - bATa)
           /Foam::max
            (
                Foam::max(mag(aAb), mag(bATa)),
                scalar(1e-30)
            );
        maxTransposeError = Foam::max(maxTransposeError, transposeError);
    }

    Info<< "Thermal matrix transpose dot-test max relative error="
        << maxTransposeError << endl;
    if (maxTransposeError > 1e-10)
    {
        FatalErrorInFunction
            << "Thermal matrix transpose dot test failed with error "
            << maxTransposeError << abort(FatalError);
    }

    {
        SolverPerformance<scalar> perfTb = thermalTransposeEqn.solve();
        adjointResidualT = perfTb.finalResidual();
    }
    Tb.correctBoundaryConditions();
}
else
{
    volScalarField thermalObjectiveAdjointSource
    (
        IOobject
        (
            "thermalObjectiveAdjointSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zeroThermalObjectiveAdjointSource",
            dimensionSet(0, -3, 0, -1, 0, 0, 0),
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(thermalObjectiveDerivative, celli)
    {
        thermalObjectiveAdjointSource[celli] =
            thermalObjectiveDerivative[celli]
           /Foam::max(mesh.V()[celli], SMALL);
    }
    thermalObjectiveAdjointSource.correctBoundaryConditions();

    // Use DTEffective (molecular + turbulent) to match the primal thermal
    // operator.  Using DT alone gives an inconsistent adjoint and wrong Tb.
    fvScalarMatrix TEqnb
    (
        fvm::div
        (
           -activeAdjointThermalFlux,
            Tb,
            "div(-phi,Tb)"
        )
      - fvm::laplacian(DTEffective, Tb)
     ==
        thermalObjectiveAdjointSource
    );
    {
        SolverPerformance<scalar> perfTb = TEqnb.solve();
        adjointResidualT = perfTb.finalResidual();
    }
}

// Compute exact face-based Laplacian sensitivity for BOTH adjoint modes.
// The continuous adjoint previously used fvc::grad(T) & fvc::grad(Tb),
// a cell-gradient approximation that introduces O(10-40%) error in FD
// gradient validation.  The face-based expression matches the discrete
// operator exactly.
//
// Tensor upgrade: DTMolecular is now a volSymmTensorField.  The face-level
// Laplacian flux depends only on the direction conductivity e_f^T D e_f.
// The derivative dD/dxh is folded in at the face level via the direction
// derivative e_f^T (dD/dxh) e_f, so thermalDiffusionDerivativeDTCell stores
// dJ/dxh directly (no separate dDTDxh multiplication needed in sensitivity.H).
{
    const labelUList& fOwner = mesh.owner();
    const labelUList& fNeighbour = mesh.neighbour();
    const scalarField& fWeights = mesh.weights().primitiveField();
    const scalarField& fDeltaCoeffs = mesh.deltaCoeffs().primitiveField();
    const scalarField& fMagSf = mesh.magSf().primitiveField();
    // mesh.delta() may return tmp<surfaceVectorField> in some OpenFOAM 7 builds
    const tmp<surfaceVectorField> tDeltaSF = mesh.delta();
    const surfaceVectorField& deltaSF = tDeltaSF();
    const vectorField& fDelta = deltaSF.primitiveField();
    const symmTensorField& dDTDxhInternal = dDTDxh.internalField();

    // --- Internal faces ---
    forAll(fOwner, facei)
    {
        const label own = fOwner[facei];
        const label nei = fNeighbour[facei];

        // Face unit direction (owner -> neighbour)
        const vector& dVec = fDelta[facei];
        const scalar dMag = mag(dVec);
        const vector e_f = dVec/Foam::max(dMag, SMALL);

        // Owner-side direction derivative: e_f^T * dDTDxh[own] * e_f
        const scalar dDnDxh_own = (e_f & dDTDxhInternal[own]) & e_f;

        // Neighbour-side direction derivative
        const scalar dDnDxh_nei = (e_f & dDTDxhInternal[nei]) & e_f;

        // Base face derivative (without D derivative)
        const scalar baseDerivative =
            (Tb[own] - Tb[nei])
           *(T[own] - T[nei])
           *fDeltaCoeffs[facei]*fMagSf[facei];

        thermalDiffusionDerivativeDTCell[own] -=
            fWeights[facei]*dDnDxh_own*baseDerivative;
        thermalDiffusionDerivativeDTCell[nei] -=
            (1.0 - fWeights[facei])*dDnDxh_nei*baseDerivative;
    }

    // --- Boundary faces ---
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const labelUList& pFaceCells = patch.faceCells();
        const scalarField& pMagSf = mesh.magSf().boundaryField()[patchi];
        const scalarField& pDeltaCoeffs =
            mesh.deltaCoeffs().boundaryField()[patchi];
        const vectorField& pDelta = deltaSF.boundaryField()[patchi];

        const scalarField& TbPatch = Tb.boundaryField()[patchi];
        const scalarField& TPatch = T.boundaryField()[patchi];
        // snGrad() is only implemented for physical (non-coupled) patches;
        // processor/cyclic patches use the face-difference form below.
        // Calling it on a coupled patch aborts parallel runs
        // ("coupledFvPatchField::snGrad() Not implemented"), so guard it.
        scalarField snGradTPatch;
        if (!patch.coupled())
        {
            snGradTPatch = T.boundaryField()[patchi].snGrad();
        }

        const symmTensorField dDdxhPatch =
            dDTDxh.boundaryField()[patchi].patchInternalField();

        // Use patch.coupled() instead of isA<processorFvPatch> to avoid
        // additional header dependencies
        const bool isCoupled = patch.coupled();

        // Precompute neighbour fields for coupled patches
        scalarField TbNei, TNei;
        if (isCoupled)
        {
            TbNei = Tb.boundaryField()[patchi].patchNeighbourField();
            TNei = T.boundaryField()[patchi].patchNeighbourField();
        }

        // For physical boundaries, get Tb cell values
        scalarField TbCellPatch;
        if (!isCoupled)
        {
            TbCellPatch = Tb.boundaryField()[patchi].patchInternalField();
        }

        forAll(patch, facei)
        {
            const label celli = pFaceCells[facei];

            // zeroGradient: snGrad=0, contribution is zero
            if (!isCoupled && mag(snGradTPatch[facei]) < SMALL)
                continue;

            // Face unit direction (cell -> face)
            const vector e_f =
                pDelta[facei]/Foam::max(mag(pDelta[facei]), SMALL);

            // Direction derivative
            const scalar dDnDxh =
                (e_f & dDdxhPatch[facei]) & e_f;

            if (isCoupled)
            {
                // processor/cyclic: same structure as internal faces.  The
                // local cell receives only its OWN interpolation-weighted
                // share (mesh.weights() on each side gives that side's
                // weight), matching the serial internal-face split
                // (fWeights / (1-fWeights)) so serial and parallel
                // sensitivities are identical.
                const scalar wLocal =
                    mesh.weights().boundaryField()[patchi][facei];
                const scalar faceDerivative =
                    dDnDxh
                   *(TbPatch[facei] - TbNei[facei])
                   *(TPatch[facei] - TNei[facei])
                   *pDeltaCoeffs[facei]
                   *pMagSf[facei];

                thermalDiffusionDerivativeDTCell[celli] -=
                    wLocal*faceDerivative;
            }
            else
            {
                // Physical boundary:
                // dJ/dD_bf = Tb_cell * snGrad(T) * magSf
                // (snGrad already includes deltaCoeffs)
                const scalar faceDerivative =
                    dDnDxh
                   *TbCellPatch[facei]
                   *snGradTPatch[facei]
                   *pMagSf[facei];

                thermalDiffusionDerivativeDTCell[celli] += faceDerivative;
            }
        }
    }
}