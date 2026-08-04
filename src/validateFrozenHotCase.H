if (freezeColdFlowForValidation && mmaUpdateEnabled)
{
    FatalErrorInFunction
        << "freezeColdFlowForValidation is a validation-only option and "
        << "cannot be combined with mmaUpdateEnabled=true."
        << abort(FatalError);
}
if (!solveFlowAdjoints && mmaUpdateEnabled)
{
    FatalErrorInFunction
        << "solveFlowAdjoints=false is only valid for forward finite-"
        << "difference cases and cannot be combined with MMA."
        << abort(FatalError);
}if (freezeTurbulenceForValidation && mmaUpdateEnabled)
{
    FatalErrorInFunction
        << "freezeTurbulenceForValidation is a finite-difference diagnostic "
        << "and cannot be combined with MMA." << abort(FatalError);
}

if (frozenHotRegionEnabled)
{
    IOdictionary turbulencePropertiesCheck
    (
        IOobject
        (
            "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word simulationType
    (
        turbulencePropertiesCheck.lookup("simulationType")
    );
    if (flowModel == "laminar" && simulationType != "laminar")
    {
        FatalErrorInFunction
            << "flowModel=laminar requires simulationType laminar. Found '"
            << simulationType << "'." << abort(FatalError);
    }
    if
    (
        flowModel == "incompressibleRANSFrozen"
     && simulationType != "RAS"
    )
    {
        FatalErrorInFunction
            << "flowModel=incompressibleRANSFrozen requires simulationType "
            << "RAS. Found '" << simulationType << "'."
            << abort(FatalError);
    }
    if (flowModel == "incompressibleRANSFrozen")
    {
        const dictionary& rasDict =
            turbulencePropertiesCheck.subDict("RAS");
        const word rasModel(rasDict.lookup("RASModel"));
        if (rasModel != "kOmegaSST")
        {
            FatalErrorInFunction
                << "Frozen-RANS mode requires RASModel kOmegaSST. Found '"
                << rasModel << "'." << abort(FatalError);
        }
    }

    Info<< "Cold-side flowModel=" << flowModel
        << ", frozenTurbulenceAdjoint=" << frozenTurbulenceAdjoint
        << ", ransDirectionValidated=" << ransDirectionValidated
        << ", freezeTurbulenceForValidation="
        << freezeTurbulenceForValidation << endl;


    if
    (
        mmaUpdateEnabled
     && !gradientValidated
     && !directionOnlyOptimizationApproved
    )
    {
        FatalErrorInFunction
            << "MMA requires either gradientValidated=true or the explicit "
            << "directionOnlyOptimizationApproved=true safeguard."
            << abort(FatalError);
    }
    if
    (
        mmaUpdateEnabled
     && !gradientValidated
     && directionOnlyOptimizationApproved
    )
    {
        WarningInFunction
            << "MMA is running under the user-approved direction-only gate. "
            << "Every update must be forward-checked; strict gradient "
            << "validation remains false." << endl;
    }

    if (!T.boundaryField()[hotInletPatchID].fixesValue())
    {
        FatalErrorInFunction
            << "Hot inlet temperature patch '" << hotInletPatchName
            << "' must fix T." << abort(FatalError);
    }

    if (!T.boundaryField()[coldInletPatchID].fixesValue())
    {
        FatalErrorInFunction
            << "Cold inlet temperature patch '" << coldInletPatchName
            << "' must fix T." << abort(FatalError);
    }

    if
    (
        T.boundaryField()[coldOutletPatchID].type()
     != zeroGradientFvPatchScalarField::typeName
    )
    {
        FatalErrorInFunction
            << "The implemented cold-outlet-temperature objective assumes "
            << "zeroGradient T on patch '" << coldOutletPatchName
            << "'. Found type '"
            << T.boundaryField()[coldOutletPatchID].type() << "'."
            << abort(FatalError);
    }

    if
    (
        T.boundaryField()[hotOutletPatchID].type()
     != zeroGradientFvPatchScalarField::typeName
    )
    {
        FatalErrorInFunction
            << "Hot outlet patch '" << hotOutletPatchName
            << "' must use zeroGradient T in frozen-hot mode."
            << abort(FatalError);
    }

    // U is the cold-stream velocity field. It must not introduce cold flow
    // through the hot-stream inlet or outlet patches.
    if
    (
       !U.boundaryField()[hotInletPatchID].fixesValue()
     || !U.boundaryField()[hotOutletPatchID].fixesValue()
    )
    {
        FatalErrorInFunction
            << "Cold velocity U must use fixedValue zero on both hot-side "
            << "inlet and outlet patches." << abort(FatalError);
    }

    scalar maxColdVelocityOnHotPatches = Foam::max
    (
        gMax(mag(U.boundaryField()[hotInletPatchID])),
        gMax(mag(U.boundaryField()[hotOutletPatchID]))
    );
    reduce(maxColdVelocityOnHotPatches, maxOp<scalar>());

    const scalar coldHotPatchVelocityTolerance =
        optProperties.lookupOrDefault<scalar>
        (
            "coldHotPatchVelocityTolerance",
            1e-12
        );
    if (maxColdVelocityOnHotPatches > coldHotPatchVelocityTolerance)
    {
        FatalErrorInFunction
            << "Cold velocity is non-zero on a hot-side patch: max|U|="
            << maxColdVelocityOnHotPatches << ", tolerance="
            << coldHotPatchVelocityTolerance << abort(FatalError);
    }

    label invalidPatchOwnerCells = 0;

    const labelUList& hotInletCells =
        mesh.boundary()[hotInletPatchID].faceCells();
    const labelUList& hotOutletCells =
        mesh.boundary()[hotOutletPatchID].faceCells();
    const labelUList& coldInletCells =
        mesh.boundary()[coldInletPatchID].faceCells();
    const labelUList& coldOutletCells =
        mesh.boundary()[coldOutletPatchID].faceCells();

    forAll(hotInletCells, facei)
    {
        invalidPatchOwnerCells +=
            hotRegionMask[hotInletCells[facei]] > 0.5 ? 0 : 1;
    }
    forAll(hotOutletCells, facei)
    {
        invalidPatchOwnerCells +=
            hotRegionMask[hotOutletCells[facei]] > 0.5 ? 0 : 1;
    }
    forAll(coldInletCells, facei)
    {
        invalidPatchOwnerCells +=
            coldRegionMask[coldInletCells[facei]] > 0.5 ? 0 : 1;
    }
    forAll(coldOutletCells, facei)
    {
        invalidPatchOwnerCells +=
            coldRegionMask[coldOutletCells[facei]] > 0.5 ? 0 : 1;
    }

    reduce(invalidPatchOwnerCells, sumOp<label>());
    if (invalidPatchOwnerCells > 0)
    {
        FatalErrorInFunction
            << invalidPatchOwnerCells
            << " inlet/outlet faces are attached to cells in the wrong "
            << "region." << abort(FatalError);
    }

    if (!Pstream::parRun())
    {
        if
        (
            p.needReference()
         && pRefCell >= 0
         && coldRegionMask[pRefCell] < 0.5
        )
        {
            FatalErrorInFunction
                << "Primal pressure reference cell " << pRefCell
                << " is outside the cold region." << abort(FatalError);
        }

        if
        (
            (
                pa.needReference()
             && paRefCell >= 0
             && coldRegionMask[paRefCell] < 0.5
            )
         || (
                pb.needReference()
             && pbRefCell >= 0
             && coldRegionMask[pbRefCell] < 0.5
            )
         || (
                pc.needReference()
             && pcRefCell >= 0
             && coldRegionMask[pcRefCell] < 0.5
            )
        )
        {
            FatalErrorInFunction
                << "At least one adjoint pressure reference cell is outside "
                << "the cold region." << abort(FatalError);
        }
    }
}
