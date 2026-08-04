Info<< "Validating discrete objective derivatives" << endl;

scalarField objectiveTestDirection(mesh.nCells(), 0.0);
forAll(objectiveTestDirection, celli)
{
    objectiveTestDirection[celli] = Foam::sin(0.371*(celli + 1));
}

const scalar thermalObjectiveTestEpsilon = 1.0;
const scalar pressureObjectiveTestEpsilon = 100.0;
scalar thermalAnalytic = 0.0;
scalar thermalPlus = 0.0;
scalar thermalMinus = 0.0;

forAll(thermalObjectiveDerivative, celli)
{
    thermalAnalytic +=
        thermalObjectiveDerivative[celli]*objectiveTestDirection[celli];
}

if (thermalObjectiveType == "maximizeColdOutletTemperature")
{
    // At initialization the forward field has not been solved yet, so phi
    // may be zero.  Validate with area weights here — the derivative
    // formula is identical to the mass-flow-weighted version used at
    // runtime (computeObjective.H), only the weight differs.
    const labelUList& outletCells =
        mesh.boundary()[coldOutletPatchID].faceCells();
    const scalarField& outletAreas =
        mesh.magSf().boundaryField()[coldOutletPatchID];

    forAll(outletCells, facei)
    {
        const label celli = outletCells[facei];
        thermalPlus -= outletAreas[facei]
           *(T[celli] + thermalObjectiveTestEpsilon*objectiveTestDirection[celli])
           /Foam::max(coldOutletArea*thermalReferenceTemperature, SMALL);
        thermalMinus -= outletAreas[facei]
           *(T[celli] - thermalObjectiveTestEpsilon*objectiveTestDirection[celli])
           /Foam::max(coldOutletArea*thermalReferenceTemperature, SMALL);
    }
}

reduce(thermalAnalytic, sumOp<scalar>());
reduce(thermalPlus, sumOp<scalar>());
reduce(thermalMinus, sumOp<scalar>());
const scalar thermalFiniteDifference =
    (thermalPlus - thermalMinus)/(2.0*thermalObjectiveTestEpsilon);
const scalar thermalDerivativeError =
    mag(thermalAnalytic - thermalFiniteDifference)
   /Foam::max
    (
        Foam::max(mag(thermalAnalytic), mag(thermalFiniteDifference)),
        scalar(1e-12)
    );

scalar pressureAnalytic = 0.0;
scalar pressurePlus = 0.0;
scalar pressureMinus = 0.0;
forAll(pressureConstraintDerivative, celli)
{
    pressureAnalytic +=
        pressureConstraintDerivative[celli]*objectiveTestDirection[celli];
}

if (nObjPatch >= 2)
{
    const label inletPatchI = conPatchList[0];
    const label outletPatchI = conPatchList[1];
    scalar inletArea = sum(mesh.magSf().boundaryField()[inletPatchI]);
    scalar outletArea = sum(mesh.magSf().boundaryField()[outletPatchI]);
    reduce(inletArea, sumOp<scalar>());
    reduce(outletArea, sumOp<scalar>());

    if (!p.boundaryField()[inletPatchI].fixesValue())
    {
        const labelUList& inletCells =
            mesh.boundary()[inletPatchI].faceCells();
        const scalarField& inletAreas =
            mesh.magSf().boundaryField()[inletPatchI];
        forAll(inletCells, facei)
        {
            const label celli = inletCells[facei];
            pressurePlus += rhoFluid*inletAreas[facei]
               *(p[celli] + pressureObjectiveTestEpsilon*objectiveTestDirection[celli])
               /Foam::max(pressureDropMaxPa*inletArea, SMALL);
            pressureMinus += rhoFluid*inletAreas[facei]
               *(p[celli] - pressureObjectiveTestEpsilon*objectiveTestDirection[celli])
               /Foam::max(pressureDropMaxPa*inletArea, SMALL);
        }
    }

    if (!p.boundaryField()[outletPatchI].fixesValue())
    {
        const labelUList& outletCells =
            mesh.boundary()[outletPatchI].faceCells();
        const scalarField& outletAreas =
            mesh.magSf().boundaryField()[outletPatchI];
        forAll(outletCells, facei)
        {
            const label celli = outletCells[facei];
            pressurePlus -= rhoFluid*outletAreas[facei]
               *(p[celli] + pressureObjectiveTestEpsilon*objectiveTestDirection[celli])
               /Foam::max(pressureDropMaxPa*outletArea, SMALL);
            pressureMinus -= rhoFluid*outletAreas[facei]
               *(p[celli] - pressureObjectiveTestEpsilon*objectiveTestDirection[celli])
               /Foam::max(pressureDropMaxPa*outletArea, SMALL);
        }
    }
}

reduce(pressureAnalytic, sumOp<scalar>());
reduce(pressurePlus, sumOp<scalar>());
reduce(pressureMinus, sumOp<scalar>());
const scalar pressureFiniteDifference =
    (pressurePlus - pressureMinus)/(2.0*pressureObjectiveTestEpsilon);
const scalar pressureDerivativeError =
    mag(pressureAnalytic - pressureFiniteDifference)
   /Foam::max
    (
        Foam::max(mag(pressureAnalytic), mag(pressureFiniteDifference)),
        scalar(1e-12)
    );

Info<< "Discrete objective derivative errors: thermal="
    << thermalDerivativeError << ", pressure="
    << pressureDerivativeError << endl;

if
(
    thermalDerivativeError > 1e-9
 || pressureDerivativeError > 1e-9
)
{
    FatalErrorInFunction
        << "Discrete objective derivative validation failed: thermal="
        << thermalDerivativeError << ", pressure="
        << pressureDerivativeError << abort(FatalError);
}