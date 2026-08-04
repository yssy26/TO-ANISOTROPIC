/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocityHeatFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"

Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}

Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
}

Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const adjointOutletVelocityHeatFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}

Foam::adjointOutletVelocityHeatFvPatchVectorField::
adjointOutletVelocityHeatFvPatchVectorField
(
    const adjointOutletVelocityHeatFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
{}

void Foam::adjointOutletVelocityHeatFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& adjointFlux =
        patch().lookupPatchField<surfaceScalarField, scalar>("phib");
    const fvPatchField<vector>& adjointVelocity =
        patch().lookupPatchField<volVectorField, vector>("Ub");
    const fvsPatchField<scalar>& primalFlux =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const dictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");
    const dimensionedScalar nu(transportProperties.lookup("nu"));

    const scalarField& deltaInv = patch().deltaCoeffs();
    const scalarField primalNormalVelocity =
        primalFlux/patch().magSf();

    const vectorField adjointNeighbour =
        adjointVelocity.patchInternalField();
    const vectorField adjointNeighbourNormal =
        (adjointNeighbour & patch().nf())*patch().nf();
    const vectorField adjointNeighbourTangential =
        adjointNeighbour - adjointNeighbourNormal;

    scalarField denominator =
        primalNormalVelocity + nu.value()*deltaInv;
    forAll(denominator, facei)
    {
        if (mag(denominator[facei]) < SMALL)
        {
            denominator[facei] =
                denominator[facei] < 0.0 ? -SMALL : SMALL;
        }
    }

    const vectorField adjointTangential =
        nu.value()*deltaInv*adjointNeighbourTangential/denominator;
    const vectorField adjointNormal =
        (adjointFlux*patch().Sf())/sqr(patch().magSf());

    operator==(adjointTangential + adjointNormal);
    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::adjointOutletVelocityHeatFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocityHeatFvPatchVectorField
    );
}

// ************************************************************************* //
