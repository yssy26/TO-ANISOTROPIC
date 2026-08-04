/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "adjointInletPressurePDFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "RASModel.H"

Foam::adjointInletPressurePDFvPatchScalarField::
adjointInletPressurePDFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::adjointInletPressurePDFvPatchScalarField::
adjointInletPressurePDFvPatchScalarField
(
    const adjointInletPressurePDFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adjointInletPressurePDFvPatchScalarField::
adjointInletPressurePDFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::adjointInletPressurePDFvPatchScalarField::
adjointInletPressurePDFvPatchScalarField
(
    const adjointInletPressurePDFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


void Foam::adjointInletPressurePDFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& optDict =
        db().lookupObject<IOdictionary>("optProperties");
    scalar pressureDropMaxPa =
        optDict.lookupOrDefault<scalar>("pressureDropMaxPa", 40.0);
    scalar rhoFluidVal =
        optDict.lookupOrDefault<scalar>("rhoFluid", 1.0);

    scalar inletArea = gSum(patch().magSf());

    scalarField pdSensInlet(patch().size(), rhoFluidVal / Foam::max(pressureDropMaxPa * inletArea, SMALL));

    operator==(pdSensInlet);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointInletPressurePDFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointInletPressurePDFvPatchScalarField
    );
}
