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

#include "adjointOutletPressurePDFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "RASModel.H"

Foam::adjointOutletPressurePDFvPatchScalarField::
adjointOutletPressurePDFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::adjointOutletPressurePDFvPatchScalarField::
adjointOutletPressurePDFvPatchScalarField
(
    const adjointOutletPressurePDFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adjointOutletPressurePDFvPatchScalarField::
adjointOutletPressurePDFvPatchScalarField
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


Foam::adjointOutletPressurePDFvPatchScalarField::
adjointOutletPressurePDFvPatchScalarField
(
    const adjointOutletPressurePDFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


void Foam::adjointOutletPressurePDFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar nu(transportProperties.lookup("nu"));

    const dictionary& optDict =
        db().lookupObject<IOdictionary>("optProperties");
    scalar pressureDropMaxPa =
        optDict.lookupOrDefault<scalar>("pressureDropMaxPa", 40.0);
    scalar rhoFluidVal =
        optDict.lookupOrDefault<scalar>("rhoFluid", 1.0);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvsPatchField<scalar>& phicp =
        patch().lookupPatchField<surfaceScalarField, scalar>("phic");

    const fvPatchField<vector>& Ucp =
        patch().lookupPatchField<volVectorField, vector>("Uc");

    scalarField Up_n = phip / patch().magSf();
    scalarField Ucp_n = phicp / patch().magSf();

    const scalarField& deltainv = patch().deltaCoeffs();

    scalarField Ucneigh_n = (Ucp.patchInternalField() & patch().nf());

    scalar outletArea = gSum(patch().magSf());

    scalarField pdSensOutlet(patch().size(), -rhoFluidVal / Foam::max(pressureDropMaxPa * outletArea, SMALL));

    operator==(pdSensOutlet + (Up_n * Ucp_n) + 2*nu.value()*deltainv*(Ucp_n - Ucneigh_n));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointOutletPressurePDFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointOutletPressurePDFvPatchScalarField
    );
}
