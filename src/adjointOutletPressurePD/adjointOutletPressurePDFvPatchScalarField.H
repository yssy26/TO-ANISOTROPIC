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

Class
    adjointOutletPressurePDFvPatchScalarField

Description
    Adjoint pressure outlet boundary condition for pressure-drop constraint.
    At the outlet (where primal p = fixedValue 0):
        pc = -rhoFluid / (pressureDropMaxPa * A_outlet)
    This encodes the sensitivity of the pressure-drop constraint g_pd = Δp/Δp_max - 1
    with respect to the outlet pressure.

SourceFiles
    adjointOutletPressurePDFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointOutletPressurePDFvPatchScalarField_H
#define adjointOutletPressurePDFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class adjointOutletPressurePDFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    TypeName("adjointOutletPressurePD");

    adjointOutletPressurePDFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    adjointOutletPressurePDFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    adjointOutletPressurePDFvPatchScalarField
    (
        const adjointOutletPressurePDFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new adjointOutletPressurePDFvPatchScalarField(*this)
        );
    }

    adjointOutletPressurePDFvPatchScalarField
    (
        const adjointOutletPressurePDFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new adjointOutletPressurePDFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif
