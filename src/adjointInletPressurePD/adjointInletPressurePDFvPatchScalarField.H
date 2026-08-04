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
    adjointInletPressurePDFvPatchScalarField

Description
    Adjoint pressure inlet boundary condition for pressure-drop constraint.
    At the inlet (where primal U = fixedValue):
        pc = rhoFluid / (pressureDropMaxPa * A_inlet)
    This encodes the sensitivity of the pressure-drop constraint g_pd = Δp/Δp_max - 1
    with respect to the inlet pressure.

SourceFiles
    adjointInletPressurePDFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointInletPressurePDFvPatchScalarField_H
#define adjointInletPressurePDFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class adjointInletPressurePDFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    TypeName("adjointInletPressurePD");

    adjointInletPressurePDFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    adjointInletPressurePDFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    adjointInletPressurePDFvPatchScalarField
    (
        const adjointInletPressurePDFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new adjointInletPressurePDFvPatchScalarField(*this)
        );
    }

    adjointInletPressurePDFvPatchScalarField
    (
        const adjointInletPressurePDFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new adjointInletPressurePDFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif
