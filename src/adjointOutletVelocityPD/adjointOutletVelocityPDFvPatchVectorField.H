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
    adjointOutletVelocityPDFvPatchVectorField

Description
    Adjoint velocity outlet boundary condition for pressure-drop constraint.
    Same analytical treatment as adjointOutletVelocityPower but for the
    pressure-drop adjoint field (Uc, pc, phic).

SourceFiles
    adjointOutletVelocityPDFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointOutletVelocityPDFvPatchVectorField_H
#define adjointOutletVelocityPDFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class adjointOutletVelocityPDFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    TypeName("adjointOutletVelocityPD");

    adjointOutletVelocityPDFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    adjointOutletVelocityPDFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    adjointOutletVelocityPDFvPatchVectorField
    (
        const adjointOutletVelocityPDFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new adjointOutletVelocityPDFvPatchVectorField(*this)
        );
    }

    adjointOutletVelocityPDFvPatchVectorField
    (
        const adjointOutletVelocityPDFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new adjointOutletVelocityPDFvPatchVectorField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

}

#endif
