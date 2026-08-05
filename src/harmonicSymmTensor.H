#ifndef harmonicSymmTensor_H
#define harmonicSymmTensor_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

// Harmonic face interpolation for symmetric positive-definite tensor fields.
//
// For an internal face with OpenFOAM linear owner weight lambda,
//
//   K_f = inv((1-lambda) inv(K_owner) + lambda inv(K_neighbour)).
//
// On an equally spaced mesh this reduces to
//
//   K_f = 2 inv(inv(K_owner) + inv(K_neighbour)).
//
// This scheme was introduced for the Stage-A three-layer conduction test.
// The current topology sensitivity in AdjHeatTransfer.H differentiates the
// standard linear tensor interpolation, so harmonicSymmTensor must remain a
// validation/forward-verification scheme until its coefficient derivative is
// added to the adjoint.
class harmonicSymmTensor
:
    public surfaceInterpolationScheme<symmTensor>
{
public:

    TypeName("harmonicSymmTensor");

    harmonicSymmTensor(const fvMesh& mesh)
    :
        surfaceInterpolationScheme<symmTensor>(mesh)
    {}

    harmonicSymmTensor(const fvMesh& mesh, Istream&)
    :
        surfaceInterpolationScheme<symmTensor>(mesh)
    {}

    harmonicSymmTensor
    (
        const fvMesh& mesh,
        const surfaceScalarField&,
        Istream&
    )
    :
        surfaceInterpolationScheme<symmTensor>(mesh)
    {}

    virtual tmp<surfaceScalarField> weights
    (
        const volSymmTensorField&
    ) const;

    virtual tmp<surfaceSymmTensorField> interpolate
    (
        const volSymmTensorField& vf
    ) const;

    void operator=(const harmonicSymmTensor&) = delete;
};

} // End namespace Foam

#endif
