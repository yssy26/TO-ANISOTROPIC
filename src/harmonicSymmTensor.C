#include "harmonicSymmTensor.H"
#include "fvMesh.H"

namespace Foam
{

defineTypeNameAndDebug(harmonicSymmTensor, 0);

surfaceInterpolationScheme<symmTensor>::
    addMeshConstructorToTable<harmonicSymmTensor>
    addHarmonicSymmTensorMeshConstructorToTable_;

surfaceInterpolationScheme<symmTensor>::
    addMeshFluxConstructorToTable<harmonicSymmTensor>
    addHarmonicSymmTensorMeshFluxConstructorToTable_;


tmp<surfaceScalarField> harmonicSymmTensor::weights
(
    const volSymmTensorField&
) const
{
    return this->mesh().surfaceInterpolation::weights();
}


tmp<surfaceSymmTensorField> harmonicSymmTensor::interpolate
(
    const volSymmTensorField& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& lambdaField =
        mesh.surfaceInterpolation::weights();
    const scalarField& lambda = lambdaField.primitiveField();
    const symmTensorField& vfi = vf.primitiveField();

    tmp<surfaceSymmTensorField> tFace
    (
        new surfaceSymmTensorField
        (
            IOobject
            (
                "harmonicInterpolate(" + vf.name() + ')',
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            vf.dimensions()
        )
    );
    surfaceSymmTensorField& faceField = tFace.ref();
    symmTensorField& faceInternal = faceField.primitiveFieldRef();

    forAll(owner, facei)
    {
        const symmTensor inverseResistance =
            (1.0 - lambda[facei])*inv(vfi[owner[facei]])
          + lambda[facei]*inv(vfi[neighbour[facei]]);
        faceInternal[facei] = inv(inverseResistance);
    }

    surfaceSymmTensorField::Boundary& faceBoundary =
        faceField.boundaryFieldRef();

    forAll(faceBoundary, patchi)
    {
        if (vf.boundaryField()[patchi].coupled())
        {
            const scalarField& patchLambda =
                lambdaField.boundaryField()[patchi];
            const symmTensorField local =
                vf.boundaryField()[patchi].patchInternalField();
            const symmTensorField neighbourField =
                vf.boundaryField()[patchi].patchNeighbourField();

            forAll(faceBoundary[patchi], facei)
            {
                const symmTensor inverseResistance =
                    (1.0 - patchLambda[facei])*inv(local[facei])
                  + patchLambda[facei]*inv(neighbourField[facei]);
                faceBoundary[patchi][facei] = inv(inverseResistance);
            }
        }
        else
        {
            faceBoundary[patchi] = vf.boundaryField()[patchi];
        }
    }

    return tFace;
}

} // End namespace Foam
