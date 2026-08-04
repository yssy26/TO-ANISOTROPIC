// Ensure MMA cannot feed altered fixed-region values into the PDE filter.
#include "enforceFixedDesignRegions.H"

// Restrict the Helmholtz filter to faces whose owner and neighbour are both
// active design cells. Outside the design region the equation reduces to
// xp=x, so prescribed solid/fluid values cannot diffuse into the design.
fvScalarMatrix Eqn4
(
    fvm::laplacian(designFilterFaceMask, xp)
  - fvm::Sp(b, xp)
  + x*b
);
Eqn4.solve();

#include "enforceFixedDesignRegions.H"

// Heaviside projection. The volume-preserving threshold is computed over the
// active design region only.
del = mmaUpdateEnabled
    ? Foam::min(0.2*opt, 60.0)
    : validationProjectionBeta;
eta0 = 0.0001;
eta1 = 0.9999;
projectionThresholdAdaptive = false;

y0 = diff(xp, mesh.V(), designMask, del, eta0, n);
scalar y1 = diff(xp, mesh.V(), designMask, del, eta1, n);
reduce(y0, sumOp<scalar>());
reduce(y1, sumOp<scalar>());

if (y0*y1 > 0.0)
{
    WarningInFunction
        << "Projection volume equation is not bracketed: y(eta0)=" << y0
        << ", y(eta1)=" << y1
        << ". Using eta=0.5 for this iteration." << endl;
    eta5 = 0.5;
}
else
{
    projectionThresholdAdaptive = del > SMALL;
    do
    {
        eta5 = 0.5*(eta0 + eta1);
        y5 = diff(xp, mesh.V(), designMask, del, eta5, n);
        reduce(y5, sumOp<scalar>());

        if (y0*y5 < 0.0)
        {
            eta1 = eta5;
            y1 = y5;
        }
        else
        {
            eta0 = eta5;
            y0 = y5;
        }
    }
    while ((eta1 - eta0) > projectionEtaTolerance);
}

for (i = 0; i < n; ++i)
{
    if (designMask[i] <= 0.5)
    {
        continue;
    }

    if (xp[i] <= eta5)
    {
        xh[i] = eta5*
        (
            Foam::exp(-del*(1.0 - xp[i]/eta5))
          - (1.0 - xp[i]/eta5)*Foam::exp(-del)
        );
    }
    else
    {
        xh[i] = eta5 + (1.0 - eta5)*
        (
            1.0
          - Foam::exp(-del*(xp[i] - eta5)/(1.0 - eta5))
          + (xp[i] - eta5)*Foam::exp(-del)/(1.0 - eta5)
        );
    }
}

#include "enforceFixedDesignRegions.H"
