scalar diff
(
    const volScalarField& gamma,
    const scalarField& V,
    const volScalarField& designMask,
    const double del,
    const double eta,
    const int n
)
{
    scalar z = 0.0;
    const label nLocal = gamma.primitiveField().size() < label(n)
      ? gamma.primitiveField().size()
      : label(n);

    for (label i = 0; i < nLocal; ++i)
    {
        if (designMask[i] <= 0.5)
        {
            continue;
        }

        scalar projected = 0.0;
        if (gamma[i] <= eta)
        {
            projected = eta*
            (
                Foam::exp(-del*(1.0 - gamma[i]/eta))
              - (1.0 - gamma[i]/eta)*Foam::exp(-del)
            );
        }
        else
        {
            projected = eta + (1.0 - eta)*
            (
                1.0
              - Foam::exp(-del*(gamma[i] - eta)/(1.0 - eta))
              + (gamma[i] - eta)*Foam::exp(-del)/(1.0 - eta)
            );
        }

        z += (gamma[i] - projected)*V[i];
    }

    return z;
}
