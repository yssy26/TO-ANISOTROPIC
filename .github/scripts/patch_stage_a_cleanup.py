#!/usr/bin/env python3
"""Apply the two post-Stage-A diagnostic cleanups exactly once."""

from pathlib import Path


def patch_report_label() -> None:
    path = Path("src/tests/stage_a_acceptance.py")
    text = path.read_text(encoding="utf-8")
    old = 'f"{axis}: flux={flux_error:.3e}, R2={linearity:.6f}"'
    new = (
        'f"{axis}: relative flux error={flux_error:.3e}, "\n'
        '                f"R2={linearity:.6f}"'
    )
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"expected one report label occurrence, found {count}")
    path.write_text(text.replace(old, new), encoding="utf-8")


def patch_partition_invariant_directions() -> None:
    path = Path("src/validateFrozenGradient.H")
    text = path.read_text(encoding="utf-8")

    old_direction = "Foam::sin(0.271*(dirI + 1)*(celli + 1) + 0.5);"
    count = text.count(old_direction)
    if count != 2:
        raise RuntimeError(
            f"expected two local-cell direction occurrences, found {count}"
        )

    marker = "    for (label dirI = 0; dirI < nDirF1; ++dirI)\n"
    marker_count = text.count(marker)
    if marker_count != 1:
        raise RuntimeError(
            f"expected one F1 direction-loop marker, found {marker_count}"
        )

    helper = """    // Partition-invariant synthetic FD directions. Local cell labels
    // restart on every MPI rank and therefore cannot define a serial/parallel
    // comparison direction. Normalize physical cell-centre coordinates by
    // the global mesh bounds so the same physical cell receives the same
    // value for every decomposition.
    const boundBox fdDirectionBounds(mesh.points(), true);
    const vector fdDirectionMin = fdDirectionBounds.min();
    const vector fdDirectionSpan =
        fdDirectionBounds.max() - fdDirectionMin;

    auto partitionInvariantFDValue =
        [&](const label celli, const label dirI) -> scalar
    {
        const vector& centre = mesh.C()[celli];
        const scalar xi =
            (centre.x() - fdDirectionMin.x())
           /Foam::max(mag(fdDirectionSpan.x()), SMALL);
        const scalar eta =
            (centre.y() - fdDirectionMin.y())
           /Foam::max(mag(fdDirectionSpan.y()), SMALL);
        const scalar zeta =
            (centre.z() - fdDirectionMin.z())
           /Foam::max(mag(fdDirectionSpan.z()), SMALL);
        const scalar mode = scalar(dirI + 1);

        return Foam::sin
        (
            mode*(2.137*xi + 3.071*eta + 4.113*zeta) + 0.5
        );
    };

"""

    text = text.replace(marker, helper + marker)
    text = text.replace(
        old_direction,
        "partitionInvariantFDValue(celli, dirI);",
    )

    if old_direction in text:
        raise RuntimeError("local-cell direction expression remains")
    replacement_count = text.count("partitionInvariantFDValue(celli, dirI);")
    if replacement_count != 2:
        raise RuntimeError(
            f"expected two partition-invariant calls, found {replacement_count}"
        )

    path.write_text(text, encoding="utf-8")


def main() -> None:
    patch_report_label()
    patch_partition_invariant_directions()


if __name__ == "__main__":
    main()
