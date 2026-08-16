#!/usr/bin/env python3
"""Algebraic regression for serial/processor-face sensitivity assembly.

This test does not replace the OpenFOAM A7 case.  It verifies that splitting an
internal face into two processor patches produces exactly the same local-cell
contributions when each side uses patchInternalField/patchNeighbourField and
complementary OpenFOAM interpolation weights.
"""

from __future__ import annotations

import math
import random


def serial_contributions(
    tb_owner: float,
    tb_neighbour: float,
    t_owner: float,
    t_neighbour: float,
    d_owner: float,
    d_neighbour: float,
    owner_weight: float,
    geometric_factor: float,
) -> tuple[float, float]:
    base = (
        (tb_owner - tb_neighbour)
        * (t_owner - t_neighbour)
        * geometric_factor
    )
    return (
        -owner_weight * d_owner * base,
        -(1.0 - owner_weight) * d_neighbour * base,
    )


def processor_contributions(
    tb_owner: float,
    tb_neighbour: float,
    t_owner: float,
    t_neighbour: float,
    d_owner: float,
    d_neighbour: float,
    owner_weight: float,
    geometric_factor: float,
) -> tuple[float, float]:
    # Owner rank: local=owner, remote=neighbour, local weight=owner_weight.
    owner_base = (
        (tb_owner - tb_neighbour)
        * (t_owner - t_neighbour)
        * geometric_factor
    )
    owner_value = -owner_weight * d_owner * owner_base

    # Neighbour rank: local/remote ordering and face direction are reversed.
    # The product of the two jumps is unchanged, and the local interpolation
    # weight is complementary.
    neighbour_base = (
        (tb_neighbour - tb_owner)
        * (t_neighbour - t_owner)
        * geometric_factor
    )
    neighbour_value = -(1.0 - owner_weight) * d_neighbour * neighbour_base
    return owner_value, neighbour_value


def old_evaluated_patch_owner_value(
    tb_owner: float,
    tb_neighbour: float,
    t_owner: float,
    t_neighbour: float,
    d_owner: float,
    owner_weight: float,
    geometric_factor: float,
) -> float:
    # Models the removed implementation: boundaryField()[patch] is an
    # evaluated/interpolated patch value, not the local cell-centre state.
    tb_patch = owner_weight * tb_owner + (1.0 - owner_weight) * tb_neighbour
    t_patch = owner_weight * t_owner + (1.0 - owner_weight) * t_neighbour
    base = (
        (tb_patch - tb_neighbour)
        * (t_patch - t_neighbour)
        * geometric_factor
    )
    return -owner_weight * d_owner * base


def relative_error(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-30)


def main() -> None:
    rng = random.Random(20260805)
    worst_corrected = 0.0
    old_defect_observed = False

    for _ in range(1000):
        tb_owner = rng.uniform(-3.0, 3.0)
        tb_neighbour = rng.uniform(-3.0, 3.0)
        t_owner = rng.uniform(250.0, 900.0)
        t_neighbour = rng.uniform(250.0, 900.0)
        d_owner = rng.uniform(0.01, 20.0)
        d_neighbour = rng.uniform(0.01, 20.0)
        owner_weight = rng.uniform(0.05, 0.95)
        geometric_factor = rng.uniform(0.01, 100.0)

        serial = serial_contributions(
            tb_owner,
            tb_neighbour,
            t_owner,
            t_neighbour,
            d_owner,
            d_neighbour,
            owner_weight,
            geometric_factor,
        )
        parallel = processor_contributions(
            tb_owner,
            tb_neighbour,
            t_owner,
            t_neighbour,
            d_owner,
            d_neighbour,
            owner_weight,
            geometric_factor,
        )

        for serial_value, parallel_value in zip(serial, parallel):
            worst_corrected = max(
                worst_corrected,
                relative_error(serial_value, parallel_value),
            )

        old_owner = old_evaluated_patch_owner_value(
            tb_owner,
            tb_neighbour,
            t_owner,
            t_neighbour,
            d_owner,
            owner_weight,
            geometric_factor,
        )
        if relative_error(old_owner, serial[0]) > 1.0e-3:
            old_defect_observed = True

    if worst_corrected > 5.0e-15:
        raise SystemExit(
            f"FAIL corrected processor split: worst relative error="
            f"{worst_corrected:.3e}"
        )
    if not old_defect_observed:
        raise SystemExit("FAIL old evaluated-patch defect was not detected")

    print(
        "PASS processor/internal sensitivity equivalence: "
        f"worst relative error={worst_corrected:.3e}"
    )
    print("PASS regression detects evaluated-patch-value defect")
    print("PROCESSOR PATCH SENSITIVITY UNIT TEST: PASS")


if __name__ == "__main__":
    main()
