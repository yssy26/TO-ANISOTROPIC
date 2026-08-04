#!/usr/bin/env python3
"""Deterministic unit tests for the anisotropic material interpolation.

This test is independent of OpenFOAM. It validates the algebra used in
readThermalProperties.H and updateMaterialProperties.H before case-level tests.
"""

from __future__ import annotations

import math
import random
from typing import Iterable, List

Matrix = List[List[float]]


def eye(scale: float) -> Matrix:
    return [[scale, 0.0, 0.0], [0.0, scale, 0.0], [0.0, 0.0, scale]]


def add(a: Matrix, b: Matrix) -> Matrix:
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def sub(a: Matrix, b: Matrix) -> Matrix:
    return [[a[i][j] - b[i][j] for j in range(3)] for i in range(3)]


def mul(a: Matrix, scalar: float) -> Matrix:
    return [[a[i][j] * scalar for j in range(3)] for i in range(3)]


def max_abs(a: Matrix) -> float:
    return max(abs(value) for row in a for value in row)


def relative_error(a: Matrix, b: Matrix) -> float:
    denominator = max(max_abs(a), max_abs(b), 1.0e-30)
    return max_abs(sub(a, b)) / denominator


def g(xh: float, qu: float) -> float:
    return xh * (1.0 + qu) / (qu + xh)


def dg_dxh(xh: float, qu: float) -> float:
    return (1.0 + qu) * qu / (qu + xh) ** 2


def diffusivity(xh: float, qu: float, ks: Matrix, kf: float, rhoc: float) -> Matrix:
    return mul(add(ks, mul(sub(eye(kf), ks), g(xh, qu))), 1.0 / rhoc)


def diffusivity_derivative(
    xh: float, qu: float, ks: Matrix, kf: float, rhoc: float
) -> Matrix:
    return mul(sub(eye(kf), ks), dg_dxh(xh, qu) / rhoc)


def directional_conductivity(k: Matrix, direction: Iterable[float]) -> float:
    e = list(direction)
    norm = math.sqrt(sum(component * component for component in e))
    if norm <= 0.0:
        raise ValueError("direction must be non-zero")
    e = [component / norm for component in e]
    return sum(e[i] * k[i][j] * e[j] for i in range(3) for j in range(3))


def test_isotropic_regression() -> None:
    rng = random.Random(20260804)
    worst = 0.0
    for _ in range(1000):
        xh = rng.uniform(0.0, 1.0)
        qu = rng.uniform(0.01, 0.5)
        ks_scalar = rng.uniform(1.0, 100.0)
        kf = rng.uniform(0.01, 1.0)
        rhoc = rng.uniform(100.0, 5000.0)
        tensor_value = diffusivity(xh, qu, eye(ks_scalar), kf, rhoc)
        scalar_value = (ks_scalar + (kf - ks_scalar) * g(xh, qu)) / rhoc
        reference = eye(scalar_value)
        worst = max(worst, relative_error(tensor_value, reference))
    assert worst < 1.0e-14, f"isotropic regression error={worst:.3e}"
    print(f"PASS isotropic regression: worst relative error={worst:.3e}")


def test_tensor_derivative() -> None:
    rng = random.Random(260804)
    worst = 0.0
    for _ in range(500):
        diagonal = [rng.uniform(5.0, 50.0) for _ in range(3)]
        xy = rng.uniform(-0.05, 0.05) * min(diagonal[0], diagonal[1])
        xz = rng.uniform(-0.05, 0.05) * min(diagonal[0], diagonal[2])
        yz = rng.uniform(-0.05, 0.05) * min(diagonal[1], diagonal[2])
        ks = [
            [diagonal[0], xy, xz],
            [xy, diagonal[1], yz],
            [xz, yz, diagonal[2]],
        ]
        xh = rng.uniform(0.05, 0.95)
        qu = rng.uniform(0.02, 0.4)
        kf = rng.uniform(0.02, 1.0)
        rhoc = rng.uniform(200.0, 3000.0)
        step = 1.0e-7
        plus = diffusivity(xh + step, qu, ks, kf, rhoc)
        minus = diffusivity(xh - step, qu, ks, kf, rhoc)
        finite_difference = mul(sub(plus, minus), 1.0 / (2.0 * step))
        analytic = diffusivity_derivative(xh, qu, ks, kf, rhoc)
        worst = max(worst, relative_error(finite_difference, analytic))
    assert worst < 2.0e-7, f"tensor derivative error={worst:.3e}"
    print(f"PASS tensor derivative: worst relative error={worst:.3e}")


def test_directional_conductivity() -> None:
    k = [[21.6, 0.0, 0.0], [0.0, 21.6, 0.0], [0.0, 0.0, 16.6]]
    values = {
        "x": directional_conductivity(k, (1.0, 0.0, 0.0)),
        "y": directional_conductivity(k, (0.0, 1.0, 0.0)),
        "z": directional_conductivity(k, (0.0, 0.0, 1.0)),
    }
    assert abs(values["x"] - 21.6) < 1.0e-14
    assert abs(values["y"] - 21.6) < 1.0e-14
    assert abs(values["z"] - 16.6) < 1.0e-14
    assert values["x"] > values["z"]
    print(f"PASS directional conductivity: {values}")


def test_material_endpoints() -> None:
    ks = [[21.6, 0.0, 0.0], [0.0, 21.6, 0.0], [0.0, 0.0, 16.6]]
    kf = 0.0469
    rhoc = 619.2492
    qu = 0.1
    solid = diffusivity(0.0, qu, ks, kf, rhoc)
    expected_solid = mul(ks, 1.0 / rhoc)
    assert relative_error(solid, expected_solid) < 1.0e-14

    fluid = diffusivity(1.0, qu, ks, kf, rhoc)
    expected_fluid = eye(kf / rhoc)
    assert relative_error(fluid, expected_fluid) < 5.0e-14
    print("PASS material endpoints")


def main() -> None:
    test_isotropic_regression()
    test_tensor_derivative()
    test_directional_conductivity()
    test_material_endpoints()
    print("ANISOTROPIC MATERIAL MODEL UNIT TESTS: PASS")


if __name__ == "__main__":
    main()
