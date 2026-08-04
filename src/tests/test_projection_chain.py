#!/usr/bin/env python3
"""Finite-difference check for the adaptive volume-preserving projection."""

import math


def project(value, beta, eta):
    decay = math.exp(-beta)
    if value <= eta:
        return eta * (
            math.exp(-beta * (1.0 - value / eta))
            - (1.0 - value / eta) * decay
        )
    return eta + (1.0 - eta) * (
        1.0
        - math.exp(-beta * (value - eta) / (1.0 - eta))
        + (value - eta) * decay / (1.0 - eta)
    )


def projection_derivatives(value, beta, eta):
    decay = math.exp(-beta)
    if value <= eta:
        exponential = math.exp(-beta * (1.0 - value / eta))
        dx = beta * exponential + decay
        deta = exponential * (1.0 - beta * value / eta) - decay
    else:
        exponential = math.exp(-beta * (value - eta) / (1.0 - eta))
        dx = beta * exponential + decay
        deta = (
            exponential
            * (1.0 - beta * (1.0 - value) / (1.0 - eta))
            - decay
        )
    return dx, deta


def solve_eta(values, volumes, beta):
    lower = 1.0e-8
    upper = 1.0 - 1.0e-8

    def residual(eta):
        return sum(
            volume * (value - project(value, beta, eta))
            for value, volume in zip(values, volumes)
        )

    lower_residual = residual(lower)
    for _ in range(100):
        middle = 0.5 * (lower + upper)
        middle_residual = residual(middle)
        if lower_residual * middle_residual < 0.0:
            upper = middle
        else:
            lower = middle
            lower_residual = middle_residual
        if upper - lower < 1.0e-13:
            break
    return 0.5 * (lower + upper)


def objective(values, volumes, weights, beta):
    eta = solve_eta(values, volumes, beta)
    return sum(
        volume * weight * project(value, beta, eta)
        for value, volume, weight in zip(values, volumes, weights)
    )


def analytic_gradient(values, volumes, weights, beta):
    eta = solve_eta(values, volumes, beta)
    derivatives = [
        projection_derivatives(value, beta, eta) for value in values
    ]
    denominator = sum(
        volume * deta
        for volume, (_, deta) in zip(volumes, derivatives)
    )
    correction = sum(
        volume * weight * deta
        for volume, weight, (_, deta) in zip(volumes, weights, derivatives)
    ) / denominator
    return [
        volume * (weight * dx + correction * (1.0 - dx))
        for volume, weight, (dx, _) in zip(volumes, weights, derivatives)
    ]


def main():
    values = [0.13, 0.29, 0.47, 0.66, 0.88]
    volumes = [0.7, 1.1, 0.9, 1.6, 0.8]
    weights = [1.2, -0.4, 0.8, 2.0, -1.1]
    beta = 8.0
    epsilon = 1.0e-6

    analytic = analytic_gradient(values, volumes, weights, beta)
    finite_difference = []
    for index in range(len(values)):
        plus = list(values)
        minus = list(values)
        plus[index] += epsilon
        minus[index] -= epsilon
        finite_difference.append(
            (
                objective(plus, volumes, weights, beta)
                - objective(minus, volumes, weights, beta)
            )
            / (2.0 * epsilon)
        )

    max_relative_error = max(
        abs(adjoint - fd) / max(1.0, abs(adjoint), abs(fd))
        for adjoint, fd in zip(analytic, finite_difference)
    )

    eta = solve_eta(values, volumes, beta)
    projected_volume = sum(
        volume * project(value, beta, eta)
        for value, volume in zip(values, volumes)
    )
    raw_volume = sum(
        volume * value for value, volume in zip(values, volumes)
    )

    if max_relative_error > 2.0e-7:
        raise SystemExit(
            "projection gradient check failed: "
            f"relative error={max_relative_error:.3e}"
        )
    if abs(projected_volume - raw_volume) > 1.0e-10:
        raise SystemExit("volume-preserving projection check failed")

    print(
        "PASS adaptive projection chain rule: "
        f"relative error={max_relative_error:.3e}"
    )


if __name__ == "__main__":
    main()
