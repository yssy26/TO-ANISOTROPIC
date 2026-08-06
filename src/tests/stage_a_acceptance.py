#!/usr/bin/env python3
"""Evaluate Stage A anisotropic-conduction validation results.

The OpenFOAM runs are performed externally. This script consumes one JSON
summary, applies fixed acceptance criteria, writes a Markdown report, and exits
with status 0 only when every mandatory check passes.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple


def finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def require_number(data: Dict[str, Any], key: str) -> float:
    value = data.get(key)
    if not finite_number(value):
        raise ValueError(f"missing or non-finite numeric field: {key}")
    return float(value)


def status(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def evaluate(payload: Dict[str, Any]) -> Tuple[bool, List[str]]:
    rows: List[str] = []
    overall = True

    iso = payload.get("isotropicRegression", {})
    iso_error = require_number(iso, "maxRelativeFieldError")
    iso_files = int(require_number(iso, "filesCompared"))
    iso_pass = iso_error <= 1.0e-12 and iso_files >= 1
    rows.append(
        f"| A1 Isotropic regression | {status(iso_pass)} | "
        f"max relative field error={iso_error:.3e}; files={iso_files} |"
    )
    overall &= iso_pass

    axis_entries = payload.get("axisConduction", [])
    axis_by_name = {str(item.get("axis")): item for item in axis_entries}
    axis_pass = set(axis_by_name) >= {"x", "y", "z"}
    axis_details: List[str] = []
    for axis in ("x", "y", "z"):
        item = axis_by_name.get(axis, {})
        try:
            flux_error = require_number(item, "relativeFluxError")
            linearity = require_number(item, "temperatureLinearityR2")
            case_pass = flux_error <= 0.01 and linearity >= 0.999
            axis_pass &= case_pass
            axis_details.append(
                f"{axis}: relative flux error={flux_error:.3e}, "
                f"R2={linearity:.6f}"
            )
        except ValueError:
            axis_pass = False
            axis_details.append(f"{axis}: missing")
    rows.append(
        f"| A2 Axis conduction x/y/z | {status(axis_pass)} | "
        + "; ".join(axis_details)
        + " |"
    )
    overall &= axis_pass

    layers = payload.get("threeLayerConduction", {})
    layer_flux_error = require_number(layers, "relativeFluxError")
    interface_jump = require_number(layers, "maxInterfaceTemperatureJump")
    layer_pass = layer_flux_error <= 0.02 and interface_jump <= 1.0e-6
    rows.append(
        f"| A3 Three-layer resistance | {status(layer_pass)} | "
        f"flux error={layer_flux_error:.3e}; interface jump={interface_jump:.3e} K |"
    )
    overall &= layer_pass

    direction = payload.get("directionality", {})
    ratio_error = require_number(direction, "relativeConductivityRatioError")
    ordering_pass = bool(direction.get("orderingPass", False))
    direction_pass = ratio_error <= 0.05 and ordering_pass
    rows.append(
        f"| A4 Physical directionality | {status(direction_pass)} | "
        f"ratio error={ratio_error:.3e}; ordering={ordering_pass} |"
    )
    overall &= direction_pass

    boundary = payload.get("boundarySensitivity", {})
    boundary_error = require_number(boundary, "maximumRelativeError")
    boundary_directions = int(require_number(boundary, "directions"))
    boundary_pass = boundary_error <= 0.05 and boundary_directions >= 3
    rows.append(
        f"| A5 Boundary sensitivity FD | {status(boundary_pass)} | "
        f"max error={boundary_error:.3e}; directions={boundary_directions} |"
    )
    overall &= boundary_pass

    adjoint = payload.get("thermalAdjointFD", {})
    adjoint_error = require_number(adjoint, "maximumRelativeError")
    adjoint_directions = int(require_number(adjoint, "directions"))
    plateau_pass = bool(adjoint.get("fdPlateauPass", False))
    adjoint_pass = (
        adjoint_error <= 0.05 and adjoint_directions >= 4 and plateau_pass
    )
    rows.append(
        f"| A6 Thermal adjoint FD | {status(adjoint_pass)} | "
        f"max error={adjoint_error:.3e}; directions={adjoint_directions}; "
        f"plateau={plateau_pass} |"
    )
    overall &= adjoint_pass

    parallel = payload.get("serialParallelConsistency", {})
    field_error = require_number(parallel, "maxRelativeFieldError")
    objective_error = require_number(parallel, "objectiveRelativeError")
    gradient_error = require_number(parallel, "gradientRelativeError")
    parallel_pass = (
        field_error <= 1.0e-10
        and objective_error <= 1.0e-10
        and gradient_error <= 1.0e-8
    )
    rows.append(
        f"| A7 Serial/parallel consistency | {status(parallel_pass)} | "
        f"field={field_error:.3e}; objective={objective_error:.3e}; "
        f"gradient={gradient_error:.3e} |"
    )
    overall &= parallel_pass

    return overall, rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("STAGE_A_VALIDATION_REPORT.md"),
    )
    args = parser.parse_args()

    payload = json.loads(args.results.read_text(encoding="utf-8"))
    try:
        passed, rows = evaluate(payload)
        error_message = ""
    except (ValueError, TypeError, KeyError) as exc:
        passed = False
        rows = []
        error_message = str(exc)

    metadata = payload.get("metadata", {})
    report_lines = [
        "# Stage A Anisotropic-Conduction Validation Report",
        "",
        f"- Repository commit: `{metadata.get('commit', 'UNKNOWN')}`",
        f"- OpenFOAM version: `{metadata.get('openfoamVersion', 'UNKNOWN')}`",
        f"- Generated by: `{metadata.get('generatedBy', 'UNKNOWN')}`",
        f"- Overall result: **{status(passed)}**",
        "",
    ]
    if error_message:
        report_lines.extend([f"Schema error: `{error_message}`", ""])
    report_lines.extend(
        [
            "| Check | Result | Evidence |",
            "|---|---:|---|",
            *rows,
            "",
            "`anisotropicConductivityValidated` may be set to `true` only "
            "when the overall result is PASS.",
            "",
        ]
    )
    args.report.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"Stage A result: {status(passed)}")
    print(f"Report: {args.report}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
