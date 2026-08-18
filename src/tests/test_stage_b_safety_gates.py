#!/usr/bin/env python3
"""Static regression checks for Stage-B pressure and MMA safety gates."""

from pathlib import Path
import re
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]


class StageBSafetyGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.production = (
            REPO_ROOT / "src" / "solveDiscreteFlowAdjointProduction.H"
        ).read_text(encoding="utf-8")
        cls.diagnostic = (
            REPO_ROOT / "src" / "solveDiscreteFlowAdjoint.H"
        ).read_text(encoding="utf-8")
        cls.mma_gate = (
            REPO_ROOT / "src" / "validateMmaUnlockGate.H"
        ).read_text(encoding="utf-8")

    def test_pressure_reference_follows_openfoam_need_reference(self):
        self.assertIn(
            "const bool prodPressureNeedsReference = p.needReference();",
            self.production,
        )
        self.assertIn(
            "const bool discretePressureNeedsReference = p.needReference();",
            self.diagnostic,
        )

    def test_no_legacy_unconditional_reference_zeroing_remains(self):
        forbidden = (
            "prodRhs[prodPIndex(pRefCell)] = 0.0;",
            "prodSolution[prodPIndex(pRefCell)] = 0.0;",
            "prodResidual[prodPIndex(pRefCell)] = 0.0;",
            "discreteRhs[discretePIndex(pRefCell)] = 0.0;",
            "discreteSolution[discretePIndex(pRefCell)] = 0.0;",
            "physRhs[discretePIndex(pRefCell)] = 0.0;",
            "physSolution[discretePIndex(pRefCell)] = 0.0;",
        )
        combined = self.production + self.diagnostic
        for statement in forbidden:
            self.assertNotIn(statement, combined)

    def test_explicit_reference_row_is_conditioned(self):
        self.assertRegex(
            self.diagnostic,
            re.compile(
                r"if \(discretePressureNeedsReference\)\s*\{\s*"
                r"const label pRefRow = discretePIndex\(pRefCell\);"
            ),
        )

    def test_mma_has_no_direction_only_bypass(self):
        self.assertNotIn("guardedDirectionOnly", self.mma_gate)
        self.assertNotIn("The only exception", self.mma_gate)
        self.assertIn(
            "if (mmaUpdateEnabled && !frozenGradientValidated)",
            self.mma_gate,
        )
        self.assertIn(
            "Direction-only approval cannot bypass this gate.",
            self.mma_gate,
        )


if __name__ == "__main__":
    unittest.main()
