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

    def test_production_preconditioner_default_is_linear_setup(self):
        self.assertIn(
            '"discreteProdPreconditionerSetup",',
            self.production,
        )
        self.assertRegex(
            self.production,
            re.compile(
                r'"discreteProdPreconditionerSetup",\s*"pressureGAMG"'
            ),
        )
        self.assertIn(
            'prodPreconditionerSetup == "pressureGAMG"',
            self.production,
        )
        self.assertIn(
            "-fvm::laplacian(prodPressurePrecMobility, prodPressurePrecPsi)",
            self.production,
        )
        self.assertIn(
            "prodPressurePrecMatrix.solve(prodPressurePrecSolverDict);",
            self.production,
        )
        self.assertIn(
            'else if (prodPreconditionerSetup == "bruteForceL1")',
            self.production,
        )
        self.assertIn(
            '"discreteProdEnableRitzPilot",',
            self.production,
        )
        self.assertRegex(
            self.production,
            re.compile(
                r"prodEnableRitzPilot\s*&&\s*\(prodIteration"
            ),
        )

    def test_pressure_gamg_equivalence_decomposition_gate(self):
        # BFINAL-010: the pressure-GAMG preconditioner must be gated by a
        # four-group matrix/operator equivalence decomposition (hard
        # thresholds, not configurable) plus a laplacianSchemes guard.
        for metric in (
            "interiorDiagRelL2",
            "boundaryDiagRelL2",
            "effectiveDiagRelL2",
            "offDiagRelL2",
        ):
            self.assertIn(metric, self.production)
        self.assertIn("PRODPRECGAMGCHECK", self.production)
        self.assertIn(
            "pressureGAMG preconditioner matrix is NOT equivalent to ",
            self.production,
        )
        self.assertIn(
            "PRODPRECGAMGSCHEME",
            self.production,
        )
        self.assertIn(
            'subDict("laplacianSchemes")',
            self.production,
        )
        self.assertIn(
            "prodPressurePrecMatrix.internalCoeffs()[patchi]",
            self.production,
        )
        self.assertIn(
            "addBoundaryDiag",
            self.production,
        )
        self.assertIn(
            "effectiveDiagRelL2 > 1e-8",
            self.production,
        )
        self.assertIn(
            "offDiagRelL2 > 1e-8",
            self.production,
        )
        # The old incompatible bare comparison (boundary-less fvMatrix
        # diagonal vs boundary-inclusive operator diagonal) must not
        # regress.
        self.assertNotIn(
            "prodPressurePrecMatrix.diag()[celli] - refDiag",
            self.production,
        )
        self.assertNotIn(
            "): diagRelL2=",
            self.production,
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
