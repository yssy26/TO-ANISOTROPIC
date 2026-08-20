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
        cls.main_solver = (
            REPO_ROOT / "src" / "MTO_HF.C"
        ).read_text(encoding="utf-8")
        cls.create_fields = (
            REPO_ROOT / "src" / "createFields.H"
        ).read_text(encoding="utf-8")
        cls.production_probe = (
            REPO_ROOT / "src" / "stageB13GradientProbe.H"
        ).read_text(encoding="utf-8")
        cls.contraction_probe = (
            REPO_ROOT / "src" / "stageB13ContractionProbe.H"
        ).read_text(encoding="utf-8")
        cls.b2_amplitude = (
            REPO_ROOT / "src" / "validateStageB2GradientAmplitude.H"
        ).read_text(encoding="utf-8")
        cls.adj_ht = (
            REPO_ROOT / "src" / "AdjNS_HT.H"
        ).read_text(encoding="utf-8")
        cls.sensitivity = (
            REPO_ROOT / "src" / "sensitivity.H"
        ).read_text(encoding="utf-8")
        cls.rx_helper = (
            REPO_ROOT / "src" / "rxPressureRowTranspose.H"
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

    def test_bfinal011_label_isolation_switches_and_inner_solver_controls(self):
        # BFINAL-011 cycle-1: per-label adjoint isolation switches (default
        # true, baseline unchanged) and inner-solver options for the
        # pressure-Poisson preconditioner (GAMG default; PCG+DIC fallback;
        # GAMG robustness coefficients overridable; fast-fail on non-finite
        # inner psi).
        for switch in (
            "solveThermalCouplingFlowAdjoint",
            "solvePressureDropFlowAdjoint",
        ):
            self.assertIn(switch, self.main_solver)
            self.assertRegex(
                self.create_fields,
                re.compile(rf'"{switch}",\s*true'),
            )
        self.assertIn(
            '"discreteProdPressurePrecSolver",',
            self.production,
        )
        self.assertRegex(
            self.production,
            re.compile(r'"discreteProdPressurePrecSolver",\s*"GAMG"'),
        )
        self.assertIn(
            'prodPressurePrecSolverDict.add("solver", "PCG");',
            self.production,
        )
        self.assertIn(
            'prodPressurePrecSolverDict.add("preconditioner", "DIC");',
            self.production,
        )
        self.assertIn(
            '"discreteProdPressurePrecScaleCorrection",',
            self.production,
        )
        self.assertIn(
            '"discreteProdPressurePrecNPreSweeps",',
            self.production,
        )
        self.assertIn(
            '"discreteProdPressurePrecDirectSolveCoarsest",',
            self.production,
        )
        self.assertIn(
            "prodPressurePrecApplyCount",
            self.production,
        )
        self.assertIn(
            "non-finite psi",
            self.production,
        )
        # The equivalence gate must read the matrix through a const
        # lduMatrix& alias: the non-const lower() materialises the
        # symmetric-storage lower array (side effect that rejects PCG and
        # flips the GAMG scaleCorrection default).
        self.assertIn(
            "const lduMatrix& prodPrecLduMatrix = prodPressurePrecMatrix;",
            self.production,
        )
        self.assertNotIn(
            "const scalarField& prodPrecLower = prodPressurePrecMatrix.lower();",
            self.production,
        )

    def test_bfinal013_localization_probes_are_gated_and_readonly(self):
        # BFINAL-013: P1 source dot tests and P2 contraction decomposition
        # are switch-gated diagnostics (default off) and must not alter the
        # formal acceptance machinery.
        self.assertIn(
            "stageB13GradientProbe.H",
            self.main_solver,
        )
        self.assertIn(
            "stageB13ContractionProbe.H",
            self.b2_amplitude,
        )
        self.assertIn(
            "B13SOURCEPROBE P1a gDP/p-source",
            self.production_probe,
        )
        self.assertIn(
            "B13SOURCEPROBE P1b J/T-source",
            self.production_probe,
        )
        self.assertIn(
            "B13SOURCEPROBE P1c J/phiOut-source",
            self.production_probe,
        )
        self.assertIn(
            '"stageB13GradientProbe"',
            self.b2_amplitude,
        )
        self.assertIn(
            "B13CONTRACT P2 decomposition",
            self.contraction_probe,
        )

    def test_bfinal015_state_export_is_gated_and_measurement_only(self):
        # BFINAL-015 T2/T3: the state export is switch-gated (default off)
        # and only writes files inside the B2 output directory (pure
        # measurement round: no operator/gradient/objective change).
        self.assertIn(
            '"stageB15StateExport"',
            self.b2_amplitude,
        )
        self.assertIn(
            "wstate_baseline_U.mtx",
            self.b2_amplitude,
        )
        self.assertIn(
            "wstate_",
            self.b2_amplitude,
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

    def test_bfinal018_tolerance_and_functional_stability_gate(self):
        # BFINAL-018 fix 2: the production/diagnostic adjoint tolerance
        # default is tightened to 1e-12 (a 1e-9 true residual leaves an
        # O(1) error in the ~1e8-amplified gradient functionals), and the
        # FGMRES loop is guarded by a cross-restart gradProxy stability
        # check (default ON, threshold 1e-6) so a lambda_679-class
        # functional false convergence cannot pass the residual criterion
        # alone.
        for text in (self.production, self.diagnostic):
            self.assertIn(
                '"discreteFlowAdjointTolerance", 1e-12',
                text,
            )
        self.assertRegex(
            self.production,
            re.compile(
                r'"discreteProdGradientStabilityCheck",\s*\n\s*true'
            ),
        )
        self.assertRegex(
            self.production,
            re.compile(
                r'"discreteProdGradientStabilityTolerance",\s*\n\s*1e-6'
            ),
        )
        self.assertIn("GRADSTABLE", self.production)
        self.assertIn("prodGradientStable", self.production)
        self.assertIn("B18RHSEXPORT", self.production)

    def test_bfinal018_source_folding_matches_operator_flux_map(self):
        # BFINAL-018 fix 3: the b_TC face-functional routing must be the
        # transpose of the SAME relaxed-SIMPLE flux tangent the J operator's
        # P rows implement (alphaRel*rAU hA path + (1-alphaRel) direct part
        # + kf with the operator's own J_PP sign).  The retired helpers
        # (plain-interpolation routing with opposite-sign interior kf) must
        # not reappear.
        self.assertIn("alphaRel*prodRAU[own]*weight*faceTranspose",
                      self.production)
        self.assertIn(
            "prodExternalTranspose[prodPIndex(own)] += kf*faceFunctional",
            self.production,
        )
        self.assertIn("alphaRel*rAUAdj[own]*weight*faceTranspose",
                      self.diagnostic)
        self.assertIn(
            "discreteExternalFluxTranspose[discretePIndex(own)]"
            " += kf*faceFunctional",
            self.diagnostic,
        )
        self.assertNotIn(
            "applyProdPressureFluxCorrectionTranspose",
            self.production,
        )
        self.assertNotIn(
            "applyPressureFluxCorrectionTranspose",
            self.diagnostic,
        )
        # the thermal face functional itself (B0.2-verified) is unchanged
        self.assertIn("adjointDownwind*temperatureJump", self.adj_ht)

    def test_bfinal019_j_assembly_pressure_row_and_flux_direct(self):
        # BFINAL-019: the thermal-objective (J) xh assembly must carry the
        # two DERIVATION_FIX3 terms that no (U,p) source can absorb:
        #   (1) -pb^T R_P,x  -> gsenshMeanTPressureRow = -rxPressureRowTb
        #       *dAlphaDxh   (second exact-transpose contraction of the SAME
        #                       BFINAL-004-validated J_P operator with the
        #                       thermal adjoint pressure pb)
        #   (2) +g^T phi_x   -> gsenshMeanTFluxDirect = +rxFluxDirectT
        #       *dAlphaDxh   (Gx direct term; phi_x mirrors the LOCKED
        #                       dHbyA/drAU/g0 channels of the operator)
        # The helper must stay token-compatible with the historical pc path
        # when no parameterization macros are defined (default include first,
        # macro-parameterized pb include second, #undef cleanup after).
        self.assertEqual(
            self.sensitivity.count('#include "rxPressureRowTranspose.H"'), 2
        )
        self.assertEqual(
            self.sensitivity.count("#define RX_ADJ_PRESSURE"), 1
        )
        # the FIRST include must be the macro-free default (pc) path: the
        # macro block appears only for the second, pb-parameterized include
        first = self.sensitivity.index('#include "rxPressureRowTranspose.H"')
        second = self.sensitivity.index('#include "rxPressureRowTranspose.H"',
                                        first + 1)
        macro_block = self.sensitivity.index("#define RX_ADJ_PRESSURE")
        self.assertLess(first, macro_block)
        self.assertGreater(macro_block, first)
        self.assertGreater(second, macro_block)
        self.assertIn("#define RX_ADJ_PRESSURE pb", self.sensitivity)
        self.assertIn("#define RX_OUTPUT_FIELD rxPressureRowTb",
                      self.sensitivity)
        self.assertIn(
            "#define RX_FLUX_DIRECT_SOURCE thermalCouplingFaceFunctional",
            self.sensitivity,
        )
        self.assertIn("#undef RX_ADJ_PRESSURE", self.sensitivity)
        # the two production terms and their wiring into fsenshMeanT
        self.assertIn("gsenshMeanTPressureRow = -rxPressureRowTb*dAlphaDxh;",
                      self.sensitivity)
        self.assertIn("gsenshMeanTFluxDirect = rxFluxDirectT*dAlphaDxh;",
                      self.sensitivity)
        self.assertIn("fsenshMeanT += gsenshMeanTPressureRow;",
                      self.sensitivity)
        self.assertIn("fsenshMeanT += gsenshMeanTFluxDirect;",
                      self.sensitivity)
        # the face functional consumed by Gx is the value-replica of the
        # AdjNS_HT.H b_TC folding (masked upwind-downwind interior + cold
        # outlet dJ/dphi_out)
        self.assertIn("adjointDownwind*temperatureJump", self.sensitivity)
        self.assertIn("thermalObjectiveFluxDerivative", self.sensitivity)
        # frozen-validation mode still zeroes the new decomposition fields
        self.assertIn("gsenshMeanTPressureRow = dimensionedScalar",
                      self.sensitivity)
        self.assertIn("gsenshMeanTFluxDirect = dimensionedScalar",
                      self.sensitivity)
        # helper contract: default macros + #undef cleanup + Gx basis
        self.assertIn("#ifndef RX_ADJ_PRESSURE", self.rx_helper)
        self.assertIn("#define RX_ADJ_PRESSURE pc", self.rx_helper)
        self.assertIn("#undef RX_FLUX_DIRECT_SOURCE", self.rx_helper)
        self.assertIn("rxFluxDirectT", self.rx_helper)
        self.assertIn("RX_FLUX_DIRECT_SOURCE[fi]", self.rx_helper)
        # the historical pc production term must remain untouched
        self.assertIn("gsenshPressureDropPressureRow = -rxPressureRowT"
                      "*dAlphaDxh;", self.sensitivity)


if __name__ == "__main__":
    unittest.main()
