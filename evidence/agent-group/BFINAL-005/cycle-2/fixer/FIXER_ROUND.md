# BFINAL-005 — Fixer round (cycle-2/fixer): un-gate the pressure-row term

- Task: BFINAL-005 PATCH_RX_PRESSURE_ROW — bounded Fixer round after
  Post-Reviewer `FAIL_IMPLEMENTATION` (single implementation defect, no
  hypothesis/scope change).
- Reviewer-identified defect (verbatim intent): the new term
  `gsenshPressureDropPressureRow = -(J_P^T pc)*dAlphaDxh` was gated behind the
  diagnostic switch `stageB6RxDesignOracle`; with the switch OFF (production
  config) `rxPressureRowT == 0` so production R_x remained momentum-only,
  violating task §6/§22.
- Reviewer-mandated fix: compute `rxPressureRowT` UNCONDITIONALLY (move the
  reduced-SIMPLE rebuild, drAU/dHbyA, g0, T1/T2 face loops, and the
  rxPressureRowT assignment outside the switch guard); keep ONLY the runtime
  dot-test self-check (forward re-assembly + FatalError) and the diagnostic
  Info/exports gated behind `stageB6RxDesignOracle`.

## Provenance

| item | value |
|---|---|
| workspace | `/home/ys/dsH/TO-ANISOTROPIC` == git root (verified) |
| HEAD | `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged) |
| branch | `agent/dsH-stage-b-validation` |
| pre-fix sensitivity.H sha256 | `d777191b…` (executor patch) → `f0ed0726…` (post comment fix) |
| pre-fix helper sha256 | `f8f93796…` → post-fix `fd648ae0…` |
| locked heads | solveDiscreteFlowAdjoint.H `f0c81497…`, solveDiscreteFlowAdjointProduction.H `8c4901bf…` (UNCHANGED) |
| stageB6RxDesignOracle.H | `09832dd6…` (BFINAL-004 locked, UNCHANGED) |
| pre-fix evidence | `cycle-2/fixer/pre-state/` (00–04) |
| build log | `cycle-2/fixer/build/build_fixer.log` |

## What changed (exactly two authorized files)

### src/rxPressureRowTranspose.H (helper, authorized)
Restructure ONLY:
- The reduced-SIMPLE rebuild (a), drAU/dHbyA (b), g0 (c), T1/T2 face loops and
  `rxT` assembly, and the `rxPressureRowT = T` assignment now run
  UNCONDITIONALLY (outside any switch guard).
- The runtime dot-test self-check (deterministic `w`, forward `J_P*w`
  re-assembly, `sum(T*w)` vs `sum(pc*(J_P*w))`, FatalError on relErr > 1e-8)
  and its diagnostic Info remain gated behind
  `optProperties.lookupOrDefault<Switch>("stageB6RxDesignOracle", false)`.
- Header comment updated to document the production semantics.
- Verified: normalized executable-statement multiset is IDENTICAL pre/post
  (only guard placement and comments changed) — the math is untouched.

### src/sensitivity.H (authorized)
- Comment-only update: the helper is no longer described as "a no-op when the
  stageB6RxDesignOracle switch is off"; now documents unconditional
  computation.  No code change (verified by reconstructing the executor's
  version from pre-state diff and diffing — the only difference is the 4-line
  comment).

Forbidden modules: zero diff (only the two locked BFINAL-003 heads appear in
`git diff --stat HEAD`, hashes unchanged).  `git diff --check` clean.

## Build

- Env exactly per approved plan: clean PATH, `source /opt/openfoam7/etc/bashrc`,
  `FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin`, `unset FOAM_SIGFPE`,
  no `set -e`; `cd src && touch MTO_HF.C && wmake`.
- `WMAKE_EXIT=0`; 0 errors; only pre-existing unused-variable warnings.
- New binary: `build/bin/MTO_HF` sha256 `744d8e7433fd1ae0f62d65a0d6667072378946adc42116be92e71d2d85c705ac`
  (pre-fix verified binary was `66a75760…`).

## Verification run A — switch ON (official case, diagnostic enabled)

- CWD `/home/ys/dsH/b2_case_smoke`, `stageB6RxDesignOracle true`; binary by
  absolute path; MTO_RC=0, ExecutionTime 475.42 s, 0 crash markers.
- Log: `run_switchON_fixer.log`.
- Dot test: `relErr=2.57680498181e-15` (identical to cycle-2/post-reviewer).
- Production xh decomposition: `|momentum|L2=0.0373267050442
  |pressureRow|L2=0.0145680567705 |total|L2=0.0516444470747` (identical).
- G1–G5 independently recomputed from raw .mtx (`verify_gates_fixer.py`, pure
  python, math.fsum-free independent code path):
  - G1 D_pressure_production vs oracle: relErr 4.37e-15 / 2.02e-15 / 4.55e-14
    (D1/D2/D3), sign identical — PASS (≪1e-4 preferred).
  - G2 D_total_production vs oracle: relErr 1.23e-15 / 2.42e-15 / 1.26e-13,
    |D_pressure/D_total| = 0.134 / 0.230 / 0.171 — PASS.
  - G3 raw-design projection: relErr 2.42e-9 / 5.38e-10 / 1.04e-9 — PASS.
  - G4 anchors (from raw log): momentum relL2=1.65311363075e-4,
    P-total relL2=1.52828911949e-4, J/J^T dot 6.42244770078e-14,
    ExplicitJToracle maxRelL_U=3.52693970502e-16 / maxRelL_P=3.98062510668e-16,
    GatePR relL2=1.03546813403 cos=0.306822860167 — all byte-identical to
    cycle-2 / BFINAL-003 — PASS.
  - G5: volume proj D1/D2/D3 = 0.1555186511144686 / -0.08449667741658891 /
    0.06396237684432686; anisotropy solidK=(21.6 0 0 21.6 0 16.6); discrete
    objective thermal=2.7025486001e-11 pressure=7.07793481308e-12 —
    unchanged — PASS.
- Fresh artifacts sha256-IDENTICAL to cycle-2 archived artifacts (9/9 key
  files: rxpr_prod_gsens*.mtx, stageB6_dirs/lambda/rxc_*/prod_gsensVol).

## Verification run B — production (switch OFF, reviewer's missing case)

- Case copy `/home/ys/dsH/b2_case_smoke_prod_fix` with
  `stageB6RxDesignOracle false` (only that line changed); same binary; MTO_RC=0,
  ExecutionTime 476 s, 0 crash markers.
- Log: `run_production_switchOFF_fixer.log`.
- `verify_production_activation.py` results:
  - C1/C2: written `1/gsenshPressureDrop` and `1/gsensPressureDrop`
    BYTE-IDENTICAL between switch-ON and switch-OFF runs — the production
    computation path is identical, pressure-row term active.
  - C3: production `dDP` range = `[-0.0277320682043, 0.0259232358427]`
    (patched total) — NOT the momentum-only BFINAL-004 range
    `[-0.0200726738411, 0.0194203165319]`.
  - C4: dot-test line and xh-fields Info ABSENT in production log (gated as
    intended), present in switch-ON log.
  - C5/G3-in-production: projection of the PRODUCTION-written raw-design
    `1/gsensPressureDrop` onto D1/D2/D3 = -0.05623563324297 /
    +0.29767243191247 / -0.03634142980197 — matches the oracle TOTAL
    (relErr 2.4e-9 / 5.4e-10 / 1.0e-9), NOT the momentum-only values
    (-0.0486867554277 / +0.229079662362 / -0.0301314016208).  This is the
    direct proof that production R_x now contains BOTH rows with the switch
    OFF.
  - Gated exports (`rxpr_prod_*.mtx`) NOT regenerated in the production run
    (mtimes 11:15 = pre-run copy), confirming exports stay diagnostic-only.
  - G4 anchors in production log byte-identical to switch-ON (10 matches:
    momentum/P-total relL2, J/J^T dot, ExplicitJToracle, GatePR, discrete
    objective).

## Conclusion

The reviewer's single BLOCKING implementation defect is fixed:
`rxPressureRowT` is now computed unconditionally, so the production
configuration (`stageB6RxDesignOracle=false`) yields a NON-ZERO
pressure/continuity-row term
(`gsenshPressureDropPressureRow = -rxPressureRowT*dAlphaDxh`), and the
production raw-design `dgdx[1]` matches the BFINAL-004-locked
`D_total = D_momentum + D_pressure` oracle to ~1e-9 for D1/D2/D3.
Only the runtime dot-test self-check and diagnostic exports remain gated
behind the diagnostic switch.  G1–G5 pass in the switch-ON run; production
activation is proven in the switch-OFF run; R_w/J and all trusted anchors are
unchanged.  No forbidden module touched; no empirical coefficient; no
D1/D2/D3 fitting; no sign/convention change.
