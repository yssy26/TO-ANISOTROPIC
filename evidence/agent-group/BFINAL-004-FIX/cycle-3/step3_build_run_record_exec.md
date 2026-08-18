# BFINAL-004-FIX cycle-3 — step3 build-and-run record (executor, independent re-run)

- Date: 2026-08-17 · Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root verified,
  `git rev-parse --show-toplevel` == workspace) · HEAD
  `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged)
- Scope of this record: **step3-build-and-run only** — a fresh, independently
  executed header-only incremental build + case run on the writable original
  `/home/ys/dsH/b2_case_smoke`, with full log captured under
  `evidence/agent-group/BFINAL-004-FIX/cycle-3/`. Prior executor text reports
  were treated as context only; every acceptance criterion below was re-verified
  from the raw artifacts of THIS run.

## 1. Probe baseline

- `src/stageB6RxDesignOracle.H` sha256 = `09832dd6ae5a3c7f37201b200bfb19dc0c48f42a3989f9593979c73ca098baec`
  — identical to cycle-3/pre-state (content untouched by the forced recompile
  `touch`; only mtime changed 04:20 → 04:58).
- B1/B2/B3 fix markers present in source: `anRUaW`/`anRPaW` deltaAlpha weighting
  (RX-A), `zeroGradientFvPatchScalarField::typeName` +
  `correctBoundaryConditions()` for `stageB6_d/y` and the `gsensVolR`/`gsensR`
  (+ `gsenshVolR`/`gsenshR`) replicas (RX-C/§9/§10), and the verified B3 chain
  `drAU=-rAU^2/alphaRel`, `dHbyA=(rAU/alphaRel)*((1-alphaRel)*U-HbyA)`.
- Fresh binary fixed-header markers: `J_P*deltaAlpha`(1), `degenerate(nB=0)`(1),
  `projectionEtaDenominator singular`(1); old `anRPa[celli]*da` string absent (0).

## 2. Build (header-only incremental wmake)

Environment exactly per approved plan:

```
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src && wmake
```

- **WMAKE_EXIT = 0** (START 2026-08-16T20:58:41Z → END 2026-08-16T21:21:59Z)
- Exactly **one** compile unit rebuilt: `-c MTO_HF.C` (grep count of
  `g++ .*-c .*\.C` = 1), followed by the relink of all 10 objects into
  `FOAM_USER_APPBIN/MTO_HF`. No other object recompiled (header-only change).
- Binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`, 4747104 B,
  mtime **2026-08-17 05:21:59** (fresh; prior 04:43:55), sha256
  `15172c3aeed061ebb3fa5f2edb5d7954f36730ecb673d04a22ea689b9c20fdab`.
- Warnings: pre-existing unused-variable warnings from production headers only;
  no errors.
- Log: `cycle-3/build_step3_exec.log`.

## 3. Run (absolute-path binary on writable case)

- CWD: `/home/ys/dsH/b2_case_smoke` (writable original; `constant/optProperties`
  line 87 `stageB6RxDesignOracle true;` verified).
- Binary invoked by absolute path `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`
  (same clean-PATH env; FOAM_USER_APPBIN exported after sourcing bashrc;
  `unset FOAM_SIGFPE`).
- **MTO_RC = 0** (EXIT=0); START 2026-08-16T21:22:06Z → END 2026-08-16T21:30:02Z
  (476 s clock).
- Crash markers in whole log: `SIGABRT`/`calculatedFvPatchField`/`FOAM FATAL`/
  `Aborted`/`signal` count = **0**.
- Log tail clean and complete (no truncation): ends with
  `mmaUpdateEnabled=false: primal/adjoint fields evaluated without advancing the
  design.` + `MTO_HF: Finished sensitivity.H` + `MTO_RC=0` + `END=...`.
- Log: `cycle-3/run_stageb6_cycle3_exec.log` (1916 lines).

## 4. Section presence

- `StageB6 RX-A:` (L1750), `StageB6 RX-B:` (L1763), `StageB6 RX-C:` (L1774),
  `StageB6 §9/§10:` (L1874), `--- StageB6 oracle complete` (L1882) — all print.
- Volume anchor: `StageB6 volume-anchor D1: projV="0.1555186511144686"` (==
  required 0.1555186511), D2 -0.08449667741658891, D3 0.06396237684432686.

## 5. Artifacts

- 21 `stageB6_*.mtx` in the case dir, all non-empty (zero-byte count 0),
  including all `rxc_*` (rxc_FD_all_eps 10.67 MB, rxc_alpha_FD_all_eps
  2.25 MB, rxc_analytic 2.22 MB, rxc_xh_FD_all_eps 2.26 MB, rxc_xh_analytic
  452 KB, rxc_z_analytic 452 KB). Fresh mtimes 05:29–05:30.
- Copies refreshed into `cycle-3/artifacts/` (21 files, 0 zero-byte).
- **Determinism:** every compared artifact byte-identical to the prior
  cycle-3 artifacts (`cmp -s` IDENTICAL for rxa/rxb/rxc analytic+FD and
  rpa_terms, lambda).

## 6. Independent recompute (executor-side, pure Python, no shared formula)

`python3 cycle-3/verify_recompute_cycle3.py <case dir>` →
`cycle-3/verify_recompute_cycle3_exec.out` (VERIFY_RC=0). Highlights:

- RX-A rpa_terms identity: max|col1+col2−col3| = 1.29e-26 (machine zero);
  max|col3 − analytic R_P block| = 0.
- RX-A (weighted): R_U,alpha relL2 9.24e-8…9.28e-6 cos=1;
  R_P,alpha relL2 6.80e-7…6.88e-5 cos>0.99999999 (5 eps plateau).
- RX-B: R_U,xh relL2 1.114e-6…1.278e-10 cos=1;
  R_P,xh relL2 5.61e-7…4.50e-10 cos=1 (5 eps plateau; FD R_P,xh ≈ 9.7126e-5).
- RX-C (per D1/D2/D3, eps≤1e-4): xh-tangent relL2 2.1e-5…2.9e-5,
  alpha-tangent relL2 1.95e-5…2.58e-5, R_U,xd relL2 2.3e-5…2.9e-5,
  R_P,xd relL2 4.1e-5…1.1e-4 — all below the 1e-3 diagnostic ceiling,
  cos>0.9999999.
- §9/§10 (lambda-weighted, pressureGradientScale=1):
  D1: |D_pressure/D_momentum| = 0.155; D2 = 0.299; D3 = 0.206;
  prodGsenDPressDrop·d matches −pcᵀR_U,x·d to |Δ|/|prod| ≈ 1e-9 (production
  contains only the momentum piece).

## 7. No unintended re-solve / turbulence frozen

- Between the StageB6 start and `oracle complete` markers, the only DICPCG
  solves are `stageB6_y` (filter tangent, zeroGradient BC), `xp` (Helmholtz
  filter tangent) and the §9/§10 replicas `stageB6_gsensVol` /
  `stageB6_gsensPressureDrop`; **zero** solves of U/p/T/k/omega/nuTilda
  (count 0) — no state variable re-solved.
- Turbulence frozen throughout: `frozenTurbulenceAdjoint=true`,
  `freezeTurbulenceForValidation=true`; StageB6 fingerprint
  nutFrozen|2=0, nuEffFrozen|2=0.00951359211546, k/omega unchanged between
  the two fingerprint prints (L1747 vs L1881 identical).
- Regression anchors (production, byte-identical to BFINAL-003):
  ExplicitJToracle maxRelL2=3.97863790105e-16 / maxRelL_U=3.52693970502e-16 /
  maxRelL_P=3.98062510668e-16; StageB5 momentum relL2=1.65312e-4 /
  cos=0.999999986345; P-total 1.52829e-4; StageB4 dPhi_J L2=2.19206493915e-05;
  RxProbe Rx-A relL2=1.61e-14 cos=1.

## 8. Scope integrity

- `git status --porcelain=v1` identical to cycle-3/pre-state: only the two
  pre-existing BFINAL-003 production heads `M`; hashes unchanged
  (solveDiscreteFlowAdjoint.H=f0c81497…, solveDiscreteFlowAdjointProduction.H=
  8c4901bf…). Probe header content unchanged (sha256 identical). No production
  file, optProperties, or case file modified by this step.

## 9. Acceptance summary (step3)

| criterion | result |
|---|---|
| wmake exit 0 + fresh binary mtime | PASS (exit 0; binary 04:43:55 → 05:21:59) |
| only MTO_HF.C recompiled/relinked | PASS (1 compile line + relink) |
| MTO_RC == 0 | PASS (0) |
| no SIGABRT / calculatedFvPatchField | PASS (0 occurrences) |
| RX-A/B/C + §9/§10 sections print | PASS |
| all stageB6_*.mtx non-empty incl. rxc_* | PASS (21 files, 0 zero-byte) |
| log saved under cycle-3/ | PASS (build_step3_exec.log, run_stageb6_cycle3_exec.log) |
