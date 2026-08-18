# BFINAL-005 cycle-2 — step3 build-and-run record (executor)

- Task: BFINAL-005 PATCH_RX_PRESSURE_ROW (mode=PATCH, maxReplans=0) — step3_build_run
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` == git root (verified by
  `git rev-parse --show-toplevel`); HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280`
  (unchanged); branch `agent/dsH-stage-b-validation`.
- Prior steps verified from real workspace/evidence (not trusted blindly):
  step1 pre-state present (`cycle-2/pre-state/` 00-06), step2 patch present
  (`src/rxPressureRowTranspose.H` + `src/sensitivity.H` diff 102+/2-),
  cycle-1 marked ABORTED_PREAUTHORIZATION, locked BFINAL-003 heads unchanged.

## 1. Entry state fingerprint

- See `00_entry_fingerprint.txt`. Key facts before this step:
  - `src/sensitivity.H` sha256 = `d777191b472360e4710c6dc78cf3e3c87e3737daf02ef2c040b3e81a7e457efd` (post-step2)
  - `src/rxPressureRowTranspose.H` sha256 = `f8f9379660a0d155284b8a6ba15003befd37873410bb29d4dc10a3fd4aee9a6d`
  - Locked heads: solveDiscreteFlowAdjoint.H `f0c81497…`, solveDiscreteFlowAdjointProduction.H `8c4901bf…`
  - Forbidden-module diff vs HEAD: ONLY the two locked BFINAL-003 heads
    (260 + 69 lines); all other forbidden modules zero diff.
  - Binary `build/bin/MTO_HF` mtime 10:16:15 (from step2 check build).
  - No `rxpr_prod_*.mtx` in case dir before run.

## 2. Build (official, header-only incremental wmake)

Environment exactly per approved plan (clean PATH, source bashrc, unset FOAM_SIGFPE, no set -e):

```
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src && touch MTO_HF.C && wmake
```

- **WMAKE_EXIT = 0** (BUILD_START 2026-08-17T02:17:55Z → BUILD_END 02:43:48Z, ~26 min).
- Exactly **one** compile unit rebuilt: `-c MTO_HF.C` (grep count of `-c .*\.C` = 1),
  followed by the relink into `FOAM_USER_APPBIN/MTO_HF` (grep count = 1). No other object rebuilt.
- `error:` count in build log = **0** (only pre-existing unused-variable warnings
  from NS.H / createFields.H / solveDiscreteFlowAdjoint.H).
- Fresh binary: `build/bin/MTO_HF`, 4816856 B, mtime **2026-08-17 10:43:48**
  (was 10:16:15), sha256 `66a75760ab9b19cd12f536c14ed95b492ac131f98573c70e91816a7c4c90fd8d`.
- Log: `build_step3.log`.

## 3. Run (absolute-path binary on writable case)

- CWD `/home/ys/dsH/b2_case_smoke`; `constant/optProperties` verified at run start:
  line 79 `stageB4JacobianProbe true;`, line 87 `stageB6RxDesignOracle true;`
  (also `mmaUpdateEnabled false`, `frozenGradientValidated false`).
- Binary by absolute path `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`; same clean env.
- **MTO_RC = 0**; RUN_START 2026-08-17T02:43:56Z → RUN_END 02:52:03Z (~8 min);
  ExecutionTime = 485.76 s.
- Crash markers in whole log: `FOAM FATAL` / `SIGABRT` / `Aborted` / `signal` /
  `calculatedFvPatchField` count = **0**.
- Log tail complete: `MTO_HF: Finished sensitivity.H` + `MTO_RC=0`.
- Log: `run_stageb6_cycle2.log` (1916 lines).

## 4. Key log markers (this patch)

- **Dot test (transpose pinned before gates):** L1895
  `RxPressureRowTranspose: sum(T*w)="-2.696485451749351e-08" sum(pc*(J_P*w))="-2.696485451749358e-08" relErr=2.57680498181e-15 alphaRel=0.4 |T|L2=0.000256656999292 |J_P*w|L2=6.22344977639e-10 |g0_bnd|max=0.0650620425279 w-support=5040`
  — relErr 2.58e-15 << 1e-8 threshold (machine precision): `rxPressureRowT == J_P^T pc` is exact.
- **Production xh decomposition:** L1896
  `RxPressureRow production xh fields: |momentum|L2=0.0373267050442 |pressureRow|L2=0.0145680567705 |total|L2=0.0516444470747`
  — both components nonzero; total = momentum + pressureRow in the exported fields.
- **New decomposition exports:** `rxpr_prod_gsensh_momentum.mtx`,
  `rxpr_prod_gsensh_pressurerow.mtx`, `rxpr_prod_gsensh_total.mtx` (xh-level pre-chain),
  `rxpr_prod_gsens.mtx` (raw-design post-chain). All non-empty, fresh mtime 10:52.

## 5. Oracle anchors (BFINAL-004 §9/§10) — reproduced EXACTLY

| dir | D_momentum | D_pressure | D_total | prod·d (momentum replica) |
|---|---|---|---|---|
| D1 | -0.0486867554277 | -0.00754887795106 | -0.0562356333787 | -0.0486867553288 (rel 2.03e-9) |
| D2 | +0.229079662362 | +0.0685927693909 | +0.297672431753 | +0.229079662485 (rel 5.40e-10) |
| D3 | -0.0301314016208 | -0.00621002821874 | -0.0363414298395 | -0.0301314015977 (rel 7.65e-10) |

All identical to the BFINAL-004 locked anchors (eps-independent analytic chain).
`prodGsenDPressDrop*d` == D_momentum to ~1e-9 confirms the stageB6 production replica is
still momentum-only; the NEW production pressure-row term lives in the `rxpr_prod_*.mtx`
decomposition (G1/G2 evaluation is the next step).

## 6. R_w/J no-regression anchors (G4) — all MATCH BFINAL-003

- momentum relL2 = **1.65312e-4** (L1290, anchor 1.653e-4)
- P-total relL2 = **1.52829e-4** (L1294, anchor 1.53e-4)
- J/J^T full dot: `Reduced cold-flow operator transpose dot-test (pressureDrop) max relative error = 6.42244770078e-14` (L1394; thermalCoupling same value L1004) — anchor **6.42e-14**
- Explicit oracle: `ExplicitJToracle: maxRelL2=3.97863790105e-16 maxCosErr=2.22044604925e-16 maxRelL_U=3.52693970502e-16 maxRelL_P=3.98062510668e-16` (L1266) — anchors 3.53e-16 / 3.98e-16
- GatePR (h=1e-3/3e-4/1e-4/3e-5) prints, values identical across both occurrences (L1233-1249, L1623-1639)
- BlockDot (thermalCoupling): PU=2.078e-13 PP=2.666e-14 UU=7.310e-14 UP=1.796e-12 boundaryU=1.804e-13 (L1201)
- RxProbe Rx-A relL2=1.6135e-14 cos=1 (L147)

## 7. Non-flow anchors (G5)

- StageB6 volume anchors: D1 `0.1555186511144686`, D2 `-0.08449667741658891`, D3 `0.06396237684432686` — identical to BFINAL-004-FIX cycle-3.
- StageB6 fingerprint (L1888): |U|2=16796.0670769 |p|2=7965118.75546 |phi|2=0.0041733366884 |alpha|2=14587858210.1 |k|2=11566.7580246 |omega|2=13527763.4515 |nutFrozen|2=0 |nuEffFrozen|2=0.00951359211546 — matches cycle-3 fingerprint; turbulence frozen.

## 8. Determinism of oracle artifacts

All compared analytic artifacts **byte-identical** to BFINAL-004-FIX cycle-3/artifacts:
`stageB6_rpa_terms.mtx`, `rxa_analytic`, `rxb_analytic`, `rxc_analytic`,
`rxc_xh_analytic`, `rxc_z_analytic`, `stageB6_lambda.mtx` — `cmp -s` IDENTICAL (7/7).

## 9. Artifacts collected

- `artifacts/`: 45 files (2 explicitJToracle/rhs, 11 stageB2_*, 8 stageB5_*, 21 stageB6_*,
  4 rxpr_prod_*) — zero-byte count **0**; manifest `artifacts_sha256.txt` (45 entries).
- Freshly regenerated by this run: all stageB6_* (10:51-10:52), rxpr_prod_* (10:52),
  explicitJT/explicitRhs/stageB2/stageB5 (10:51).

## 10. Scope integrity

- `git status --porcelain=v1` at exit identical to entry: only `src/sensitivity.H`
  (the authorized patch) + the two pre-existing LOCKED BFINAL-003 heads `M`;
  `src/rxPressureRowTranspose.H` untracked (authorized new helper); no forbidden module diff.
- sensitivity.H sha256 unchanged during step3 (d777191b…); helper sha256 f8f93796….
- State fingerprint entry==exit: HEAD/branch/porcelain/sha256/forbidden-diff-stat IDENTICAL
  (modulo binary mtime 10:16:15 → 10:43:48 from the rebuild, expected).
- No case or optProperties file modified by this step.

## 11. Acceptance summary (step3)

| criterion | result |
|---|---|
| wmake exit 0 + fresh binary mtime | PASS (exit 0; binary 10:16:15 → 10:43:48, sha 66a75760…) |
| exactly the MTO_HF unit + relink | PASS (1 compile line + 1 relink, 0 errors) |
| MTO_RC == 0 | PASS (0) |
| no FOAM FATAL / SIGABRT / signal | PASS (0 occurrences) |
| oracle artifacts regenerated | PASS (21 stageB6 + stageB2/stageB5/explicit regenerated, analytic byte-identical to cycle-3) |
| new decomposition .mtx non-empty | PASS (4 rxpr_prod_*.mtx, 0 zero-byte) |
| state fingerprint unchanged entry==exit | PASS (IDENTICAL modulo expected binary mtime) |
| dot test pins transpose | PASS (relErr 2.58e-15 ≤ 1e-8) |
| R_w/J anchors | PASS (momentum 1.653e-4, P-total 1.53e-4, J/J^T 6.42e-14, explicit 3.53e-16/3.98e-16) |
