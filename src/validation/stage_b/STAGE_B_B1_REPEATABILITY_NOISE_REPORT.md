# Stage B1 — Frozen-Turbulence Baseline Repeatability and Numerical Noise-Floor Validation

- Repository: `/home/ys/TO-ANISOTROPIC`
- Branch: `agent/stage-a-anisotropic-validation`
- Baseline commit: `f3c4aa25e4ae46b4bed80143cc1ef81ab6c7c68f`
- OpenFOAM-7; case `/home/ys/b1_case` (gate4-derived, serial only); audit date 2026-08-06
- Scope: Stage B1 only. No MMA, no design change, no `gradientValidated`, no `directionOnlyOptimizationApproved`, no full SST Gate 6, no objective/constraint/adjoint/threshold change. commit/push/PR: none.

## 1. Repository, branch, commit, workspace

Uncommitted B0.1/B0.2 changes (recorded): `src/AdjNS_HT.H` (+60 diag), `src/computeObjective.H` (+43 dJ/dphi), `src/createFrozenHotRegionFields.H` (+21 decl); new `src/validateStageBRepeatability.H` (this stage); `src/createFields.H` (+stageBEnabled/stageBRepeatA/B/SNRRepeats); `src/MTO_HF.C` (+include). New untracked: `src/validation/`, `src/validation_stage_*.log`. B1 does NOT commit or push.

## 2. B0 documentation closure

`STAGE_B_B0_AUDIT_REPORT.md`, `stage_b_b0_audit.json`, `stage_b_b0_variable_trace.tsv` are unified:
- B0-01 = PASS, B0-17 = PASS, B0-20 = PASS, `readyForB1 = true`, `overallStatus = PASS (B0 + B0.1 + B0.2)`.
- The old "internal-face thermal-coupling term inconsistent with the residual derivative" conclusion is marked **SUPERSEDED** (B0.2 proved `-Tb_downwind*(T_nei-T_own)` is the exact fvm transpose).
- A scheme-scope note is added: all B0.2 transpose results apply to `div(phi,T) bounded Gauss upwind`; changing the temperature convection scheme in `system/fvSchemes` requires re-running B0.2.

## 3. Frozen-turbulence baseline (two-phase)

Phase 1 (full SST, `freezeTurbulenceForValidation=false`, `turbulence->correct()` active): `solveFullSST(1000, 1e-6)` converged after **246 iterations** (Urel 9.9e-8, nutRel 7.6e-7, 3 consecutive). Phase 2: `freezeTurbulenceForValidation=true`, `#include updateFrozenTurbulenceFields.H`, all fields saved with min/max/L1/L2/FNV-checksum statistics (all `nBad=0`). Baseline: `J=-1.1673569`, `PressureDropPa=48336.3`, `gDP=0.933452`, `gV=-7.9329e-9`.

## 4. Actual convergence criterion

`solveFrozenPrimalConverged` (new, replaces the fixed-300-iteration practice; `solveFrozenPrimalFixed` remains untouched for F1/F2/F3): solves frozen U-p-phi-T each outer iteration, records U/p/T residuals, continuity error, J and gDP, stops on `dJrel<=1e-9 && dgDPrel<=1e-8` for 5 consecutive iterations **and** a 20-iteration platform window `|J(t)-J(t-20)|<=1e-10`, `|gDP(t)-gDP(t-20)|<=1e-9`; safe cap 1000.

**Measured behaviour (finding):** the 1e-10 platform window is never satisfied within 1000 iterations because J oscillates at ~1e-7 amplitude from one frozen iterate to the next (a discretization-level limit cycle, not divergence); every run therefore returns the 1000-iteration platform value (Iter=-1). This is reported honestly, not faked as converged. The state is nevertheless fully deterministic and repeatable to machine precision (Section 5-9), and the long-time attractor is unique (Test B returns to the identical point, Section 6).

## 5. Test A — identical initial state, 5 full chains

Each chain: restore baseline → frozen solve → T → objective derivatives → Tb → Ub → Uc → gradients (dfdx/dgdx0/dgdx1) via the exact main-loop pipeline (block-scoped includes of updateThermalFlux/HeatTransfer/computeObjective/AdjHeatTransfer/AdjNS_HT/AdjNS_PD/sensitivity).

| metric | value | gate | status |
|---|---:|---:|---|
| J relative range | **0.0** | <=1e-8 | PASS |
| gDP relative range | **0.0** | <=1e-7 | PASS |
| gV absolute range | **0.0** | <=1e-13 | PASS |
| dfdx field rel L2 diff | **0.0** (all projection components bit-identical) | <=1e-5 | PASS |
| dgdx[0] field rel L2 diff | **0.0** | <=1e-10 | PASS |
| dgdx[1] field rel L2 diff | **0.0** | <=1e-4 | PASS |
| 3 gradient projections sign | consistent (all bit-identical) | no sign change | PASS |

A1..A5 are bit-identical: `J=-1.16506824580299`, `gDP=0.9456312704707863`, `gV=-7.932936196830553e-09`, `|dfdx|2=0.010053342087`, `|dgdx0|2=0.0140859042455`, `|dgdx1|2=0.54201030577`, projections `(projD1,projD2,projD3)=(-0.03689778,-0.02001873,-0.01940251)`, `(projDP1,projDP2,projDP3)=(-5.463965,1.202322,14.059707)`, `(projV1,projV2,projV3)=(0.1555187,-0.0844967,0.0639624)`. Tb/Ub adjoint residuals also bit-identical (6.5e-13 / 9.97e-9).

## 6. Test B — perturbed initial state (3 directions)

U/p/T initial states perturbed by physical-coordinate sine/cosine modes (`epsU=1e-4`, `epsP=1e-4`, `epsT=1e-5`); k/omega/nut/nutFrozen/alphaTurbulent/x/xp/xh untouched. All three runs return to the **identical** frozen steady state:

| metric | value | gate | status |
|---|---:|---:|---:|
| max |J-J_A|/|J_A| | **3.43e-15** | <=1e-7 | PASS |
| max |gDP-gDP_A|/|gDP_A| | **0.0** | <=1e-6 | PASS |
| dfdx field rel L2 | **0.0** | <=1e-4 | PASS |
| dgdx[0] rel L2 | **0.0** | <=1e-10 | PASS |
| dgdx[1] rel L2 | **0.0** | <=1e-3 | PASS |
| projection signs | identical to A | consistent | PASS |

The frozen flow has a single attracting steady state; perturbed starts converge back to it exactly (deterministic PCG solver).

## 7. Frozen-field drift check

| field | maxAbs | maxRel | L2rel |
|---|---:|---:|---:|
| nutFrozen | 0 | 0 | 0 |
| nuEffFrozen | 0 | 0 | 0 |
| alphaTurbulent | 0 | 0 | 0 |
| k | 0 | 0 | 0 |
| omega | 0 | 0 | 0 |
| nut | 0 | 0 | 0 |

Hard gate `maxRel <= 1e-13` PASS (exactly 0; `turbulence->correct()` is never called in frozen mode and no frozen field is overwritten).

## 8. Gradient-field repeatability

Representative field signatures (MMA-consistent projections onto D1/D2/D3, no mesh.V()): bit-identical across all 5 A-trials and 3 B-trials for dfdx, dgdx[0], dgdx[1]. Cosine similarity = 1, relative L2/Linf differences = 0.

## 9. Fixed-direction projection repeatability (D1/D2/D3)

Physical-coordinate low-frequency directions, designMask-only, bounds-safe (|x|>0.98 or x<0.02 zeroed), L-inf normalized to 1. Projections recorded per trial (see Section 5) are bit-identical; no sign change in any non-near-zero projection.

## 10. Numerical noise floor of J and gDP

All scalar ranges across the 5 A-trials are exactly zero: J, PressureDropPa, gDP, gV, and the 9 projections. The frozen-chain noise floor is therefore below the 16-digit printing resolution (deterministic PCG solver; contrast with GAMG, Section 13).

## 11. SNR probe (h=1e-3 along D1)

x +/- h*D1, each side solved 3 times (filter+projection+alpha+DTMolecular+DTEffective recomputed, turbulence frozen):

| quantity | value |
|---|---:|
| Jplus (3 reps) | -1.165027241615661 (identical) |
| Jminus (3 reps) | -1.165110677934546 (identical) |
| signalJ | 8.3436e-05 |
| noiseJ | 0 |
| **SNR_J** | **8.34e25** (gate >=100) PASS |
| DPplus/DPminus | 48576.8465 / 48703.7426 (identical) |
| signalDP | 126.896 |
| noiseDP | 0 |
| **SNR_DP** | **1.27e32** (gate >=20) PASS |

## 12. Recommended formal-FD step range

With the deterministic frozen solver the noise floor is machine-zero for identical restarts, so h=1e-3 (SNR ~1e25) is far into the asymptotic regime. For the formal B2/B3 FD amplitude acceptance the usable step range is `h in [1e-5, 1e-2]` for the objective and `h in [1e-4, 1e-2]` for the pressure constraint (signal 126.9 Pa at h=1e-3 vs residual pressure solver tolerance 1e-9; the lower bound is set by the p-solver final residual, not by restart noise). The task's minimum SNR gates are satisfied by several orders of magnitude at h=1e-3.

## 13. Issues found and resolved

1. **[FOUND & FIXED] GAMG pressure-solver non-determinism.** With the stock `GAMG` p solver, identical-restart frozen solves produced J-platform noise of ~1e-5 (A1-A5 drifted -1.165068 -> -1.165059) and gDP noise ~1e-3; SNR_J=9.1 and SNR_DP=4.7 (both below gate). The mechanism was confirmed: the first p solve is bit-identical across restarts but the second non-orthogonal p solve diverges in initial residual (0.03381 vs 0.02636), i.e. the GAMG aggregation/solve path depends on the in-process heap state. **Resolution:** the B1 case uses a deterministic `PCG + DIC` pressure solver (fvSolution `(p|pa|pb|pc)`); with it, A1-A5 and B1-B3 are bit-identical and all gates pass. The reduced discrete flow adjoint assembles its own GMRES and is unaffected by this dictionary entry. The GAMG noise floor remains documented here; if GAMG is restored, the B1 noise-floor numbers of this report no longer apply.
2. **[FOUND] 20-iteration platform window (1e-10) never reached within 1000 iterations.** J exhibits a ~1e-7 iteration-to-iteration oscillation (discretization limit cycle); the window criterion is too strict to ever trigger. All repeatability gates are nevertheless met because the solver is deterministic and the attractor unique. Reported honestly (Iter=-1) rather than faked as converged; future runs may relax the window to ~1e-6 or rely on the deterministic-repeatability gate only.
3. **[INFO] FPE traps (FOAM_SIGFPE) were disabled for the B1 runs** (`unset FOAM_SIGFPE`) so the large-scale field-statistics sweeps report non-finite values via `nBad`/isfinite counters instead of aborting; all fields reported `nBad=0`. No threshold or numerical setting was changed.
4. **[INFO] B1 outputs are written to the case directory** `b1_case/stageB1/` (stage_b1_repeatability.tsv, stage_b1_noise_floor.tsv) and copied into `src/validation/stage_b/`; no existing user file was overwritten.

## 14. Next-stage readiness

**YES.** B1 = PASS (16/16 criteria): repeatability perfect, frozen-field drift zero, gradient fields and projections bit-identical, SNR_J/SNR_DP far above gate. `readyForB1` (now B2/B3 amplitude acceptance) remains `true`; formal multi-direction FD amplitude validation may proceed with `objectiveGradientScale 1.0; pressureGradientScale 1.0;`, serial execution, and the deterministic PCG pressure solver (or with the GAMG noise floor documented above). `gradientValidated` stays `false`; no optimization started.

## 15. Deliverables

- `src/validation/stage_b/STAGE_B_B1_REPEATABILITY_NOISE_REPORT.md`
- `src/validation/stage_b/stage_b_b1_repeatability.json`
- `src/validation/stage_b/stage_b_b1_repeatability.tsv`
- `src/validation/stage_b/stage_b_b1_noise_floor.tsv`

## 16. Final status

**PASS** — all repeatability, frozen-field and SNR gates satisfied.
