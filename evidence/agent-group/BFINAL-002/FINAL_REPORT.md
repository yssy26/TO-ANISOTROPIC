# BFINAL-002 FINAL_REPORT — Boundary-consistent P-row and dphi/dU relaxation discrimination

- Stage: B-final · Mode: **DIAGNOSTIC_ONLY** · Hypothesis: `BFINAL-002-PROW-BOUNDARY-RELAXATION`
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) · HEAD `ca8a772b8339a36e7906ea32a7125061c47aa280` (== task baseline) · branch `agent/dsH-stage-b-validation`
- Date: 2026-08-16 (S4 report step; S1 audit + S2 probe + S3 run completed earlier today)
- **No production patch was made in this round** (task §11; Candidate B is diagnostic-only). The only tracked modification is the S2 guarded probe invocation (`src/solveDiscreteFlowAdjoint.H`, +20 lines, switch-guarded; `git diff --stat` = 1 file, 20 insertions).

---

## 0. Executive summary

**Q-A (Gate A).** The P-row O(1) error is **not** a boundary-flux probe artifact. Under a boundary-consistent FD that rebuilds phi through the *exact* production SIMPLE path (UEqn.relax → rAU=1/A → constrainHbyA → fvc::flux → adjustPhi → rAtU → constrainPressure → pEqn=laplacian(rAtU,p)+setReference(pRefCell) → phi=phiHbyA−pEqn.flux, boundary flux rebuilt, dp=0), the mismatch **persists at O(1)**: P-total relL2 = **1.036** (vs 0.824 before), P-internal relL2 = **0.858** (vs 1.09e-7 before — the error now also covers internal rows), P-boundary relL2 = **1.145** (vs 1.41). The old 0.824 measurement *was* partly an artifact (the old F1 reference pinned boundary fluxes at baseline → boundary dphi ≡ 0 in FD while J has a boundary dphi term → the error appeared 100 % boundary-concentrated), but fixing the boundary semantics does **not** collapse the error: J's P-row is genuinely not the derivative of the production relaxed residual.

**Q-B (Gate B).** **CONFIRMED** — the unrelaxed rAUAdj in J's dphi/dU tangent makes the production reduced Jacobian's dphi tangent wrong at O(1). Isolated oracle (p/design/alpha/k/omega/nutFrozen/nuEffFrozen fixed, deterministic dU, phi rebuilt with production semantics):
- **Candidate A** (current J formulation, unrelaxed rAUAdj): normRatio = 0.591, relL2 = **0.795** (internal 0.790, boundary 1.388), cos = 0.607 — O(1) mismatch at every eps.
- **Candidate B** (diagnostic-only production relaxed mobility: relaxed rAU × relaxed-H tangent incl. relax source + constrainHbyA pinning): normRatio = 1.000, relL2 = **1.06e-11** (roundoff), cos = 1.000 — matches the true production dphi/dU at roundoff for all 5 eps, internal *and* boundary faces.

All five task §6 patch-authorization criteria are met → a relaxation patch is *authorized* for the next Agent Group round (not executed here).

**Recommended next action: A — PATCH_J_RELAXATION** (exactly one; see §8 for the exact patch spec = Candidate B semantics, *not* a naive 0.4× rescale).

---

## 1. Q-A (Gate A): boundary-consistent P-row B-F2

### 1.1 The three phi-boundary conventions (S1 audit, confirmed at ca8a772b)

| Convention | Boundary-flux tangent dphi_b/dU (dp=0) | Who |
|---|---|---|
| **P** (production state map) | 0 at the 7 fixedValue/noSlip U patches (constrainHbyA pins HbyA_b = U_b); `Sf_b·d(rAU_rel·H_rel)_b/dU` at the zeroGradient-U outlet | NS.H L45, L113–116, L1020–1032 |
| **J** (reduced Jacobian) | `Sf_b·(rAUAdj·dH_b)` at **all 8** patches (extrapolated, unrelaxed; no constrainHbyA) | applyDiscreteFlowJ L689–697 |
| **F1** (old BFINAL-001 FD reference) | **0 everywhere** (boundary flux pinned to stored base phi) | NS.H L676–681 |

adjustPhi / constrainPressure / pEqn.setReference are **no-ops** in this case (p.needReference()==false, no fixedFluxPressure patches) — verified from `/opt/openfoam7` sources (S1). Production boundary semantics reduce to **constrainHbyA pinning (7 patches) + outlet extrapolation**.

### 1.2 Oracle construction (S2 probe, `src/stageB5BoundaryRelaxOracle.H`, switch `stageB5BoundaryRelaxOracle`, default false)

- Deterministic direction dUdir = sin(0.271·(celli+1)) pattern (identical to BFINAL-001 `stageB2_dUdirPhys.mtx`); v = [dU; dp=0].
- Jv = `applyDiscreteFlowJ(v)` (existing matrix-free reduced Jacobian; `collectDPhiJ=true` captures J's internal-face dphi = Candidate A internal).
- FD reference = central difference of the **full production reduced residual** R = (R_U, R_P), R_U = −UEqn.residual() (relaxation-invariant), R_P = div(rebuilt phi) with the phi rebuilt through the exact NS.H SIMPLE sequence listed above and the **boundary flux rebuilt** (not pinned).
- 5 eps {1e-3, 3e-4, 1e-4, 3e-5, 1e-5}; per block (momentum 3N / P-internal 26208 cells / P-boundary 7392 cells / P-total N) report |Jv|, |FD|, |err|, relL2, cos, max-abs-diff+location.

### 1.3 Results (eps = 1e-3; the table is eps-stable to ≥7 significant digits at all 5 eps — converged FD plateau, no truncation/roundoff contamination of the systematic mismatch)

| Block | |Jv| L2 | |FD| L2 | |err| L2 | relL2 | cos | max \|diff\| @ |
|---|---|---|---|---|---|---|
| momentum (3N) | 1.47703641607 | 1.47704274460 | 2.44171949431e-4 | **1.65311e-4** | 0.999999986 | 1.4286e-5 @ idx 18693 |
| P-internal (26208) | 2.65039396759e-6 | 2.79956523157e-6 | 2.40303538019e-6 | **0.858360** | 0.612374 | 8.4864e-8 @ cell 6161 |
| P-boundary (7392) | 2.42451109717e-6 | 3.31295910394e-6 | 3.79456365009e-6 | **1.145370** | 0.152835 | 1.9695e-7 @ cell 7761 (sideWalls) |
| P-total (33600) | 3.59205265047e-6 | 4.33742590832e-6 | 4.49146883915e-6 | **1.035515** | 0.370431 | 1.9695e-7 @ cell 7761 (sideWalls) |

Epsilon table (relL2 / cos; identical at all eps to the digits shown):

| eps | momentum | P-internal | P-boundary | P-total |
|---|---|---|---|---|
| 1e-3 | 1.653e-4 / 0.999999986 | 0.85836 / 0.612374 | 1.14537 / 0.152835 | 1.03551 / 0.370431 |
| 3e-4 | 1.653e-4 / 0.999999986 | 0.85836 / 0.612374 | 1.14537 / 0.152835 | 1.03551 / 0.370431 |
| 1e-4 | 1.653e-4 / 0.999999986 | 0.85836 / 0.612374 | 1.14537 / 0.152835 | 1.03551 / 0.370431 |
| 3e-5 | 1.653e-4 / 0.999999986 | 0.85836 / 0.612374 | 1.14537 / 0.152835 | 1.03551 / 0.370431 |
| 1e-5 | 1.653e-4 / 0.999999986 | 0.85836 / 0.612374 | 1.14537 / 0.152835 | 1.03551 / 0.370431 |

Comparison vs the BFINAL-001 reference (F1: unrelaxed fresh matrix, boundary pinned):

| | F1 (old) | P-consistent (new) |
|---|---|---|
| P-internal relL2 | 1.09e-7 (J == unrelaxed convention internally) | **0.858** (FD now uses the *relaxed* production tangent) |
| P-boundary relL2 | 1.41 (FD boundary dphi ≡ 0 by pinning) | **1.145** |
| P-total relL2 | 0.824 | **1.036** |

### 1.4 Per-patch boundary decomposition (independent cross-check X1, from the raw artifacts)

| patch (U BC) | nFaces | \|dphi_FD\| L2 (production) | \|dphi_A\| L2 (J) | A-vs-FD relL2 / cos |
|---|---|---|---|---|
| inlet (fixedValue) | 84 | **0.0** | 2.559e-9 | — |
| outlet (zeroGradient) | 84 | 1.5781e-6 | 1.5773e-6 | 0.169 / 0.986 |
| hotInlet (fixedValue) | 196 | **0.0** | 8.12e-12 | — |
| hotOutlet (fixedValue) | 196 | **0.0** | 8.13e-12 | — |
| solidEndWalls (noSlip) | 280 | **0.0** | 1.21e-9 | — |
| bottomWall (noSlip) | 1120 | **0.0** | 1.41e-8 | — |
| topWall (noSlip) | 1120 | **0.0** | 4.06e-12 | — |
| sideWalls (noSlip) | 4800 | **0.0** | **2.174e-6** | — |

- Production boundary dphi ≡ **exactly 0** at the 7 constrainHbyA-pinned patches; nonzero only at the outlet (1.5781e-6).
- J's Candidate A boundary dphi is nonzero at **all 8** patches — dominated by sideWalls (2.174e-6, a noSlip wall where production has 0) → the missing constrainHbyA pinning is part of the boundary gap; at the outlet the residual gap (relL2 0.169, cos 0.986) is the relaxation convention.

### 1.5 Oracle trustworthiness (independent cross-check X2)

Reconstructing `div(dphi_FD)` per cell from the mesh (my own owner/neighbour + Sf from faces/points + cell volumes) reproduces the Gate-A per-cell P-row FD **exactly**: relL2 = 3.19e-11, cos = 1.000000000 (internal 4.7e-11, boundary 1.4e-11). Both oracles therefore measure the identical derivative object (div of the production phi tangent). (Note: the probe's P-row residual convention is the raw face-flux sum Σ_f φ_f — no 1/V — shared by J and FD; the probe's own selfcheck of J's P rows vs its divergence assembly is 1.25e-17.) Further trust anchors reproduced in the same runs: momentum 1.6531e-4 / cos 0.999999986; GatePR 1.12295958888e-10 (J vs the extrapolated-unrelaxed F2 FD); |phiRb0−phi| = 3.49247878283e-11; NS UEqn A() avg 158965329.116; frozen-RANS 44 correctors.

### 1.6 Q-A verdict

The O(1) P-row error is **real** (survives a boundary-consistent, production-semantics oracle) and its root cause is **not** the probe's boundary-flux handling: it is the **relaxation convention in J's dphi tangent** (unrelaxed rAUAdj, no relax-source term), i.e. the same root cause as Q-B (row 2, CONFIRMED), plus the missing constrainHbyA pinning at the 7 fixedValue-U patches on the boundary part. The old 0.824's *boundary concentration* was partly a pinning artifact (F1), but the O(1) *magnitude* was not. → **GENERAL_P_ROW_JACOBIAN_DEFECT = NOT_SUPPORTED as a blanket statement; a localized, fully explained J-vs-production gap exists (relaxation convention), quantified at roundoff by Gate B.**

---

## 2. Q-B (Gate B): isolated dphi/dU relaxation oracle

### 2.1 Construction

- dU = deterministic dUdir; p / design / alpha / k / omega / nutFrozen / nuEffFrozen **fixed**.
- phi(U±eps·dU) rebuilt through the exact NS.H production sequence (fresh UEqn → .relax() → fvOptions.constrain → rAU=1/A → constrainHbyA → fvc::flux → adjustPhi → rAtU → constrainPressure → pEqn=laplacian(rAtU,p)+setReference → phi=phiHbyA−pEqn.flux()).
- FD_dphi = [phi(U+)−phi(U−)]/(2eps) per face.
- **Candidate A** = current production/J formulation: internal = collectedDPhiJ (J's dphi), boundary = Sf_b·(rAUAdj·dH_J) at all patches (unrelaxed).
- **Candidate B** (diagnostic-only, NOT installed in production) = production relaxed mobility: internal Sf·(w·rAU_rel·dH_rel_o + (1−w)·rAU_rel·dH_rel_n), boundary 0 at the 7 non-assignable U patches (constrainHbyA) + extrapolated Sf_b·(rAU_rel·dH_rel) at the outlet; rAU_rel = production relaxed rAU (=primalPressureMobility); dH_rel = exact tangent of the **production relaxed H** (FD of the reassembled relaxed UEqn at U±eps·dU, affine ⇒ exact; includes the relax source and dev2 source). **Not** a naive "0.4×rAUAdj × unrelaxed dH".

### 2.2 Results (all 5 eps; Candidate-A numbers eps-independent = systematic O(1); Candidate-B numbers grow 1e-11→1e-9 as eps→1e-5 = pure roundoff 1/eps scaling)

| eps | Candidate | \|analytic\| L2 | \|FD\| L2 | norm ratio | relL2 (all) | relL2 int / bnd | cos | max \|diff\| @ |
|---|---|---|---|---|---|---|---|---|
| 1e-3 | A (J, unrelaxed) | 1.298506e-5 | 2.197735e-5 | **0.590838** | **0.794617** | 0.790385 / 1.387680 | **0.607335** | 1.486e-7 @ face 83557 |
| 1e-3 | B (production relaxed) | 2.197735e-5 | 2.197735e-5 | **1.000000** | **1.0592e-11** | 1.0576e-11 / 1.3307e-11 | **1.000000** | 1.55e-17 @ face 24222 |
| 3e-4 | A | 1.298506e-5 | 2.197735e-5 | 0.590838 | 0.794617 | 0.790385 / 1.387680 | 0.607335 | 1.486e-7 @ 83557 |
| 3e-4 | B | 2.197735e-5 | 2.197735e-5 | 1.000000 | 3.4705e-11 | 3.4656e-11 / 4.3158e-11 | 1.000000 | 5.49e-17 @ 21583 |
| 1e-4 | A | 1.298506e-5 | 2.197735e-5 | 0.590838 | 0.794617 | 0.790385 / 1.387679 | 0.607335 | 1.486e-7 @ 83557 |
| 1e-4 | B | 2.197735e-5 | 2.197735e-5 | 1.000000 | 1.0709e-10 | 1.0683e-10 / 1.4910e-10 | 1.000000 | 1.72e-16 @ 24882 |
| 3e-5 | A | 1.298506e-5 | 2.197735e-5 | 0.590838 | 0.794617 | 0.790385 / 1.387680 | 0.607335 | 1.486e-7 @ 83557 |
| 3e-5 | B | 2.197735e-5 | 2.197735e-5 | 1.000000 | 3.5385e-10 | 3.5277e-10 / 5.2216e-10 | 1.000000 | 5.70e-16 @ 24383 |
| 1e-5 | A | 1.298506e-5 | 2.197735e-5 | 0.590838 | 0.794617 | 0.790385 / 1.387679 | 0.607335 | 1.486e-7 @ 83557 |
| 1e-5 | B | 2.197735e-5 | 2.197735e-5 | 1.000000 | 1.0551e-9 | 1.0549e-9 / 1.0976e-9 | 1.000000 | 1.65e-15 @ 21111 |

rAU statistics (log line `StageB5 rAU:`, identical in both inclusion blocks and all runs):

| quantity | avg | min | max |
|---|---|---|---|
| rAUAdj (unrelaxed = 1/discretePrimalMomentum.A()) | 7.05844426338e-7 | 9.99981274586e-9 | 5.1321646829e-6 |
| rAUrel (production relaxed = primalPressureMobility) | 2.82337770533e-7 | 3.99992509834e-9 | 2.05286587567e-6 |

relaxFactorU (fvSolution relaxationFactors/equations/U) = **0.4** · avgRatio(rAUAdj/rAUrel) = **2.50000000002** · normRatio(|rAUAdj|/|rAUrel|) = **2.49999999999**.

### 2.3 Why the norm ratio is 0.591 and not 0.4 (independent cross-check X3)

From the OpenFOAM-7 source (fvMatrix::relax: `D /= alpha` with α=0.4, `S += (D−D0)·ψ`), the production HbyA tangent is

**dphi_prod = flux(0.4·rAU_u·dH_u + 0.6·dU)**   (dH_u = unrelaxed H tangent; the 0.6·dU term is the relax-source contribution (D_rel−D0)·dU entering H)

while J's Candidate A is dphi_A = flux(rAU_u·dH_u). Verified numerically from the artifacts: |FD − 0.4·A| = 1.9270e-5 vs |flux(0.6·dU)| = 1.9251e-5 (0.1 % magnitude match, cos = 0.9989). Hence the J formulation is missing **both** the 0.4 mobility factor **and** the 0.6·dU relax-source term — the observed norm ratio 0.591 (not 0.4 or 2.5) is the fingerprint of the partial compensation. **A naive "0.4 × rAUAdj" rescale would be wrong**; the correct patch must reproduce Candidate B (relaxed rAU × relaxed H incl. the relax source).

### 2.4 Task §6 PATCH-authorization criteria — all met

1. Candidate A shows a clear O(1) mismatch: relL2 = 0.795, cos = 0.607, norm ratio 0.591. ✓
2. Candidate B reduces relL2 to **1.06e-11** (preference <1e-3 exceeded by 8 orders). ✓
3. Candidate B cosine = 1.000 > 0.999. ✓
4. Multiple eps with a stable FD plateau (5 eps, Candidate-B roundoff scaling 1e-11→1e-9). ✓
5. The improvement is **not** from excluding boundary rows/faces: internal 0.790 → 1.06e-11 **and** boundary 1.388 → 1.33e-11 both improve to roundoff. ✓

### 2.5 Q-B verdict

**CONFIRMED** — the unrelaxed rAUAdj in J's dphi/dU tangent makes the production reduced Jacobian's dphi/dU incorrect at O(1), and the correct formulation is the relaxed primal-equivalent mobility with the full relaxed-H tangent (Candidate B). The observed 2.5× rAU ratio does **not** appear 1:1 in dphi (norm ratio 0.591) because the relax-source term partially compensates — so the "2.5 factor" is a clue, and Candidate B (not 0.4×rAUAdj) is the proof and the fix spec.

---

## 3. Post-Reviewer independent verification (task §10)

Performed independently for this report (not relying on the Executor's S2/S3 summaries):

1. **Fresh full-case rerun** (independent background job, same environment as `run_mtohf.sh` but invoking the workspace binary `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`, switch ON): EXIT = 0, ExecutionTime 495.6 s. Raw log `evidence/agent-group/BFINAL-002/stageB5_run_s4_independent.log`. Every StageB5 metric and every anchor (momentum 1.6531e-4 / cos 0.999999986, GatePR 1.12295958888e-10, |phiRb0−phi| 3.49247878283e-11, 44 correctors, rAU stats) is **byte-identical** to the S3 run; all 8 exported artifacts are **md5-identical** to `artifacts_s3/` ⇒ deterministic across 3 independent runs (S2 on2, S3, this rerun).
2. **Independent recomputation from raw artifacts** (`recompute_s4_independent.py`, numpy, deliberately different code path from `recompute_stageB5.py`; output `recompute_s4.out`): reproduces every Gate-A and Gate-B metric **exactly** (all relL2/cos/norm ratios to printed digits), plus three cross-checks the probe does not print:
   - **X1** per-patch boundary decomposition (§1.4): production dphi ≡ 0 on the 7 pinned patches; J nonzero at all patches (sideWalls-dominated).
   - **X2** divergence consistency (§1.5): div(dphi_FD) from my own mesh reconstruction == Gate-A P-row FD at relL2 3.2e-11, cos 1.0.
   - **X3** relaxation algebra (§2.3): production tangent = flux(0.4·rAU_u·dH_u + 0.6·dU), 0.1 % magnitude match, cos 0.9989.
3. **Source verification**: `git diff` shows only the guarded 20-line invocation in `src/solveDiscreteFlowAdjoint.H`; `src/stageB5BoundaryRelaxOracle.H` is new and switch-guarded; no production math line changed; switch-OFF runs are byte-identical to baseline (S2 `stageB5_run_off.log`, anchors identical to BFINAL-001 `rerun.log`).

---

## 4. Decision table (required)

| Hypothesis | Status | Evidence |
|---|---|---|
| P-row O(1) error caused by boundary-oracle mismatch | **REJECTED** | Boundary-consistent production FD does **not** collapse the error: P-total relL2 0.824 → 1.036 (not −2 orders); P-internal 1.09e-7 → 0.858; P-boundary stays O(1) at 1.145, all stable over 5 eps. The old 0.824's boundary *concentration* was partly a pinning artifact (F1 pinned boundary fluxes → boundary dphi ≡ 0), but the O(1) magnitude is real: J's P-row is not the derivative of the production relaxed residual. Localized gap: (a) unrelaxed rAUAdj + missing relax-source term in J's dphi tangent (L628–629 internal, L689–697 boundary) — the same root cause as row 2; (b) missing constrainHbyA pinning at the 7 fixedValue-U patches (production dphi_b ≡ 0 there; J extrapolates nonzero, sideWalls-dominated 2.17e-6). Both are fixed by making J's dphi = Candidate B semantics (X2: oracle self-consistent at 3e-11). |
| unrelaxed rAU causes incorrect dphi/dU | **CONFIRMED** | All 5 task §6 criteria met: Candidate A O(1) (relL2 0.795, cos 0.607, normRatio 0.591 at every eps); Candidate B relL2 1.06e-11 (≪ 1e-3), cos 1.000; 5-eps stable FD plateau; improvement covers internal *and* boundary faces (0.790→1.06e-11, 1.388→1.33e-11). Production tangent = flux(0.4·rAU_u·dH_u + 0.6·dU) verified algebraically (X3, cos 0.9989) and by Candidate B at roundoff. Note: the 2.5× rAU ratio does not appear 1:1 in dphi (normRatio 0.591) — the relax-source term 0.6·dU partially compensates; a naive 0.4×rAUAdj patch would be wrong. |
| missing R_x pressure-row derivative | **DEFERRED** | B-F3 not yet authorized (task §7). The residual semantics B-F3 must use are now settled: the production relaxed-rAU phi reconstruction (Candidate B semantics) survives Gate A as the reference. |

---

## 5. Exactly one recommended next action

### **A — PATCH_J_RELAXATION** (next Agent Group round; NOT executed in this DIAGNOSTIC_ONLY round)

- **What:** change J's P-row dphi/dU tangent (solveDiscreteFlowAdjoint.H SIMPLE dphi tangent, L625–644 internal and L689–697 boundary; rAUAdj usage L76–80/L628–629) to the production relaxed semantics **exactly as Candidate B**: use the production relaxed mobility (rAU_rel = primalPressureMobility = 1/A of the *relaxed* UEqn) with the *relaxed* H tangent **including the relax-source term** (i.e., the derivative of H as computed from the relaxed matrix, `flux(0.4·rAU_u·dH_u + 0.6·dU)` equivalent), and apply constrainHbyA boundary pinning (dphi_b = 0 at the 7 non-assignable U patches; extrapolated at the outlet). **Do NOT** implement as a naive `0.4·rAUAdj·dH` rescale (would miss the 0.6·dU relax-source term and stay wrong at relL2 ≈ 0.4–0.8 level).
- **Why A and not B/C/D:**
  - *Not B (PATCH_P_BOUNDARY_JACOBIAN)*: the O(1) error is not boundary-localized under production semantics — P-internal relL2 = 0.858. A boundary-only patch would leave the internal rows wrong; the boundary gap is part of the same Candidate-B semantics (pinning + relaxed mobility).
  - *Not C (START_BF3_RX_CLOSURE)*: the R_x design-derivative test (deferred by task §7) must use the residual semantics that survive Gate A — which is the *production* (Candidate B) phi tangent. Running B-F3 against the current (unrelaxed) J semantics would test the wrong Jacobian; the relaxation patch must land first.
  - *Not D (REPLAN_DIAGNOSTICS)*: both questions are answered decisively with a trustworthy, deterministic, cross-validated oracle (X1/X2/X3 + 3 identical runs); no diagnostic redesign is needed.
- **Verification for the patch round:** after patching, Gate A must give P-internal/P-boundary/P-total relL2 < 1e-3 with cos > 0.999 (i.e., J*v vs the production-consistent FD collapses), and the momentum row must stay at 1.65e-4 (R_U is relaxation-invariant, so the momentum block is untouched by this patch). Then re-validate the pressure-drop directional-derivative projection (D1–D3) before B-F3.

---

## 6. Compliance

- **DIAGNOSTIC_ONLY:** no production patch (no P-row J, no dphi/dU, no R_x, no sensitivity, no SIMPLE-transpose, no MMA). Candidate B exists only inside the guarded probe.
- Only files_allowed touched: `evidence/agent-group/BFINAL-002/**` (audit_boundary.md, stageB5_probe_design.md, recompute_stageB5.py, recompute_s4_independent.py, recompute_s3.out, recompute_s4.out, stageB5_build*.log, stageB5_run_{off,on,on2,s3,s4_independent}.log, artifacts/, artifacts_s3/, case_optProperties.*, S2/S3 summaries), `src/stageB5BoundaryRelaxOracle.H` (new, switch-guarded, no-op when disabled), guarded invocation in `src/solveDiscreteFlowAdjoint.H` (+20 lines).
- `git status` at report time: only `M src/solveDiscreteFlowAdjoint.H` (the 20 guarded lines) plus untracked new probe/evidence/build files; `git diff --stat` = 1 file / 20 insertions.

## 7. Artifacts (evidence/agent-group/BFINAL-002/)

- Raw logs: `stageB5_run_s3.log` (S3), `stageB5_run_s4_independent.log` (Post-Reviewer fresh rerun), `stageB5_run_off.log` (switch OFF, byte-identical baseline)
- Exported matrices: `artifacts_s3/stageB5_*.mtx` (8 files; md5-identical across S2/S3/fresh rerun and the case-dir exports)
- Recomputation: `recompute_s3.out` (recompute_stageB5.py), `recompute_s4.out` + `recompute_s4_independent.py` (independent numpy recomputation incl. X1/X2/X3)
- Documentation: `audit_boundary.md` (S1), `stageB5_probe_design.md` (S2), `stageB5_s2_summary.md`, `S3_run_summary.md`, `S3_build_record.md`, `case_optProperties.{baseline,run1,run2,s3run}`
