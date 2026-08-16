# D4 — Gap analysis + next-Gate specification (B-F2 / B-F3 single-convention spec)

Stage: B-final (DIAGNOSTIC_ONLY, cycle 2)
Step: D4 of plan-002-replan (hypothesis BFINAL-INDETERMINATE-PROBE-CONVENTION)
Author: Executor (flash), evidence-only; no src/ or case modification.
Date/HEAD: workspace /home/ys/dsH/TO-ANISOTROPIC @ df695c7 (uncommitted B4 edits
preserved); reference repo /home/ys/TO-ANISOTROPIC @ 292c711.
Synthesizes D2 (audit.md), D3 (rx_audit.md), D1 (rerun.log), and the NEW
read-only Python cross-checks (crosscheck_bf2.out, localize_prow.out) run in
this step against the exported artifacts in /home/ys/b2_case_smoke.

---

## 0. Verdict (one paragraph)

**No single survivor is named: (a) J != R_w, (b) incomplete R_x, and
(c) RHS/sign/final-assembly all remain OPEN**, consistent with the approved
hypothesis. The new read-only evidence *narrows* the search in one important
way and *identifies one concrete, quantified, previously-untested candidate
mechanism* for (a):

1. **Momentum row J*v is verified against the actual frozen residual FD to
   relL2 = 1.65e-4, cos = 0.99999999** (crosscheck_bf2.out: J·v from the
   exported complete J^T vs NS.H `stageB2_ruDudU.mtx`, the central FD of the
   real frozen momentum residual `A*U - source`, same sine direction). This is
   the strongest momentum-row verification yet — it is cleaner than T1 (0.85%)
   because the reference carries no `+2*H_ldu` formula error. The momentum
   row is therefore *not* where the end-to-end 52x/2.3x/sign-flip lives.
2. **Continuity-row J*v vs the exported rebuilt-phi FD is 0.824 overall but
   100% boundary-concentrated**: internal cells relL2 = 1.09e-7 (cos = 1.0),
   boundary cells relL2 = 1.41 (cos = 0.246), and the in-code same-convention
   GatePR test (boundary flux rebuilt, not pinned) passes at 1e-10..4e-9
   (rerun.log L1230-1245). The 0.824 is a reference-side construction
   artifact of the NS.H rebuilt-phi FD (it pins boundary phi to base,
   NS.H L676-681), NOT a J defect — **but it is another demonstration of the
   probe-convention chaos the hypothesis predicts**.
3. **NEW candidate mechanism for (a), quantified: J's SIMPLE dphi/dU dH-part
   uses the UNRELAXED rAUAdj = 1/discretePrimalMomentum.A() (= 7.0584e-07 avg),
   while the primal phi the state actually satisfies is rebuilt with the
   RELAXED rAU = 1/UEqn.A() after UEqn.relax() (= 2.8234e-07 avg) — exactly
   factor 2.5 = 1/U-relax** (verified numerically from rerun.log L145/L165/
   L965 and from the OpenFOAM-7 source `fvMatrix::relax`: `D /= alpha`,
   alpha = 0.4). Every existing phi-rebuild FD reference (GateOracle L1467/
   L1486, GatePR L1477-1505, NS.H rebuilt L634-668) rebuilds phi with a fresh
   UNRELAXED matrix, so J matches them at 1e-10 **by construction** — none of
   them can see that the actual frozen residual's phi tangent carries the
   relaxed rAU. First principles (dA/dU = 0 verified "deltaA |dA_FD|L2=0",
   rerun.log L163) give dphi/dU = flux(rAU·dH/dU), so a rebuilt-with-relax
   FD predicts J's dH part is ~2.5x too large — a testable factor in the
   continuity row and the velocityJump momentum term. This is the minimal
   single-convention B-F2 discriminator below.

No gate is passed by this round's evidence; B-F2/B-F3 below are the next
numerical Gates (handoff §8), specified in single-convention form with exact
formulas, ordering, epsilon sweeps and thresholds.

---

## 1. Status of (a)/(b)/(c) after D4 (all OPEN — no survivor)

| Candidate | Status after D4 | Constraint / new evidence | Closer needed |
|---|---|---|---|
| (a) J != R_w | **OPEN** | Momentum row clean to 1.65e-4 (NEW, crosscheck_bf2.out). Continuity row self-consistent with the UNRELAXED phi rebuild to 1e-10 (GatePR) but the ACTUAL phi is RELAXED-rAU built (NS.H L45+L113+L129; rebuilt-baseline |phiRb0-phi|=3.5e-11, rerun.log L151) — so the true R_P row dH part is a ~2.5x candidate error (NEW). VelocityJump uses the same dphi. Boundary: constrainHbyA/adjustPhi variations unmodeled on both J and GatePR. | **B-F2** with the RELAXED primal phi rebuild |
| (b) R_x incomplete | **OPEN (leading)** | Rx-B self-referential, sine direction, relaxed-diag mix (D3). Code R_x = momentum-Brinkman only (sensitivity.H L65-68); first principles dRp/dalpha != 0; Rx-B itself measures NONZERO (8.1e-10). | **B-F3** design-direction fixed-state FD, split R_U,x / R_P,x |
| (c) RHS/sign/assembly | **OPEN** | g_w = 7.08e-12 anchored; transpose self-consistent 3.3e-16; direct solve 3.48e-9. End-to-end chain (sensitivity.H L65-68 = -dAlphaDxh*(U&Uc)*V; AdjNS_PD.H L3-4 = pressureConstraintDerivative) has no independent FD. | B-F4 tangent oracle, B-F5 tangent-adjoint identity (after (a)/(b)) |

Ruled out (unchanged): linear solver, transpose self-consistency, projection
chain, volume chain, g_w RHS.

---

## 2. New read-only evidence produced in D4

Scripts: `crosscheck_bf2.py`, `localize_prow.py` (evidence/agent-group/
BFINAL-001/cycle-2/). Outputs: `crosscheck_bf2.out`, `localize_prow.out`.
All reads confined to /home/ys/b2_case_smoke/*.mtx and the mesh
owner/neighbour/boundary dictionaries; nothing written into the case.

### 2.1 Method (read-only scatter matvec)

explicitJT.mtx is the COMPLETE physical J^T (134400x134400, 13210746 nnz,
pRef row pinned to identity, ExplicitJToracle maxRelL2 = 3.9e-16, rerun.log
L1263). J*v is computed by scattering J^T entries: for each (r,c,val),
out[c] += val*v[r]  (J[c,r] = J^T[r,c]). Sanity: J*v P-row == matrix-free
P-row export (stageB2_Jdq_P_forward) to relL2 = 1.24e-15.

### 2.2 Momentum row (crosscheck_bf2.out)

```
direction dUdir = sin(0.271*(celli+1)) * (1, 0.5, -0.3)   (NS.H L480-484,
same as solveDiscreteFlowAdjoint.H L1046-1049)
v = [dUdir; 0]
J*v (momentum rows)  vs  stageB2_ruDudU.mtx
   (central FD of R_U = A*U - source, fvm::div(phi,U)-laplacian+Sp,
    fvOptions.constrain, NO dev2, NO grad(p), frozen phi, hU=0.388):
relL2 = 1.653123e-04   cos = 0.999999986
|Jv_U|L2 = 1.477036   |ruDudU|L2 = 1.477043
```
Caveat: the ruDudU reference is FROZEN-phi (phi not perturbed) and excludes
dev2, while J's momentum row includes the deviatoric term and the
velocityJump(dphi) term; the 1.65e-4 agreement therefore also confirms that
dev2 + velocityJump are small for this direction (consistent with T-dev =
5.4e-5). It does NOT validate the dphi convention itself (that is B-F2's job
with rebuilt phi).

### 2.3 Continuity row (crosscheck_bf2.out + localize_prow.out)

```
J*v (P rows)  vs  stageB2_rpDudU_rebuilt.mtx (NS.H rebuilt-phi P-row FD,
                boundary phi PINNED to base, L676-681):
ALL cells:      relL2 = 8.240837e-01  cos = 0.695445  |diff|L2 = 2.69e-06
internal cells: relL2 = 1.090701e-07  cos = 1.000000  |diff|L2 = 2.89e-13
boundary cells: relL2 = 1.407930e+00  cos = 0.246014  |diff|L2 = 2.69e-06
(7392 boundary cells of 33600; internal 26208)

J*v (P rows)  vs  stageB2_Jdq_P_forward.mtx (J matrix-free P-row, same J):
relL2 = 1.239622e-15  cos = 1.0   (explicitJT scatter == matrix-free J)
```
Interpretation: the P-row discrepancy is ENTIRELY at boundary cells, where
the NS.H rebuilt-phi FD has zero boundary dphi by construction (base-fixed
boundary, admitted in NS.H L388-391), while J's P-row includes the boundary
dphi term (solveDiscreteFlowAdjoint.H L696-697). The in-code GatePR test —
whose FD uses the REBUILT boundary flux (L1497-1502) — passes including
boundary at relL2 = 1.1e-10 (h=1e-3) .. 3.8e-9 (h=3e-5) with a plateau
(rerun.log L1230-1245). **Conclusion: J's continuity row is consistent with
the UNRELAXED-rAU rebuilt phi at internal cells and with the rebuilt-boundary
same-convention FD everywhere; the open question is the rAU convention
(§2.4) and the constrainHbyA/adjustPhi boundary variations.**

### 2.4 Face-level dphi cross-check (localize_prow.out)

```
dphiJ (J's collectedDPhiJ)  vs  stageB2_dphiDudU_FD (GateOracle FD, rAUc
   = 1/discretePrimalMomentum.A() UNRELAXED, L1467/1486):   relL2 = 9.0e-10
dphiJ                      vs  stageB2_dphiDudU_face (NS.H rebuilt, rAUrb =
   1/tUEqnRb.A() UNRELAXED, L640):                          relL2 = 4.7e-8
```
All three references agree at the face level because they share the SAME
UNRELAXED rAU convention — none of them is a relaxed-rAU rebuild, so none
tests the actual primal phi tangent (§2.5).

### 2.5 The relaxation-convention gap (NEW, quantified)

- NS.H primal: `UEqn.relax()` (L45, U-relax = 0.4 in system/fvSolution),
  then `rAU = 1/UEqn.A()` (L113) — **relaxed** rAU (avg 2.8234e-07,
  rerun.log L145/L165); HbyA = rAU*UEqn.H(); phiHbyA = flux(HbyA);
  adjustPhi; pEqn = laplacian(rAtU,p); phi = phiHbyA - pEqn.flux()
  (L1016-1041). The stored phi is the RELAXED-rAU SIMPLE reconstruction:
  the rebuilt baseline with rAUrb0 = 1/UEqn.A() (relaxed) matches phi to
  |phiRb0-phi|L2 = 3.49e-11 (rerun.log L151).
- J's phi tangent (solveDiscreteFlowAdjoint.H L625-641): dH part uses
  `rAUAdj = 1/discretePrimalMomentum.A()` — **unrelaxed** (discretePrimal-
  Momentum built at L22-27 with NO relax call; avg 7.0584e-07, rerun.log
  L965); dp part uses `mobF = primalPressureMobility = rAtU` — **relaxed**
  (NS.H L129). So J's dphi is internally MIXED: dH part unrelaxed (2.5x too
  large vs the actual phi), dp part relaxed (correct, T2 = 7.5e-21).
- Verified from /opt/openfoam7 fvMatrix.C: `relax(): D /= alpha` with
  alpha = 0.4 → A_relaxed = A_unrelaxed/0.4 = 2.5·A_unrelaxed →
  rAU_relaxed = 0.4·rAU_unrelaxed. Logs confirm the ratio exactly:
  rAUAdj/rAtU = 2.500000009, A()_adj/A()_NS = 2.500000000 (avg).
- With dA/dU = 0 ("deltaA |dA_FD|L2=0", rerun.log L163), rAU is U-independent
  and dphi/dU = flux(rAU·dH/dU): the true (relaxed) dphi/dU dH part is
  **0.4x** J's current value. Every existing probe uses unrelaxed rAU on both
  sides, so none can see this factor.

This is the concrete, quantified, untested link that B-F2 must resolve.
It is NOT yet proven (the relaxation factor could in principle cancel through
H()/constrainHbyA/adjustPhi details); B-F2 with the exact primal rebuild
measures it directly. If confirmed, it is a complete (a)-mechanism candidate
for the observed 52x / 2.3x / sign-flip (a 2.5x error in a block of J^T,
inverted through the Schur complement, can produce exactly such end-to-end
amplitude/ sign distortions).

---

## 3. B-F2 spec — single-convention J*v vs central FD of the ACTUAL residual

Goal (handoff §8): close (a) with a reference that uses the ACTUAL frozen
primal residual definition — phi rebuilt EXACTLY as the primal does (relaxed
rAU, constrainHbyA, adjustPhi, pEqn.flux()), not a second implementation
sharing J's analytic formulas, and not an unrelaxed rebuild.

### 3.1 Exact reference residual (the FD side)

For a candidate state (U,p) at the converged design point, with frozen
turbulence (nuEffFrozen, alpha fixed):

```
R_U(U,p;phi)  = A(phi)·U - source          (momentum; unrelaxed matrix,
               fvm::div(phi,U) - fvm::laplacian(nuEffFrozen,U)
               + fvm::Sp(alpha,U), + fvOptions.constrain; dev2 folded into
               source when ransFlowModel — match NS.H L419-426/L501-526;
               NO grad(p) — dp=0 in the U-column test so its variation is 0)
R_P(U,p)      = div(phi_SIMPLE(U,p))        (continuity; phi REBUILT through
               the EXACT primal path, RELAXED matrix:
                 UEqn  = div(phi,U) - laplacian(nuEffFrozen,U) + Sp(alpha,U)
                 UEqn.relax()                       // U-relax = 0.4
                 rAU   = 1/UEqn.A()
                 HbyA  = constrainHbyA(rAU*UEqn.H(), U, p)
                 phiHbyA = fvc::flux(HbyA);  adjustPhi(phiHbyA, U, p)
                 rAtU  = rAU (consistent=false)
                 pEqn  = fvm::laplacian(rAtU, p)    // p FROZEN (state)
                 phi   = phiHbyA - pEqn.flux()      // boundary flux REBUILT,
                                                     // NOT pinned to base
               R_P   = div(phi) over internal + boundary faces)
```

Central FD along direction v (real design-direction or deterministic
non-neutral direction, NOT sine-only; see §3.3):

```
FD_U = [R_U(U+eps v, p; phi(U+eps v,p)) - R_U(U-eps v, p; phi(U-eps v,p))]/(2 eps)
FD_P = [R_P(U+eps v, p) - R_P(U-eps v, p)]/(2 eps)
```

The momentum residual is evaluated at the REBUILT phi (so the phi-mediated
convection variation is captured on BOTH rows, as in the reduced system), not
at the frozen base phi. This is the crucial single-convention choice: J's
velocityJump + continuity row both use dphi_SIMPLE, so the FD must too.

### 3.2 J side

```
J*v  with  v = [dU; dp]   (dU in physical units, dp = 0 for the U-column test)
via applyDiscreteFlowJ (solveDiscreteFlowAdjoint.H L414-711) — the SAME
operator used for the adjoint solve; do NOT use the exported explicitJT
transpose for the acceptance run (the scatter matvec is exact to 1e-15 but
the gate must test the actual code path).
```

### 3.3 Ordering / indexing (fixed convention)

```
N        = 33600 cells;  4N = 134400 unknowns
U index  = 3*cell + cmpt      (cmpt 0,1,2 = x,y,z)   [0 .. 3N-1]
P index  = 3N + cell          [3N .. 4N-1]
face dphi: internal faces fi = 0..96859 (owner = lowerAddr[fi],
           neighbour = upperAddr[fi]); boundary faces start at 96860
           (inlet 96860-96943, outlet 96944-97027, hotInlet 97028-97223,
           hotOutlet 97224-97419, solidEndWalls 97420-97699, bottomWall
           97700-98819, topWall 98820-99939, sideWalls 99940-104739).
```
Comparisons are component-wise on the physical (unscaled) rows, exactly as
the exports are written (U block 3N, P block N).

### 3.4 Directions

Use at least two deterministic non-neutral directions with strong signal:
- D1: `dU = dUdir = sin(0.271*(celli+1))*(1, 0.5, -0.3)` (existing;
  signal |Jv_U| = 1.477 — strong).
- D2: design-direction physical: map dAlphaDxh (sensitivity.H L13-15) into a
  velocity-like perturbation is NOT physical; instead use a second sine with
  a different seed, e.g. `cos(0.271*(celli+1))*(0.7, -0.4, 0.2)`, and one
  with dp != 0 (see §3.6). All directions must be non-neutral (avoid
  alternating-sign cancellations; check |FD| > noise floor per block).

### 3.5 Epsilon sweep / plateau

```
eps in {1e-3, 3e-4, 1e-4, 3e-5, 1e-5}   (relative to |U| scale, matching the
GatePR sweep h = 1e-3..3e-5 that already plateaued)
Report relL2 and cos at each eps; require a plateau: the last three eps
agree within a factor 2 in relL2, and relL2 does not grow as eps shrinks
(except the expected round-off rise below ~1e-5 on the small P rows).
Threshold per block (from the approved plan):
   momentum block relL2 < 1e-3  AND  continuity block relL2 < 1e-3
   (cos >= 0.999 as an additional sanity; |Jv|/|FD| in [0.99, 1.01] if the
   block is well-conditioned — note the P row is tiny, |Jv_P| ~ 3.6e-6, so
   use the relative threshold, not absolute).
```

### 3.6 Must ALSO test (the relaxation factor + dp part)

1. **rAU convention sub-test**: rebuild phi twice — (i) with UEqn.relax()
   applied (the ACTUAL primal convention), (ii) with the fresh unrelaxed
   matrix (current probe convention). Report both FD_P vs J*v. If (i) shows
   relL2 ~ 1.5 (i.e. J's dH part ~2.5x too large) while (ii) shows ~1e-10,
   the relaxation-convention candidate is CONFIRMED as the (a)-mechanism
   (§2.5) and the fix is a src change (use the relaxed rAU/rAtU consistently
   in the dphi dH part) — to be approved as a new hypothesis after this gate.
2. **dp part**: direction with dU = 0, dp = eps·cos(0.043*(celli+1)) (as T2);
   J_P vs div(rAtU·snGrad(dp)·|Sf|) — already 7.5e-21 (T2), keep as
   regression; also J_U(dp) momentum block vs FD (the -grad(p) coupling) —
   currently NOT separately FD-tested (T2 only checks the P row).
3. **boundary sub-test**: split the P-row comparison internal vs boundary
   cells (as in localize_prow.out); the boundary rows must use the REBUILT
   boundary flux on the FD side (GatePR convention), NOT the base-pinned
   boundary, and J's boundary dphi term (L696-697) must be included. Report
   the boundary-cell relL2 separately; if it stays O(1) with the correct
   convention, that localizes constrainHbyA/adjustPhi unmodeled variations —
   a distinct (a)-subcandidate for the next diagnostic.

Acceptance for B-F2: momentum < 1e-3 AND continuity < 1e-3 at the plateau
with the RELAXED-rAU rebuilt phi ⇒ J == R_w (a) closed; otherwise (a)
confirmed open with the failing block identified.

---

## 4. B-F3 spec — fixed-state design residual derivative R_x*d

Goal (handoff §8): close (b) with a fixed-(U,p) design-direction FD of the
ACTUAL reduced-system residual, split R_U,x and R_P,x, using the REAL design
direction (not sine), and answer "is R_P,x*d actually zero?".

### 4.1 Exact reference (FD side)

Keep U, p, and all frozen turbulence fields fixed at the converged state.
Perturb ONLY the design/material field:

```
x +/- eps*d_hat     with   d_hat = dAlphaDxh/|dAlphaDxh|  (real chain
derivative, sensitivity.H L13-15, masked+clipped as L25-35)
alpha_pm = alpha(x +/- eps d_hat)     (re-evaluate alpha from rho_bar/design)
```

Rebuild the reduced residual at the SAME (U,p) but the perturbed alpha
(reuse the B-F2 rebuild machinery, RELAXED rAU, constrainHbyA, adjustPhi,
pEqn.flux(), boundary REBUILT):

```
R_U,x*d  = [R_U(U,p; alpha+, phi+) - R_U(U,p; alpha-, phi-)]/(2 eps)
           (momentum row: dA(alpha)·U*V   +   phi-mediated div(dphi/dalpha (*) U))
R_P,x*d  = [R_P(U,p; alpha+, phi+) - R_P(U,p; alpha-, phi-)]/(2 eps)
           (continuity row: div(dphi/dalpha))
```

Analytic side (code, sensitivity.H L65-68 + NS.H RxProbe):

```
R_x,code = ( dAlphaDxh*(U&Uc_side)  ->  gsenshPressureDrop = -dAlphaDxh*(U&Uc)*V,
             R_P,x = 0 assumed )
```
Split and compare per row:
```
|R_U,x_FD - R_U,x_code| / |R_U,x_FD|     (momentum-Brinkman term = dAlphaDxh*U*V,
                                         Rx-A already relL2=1.6e-14 at fixed phi;
                                         the phi-mediated part div(dphi/dalpha*U)
                                         is the MISSING term to quantify)
|R_P,x_FD - 0| / |R_P,x_FD|              (continuity; code assumes 0)
```
Also report the three candidate missing terms of D3 §6, each relative to the
momentum-Brinkman term `|Uc*(dAlphaDxh*U*V)*d_hat|`:
```
T1 = |pc * R_P,x * d_hat|                (absent: no pc in gsenshPressureDrop)
T2 = |Uc * div(dphi/dalpha * U) * d_hat| (absent: phi-mediated momentum)
T3 = |Uc * d(boundary flux)/dalpha * d_hat| (absent: boundary flux derivative)
```

### 4.2 Epsilon sweep / thresholds

```
eps in {1e-6, 3e-7, 1e-7} x alphaMax   (hAlpha = 100 was 1e-6*1e8; sweep down)
Plateau required over the last two steps.
Thresholds (approved plan):
   |R_U,x*FD - R_U,x*code|/|R_U,x*FD| < 1e-3
   |R_P,x*FD|/|R_U,x*FD| < 1e-3   before declaring R_x complete
If R_P,x*FD is significant (>1e-3 of the momentum term), hypothesis (b) is
CONFIRMED (the continuity row and/or the phi-mediated momentum row are
missing from the assembly) and the code's momentum-Brinkman-only R_x is
incomplete.
```

### 4.3 Ordering and boundary handling

Same indexing as §3.3. Boundary phi REBUILT (not pinned) on the FD side;
report inlet/outlet-cell rows separately (the pressure-drop functional lives
there). d_hat must be the ACTUAL masked/clipped design direction, evaluated
once at the base state and reused for both +/- perturbations (do not re-mask
per perturbation).

---

## 5. Ordering of the next numerical Gates

1. **B-F2 (single-convention J*v vs actual-residual FD, §3)** — closes (a)
   and resolves the relaxation-convention candidate (§2.5). If J != R_w with
   the relaxed rebuild, stop adjoint work; the failing block names the next
   hypothesis.
2. **B-F3 (design-direction R_x*d, §4)** — closes (b); explicitly answers
   "is R_P,x*d zero?" per handoff §8.
3. **B-F4 (tangent oracle: J*w' = -R_x*d vs full design FD)** and
   **B-F5 (tangent-adjoint identity)** — close (c), only after (a)/(b)
   (handoff §8 B-F4/B-F5).

---

## 6. Feasibility of a read-only Python B-F2 (cross-check done in D4)

**Partially feasible with existing artifacts — with these exact caveats
(already exercised in §2):**

- J*v for any v: FEASIBLE read-only via the exported explicitJT.mtx (complete
  physical J^T; scatter matvec reproduces the matrix-free J P-row to 1.2e-15
  and the momentum row matches the frozen-phi FD to 1.65e-4).
- FD of the ACTUAL residual with phi REBUILT (relaxed rAU, constrainHbyA,
  adjustPhi, pEqn.flux()): NOT feasible read-only — the exported FD files
  (ruDudU, rpDudU_rebuilt, dphiDudU_FD) were all produced with UNRELAXED
  fresh matrices and/or base-pinned boundaries; a relaxed-rAU rebuilt-phi FD
  requires a solver-side probe (src change, forbidden this round). This is
  the exact gap B-F2 must fill with a new probe.
- Caveats to record: (i) J's momentum row includes deviatoric + velocityJump,
  the frozen-phi FD (ruDudU) excludes both — the 1.65e-4 agreement holds only
  because they are small for this direction; (ii) the continuity FD
  (rpDudU_rebuilt) is a hand-rolled rebuilt-phi SIMPLE with base-pinned
  boundary — use GatePR's rebuilt-boundary convention instead; (iii) all
  existing phi-rebuild FDs share J's unrelaxed rAU — none tests the actual
  (relaxed) phi tangent.
- Pure-Python sparse handling was adequate (13.2M nnz scatter matvec in 18.6 s,
  62 GB RAM); no scipy/numpy available in the environment (pure stdlib used).

---

## 7. Feed-forward / decisions for the Replanner

1. No hypothesis is closed by this round; the primary hypothesis
   BFINAL-INDETERMINATE-PROBE-CONVENTION is SUPPORTED and now quantified.
2. The relaxation-convention candidate (J dH part uses unrelaxed rAUAdj =
   2.5x the primal's relaxed rAU, untested by every existing probe) is the
   strongest concrete (a)-mechanism to date — B-F2 with the relaxed rebuild
   will confirm or refute it in one run.
3. B-F2/B-F3 require NEW solver-side probes (rebuild phi with the exact
   relaxed primal path; design-direction FD split R_U,x/R_P,x) — i.e., the
   NEXT round may be a minimal src/ probe addition under a newly approved
   hypothesis; it is NOT a J/R_x/final-assembly patch.
4. Keep: g_w anchor (2.7e-11 / 7.08e-12), Rx-A (1.6e-14), momentum-row
   1.65e-4 cross-check, GatePR 1e-10 as regression anchors.
5. Do not start MMA / SIMPLE-transpose / linear-solver work (handoff §5.4,
   §12).

---

## 8. Files/lines consulted (read-only)

- src/solveDiscreteFlowAdjoint.H: L22-27 (discretePrimalMomentum, no relax),
  L76-81 (rAUAdj = 1/A()), L101-108 (U/P indexing), L414-711
  (applyDiscreteFlowJ: dev L420-474, diag L476-484, dH offdiag L515-531,
  dH boundaryDiag L554-574, dH/V L576-582, matrix action L597-623, SIMPLE
  dphi L625-641, velocityJump L646-653, continuity L656-657, boundary
  L675-710), L1029-1128 (dPhi_J export, GateH1), L1444-1558 (Gate PR:
  rAUc = 1/discretePrimalMomentum.A() L1467, tEqR fresh unrelaxed L1477-1488,
  dphiDudU_FD export L1508-1513, GateOracle L1515-1558), L1749-1982 (T1/T-cont/
  T2), L2726-2820 (CSR assembly, pRef pin), L2896-2919 (explicitJT/rhs export).
- src/NS.H: L26-45 (UEqn + dev2 + relax), L113-129 (rAU/HbyA/phiHbyA/rAtU/
  primalPressureMobility), L143-167 (Rx-B sine direction, hAlpha=100),
  L169-190 (Rx-A), L192-410 (Rx-B), L414-463 (baseline resBase = -residual,
  UNRELAXED fresh matrix), L465-577 (frozen-phi U-column FD export), L579-716
  (rebuilt-phi P-row FD: tUEqnRb FRESH UNRELAXED L634-639, rAUrb L640,
  boundary PINNED L676-681), L1016-1041 (primal pEqn, phi = phiHbyA -
  pEqn.flux()).
- src/sensitivity.H L13-15 (dAlphaDxh), L65-68 (gsenshPressureDrop =
  -dAlphaDxh*(U&Uc)*V).
- src/AdjNS_PD.H L1-29 (discrete RHS: 0 / pressureConstraintDerivative).
- /opt/openfoam7 fvMatrix.C: relax() (`D /= alpha; S += (D-D0)*psi`),
  A(), H(), operator&, residual; lduMatrixTemplates.C H();
  fvmSup.C (fvm::Sp: diag += V*sp).
- case: system/fvSolution (U-relax 0.4, p-relax 0.2, SIMPLE consistent=false,
  pRefCell 5600), constant/optProperties (pressureDropMaxPa 25000),
  constant/polyMesh (33600 cells, 104740 faces, 96860 internal, boundary
  ranges §3.3).
- Evidence: rerun.log L119 (g_w), L144-164 (Rx-A/Rx-B/baseline/deltaH),
  L1203-1204 (dPhi_J stats), L1227-1245 (GateOracle/GatePR sweeps),
  L1250-1263 (T1/T-cont/T2/T-fvm/T-dev/T-conv, ExplicitJToracle);
  crosscheck_bf2.out; localize_prow.out; audit.md (D2); rx_audit.md (D3);
  AGENT_GROUP_HANDOFF.md §8 (B-F1..B-F5), §9, §19.
