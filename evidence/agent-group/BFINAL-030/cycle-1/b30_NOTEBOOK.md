# BFINAL-030 NOTEBOOK — T3/(1-alphaRel) 语义审判轮

Task authority: `/home/ys/dsH/TO-ANISOTROPIC/NEXT_TASK_BFINAL030.md`
Cycle directory: `/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-030/cycle-1/`
Date: 2026-08-21
Executor: executor-working-time (BFINAL-030)

Real repo: `/home/ys/dsH/TO-ANISOTROPIC`, branch `agent/dsH-stage-b-validation`,
HEAD `cd9bd7a` (BFINAL-029 commit).  `/home/ys/TO-ANISOTROPIC` is an unrelated
empty shell and is NOT used.

Unified root-cause hypothesis under judgment:

> A (1-alphaRel) flux-tangent semantics error in the T3 slot simultaneously
> produces fingerprints  (1) gDP direction-uniform 2.05x factor,
> (2) b_TC T3-slot mismatch (93.5% of route-export mismatch, U row 1.82x
> amplitude gap), (3) J direction-dependent sign flips.

Known anchors: alphaRel = 0.4; production T3 formula `U += (1-alphaRel)*w*(Sf*g)`
(Production:363-369); U gap 1.82x with candidate explanations 1/(1-alphaRel)
= 1.667 or others — DECIDED BY DERIVATION, NOT PRESET.

---

## Stage 0 — Derivation (complete; written before any fingerprint data)

### 0.1 OF7 fvMatrix relax()/H()/flux() semantics (source-verified)

`/opt/openfoam7/src/finiteVolume/fvMatrices/fvMatrix/fvMatrix.C`

- `relax(const scalar alpha)`, body 521-668:
  - ~537: `scalarField D0(D);` — the unrelaxed diagonal is saved.
  - ~631-634: diagonal-dominance clamp `D[celli] = max(mag(D[celli]), sumOff[celli]);`
    (inactive for this case: `alphaRel*prodRAU == primalPressureMobility` to
    2.19e-15 in ALL cells, machine-verified B23, so the code path is the pure
    relaxed-diagonal branch).
  - ~637: `D /= alpha;` — relaxed diagonal D_rel = D0/alphaRel.
  - ~667: `S += (D - D0)*psi_.primitiveField();` — relax-source
    S_rel = D0*(1/alphaRel - 1)*U.
- `A()` at 738: `A = D()/V`.
- `H()` at 760: `H() = (boundaryDiag contribution + lduMatrix::H(psi) + source_)/V`,
  so the relax-source enters H (and only via source_).
- `flux()` at 863-946: `fieldFlux = lduMatrix::faceH(psi_cmpt)
  + internalCoeffs*psi_internal - boundaryCoeffs*psi_nei (+ faceFluxCorrectionPtr_)`.
- `operator==`:
  - `(matrix, matrix)` at 1364: returns `(A - B)` (matrix copy/assembly).
  - `(matrix, tmp<GeometricField>)` at 1427-1439: `source() += V*field`.

### 0.2 NS.H fixed-point construction (line-verified)

`/home/ys/dsH/TO-ANISOTROPIC/src/NS.H` (1213 lines):

- 28-36: `UEqn` construction `fvm::div(phi,U) - fvm::laplacian(nuEffFrozen,U)
  + fvm::Sp(alpha,U) == -fvc::grad(p) + fvOptions(U)`.  The construction-time
  `==` invokes `operator==(matrix, tmp<field>)` => `UEqn.source() += -V*grad(p)`.
- 47: `UEqn.relax();` (alphaRel = 0.4).
- 54: `solve(UEqn == -fvc::grad(p));` — the solve-time `==` again invokes
  `operator==(matrix, tmp<field>)`, adding `-V*grad(p)` to the TEMPORARY
  solve matrix only, so the PERSISTENT source retains exactly ONE `-V*grad(p)`
  copy (B24 T7 note, verified).
- 115: `rAU = 1.0/UEqn.A();`
- 116: `HbyA = constrainHbyA(rAU*UEqn.H(), U, p);`
- 117: `phiHbyA = fvc::flux(HbyA);`
- 120-129: `rAtU(rAU)` — consistent-check; SIMPLE dict has NO "consistent"
  entry => `rAtU() == rAU`.
- 131: `primalPressureMobility = rAtU();`
- 1022-1025: `pEqn fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA);`
- 1034: `phi = phiHbyA - pEqn.flux();`
- 1040-1041: `p.relax(); U = HbyA - rAtU()*fvc::grad(p);`

fvOptions: no fvOptions file in the case => chain CLOSED, ZERO contribution.
`alpha` in NS.H:32 is the Brinkman penalty (max 1e8), NOT the relaxation factor.

### 0.3 The relaxed-SIMPLE flux tangent (derivation core)

Define the code-path single-iteration map with `dU` = linearization increment
about the converged primal state (so `dU_old = dU`):

    HbyA  = rAU_rel * H_rel(U)
    rAU_rel = 1/A_rel = 1/(A_u/alphaRel) = alphaRel*rAU_u          [relax: D/=alpha]
    dH_rel/dU  = dH_u/dU  +  (dS_rel/dU)/V
               = dH_u     +  D0*(1/alphaRel-1)*dU/V
               = dH_u     +  A_u*(1/alphaRel-1)*dU                [A_u = D0/V]
    dHbyA = alphaRel*rAU_u*dH_u  +  alphaRel*rAU_u*A_u*(1/alphaRel-1)*dU
          = alphaRel*rAU_u*dH_u  +  (1-alphaRel)*dU                [rAU_u*A_u = 1]
    dphi_f = Sf & ( alphaRel*( w*rAU_o*dH_o + (1-w)*rAU_n*dH_n )
                  + (1-alphaRel)*( w*dU_o + (1-w)*dU_n ) )        (*)
             - kf_f*(dp_n - dp_o)        [kf = mob_f*deltaCoeffs*magSf]

Term count for the writeup: grad(p) appears in the UEqn.H() source ONCE
(persistent source copy at construction).  alphaRel appears in the combination
`alphaRel` on the H-mobility part and `(1-alphaRel)` on the direct-U part
(inverse-scaled relax-source derivative); fvOptions contributes 0; relax()
replaces the diagonal by D0/alphaRel and adds the source term
(D0/alphaRel-D0)*U; the B24 H7 channel `-interp(mob*grad dp).Sf` is the
second P-row channel and carries mob = primalPressureMobility = rAtU = rAU
(no alphaRel beyond the mobility itself).

(*) is EXACTLY the tangent implemented in the production operator:
`solveDiscreteFlowAdjointProduction.H` comment 334-345 and slots hA (359-362,
`alphaRel*prodRAU*w*faceTranspose`), T3 (363-369, `(1.0-alphaRel)*w*faceTranspose`
and `(1.0-alphaRel)*(1-w)*faceTranspose`), kf (370-378), H7 (381-385), and the
`applyProdFlowJT` directRel slot (617-619, `(1.0-alphaRel)*lambdaPdiff*prodSf`).
B18 `DERIVATION_FIX3.md` Sec.2 already asserts the SAME (1-alphaRel) direct-U
coefficient — so both the production operator and the B18 derivation agree on
what the CODE PATH implements.

### 0.4 The adjudication: code-path vs fixed-point semantics

Two candidate derivations of the coefficient that the ADJOINT operator must use:

- (a) CODE-PATH: the operator must linearize exactly what the code computes in
  ONE SIMPLE momentum step about the converged state => direct coefficient
  `(1-alphaRel) = 0.6`.  Under this semantics the (1-alphaRel) factor is CORRECT
  and the hypothesis would be falsified.
- (b) FIXED-POINT: the map (design -> converged (U,p)) solves the steady
  residual R(U,p)=0, and the SIMPLE relax factor is a per-iteration
  preconditioner that disappears at the fixed point (the relax-source
  derivative `(1/alphaRel-1)*A_u*dU` is a single-iteration artifact).  The
  correct steady-state flux tangent is `dphi = Sf&dU` (direct coefficient
  1.0), i.e. c* = 1.0.  This is the unified hypothesis.

The task explicitly requires the decision by DERIVATION, not preset.  The
derivation is essentially settled toward (b) FIXED-POINT: B24 GateOracle B
measures `dPhiJ` (code-path operator tangent) vs `dPhiHfd` (hand-FD of the
TRUE primal) at relL2 = 1.36 / cos 1.06, while `dPhiHfd` vs full-FD closes at
relL2 = 1.65e-11.  I.e. the code-path operator tangent does NOT match the true
primal's tangent; the hand-FD of the true primal does.  The exported
`explicitJT_H7.mtx` is the code-path operator (ExplicitJToracle relL2 =
5.48e-16 = exported matrix == implicit operator).  Hence the mismatch lives in
the operator semantics itself — consistent with the (1-alphaRel) coefficient
being a single-iteration artifact.

PRIMARY preregistered correction coefficient: **c* = 1.0** (fixed-point).
Alternates (kept for the verdict's sensitivity discussion):
- c* = 1.09 (empirical 0.6*1.82, U-row amplitude gap fit),
- c* = 1.23 (0.6*2.05, gDP-only full-explanation fit).

### 0.5 Exported-operator layout (prerequisite for the fingerprint-1 patch)

Source: `/home/ys/dsH/TO-ANISOTROPIC/src/solveDiscreteFlowAdjoint.H`
- Index helpers: `discreteUIndex(celli,cmpt) = 3*celli+cmpt` (113-116),
  `discretePIndex(celli) = discreteNVelocity + celli` (117-120).
- Export: 3451-3482 writes `explicitJT_H7.mtx` (lookupOrDefault
  "discreteExplicitMatrixFile") row by row from the CSR of the EXPLICIT J^T
  assembly (2626-2629: "Explicit assembly of J^T ... replicates
  applyDiscreteFlowJT() exactly").  Matrix values are the raw coefficients
  (no momentum/area scaling in the matrix; the RHS file applies
  discreteMomentumScale / discreteAreaScale, 3478-3480).
- Internal-face U<-P direct-relax entries (2710-2716):
  `addCoo(U(own,c),P(own), (1-alphaRel)*w*sf_c)`,
  `addCoo(U(own,c),P(nei), -(1-alphaRel)*w*sf_c)`,
  `addCoo(U(nei,c),P(own), (1-alphaRel)*(1-w)*sf_c)`,
  `addCoo(U(nei,c),P(nei), -(1-alphaRel)*(1-w)*sf_c)`.
- Assignable-boundary U<-P direct entries (~2996-2999):
  `addCoo(U(celli,c), P(celli), (1-alphaRel)*patchSf_c)`.
- P<-U pressure-gradient-transpose entries (2683-2690): `w*sf_c` / `-w*sf_c` /
  `(1-w)*sf_c` / `-(1-w)*sf_c` — NO alphaRel.
- The earlier `/tmp/b30_pu_verify.py` failure is explained: it assumed the
  matrix was J (P-rows/U-cols with (1-alphaRel)); the export is J^T, so the
  (1-alphaRel) entries are at U-rows/P-cols and the P-rows/U-cols carry the
  plain w*sf pressure-gradient transpose.

Patch recipe (fingerprint 1, coefficient-level, ZERO production changes):
  M_corrected = M + (c* - (1-alphaRel)) * T3_base
where T3_base is the sparse pattern with coefficient 1:
  internal:  T3_base[U(own,c),P(own)] += w*sf_c ; [U(own,c),P(nei)] -= w*sf_c
             [U(nei,c),P(own)] += (1-w)*sf_c ; [U(nei,c),P(nei)] -= (1-w)*sf_c
  assignable boundary: T3_base[U(cell,c),P(cell)] += patchSf_c
This is exact regardless of hA/hB overlap at the same (row,col) because the
T3 contribution there is exactly (1-alphaRel)*T3_base and the other blocks are
untouched.

### 0.6 Three-fingerprint expected values (preregistered BEFORE data)

PREREGISTRATION below (b30_PREREGISTRATION.md) — written before any
fingerprint computation.

---

## Stage 1 — execution log (chronological)

All fingerprints computed OFFLINE (route-instrument / matrix-patch / RHS text
files only).  `src/` NOT touched.  Branch HEAD `cd9bd7a` unchanged.

### 1.1 Fingerprint 2 — corrected-T3 route vs exported b_TC (relL2)

- Instrument: `b30_fp2_route_rebuild.py` (route rebuild with T3-slot
  coefficient parameterized 0.6 -> c*; T6 boundary slot toggled keep06/corr;
  self-gate = exported b_TC relL2 0.3970 cos 0.9583, bitwise b26a route_V1).
- Results (relL2 total / relL2_Ublock / amp_ratio_U):

  | variant          | relL2      | relL2_Ublk | amp_ratio_U |
  |------------------|-----------:|-----------:|------------:|
  | baseline (0.6)   | 0.38617    | —          | 0.69558      |
  | T3=1.00 keep06   | 0.29615    | 0.29615    | 1.11821      |
  | T3=1.00 corr     | 0.33833    | 0.33833    | 1.14761      |
  | T3=1.09 keep06   | 0.35457    | 0.35457    | 1.21438      |
  | T3=1.23 keep06   | 0.47191    | 0.47191    | 1.36436      |
  | T3=1.23 corr     | 0.52882    | 0.52882    | 1.40755      |

- Fine sweep (`b30_fp2b_fine_sweep.py`, c in [0.6,1.4]): **min relL2 =
  0.25819 at c*=0.85, Uamp = 0.95857** (U-row amplitude gap 0.6956 -> 0.9586,
  i.e. ~1.0 achieved).  g1=TIn/Tmix and T8=on/off identical within sweep
  (T8 slot contributes zero; g1 choice irrelevant to T3 attribution).
- P-block relL2 stays 0.363-0.367 in all variants (T3 patch touches U-row only).
- Verdict: **moves in preregistered direction** (0.397 -> 0.258 min; U gap
  0.696 -> 0.959) but **fails the <5% threshold** (best 25.8%).  PARTIAL.

### 1.2 Fingerprint 3 — corrected b_TC x existing w_true (M1 three-direction ratios)

- Instrument: `b30_fp3_m1_ratios.py`.  RHS pinned from b26a (preregistered):
  D1 0.074436786, D2 0.000638758 (known large-cancellation small residue,
  preregistered non-blaming), D3 0.010940299.
- Production baseline ratios: D1 -0.656, D2 -10.08, D3 +1.057.

  | variant        | D1      | D2      | D3      |
  |----------------|--------:|--------:|--------:|
  | T3=1.00 keep06 | -2.175  | +104.9  | +10.01  |
  | **T3=1.00 corr** | **+1.1265** | **-32.80** | **-1.7630** |
  | T3=1.09 corr   | +1.2749 | -37.76  | -2.206  |
  | T3=1.23 corr   | +1.5059 | -45.48  | -2.894  |

- Only D1 at T3=1.00 T6=corr passes |ratio-1|<=0.15 (+1.1265).  D2 remains
  unstable (non-blaming per preregistration).  D3 rhs=0.01094 is NOT a small
  residue, yet D3 flips from production +1.057 (already PASS) to -1.763 under
  T3 correction: the T3-only correction **harms D3** — a falsification signal
  for the T3-only attribution of fingerprint 3's D3.
- Verdict: **partial movement on D1 only; all-three threshold FAILED.**  PARTIAL.

### 1.3 Fingerprint 1 — coefficient-patched explicitJT_H7 re-solve (c*=1.00, primary)

- Patch design (`b30_fp1_patch_resolve.py`): M_corrected = M + (c*-0.6)*T3_base
  restricted to **internal OFF-DIAGONAL U-row/P-column slots** [U(own),P(nei)]
  and [U(nei),P(own)] ONLY.
  - Diagonal [U(own),P(own)]/[U(nei),P(nei)] positions are NOT pure T3: the
    hA/hB Stage-2 path (solveDiscreteFlowAdjoint.H 2722-2807) writes the same
    (row,col) cells (`-upperF*Vinvo*hBx`, `-lowerF*Vinn*hBx`, hA loop over
    cellFacesOf).  Measured: 'nn' orig ~ -0.6*vv (opposite sign, dev up to
    vmax); 'oo' orig/vv ~ 1.66.  EXCLUDED from patch.
  - Boundary U-row/P-col block is COMPOSITE (hB `Bdiag*Vinvc*(alphaRel*rAUb)*
    patchSf` ratio ~0.3 + direct (1-alphaRel)*patchSf + internal-face
    accumulation), NOT a pure T3 slot.  EXCLUDED from patch.
- Gate (position-level purity mask, |rel|<1e-6): 193,720 significant
  off-diagonal positions; **183,812 (94.9%) exactly 0.6*vv** (median rel
  1.048e-14, frac_exact 0.9489); **9,908 (5.11%) composite** (boundary-adjacent
  cells, maxrel 0.398) EXCLUDED and left at production composite values
  (documented in resolve JSON meta).  Gate PASSED.
- Solve: scipy splu default COLAMD (MMD_AT_PLUS_A pathological on this
  near-singular matrix, b26-proven), t=1339.8s background+poll.
  Closure PASS (|r_P|/|rxdP| = 3.45e-11 / 2.24e-11 / 5.23e-11);
  cos(w_corr,w_true)=0.9981/0.9900/0.9986; relL2(w)=0.675/0.876/0.675
  -> w_corr ~= w_true (T3 off-diag patch is a small perturbation of the full
  operator; expect the b_TC^T signal to move mostly via the corrected b_TC).
- Main criterion `b_TC_corr^T w_corr / RHS`:
  - T6=corr: D1 **-0.502**, D2 40.5, D3 13.5  (all FAIL)
  - T6=keep: D1 -3.73, D2 155, D3 17.5      (all FAIL)
  Production lhs was -0.049/-0.0064/+0.0116; c*=1.0 does NOT close the closure.
- gDP prediction (B13 decomposition, lambda_PD solved against b_PD):
  gdp_factor = (bPD@w_corr)/(bPD@w_true) = D1 **1.0134** (PASS 1.0+-0.15),
  D2 0.1606 (FAIL), D3 0.7519 (FAIL).  Direction-uniform 2.05x factor NOT
  reproduced (only D1).
  - **defcheck caveat (must be reported)**: b_PD^T w_true / ADJ_gDP =
    0.8178/0.3660/1.0362 — b_PD does NOT reproduce ADJ_gDP even with w_true
    (B13 decomposition has an open gap; the gDP factor prediction carries this
    model error).  ADJ_gDP anchors: D1 -5.2303, D2 1.0505, D3 13.6034.
- Verdict: **main criterion 否证 (all three FAIL); gDP prediction PARTIAL
  (D1 only).**

### 1.4 Cross-cutting findings (new this round)

1. Diagonal U-row/P-col T3 positions are contaminated by the hA/hB Stage-2
   path — an offline T3 patch on them is impossible; only off-diagonal slots
   are pure (median rel 1e-14).
2. 5.11% of even the off-diagonal slots are composite at boundary-adjacent
   cells (closed-face hA sum does not cancel there) -> position-level mask
   |rel|<1e-6 is the correct patch gate, not a structural mask.
3. MMD_AT_PLUS_A never finishes on this near-singular matrix; default COLAMD
   splu is the only proven solver (b26 line).
4. D2 M1 ratio is dominated by large cancellation (|C| ~ 13|rhs|),
   preregistered as non-blaming instability.
5. b_PD defcheck gap (0.818/0.366/1.036) means the gDP factor prediction is a
   projection onto a decomposition that itself does not reproduce ADJ_gDP —
   the 2.05x direction-uniform factor is NOT fully explained by the T3 slot
   (B13 already located the bulk of the 2.05x in the momentum-row channel).

---

## Stage 2 — verdict (per preregistered mapping)

Thresholds (all at primary c*=1.00):
- fp2 relL2 < 5%: **NO** (best 25.82% at c*=0.85; 29.6% at c*=1.00).
- fp3 all-three |ratio-1|<=15%: **NO** (only D1=+1.1265 passes at T6=corr;
  D2 non-blaming large-cancellation; D3 flips 1.057 -> -1.763).
- fp1 main criterion |ratio-1|<=15% all three: **NO** (D1 -0.502, D2 40.5,
  D3 13.5).
- gDP within 1.0+-0.15: **PARTIAL** (D1 1.0134 PASS; D2 0.161, D3 0.752 FAIL;
  direction-uniform 2.05x not reproduced).

Movement in preregistered direction: YES on multiple independent channels —
fp2 relL2 0.397->0.258 min and U amplitude 0.696->0.959 (~1.0); fp3 D1
-0.656->+1.126; fp1 gDP D1 factor 2.05->1.013.  This is not 无罪.

Conviction criteria (ALL four with same c*=1.00) NOT met.

### Verdict: 部分 (partial)

The (1-alphaRel)=0.6 direct-relax coefficient in the T3 slot is IMPLICATED in
the U-row amplitude and the D1/gDP direction: correcting 0.6->1.0 moves the
U amplitude to ~1.0, fp3 D1 from -0.66 to +1.13, and the gDP D1 factor from
2.05 to 1.01.  But the single-slot T3 correction does NOT close fingerprints
2/3 (relL2 floor 25.8% vs 5% target; fp1 main closure fails all three; D3
ratio HURT from +1.057 to -1.763).  The unified single-slot hypothesis is
NOT convicted; the T3 slot is one contributor among at least two.

### Open items (residual; NOT blamed on T3 by this data)

1. hA/hB Stage-2 path (solveDiscreteFlowAdjoint.H 2722-2807) contaminates the
   diagonal U-row/P-col slots — its semantics (`Bdiag*Vinvc*(alphaRel*rAUb)*
   patchSf` on the boundary; `-upperF*Vinvo*hBx` / hA cellFacesOf loop
   internally) need a separate adjudication round (candidate: contributes the
   P-block 0.363 floor in fp2).
2. Boundary U-row/P-col block is composite (hB + direct + internal-face
   accumulation) — not a pure T3 slot; excluded from the patch by design.
3. D3 ratio flips sign under T3 correction while production was already
   +1.057 (passing): fingerprint-3's D3 has another source than the T3 slot.
4. D2 M1 ratio unstable (large-cancellation small residue, preregistered
   non-blaming).
5. b_PD defcheck gap (0.818/0.366/1.036): the gDP 2.05x direction-uniform
   factor is not fully explained by T3; B13 locates the bulk in the
   momentum-row channel.
6. MMD_AT_PLUS_A pathological on this matrix (COLAMD-only proven).

### B31 production fix spec

NOT written.  The preregistered contingency for the B31 fix spec is
conviction, which was NOT reached (verdict: partial).  Production `src/`
remains untouched; the evidence above is the hand-off material for the user's
pause-period review.
