# D2 audit — residual-convention reconciliation (B4 J/JT probes vs NS.H primal)

Stage: B-final (DIAGNOSTIC_ONLY, cycle 2)
Step: D2 of plan-002-replan (BFINAL-INDETERMINATE-PROBE-CONVENTION)
Author: Executor (flash), evidence-only; no src/ or case modification.
Date/HEAD: workspace /home/ys/dsH/TO-ANISOTROPIC @ df695c7 (uncommitted B4 edits
preserved); reference repo /home/ys/TO-ANISOTROPIC @ 292c711.
Evidence base: rerun.log (cycle-2, reproduces log.b4export), src/NS.H,
src/solveDiscreteFlowAdjoint.H, and the installed OpenFOAM-7 sources under
/opt/openfoam7 (semantics verified directly, not from memory).

---

## 0. Verdict (one paragraph)

Every B4 block-decomposition probe recomputes its reference with a *different
hand-rolled residual formula*, and none of the formulas equals the actual
frozen primal residual `R = A·U − source` in the way its comment claims. In
particular, the shared reference formula used by T1, T-fvm, T-conv and the
NS.H deltaH probe,

```
V*(M&U) + 2*H_ldu(U)        with H_ldu(U) = M.lduMatrix::H(U)
```

is **not** `A·U − source`; per the OpenFOAM-7 operator&/lduMatrix::H
semantics it evaluates to

```
boundaryDiag·U − offdiag_action(U) − source − boundarySource
                     = R_true − 2*offdiag_action(U)
```

The `+2*H_ldu` term double-counts the off-diagonal action with the wrong sign
(lduMatrix::H returns the *negative* of the off-diagonal action). T-fvm's
`relL2 ≈ 1.0` is the expected, mechanical signature of this reference-formula
error (plus the source/boundarySource terms that are in the reference but not
in J), **not** evidence that J's raw matrix action is wrong. T1 passes at
0.85% because (i) it compares *differences* `R(U1)−R(U0)` so the fixed
source/boundarySource drop out, (ii) the `2*offdiag` error is O(1%) of the
diagonal action in this fine-mesh, small-flux case (|lower|L2=4.2e-3 vs
|diag|L2=1.8), and (iii) the deviatoric term is V-scaled and tiny. T-conv's
42% compares J's SIMPLE-phi `velocityJump` against a reference that perturbs
phi with `fvc::flux(dU)` (a *different* phi convention) and again carries the
`+2*H_ldu` formula error, so it localizes nothing. The former claim that the
gap "localizes to convection" is **RETRACTED** (§7). T1 remains the only
complete full-residual comparison, but even T1's discriminating power on the
small blocks (offdiag, SIMPLE-dphi, deviatoric, boundarySource) is limited by
the diagonal dominance of this case, so (a) J≠R_w, (b) incomplete R_x,
(c) RHS/sign/assembly all stay OPEN after D2.

---

## 1. Authoritative OpenFOAM-7 semantics used below (verified from source)

File: /opt/openfoam7/src/finiteVolume/fvMatrices/fvMatrix/fvMatrix.C
File: /opt/openfoam7/src/OpenFOAM/matrices/lduMatrix/lduMatrix/lduMatrixATmul.C
File: /opt/openfoam7/src/OpenFOAM/matrices/lduMatrix/lduMatrix/lduMatrixTemplates.C
File: /opt/openfoam7/src/finiteVolume/fvMesh/fvMeshLduAddressing.H

1. **Addressing** (fvMeshLduAddressing.H): `lowerAddr()` = subList of
   `mesh.faceOwner()`, `upperAddr()` = `mesh.faceNeighbour()`. I.e.
   `upper[face]` is the coefficient of `psi[neighbour]` in the **owner** row;
   `lower[face]` is the coefficient of `psi[owner]` in the **neighbour** row.
   (Confirmed by lduMatrix::Amul and lduMatrix::H below; also
   lduAddressing.C: `owner = lowerAddr(); neighbour = upperAddr();`)

2. **lduMatrix::Amul** (lduMatrixATmul.C):
   `(A*psi)[owner] += upper[f]*psi[nei]; (A*psi)[nei] += lower[f]*psi[own]`.
   So offdiag_action(c) = Σ_faces coeff·psi_other with the mapping of (1).

3. **lduMatrix::H(psi)** (lduMatrixTemplates.C):
   `H[owner] −= upper[f]*psi[nei]; H[nei] −= lower[f]*psi[own]`
   ⇒ **H_ldu(psi) = − offdiag_action(psi)** (no V, no source, no boundary).

4. **fvMatrix::operator&** (fvMatrix.C ≈L2201):
   per component: `Mphi = −boundaryDiag·psi + H_ldu(psi) + source +
   boundarySource; Mphi /= −V`, with `boundaryDiag = diag + internalCoeffs`
   (addBoundaryDiag). Hence
   `V*(M&psi) = boundaryDiag·psi + offdiag_action(psi) − source −
   boundarySource = R(psi)` — i.e. **operator& returns the residual
   (A·psi − source with boundary terms), per volume**.

5. **fvMatrix::A()** = `D()/V`, `D() = diag + cmptAv(internalCoeffs)`
   (per-volume diagonal, averaged boundary diag).

6. **fvMatrix::H()** (per-volume):
   `H(psi) = [ boundaryDiagH·psi + H_ldu(psi) + source + boundarySource ]/V`
   with `boundaryDiagH = −internalCoeffs_cmpt + cmptAv(internalCoeffs)`
   applied per component to the cell value (all patches, including
   fixedValue). No raw `diag()` enters H().

7. **fvMatrix::residual()** (fvMatrixSolve.C ≈L351):
   `residual = source + boundarySource − boundaryDiag·psi − offdiag(psi)`
   = `− R(psi)`.

8. **fvm::Sp(sp, vf)** (finiteVolume/finiteVolume/fvm/fvmSup.C):
   `fvm.diag() += mesh.V()*sp.field()` — the Brinkman term enters the **raw
   diagonal with a V factor**. Because A() divides by V, `A() = A_base +
   alpha` exactly (V·alpha/V), so `A0 := A() − alpha = A_base` is a *correct*
   non-Brinkman per-volume diagonal; the NS.H comment "no V factor on the
   fvMatrix diagonal" (NS.H L240-242) is misleading in its justification but
   the resulting `A0 = UEqn.A() − alpha` formula is dimensionally right.

9. **fvMatrix::operator-=(field)** (fvMatrix.C ≈L1108):
   `source() += V*field`. Therefore the RANS deviatoric line
   `UEqn −= fvc::div(nuEffFrozen*dev2(T(grad(U))))` **adds**
   `+V·div(nuEffFrozen·dev2(T(grad U)))` to the source, i.e. the explicit
   deviatoric term sits on the RHS with a plus sign, as in the standard
   simpleFoam `divDevReff` split. Its linearized variation is
   `d(source) = +V·div(dev2(grad dU))`.

10. **fvMatrix::relax(α)** (fvMatrix.C ≈L672): scales the whole diagonal
    (`diag' = α·diag`) and adjusts the source so the fixed point is
    unchanged. Hence the **relaxed** `UEqn.A()` used by NS.H is
    `α_relax·A_base + α_relax·alpha (+ cmptAv(ic)/V)`, i.e. ~`α_relax·(A_base
    + alpha)` — a factor `1/α_relax` (≈1.4 for α_relax≈0.7) different from
    the unrelaxed per-volume diagonal.

---

## 2. NS.H primal residual conventions (what the state actually satisfies)

NS.H L26-45 (per SIMPLE corrector):
```
UEqn = fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U)
       == −fvc::grad(p) + fvOptions(U)
if ransFlowModel: UEqn −= fvc::div(nuEffFrozen*dev2(T(grad U)))   // source += V·div(dev2)
UEqn.relax();
solve(UEqn == −fvc::grad(p));
rAU = 1/UEqn.A()          // RELAXED per-volume diagonal (semantics 5,10)
HbyA = constrainHbyA(rAU*UEqn.H(), U, p); phiHbyA = fvc::flux(HbyA); adjustPhi
if simple.consistent(): rAtU = 1/(1/rAU − UEqn.H1()); phiHbyA += ...; HbyA −= ...
primalPressureMobility = rAtU
pEqn: fvm::laplacian(rAtU,p) == fvc::div(phiHbyA); solve; phi = phiHbyA − pEqn.flux()
p.relax(); U = HbyA − rAtU*grad(p)
```
Key points for the audit:
- The converged `(U,p,phi)` satisfies div(phi)=0 to SIMPLE tolerance and the
  momentum equation with the **relaxed** matrix; the residual `R_U = A·U −
  source` of the *unrelaxed, no-grad(p)* matrix is ~0.09 (RxProbe baseline
  `||RU||L2=0.0909`) — not machine zero, because the probe omits grad(p),
  fvOptions and uses an unrelaxed fresh matrix (NS.H L416-463).
- `rAU`/`rAtU` used for the final phi rebuild are **relaxation-inflated**
  (factor ~1/α_relax vs the unrelaxed A()).
- NS.H RxProbe baseline residual convention: `resBase = −residual() = A·U −
  source` for the matrix `div − laplacian + Sp` with `fvOptions.constrain`
  only; **no dev2 term, no −grad(p)** (NS.H L419-427 vs the real UEqn which
  has both) — a *partial* residual.

---

## 3. B4 J operator conventions (solveDiscreteFlowAdjoint.H)

Assembly (L22-99):
```
discretePrimalMomentum = fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U)
                         // NO relax, NO dev2, NO grad(p)
discreteDiagX/Y/Z = diag + internalCoeffs  (per component, ALL patches)  // = boundaryDiag
discreteUpper/Lower = raw upper/lower
rAUAdj = 1/discretePrimalMomentum.A()      // UNRELAXED per-volume diag (semantics 5)
```
J (applyDiscreteFlowJ, L414-711), for `[dU; dp]`:
- **deviatoric** (L420-474): `output[U] += V·(−div(nuEffFrozen·dev2(T(grad dU))))`
  = `−V·div(dev2(dU))`. Per semantics 9 this **matches** `dR_U/dU` (the
  residual's source variation is `+V·div(dev2(dU))`, subtracted: −dsource).
  ✓ sign correct; T-dev = 5.4e-5 confirms.
- **diagonal** (L476-484): `output[U] += discreteDiagX·dU` (= boundaryDiag·dU).
- **off-diagonal** (L597-623): `output[own] += upper·dU_nei; output[nei] +=
  lower·dU_own` — **correct per semantics 1/2** (upper = owner-row/neighbour-col).
- **dH tangent** (L502-582): `dH[own] −= upper·dU_nei; dH[nei] −=
  lower·dU_own`; boundary `dH += (−ic_cmpt + cmptAv(ic))·dU_cell` (all
  patches); `dH /= V`. This matches fvMatrix::H()'s boundary-diag convention
  (semantics 6) for the frozen-coefficient variation of H(U).
- **SIMPLE phi tangent** (L625-641): `dphi_f = Sf·(wf·rAUAdj_o·dH_o +
  (1−wf)·rAUAdj_n·dH_n) + mobF·(dp_nei−dp_own)·dcf·maf` with `mobF =
  wf·primalPressureMobility[own] + (1−wf)·primalPressureMobility[nei]`
  (= rAtU, **relaxed** convention) — **mixes unrelaxed rAUAdj (dH part) with
  relaxed primalPressureMobility (dp part)** (semantics 10). T2=7.5e-21
  verifies the dp part against div(rAtU·snGrad(dp)·magSf) exactly.
- **velocityJump** (L646-653): `output[U_downwind] += dphi·(U_nei−U_own)` —
  the convection non-linearity through phi (upwind lumping).
- **continuity row** (L656-657, 696-697): `output[P_own] += dphi;
  output[P_nei] −= dphi; boundary cells += Sf_b·(rAUAdj·dH_b)` — i.e.
  `J_P = div(dphi_SIMPLE)`.
- **pressure gradient** (L659-672, 699-708): `output[U] ±= Sf·interp(dp)`;
  boundary `+= Sf_b·dp` if !pressureFixed.
- **absent from J's U-row**: the `−boundaryCoeffs·dU_b` term of the true
  residual Jacobian (semantics 4: dR includes −d(boundarySource)). At the
  outlet (zeroGradient, bc≈ic) this would cancel part of the ic in
  discreteDiagX; in this case |ic_outlet| is small (T1=0.85% consistent), but
  it is a genuine unvalidated boundary term (candidate for the B-F2 list).

Convention summary (verified): J's raw momentum action equals the true
frozen-coefficient Jacobian for the diagonal+off-diagonal+deviatoric parts;
the unvalidated links are (i) the SIMPLE dphi tangent (rAU convention mix,
relaxed vs unrelaxed, and the upwind velocityJump lumping), (ii) the
boundarySource/`−bc·dU_b` term, (iii) whether dphi_SIMPLE equals the true
phi(alpha,U) tangent — the exact thing no current probe tests.

---

## 4. Probe-by-probe exact reference formulas and conventions

### 4.1 T1 — momentum U-U block, dp=0 (solveDiscreteFlowAdjoint.H L1783-1952)
Reference formula (momentumResidualB4, L1783-1818):
```
M(U,phi) = fvm::div(phi,UU) − fvm::laplacian(nuEffFrozen,UU) + fvm::Sp(alpha,UU)
if ransFlowModel: M −= fvc::div(nuEffFrozen*dev2(T(grad UU)))     // source += V·div(dev2)
hldu = M.lduMatrix::H(UU)                                          // = −offdiag_action(UU)
RU_b4(UU,phi) = V*(M & UU) + 2*hldu
```
By semantics 4 + 3:
```
V*(M&UU)          = boundaryDiag·UU + offdiag(UU) − source − boundarySource
RU_b4(UU,phi)     = boundaryDiag·UU − offdiag(UU) − source − boundarySource
R_true(UU,phi)    = boundaryDiag·UU + offdiag(UU) − source − boundarySource
⇒ RU_b4 = R_true − 2·offdiag_action(UU)
```
Comparison: `Jdelta = J·[dU; dp=0]` vs `RU_b4(U1,phiper) − RU_b4(U0,phibase)`
with `dU = epsB4·sine(0.037 i)`, epsB4=1e-6 (L1749-1766), `phiper =
phibase + fvc::flux(dU)` (L1820,1838) — **phi is perturbed by the geometric
flux fvc::flux(dU), NOT by the SIMPLE tangent dphi_SIMPLE that J's
velocityJump term uses**.
Why it is clean (relL2=0.008523, cos=0.99996, |J|/|FD|=0.99996): the
difference comparison cancels source/boundarySource; the residual formula
error is `2·[offdiag(U1)−offdiag(U0)] ≈ 2·offdiag(dU)` which is O(1%) of the
boundaryDiag·dU part because the matrix is strongly diagonal-dominant in this
case (matrix norms L1251: |upper|=9.03e-5, |lower|=4.17e-3, |diag|=1.80;
per-cell diag ≈ 0.008 vs per-face offdiag ≈ 1.3e-5); the deviatoric term is
V-scaled and negligible for this direction; and the SIMPLE-dphi vs
flux(dU) difference on the velocityJump is small here (T-cont=6.5e-4 shows
div(dphi_SIMPLE) ≈ div(flux(dU)) to 0.065% for this direction).
Caveats: (i) reference formula carries the `2*H_ldu` sign error — masked by
diagonal dominance; (ii) phi convention on the reference side is flux(dU),
not SIMPLE; (iii) J's `−bc·dU_b` boundary term is absent on both sides;
(iv) with maxCellAbsFD=1.25e-8 the FD signal is small, so the test is a
*diagonal-part* check with weak sensitivity to the small blocks.

Note on the log line "R0 components: ... |2H|=17479873079": that printed norm
is `2·|fvMatrix::H()|` (the *per-volume* H() of fvMatrix, semantics 6, which
includes boundaryDiag and is inflated by the `/V` with cell volumes
V≈1e-9 m³ in this case) — it is **not** the magnitude of the `2*hldu` term
actually used inside momentumResidualB4 (raw lduMatrix::H, no /V, ~O(1e-2)
per cell). It is a display artifact, not a residual magnitude, and does not
change the T1 analysis.

### 4.2 T-fvm — frozen-phi matrix part (L2041-2093)
```
Jfvm = discreteDiagX·dU + upper·dU_nei + lower·dU_own          // raw action = true A·dU
hlduFvm = discretePrimalMomentum.lduMatrix::H(dU)              // = −offdiag(dU)
Afd    = V*(discretePrimalMomentum & dU) + 2*hlduFvm           // = boundaryDiag·dU − offdiag(dU) − source − boundarySource
relFvm = |Jfvm − Afd| / |Afd| = 0.999999999866
```
Difference by construction:
```
Jfvm − Afd = 2·offdiag_action(dU) + source + boundarySource
```
with `boundarySource = boundaryCoeffs·dU_b` ≈ `ic·dU_cell` at the zeroGradient
outlet (bc≈ic there), which is O(diagonal) at boundary cells — so the
mismatch is O(1) and `relFvm ≈ 1.0` is the **expected signature of the
reference-formula convention error**, not a defect of Jfvm. In fact Jfvm IS
the true frozen-coefficient Jacobian action `A·dU` (semantics 4), while Afd
is the same `+2*H_ldu` mis-formula applied to dU. Neither side involves the
SIMPLE dphi (both are frozen-phi), so T-fvm says nothing about the phi
tangent. **T-fvm ≈ 1.0 does NOT localize a broken block.**

### 4.3 T-dev — deviatoric part (L2095-2164)
`Jdev = V·(−div(dev2(grad dU)))` vs reference `−V·dRdev` with `dRdev =
div(dev2(grad U1)) − div(dev2(grad U0))` ⇒ reference = `−V·div(dev2(dU))`.
Per semantics 9 this equals the true `dR_U/dU` deviatoric contribution.
relL2=5.36e-5 ✓ — the deviatoric block is consistent.

### 4.4 T-conv — convection phi-variation part (L2166-2205)
```
jConvCell = Jdelta − Jfvm − Jdev            // = velocityJump(dphi_SIMPLE) (dp=0)
convDelta = resWithHldu(M(phiper),U0) − resWithHldu(M(phibase),U0)   // M = fvm::div only
relConv = |jConvCell − convDelta| / |convDelta| = 0.418
```
Reference uses `phiper = phibase + fvc::flux(dU)` (geometric flux) while J's
term uses `dphi_SIMPLE`; the reference again uses the `resWithHldu =
V*(M&U)+2*hldu` formula (same sign error as §4.1); and `convDelta` is the
phi-variation of the *convection-only* residual at fixed U0, whereas J's
velocityJump acts at the upwind cell with `U_nei−U_own`. The 42% is the
combined signature of (a) different phi conventions (SIMPLE vs flux), (b) the
`2*H_ldu` formula error, (c) the upwind lumping — **it cannot discriminate
any single broken block** and is not evidence that the convection block is
wrong (T1, which includes the velocityJump term against a complete residual
difference, is clean at 0.85%).

### 4.5 T-cont — continuity row (L1954-1969)
`J_P = div(dphi_SIMPLE)` vs `V·div(fvc::flux(dU))`. relL2=6.49e-4. Two
*different* phi conventions nearly agree for this direction; consistent with
T1 passing. It verifies that the J_PU row is close to div of the geometric
flux for this direction — it does **not** verify dphi_SIMPLE against the true
SIMPLE-rebuilt phi (which is what B-F2/B-F3 need).

### 4.6 T2 — pressure coupling (L1971-2039)
`J_PP = div(mobF·(dp_nei−dp_own)·dcf·maf)` vs `div(rAtU·snGrad(dp)·magSf)`.
relL2=7.5e-21 — exact; the dphi/dp SIMPLE tangent with mobF =
primalPressureMobility (= rAtU, relaxed convention) is verified.

### 4.7 GateH1 — same-basis dH (L1082-1173)
Analytic: `dH = (−upper·dU_nei [own], −lower·dU_own [nei]) + (−ic·dU_cell,
non-fixedValue patches only)`, then `/V`; FD: `[H(U+hG·dUG) − H(U−hG·dUG)]/(2
hG)` of fvMatrix::H() (per-volume, boundaryDiag = −ic_cmpt + cmptAv(ic) on
ALL patches + H_ldu + source + boundarySource, semantics 6).
relL2=0.117, cos=0.99943. The 11.7% is the boundary-convention residual of
the *analytic*: it omits the fixedValue-patch ic contribution, the
`+cmptAv(ic)` part, and the `boundaryCoeffs·dU_b` (extrapolated boundary
value) variation, all of which fvMatrix::H() includes. Note J's actual dH
(L554-574) uses `(−ic_cmpt + cmptAv(ic))` for ALL patches — closer to
fvMatrix::H() than the GateH1 analytic is; GateH1 is a cruder check and its
11.7% is a reference-side boundary-convention mismatch, not a J defect.

### 4.8 NS.H RxProbe deltaH (NS.H L729-778) — "frozen-phi dH relL2=4.67e7"
Analytic (raw, per-cell): `dHx0[own] −= upper·dU_nei; dHx0[nei] −=
lower·dU_own; dHx0 += (−alpha·V)·dU_cell` — **no `/V`**, and a spurious
`−alpha·V·dU` term (alpha sits on the diagonal; it enters fvMatrix::H() only
through boundaryDiag at boundary cells, and with the *opposite* sign of the
interior term the analytic adds). FD side is the per-volume fvMatrix::H()
(`HPp/HPm = tUEqnRb.ref().H()`, L688/695) whose magnitude is inflated by the
`/V` (V≈1e-9 m³ in this case). relL2=46673795.77 (4.7e7) is a pure
per-volume / spurious-alpha convention mismatch — **no physics content**; it
cannot localize anything.

### 4.9 Rx-B (NS.H L192-410) — see D3 audit for the full analysis; D2 note:
`rAUmod = 1/(A0 + alphaMod)` with `A0 = UEqn.A() − alpha` where UEqn is the
**relaxed** matrix (semantics 10): in solid cells `alpha≈1e8`, so
`A0 ≈ α_relax·A_base + (α_relax−1)·alpha ≈ −3e7` and `A0+alphaMod ≈ 7e7`
while the probe's perturbation is `±hAlpha·dA = ±100` — a relative rAU change
of ~1e-6, hence `|dRp/dalpha| ≈ 8e-10 ≈ 0` **by construction**. The relaxed-
diag offset (not a physical R_P,x=0) explains the ≈0 result; combined with
the sine (non-design) direction and the self-referential phi rebuild (same
analytic SIMPLE formula being tested), Rx-B cannot settle hypothesis (b).

---

## 5. Convention mismatch table (the D2 core)

| Probe | J side | Reference side | Reference formula error vs R_true | Verdict |
|---|---|---|---|---|
| T1 (0.0085) | full J momentum row | RU_b4(U1,φ+flux(dU))−RU_b4(U0,φ) | `−2·offdiag` (small, ~1%), φ convention = flux(dU) not SIMPLE | clean; complete over all 4 operators; weak on small blocks |
| T-fvm (0.9999…) | raw diag+upper/lower action (= true A·dU) | V*(M&dU)+2·H_ldu(dU) | `+2·offdiag` sign error + source + boundarySource (O(diagonal) at outlet) | ≈1.0 is the expected reference-formula signature; NOT a J defect; no SIMPLE dphi on either side |
| T-dev (5.4e-5) | V·(−div(dev2 dU)) | −V·(div(dev2(U1))−div(dev2(U0))) | none (semantics 9) | consistent |
| T-conv (0.418) | velocityJump(dphi_SIMPLE) | [V*(M(φ+flux(dU))&U0)+2H]−[V*(M(φ)&U0)+2H] | φ convention mismatch (flux vs SIMPLE) + `2*H_ldu` + upwind lumping | non-localizing |
| T-cont (6.5e-4) | div(dphi_SIMPLE) | V·div(flux(dU)) | φ convention differs, small for this direction | verifies J_PU≈div(flux(dU)) for this direction only |
| T2 (7.5e-21) | dphi/dp = mobF·(dp_nei−dp_own)·dcf·maf | div(rAtU·snGrad(dp)·magSf) | none | exact |
| GateH1 (0.117) | (analytic) offdiag + −ic(non-fixedValue), /V | fvMatrix::H() FD | missing fixedValue ic, cmptAv(ic), bc·dU_b | boundary-convention residual of the analytic |
| NS.H deltaH (4.7e7) | (analytic) raw, −alpha·V·dU, no /V | fvMatrix::H() FD (per-volume) | no /V + spurious alpha term | units mismatch, no content |
| Rx-B (8e-10) | (analytic) SIMPLE phi rebuild, rAU=1/(A0+α) | same formula, perturbed α | relaxed-diag offset ~(α_relax−1)·α swamps ±100 perturbation | ≈0 by construction; cannot close (b) |

Net: no single probe localizes (a), (b), or (c). The only complete
full-residual comparison is T1, and it is clean — which constrains (a) to
"J's dominant (diagonal) momentum action matches the residual difference for
this direction", leaving the small-block links (SIMPLE dphi/rAU convention,
boundarySource, completeness of R_x, final assembly) unvalidated.

---

## 6. Why T1 is (still) the only complete full-residual comparison

- It is the only probe whose reference evaluates all four operators of the
  frozen momentum residual (div, laplacian, Sp, dev2) through the actual
  fvMatrix machinery on the perturbed state (L1787-1817), and it uses the
  same J operator that the B4 path assembles.
- Its limitations (documented in §4.1): the `2*H_ldu` sign quirk of the
  reference formula (masked by diagonal dominance), the φ perturbation
  convention (flux(dU) rather than SIMPLE), the absence of the
  `−bc·dU_b` boundary term on both sides, and the weak signal of the small
  blocks.
- T-fvm/T-conv/GateH1/NS.H-deltaH references are each *different*
  hand-rolled formulas that are not equal to the true residual Jacobian, so
  their large numbers (1.0 / 0.42 / 0.117 / 4.7e7) are mutually
  contradictory only because their conventions are mutually inconsistent —
  the exact state of affairs the primary hypothesis BFINAL-
  INDETERMINATE-PROBE-CONVENTION asserts.

---

## 7. RETRACTION

The prior claim that the B4 pressure-gradient gap "localizes to the
convection block" (based on T-conv=42%) is **RETRACTED**. T-conv's reference
perturbs phi with `fvc::flux(dU)` while J's compared term uses the SIMPLE
tangent `dphi_SIMPLE`, and both carry the `V*(M&U)+2*H_ldu(U)` formula error;
the 42% is a convention-mismatch artifact that cannot attribute the gap to
convection. All three candidates (a) J≠R_w, (b) incomplete R_x, (c)
RHS/sign/final-assembly remain OPEN after D2.

---

## 8. Implications for (a)/(b)/(c) after D2

- (a) J≠R_w: constrained but not excluded. The dominant diagonal/offdiag/
  deviatoric parts match the residual difference (T1 0.85%, T-dev 5.4e-5,
  T-cont 6.5e-4); the unvalidated J links are the SIMPLE dphi tangent
  (rAUAdj-unrelaxed vs primalPressureMobility-relaxed mix; upwind
  velocityJump) and the missing `−bc·dU_b` boundary term. A single-convention
  block-wise B-F2 (J·v vs central FD of the *actual* frozen residual with
  phi rebuilt the SIMPLE way) is the minimal discriminator.
- (b) R_x incomplete: still the leading suspect, and still unproven: Rx-B is
  self-referential (same analytic SIMPLE formula), sine-directional, and its
  ≈0 is explained by the relaxed-diag offset (§4.9); first-principles
  `d(rAU)/dalpha = −rAU² ≠ 0` at fixed (U,p) implies `dphi/dalpha ≠ 0` and
  hence `R_P,x·d ≠ 0` generically. Only a design-direction, fixed-state,
  split R_U,x / R_P,x FD (B-F3) can close it.
- (c) RHS/sign/assembly: g_w=7.08e-12 and the `−λᵀR_x` convention are
  self-consistent, but the end-to-end gsenshPressureDrop chain
  (sensitivity.H L65-68 = −dAlphaDxh·(U&Uc)·V; AdjNS_PD.H L3-4 =
  pressureConstraintDerivative) has no independent FD; it is tested only
  after (a) and (b) are closed (B-F4/B-F5 of handoff §8).

The minimal missing diagnostics are therefore exactly handoff §8 B-F2
(single-convention J·v vs central FD of the true frozen residual, momentum
and continuity rows separately, non-neutral directions, epsilon sweep) and
B-F3 (fixed-state R_x·d vs design-direction central FD, split R_U,x and
R_P,x, real design direction). Their exact formulas, ordering, index
conventions (U = 3·cell+cmpt, P = 3N+cell), epsilon sweep and block-wise
thresholds are specified in gate_spec.md (D4).

---

## 9. Files/lines consulted (read-only)

- src/NS.H: L26-45 (UEqn+dev2+relax), L113-129 (rAU/HbyA/phiHbyA/rAtU/
  primalPressureMobility), L136-190 (Rx-A), L192-410 (Rx-B, A0=UEqn.A()−alpha
  at L243-244), L416-463 (baseline resBase=−residual), L465-577 (frozen-phi
  U-column FD export), L579-778 (rebuilt-phi P-row FD + deltaH checks,
  spurious −alpha·V at L750), L1016-1041 (pEqn, phi=phiHbyA−pEqn.flux()).
- src/solveDiscreteFlowAdjoint.H: L22-99 (discretePrimalMomentum,
  discreteDiagX=diag+internalCoeffs, discreteUpper/Lower, rAUAdj), L414-711
  (applyDiscreteFlowJ: dev L420-474, diag L476-484, dH offdiag L515-531,
  dH boundaryDiag L554-574, dH/V L576-582, matrix action L597-623, SIMPLE
  dphi L625-641, velocityJump L646-653, continuity L656-657, pressure
  L659-672, boundary L675-710), L1082-1173 (GateH1), L1749-1766 (dUField,
  epsB4=1e-6), L1783-1818 (momentumResidualB4), L1820-1840 (dphiField,
  phiper, RU0/RU1), L1949-1952 (T1), L1954-1969 (T-cont), L1971-2039 (T2),
  L2041-2093 (T-fvm), L2095-2164 (T-dev), L2166-2205 (T-conv).
- /opt/openfoam7 sources: fvMatrix.C (operator& ≈L2201; addBoundaryDiag
  L110; A() L738; H() L760; relax L672; operator-= L1108; residual in
  fvMatrixSolve.C L351), lduMatrixTemplates.C (H), lduMatrixATmul.C (Amul),
  fvMeshLduAddressing.H (lowerAddr=faceOwner, upperAddr=faceNeighbour),
  fvmSup.C (fvm::Sp: diag += V·sp), fvmDiv/gaussConvectionScheme.C.
- Evidence: evidence/agent-group/BFINAL-001/cycle-2/rerun.log L144-164
  (Rx-A/Rx-B/baseline/deltaH), L1205 (GateH1), L1250-1262 (R0 components,
  matrix norms, T1/T-cont/T2/T-fvm/T-dev/T-conv).
