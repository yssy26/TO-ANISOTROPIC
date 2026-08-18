# BFINAL-004-FIX step2/step3 record — probe fix + header-only rebuild

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step2-fix-probe + step3-build
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @ `ca8a772b8339a36e7906ea32a7125061c47aa280`

## step2 — src/stageB6RxDesignOracle.H only (untracked diagnostic header)

Pre-state baseline (cycle-1/pre-state, captured before any edit): probe sha256
`f0b88ba9c0c885d01e140ff3b958cdcc990a4792947f08eaa36e234481f35fe1`.
Post-edit sha256: `51a02fc6288da612fc7e9b775cfd05ac757c137b492d5e863185f53c74c4b31d`.

Three fixes, exactly per the approved plan:

- **B1 (RX-A weighted comparison):** added `anRUaW[3*celli+c]=anRUa[3*celli+c]*deltaAlpha[celli]`
  and `anRPaW[celli]=anRPa[celli]*deltaAlpha[celli]`; the RX-A metrics now compare
  `stageB6Metrics(anRUaW, fdRUa, ...)` / `stageB6Metrics(anRPaW, fdRPa, ...)` (the FD is the
  deltaAlpha-weighted directional derivative `[R(alpha+eps*deltaAlpha)-R(alpha-eps*deltaAlpha)]/(2eps)`).
  Unweighted `anRUa`/`anRPa` retained for RX-B/RX-C chain reuse.
- **B2 (filter-tangent + gsens solves crash):** `stageB6_d`, `stageB6_y` (filterTangent),
  `gsenshVolR`/`gsensVolR` (§9) and `gsenshR`/`gsensR` (§10) are now constructed with
  `zeroGradientFvPatchScalarField::typeName` (mirroring production `xp`/`gsensVol`/`gsenshVol`
  in createFields.H and `dProjectionDeta` in filter_chainrule.H); `correctBoundaryConditions()`
  is called before each solve (`yfield`, `gsensVolR`, `gsensR`; `dfield`/`gsenshVolR`/`gsenshR`
  already had it). This removes the `calculatedFvPatchField::gradientInternalCoeffs` SIGABRT.
- **B3 (analytic R_P chain, no-V):** replaced the wrong
  `drAU=-alphaRel*rAUc*rAUc` and `dHbyA=rAUc*((alphaRel-1.0)*U - alphaRel*HbyABase)` with the
  derived no-V forms `drAU=-rAUc*rAUc/alphaRel` and
  `dHbyA=(rAUc/alphaRel)*((1.0-alphaRel)*U[celli] - HbyABase[celli])`; header comment (L34-36)
  and body comment (L380-382) updated to the no-V formulas. `div(dphiHbyA)`, `-div(dflux)`,
  `total` per-term decomposition export retained.
- **OFstream durability:** RX-A (`rxaFD`/`rxaAn`/`rpaTerms`), RX-B (`rxbFD`/`rxbAn`),
  RX-C (`rxcFD`/`rxcAn`/`rxcZ`/`rxcXhFD`/`rxcXhAn`/`rxcAlphFD`), §9 (`vOut`), §10
  (`ghOut`/`gOut`/`lamOut`/`daOut`/`dxOut`/`dadxOut`/`dirOut`/`ctOut`/`prOut`) are scoped and
  flushed after each write batch (`.flush()`), closing on scope exit — no tail truncation even
  if a later gate aborts.

## Verification of change scope

- `git diff --stat`: only the two pre-existing BFINAL-003 production files
  (`solveDiscreteFlowAdjoint.H`, `solveDiscreteFlowAdjointProduction.H`), which are FORBIDDEN
  and untouched this round.
- `git diff` vs `cycle-1/pre-state/git_diff_unstaged.patch`: **byte-identical** (verified with
  `diff`), i.e. the unstaged production diff was not altered.
- sha256 of `solveDiscreteFlowAdjoint.H`, `solveDiscreteFlowAdjointProduction.H`,
  `stageB5BoundaryRelaxOracle.H` unchanged vs pre-state baseline.
- Old wrong formulas cleared: `grep` for `alphaRel*rAUc*rAUc`, `(alphaRel - 1.0)*U`,
  `mesh.V()[celli]*rAUc`, `-alphaRel*rAU^2`, `rAU*((alphaRel-1)` → empty (rc=1).
- Brace balance: 148 open / 148 close.

## step3 — header-only incremental wmake

Command (clean PATH, source bashrc, OpenFOAM bin on PATH, FOAM_USER_APPBIN exported AFTER
source, unset FOAM_SIGFPE, cd src, wmake):

```
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc >/dev/null 2>&1
export PATH="/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin:$PATH"
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd /home/ys/dsH/TO-ANISOTROPIC/src && wmake
```

- `WMAKE_EXIT=0` (log: `evidence/agent-group/BFINAL-004-FIX/cycle-1/build_step2.log`).
- Only `MTO_HF.C` recompiled (header-only change; Make/linux64GccDPInt32Opt/MTO_HF.o updated),
  then relinked.
- New binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (4696344 bytes, mtime
  2026-08-17 00:44:47; old cycle-1 binary was 4691928 bytes at 22:42).
- Confirmed the run invokes the ABSOLUTE path (stale `/home/ys/OpenFOAM/ys-7/.../MTO_HF`
  has 0 stageB6 symbols and is not on PATH in the run script).

## Independent source verification of the B3 no-V derivation (this executor)

- `fvMatrix.C` `A()`: `tAphi.primitiveFieldRef() = D()/psi_.mesh().V();` → `A()=D/V`.
- `fvMatrix.C` `H()`: `Hphi.primitiveFieldRef() /= psi_.mesh().V();` → `H()=H0/V` (per-volume).
- `fvmSup.C` `Sp`: `fvm.diag() += mesh.V()*sp.field();` → `D0 = D_rest + alpha*V`.
- `fvMatrix.C` `relax()`: `D /= alpha;` then `S += (D - D0)*psi_.primitiveField();`
  → `D_rel = D0/alphaRel`, `S_rel = S0 + (D_rel-D0)*U`.
- Hence `rAU = 1/A() = V/D_rel` (carries V), `HbyA = H0/D_rel + (1-alphaRel)*U`,
  `drAU/dalpha = -rAU^2/alphaRel`, `dHbyA/dalpha = (rAU/alphaRel)*((1-alphaRel)*U - HbyA)`
  — no extra `mesh.V()`. v1's explicit-V form (REJECTED_PHYSICS) is not in the code.
