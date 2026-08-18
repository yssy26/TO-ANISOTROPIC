# BFINAL-005 cycle-2 — step2_patch evidence (executor)

- Task: BFINAL-005 PATCH_RX_PRESSURE_ROW (mode=PATCH, maxReplans=0)
- Step: step2_patch (create src/rxPressureRowTranspose.H + edit src/sensitivity.H)
- Workspace: /home/ys/dsH/TO-ANISOTROPIC == git root (verified); HEAD ca8a772b8339a36e7906ea32a7125061c47aa280; branch agent/dsH-stage-b-validation.

## Files changed (this step)

| file | change |
|---|---|
| `src/rxPressureRowTranspose.H` | NEW (untracked). Computes rxPressureRowT = J_P^T pc: the EXACT transpose of the BFINAL-004-validated forward pressure/continuity-row operator R_P,x = J_P*(dAlphaDxh*z), with the locked reduced-SIMPLE semantics (task §7): read-only rebuild of UEqn (div - laplacian(nuEffFrozen) + Sp(alpha) == -grad(p) + fvOptions, minus dev2 term if ransFlowModel; relax; fvOptions.constrain), alphaRel = relaxEquation("U") ? equationRelaxationFactor("U") : 1.0, rAU = 1/A(), HbyA = constrainHbyA(rAU*H(),U,p); drAU = -rAU^2/alphaRel, dHbyA = (rAU/alphaRel)*((1-alphaRel)*U - HbyA) [BFINAL-004 locked, NO extra V]. Transpose: T1 HbyA face loop (own/nei + assignable-U boundary pinning), T2 flux part via g0 = flux(fvm::laplacian(one,p)) (scheme-exact face gradient; fvm::laplacian and fvMatrix::flux() are linear in the coefficient field, verified in OpenFOAM-7 source), T2 boundary = 0 (forward drAU*w coefficient field carries calculated BC with boundary 0 -> interpolated gamma_b = 0 -> dflux_b = 0; runtime dot test arbitrates). Runtime dot-test self-check: sum(T*w) vs sum(pc*(J_P*w)) with deterministic w = sin(0.271*(celli+1)) on active unclipped cells (BFINAL-004 deltaAlpha support pattern), forward J_P*w re-assembled locally as an exact stageB6 assembleWeightedRPa replica; FatalError if relErr > 1e-8 (expected ~1e-14). No-op (T = 0) when optProperties "stageB6RxDesignOracle" is false -> production byte-identical to pre-BFINAL-005 momentum-only form. Include-safe: no project-header includes, no global state, no forbidden-file references. |
| `src/sensitivity.H` | DISCRETE branch only (+ one guarded export after the chain). Momentum term UNCHANGED (now named gsenshPressureDropMomentum = -dAlphaDxh*(U&Uc)*V). New gsenshPressureDropPressureRow = -rxPressureRowT*dAlphaDxh (NO extra V; dims match gsenshPressureDrop [0 2 -3 0 0 0 0]). gsenshPressureDrop = momentum + pressureRow (single field still feeds the existing filter_chainrule.H path). freezeColdFlowForValidation zeroes BOTH components + total. Guarded (stageB6RxDesignOracle) exports: rxpr_prod_gsensh_momentum.mtx / rxpr_prod_gsensh_pressurerow.mtx / rxpr_prod_gsensh_total.mtx (xh-level, pre-chain) and rxpr_prod_gsens.mtx (raw-design, post-chain) for gates G1-G3. Continuous branch byte-identical (not in diff). |

No other production file touched. `git diff --check` CLEAN. Locked BFINAL-003 heads byte-identical to cycle-2/pre-state baseline:
- src/solveDiscreteFlowAdjoint.H f0c814975cf8bc7f98b9f69efa2f0bc01bd251579920bbedda7d266b72927507
- src/solveDiscreteFlowAdjointProduction.H 8c4901bf8da4edbc7cc1efe73a7583963e4688bba93c8c58539c4c39cae1ad91

sensitivity.H sha256 (post-patch): d777191b472360e4710c6dc78cf3e3c87e3737daf02ef2c040b3e81a7e457efd
(pre-state baseline 1a3031d90a11d780ef163cd70c29e7b356a9c351045b9affb8ca2a82580568c5 — clean at HEAD).

## Math pinned by the runtime dot test (in-code)

For every per-cell weight w:
  sum(T*w) = sum(pc*(J_P*w))  with J_P*w = div(dphiHbyA_w - dflux_w), dphiHbyA_w[f] = Sf&(wf*(dHbyA*w)_own + (1-wf)*(dHbyA*w)_nei), dflux_w = flux(fvm::laplacian(drAU*w,p)).
Then the production contraction is:
  gsenshPressureDropPressureRow = -T*dAlphaDxh
  => sum(gsenshPressureDropPressureRow * z) = -sum(T*(dAlphaDxh*z)) = -pc^T J_P*(dAlphaDxh*z) = -pc^T R_P,x d = D_pressure
with z = dxh/dx*d through the EXISTING filter_chainrule.H (linear chain), and
  sum(gsenshPressureDrop * z) = D_momentum + D_pressure = D_total (G2 target).

NOT implemented (per plan; cycle-1's rejected form): forward operator + pointwise -pc multiply, which does NOT contract to -pc^T R_P,x d (J_P non-diagonal in alpha).

## Compile check (done in step2 to catch implementation defects early)

Exact cycle-3 env: clean PATH -> source /opt/openfoam7/etc/bashrc -> PATH prepend /home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin -> FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin -> unset FOAM_SIGFPE -> cd src && wmake (no set -e).
Result: WMAKE_EXIT=0, binary /home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF refreshed (Aug 17 10:16). Only pre-existing unused-variable warnings (NS.H/createFields.H/solveDiscreteFlowAdjoint.H); no new warnings/errors from rxPressureRowTranspose.H or the sensitivity.H additions.
Log: evidence/agent-group/BFINAL-005/cycle-2/build_step2_check.log

## Acceptance checklist (step2)

- [x] git diff shows ONLY src/sensitivity.H (+ new untracked src/rxPressureRowTranspose.H) changed this step (tracked diff: sensitivity.H only; locked heads unchanged).
- [x] Discrete branch only; continuous branch byte-identical (not in diff).
- [x] Momentum term unchanged (same formula -dAlphaDxh*(U&Uc)*V, same values, renamed for transparency).
- [x] Dot-test present (runtime FatalError if |sum(T*w)-sum(pc*(J_P*w))|/scale > 1e-8).
- [x] No #include of forbidden headers (helper includes only <sstream>).
- [x] git diff --check clean.
- [x] wmake EXIT=0 (compile check; official build+run is step3).

## Gate-relevant note for the reviewer

- G1: D_pressure_production = sum_active(gsensPressureDrop_patched*d) - sum_active(gsensPressureDrop_momentumOnly*d); the momentum-only chain output is stageB6's stageB6_prod_gsens.mtx (BFINAL-004 §10 production replica, momentum UNCHANGED by this patch, chain linear). Compare vs oracle D_pressure = -pc^T R_P,x d.
- G2: D_total_production = sum_active(gsensPressureDrop_patched*d) from rxpr_prod_gsens.mtx; oracle D_total = D_momentum + D_pressure (BFINAL-004 anchors D1 -0.0562356333787 / D2 +0.297672431753 / D3 -0.0363414298395; D_pressure -0.00754887795106 / +0.0685927693909 / -0.00621002821874).
- The dot-test line in the run log pins the transpose before any gate value is read.
