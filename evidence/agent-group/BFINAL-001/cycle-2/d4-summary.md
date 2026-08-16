# D4 summary — gap analysis + next-Gate spec (B-F2/B-F3)

Stage: B-final (DIAGNOSTIC_ONLY, cycle 2)
Step: D4 of plan-002-replan (BFINAL-INDETERMINATE-PROBE-CONVENTION)
Author: Executor (flash), evidence-only; no src/ or case modification.
Date/HEAD: workspace /home/ys/dsH/TO-ANISOTROPIC @ df695c7 (uncommitted B4
edits preserved). Synthesizes D2/D3 audits, rerun.log, and NEW read-only
Python cross-checks of the exported artifacts.

## Deliverable
`gate_spec.md` in evidence/agent-group/BFINAL-001/cycle-2/ (this directory).

## Acceptance check (plan D4)
- [x] names NO single survivor: (a) J!=R_w, (b) incomplete R_x,
      (c) RHS/sign/assembly all OPEN (§0, §1 of gate_spec.md).
- [x] exact B-F2 formulas: R_U = A(phi)·U - source (unrelaxed matrix,
      dev2 in source) and R_P = div(phi_SIMPLE) rebuilt through the EXACT
      relaxed primal path (UEqn.relax, rAU=1/A(), constrainHbyA, adjustPhi,
      pEqn.flux(), boundary REBUILT not pinned); J*v from applyDiscreteFlowJ;
      central FD, momentum + continuity rows separately (§3.1-3.2).
- [x] ordering: U index = 3*cell+cmpt [0..3N-1], P index = 3N+cell
      [3N..4N-1], face dphi on internal faces 0..96859 with owner/neighbour,
      boundary ranges listed (§3.3).
- [x] epsilon sweep: eps in {1e-3,3e-4,1e-4,3e-5,1e-5} (matches GatePR
      sweep), plateau rule (§3.5); B-F3 eps in {1e-6,3e-7,1e-7}·alphaMax
      (§4.2).
- [x] block-wise thresholds: momentum <1e-3 AND continuity <1e-3 for
      J=R_w (B-F2); |R_P,x·d|/|R_U,x·d| < 1e-3 before declaring R_x
      complete (B-F3) (§3.5, §4.2).
- [x] optional read-only Python B-F2 cross-check performed and caveats
      noted (§2, §6).

## New evidence from the D4 cross-checks (read-only, exported artifacts)
1. Momentum row: J·v vs NS.H frozen-phi residual FD (stageB2_ruDudU) →
   relL2 = 1.653e-4, cos = 0.99999999 (J side from explicitJT.mtx scatter
   matvec; equals matrix-free J to 1.2e-15). Strongest momentum verification
   yet; momentum row is NOT the locus of the 52x/2.3x/sign-flip.
2. Continuity row: J·v vs NS.H rebuilt-phi FD (stageB2_rpDudU_rebuilt) →
   0.824 overall but 100% boundary-concentrated (internal 1.09e-7/cos=1.0,
   boundary 1.41/cos=0.246) — the 0.824 is the NS.H reference's base-pinned
   boundary, not a J defect; the in-code GatePR (rebuilt boundary) passes at
   1e-10..4e-9 over the sweep.
3. NEW quantified candidate mechanism for (a): J's SIMPLE dphi dH-part uses
   UNRELAXED rAUAdj = 1/discretePrimalMomentum.A() (avg 7.0584e-07) while
   the actual primal phi is rebuilt with RELAXED rAU = 1/UEqn.A() after
   UEqn.relax() (avg 2.8234e-07) — exact factor 2.5 = 1/U-relax (D /= alpha,
   alpha=0.4, verified from /opt/openfoam7 source and the run logs). Every
   existing phi-rebuild FD (GateOracle/GatePR/NS.H rebuilt) uses unrelaxed
   fresh matrices, so they agree with J by construction and cannot see the
   factor. B-F2 with the relaxed rebuild settles it.
4. Ruled out/anchored: g_w 2.7e-11/7.08e-12, Rx-A 1.6e-14, momentum 1.65e-4,
   GatePR 1e-10, ExplicitJToracle 3.9e-16.

## Files written (all under allowed evidence path)
- gate_spec.md (main deliverable)
- crosscheck_bf2.py / crosscheck_bf2.out (read-only J*v vs FD cross-check)
- localize_prow.py / localize_prow.out (boundary/internal P-row localization)

No src/ or case-input file touched; git status after the step identical to
pre-state (verified).
