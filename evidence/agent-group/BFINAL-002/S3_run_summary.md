# BFINAL-002 S3 — case run summary (Executor's own clean run)

Date: 2026-08-16, 17:12:18–17:20:36 CST. Binary:
/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF (fresh rebuild 17:12:13,
see S3_build_record.md). Case: /home/ys/dsH/b2_case_smoke with
`stageB4JacobianProbe true; stageB5BoundaryRelaxOracle true;`
(optProperties snapshot: case_optProperties.s3run).

## Deliverables (this S3 execution)

- Raw log: `stageB5_run_s3.log` (1743 lines; oracle block runs once per
  inclusion site — AdjNS_HT.H and AdjNS_PD.H — both blocks byte-identical:
  deterministic, no state drift; same behavior as S2 on2 run).
- Exported artifacts: `artifacts_s3/stageB5_*.mtx` (8 files, md5-identical to
  the case-dir exports of this run).
- Independent recomputation from raw artifacts: `recompute_s3.out`
  (recompute_stageB5.py; reproduces every log number exactly).
- Build record: `S3_build_record.md`.

## Acceptance verification (all pass)

1. Run completes: EXIT=0, ExecutionTime=495.6 s, no FOAM errors.
2. Gate-A epsilon sweep present: {1e-3, 3e-4, 1e-4, 3e-5, 1e-5} (5 eps) with
   per-block |Jv|, |FD|, |err|, relL2, cosine, max-abs-diff+location.
3. Gate-B: 5 eps, per candidate (A and B) |analytic|, |FD|, norm ratio,
   relL2, cosine, internal/boundary-face rel error; relax factor 0.4; rAUAdj
   vs relaxed-rAU avg/min/max/norm ratio.
4. Anchors reproduced: momentum relL2=1.653e-4 cos=0.999999986 (frozen-phi
   FD, hU=0.388); GatePR h=0.001 relL2=1.12295958888e-10; |phiRb0-phi|
   L2=3.49247878283e-11. P-internal/P-boundary/P-total FD plateau stable
   across all 5 eps (P-total relL2 spread 0.000%).

## Key numbers from this run (identical to S2; interpretation is S4's job)

Gate A (J*v vs production-consistent FD with relaxed-rAU rebuild, dp=0):
  momentum  relL2=1.653e-4 cos=0.999999986
  P-internal relL2=0.858 cos=0.612  (max @ P cell 6161)
  P-boundary relL2=1.145 cos=0.153  (max @ P cell 7761 patch sideWalls)
  P-total    relL2=1.036 cos=0.370
  => the O(1) P-row mismatch PERSISTS under boundary-consistent production
     FD; it is NOT removed by rebuilding boundary fluxes through the exact
     NS.H SIMPLE path.

Gate B (isolated dphi/dU, p/design/alpha/k/omega/nutFrozen/nuEffFrozen fixed):
  Candidate A (current J, unrelaxed rAUAdj): normRatio=0.591, relL2=0.795
    (int 0.790, bnd 1.388), cos=0.607 — O(1) mismatch at every eps.
  Candidate B (production relaxed mobility x relaxed-H tangent incl. relax
    source + constrainHbyA pinning): normRatio=1.000, relL2=1.06e-11..1.06e-9
    (int ~ same, bnd ~ same), cos=1.000 — roundoff-level match at every eps.
  rAU: relaxFactorU=0.4; rAUAdj avg=7.05844426338e-07 min=9.99981274586e-09
    max=5.1321646829e-06; relaxed rAU avg=2.82337770533e-07
    min=3.99992509834e-09 max=2.05286587567e-06; avgRatio=2.50000000002;
    normRatio=2.49999999999.
  Note: the 2.5x rAU ratio does NOT appear 1:1 in dphi (Candidate-A normRatio
  is 0.59, not 0.4 or 2.5): the relax-source term in the relaxed H()
  (1.5*D0*psi/V) partially compensates the 0.4 mobility factor, exactly as
  audit_boundary.md §6 predicted.
