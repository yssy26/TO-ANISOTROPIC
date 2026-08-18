# BFINAL-004-FIX cycle-3 — step4-run-record (executor)

- Stage: B-final · Mode: DIAGNOSTIC_ONLY · Step: step4-run
- Workspace: `/home/ys/dsH/TO-ANISOTROPIC` (git root == workspace) @
  `ca8a772b8339a36e7906ea32a7125061c47aa280` (unchanged)
- Run: `/home/ys/dsH/b2_case_smoke` (writable original; `constant/optProperties`
  line 87 `stageB6RxDesignOracle true;` — verified; nothing else changed)
- Binary: `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF` (absolute path), size
  4747104 B, mtime 2026-08-17 04:43 (fresh, built by cycle-3 step3 from the
  audited fixed header; sha256 15172c3a...)
- Log: `cycle-3/run_stageb6_cycle3.log` · **MTO_RC=0 (EXIT=0)**; START
  2026-08-16T20:44:18Z → END 2026-08-16T20:52:16Z (478 s)
- No SIGABRT / no `calculatedFvPatchField` / no `FOAM FATAL` (grep count = 0);
  21 `stageB6_*.mtx` artifacts fresh (04:52), copied to `cycle-3/artifacts/`
  (sha256 byte-identical to the case copies; `stageB6_rxa_FD_all_eps.mtx`
  verified identical).
- Independent recompute: `cycle-3/verify_recompute_cycle3.py` (pure Python, no
  numpy, no shared formula with the C++ probe) → `verify_recompute_cycle3.out`
  reproduces every log metric from the raw artifacts (rpa_terms identity
  max|col1+col2-col3|=1.29e-26; col3 == analytic R_P block exactly).

## 1. Regression anchors — byte-identical to BFINAL-003 (R_w/J CLOSED)

| anchor | BFINAL-003 locked | this run (run_stageb6_cycle3.log) |
|---|---|---|
| ExplicitJToracle maxRelL2 | 3.97863790105e-16 | 3.97863790105e-16 (L1267, L1657) |
| ExplicitJToracle maxRelL_U | 3.52693970502e-16 | 3.52693970502e-16 |
| ExplicitJToracle maxRelL_P | 3.98062510668e-16 | 3.98062510668e-16 |
| StageB5 momentum (eps 3e-5) | 1.65311e-4 / 0.999999986345 | 0.000165311363075 / 0.999999986345 |
| StageB5 P-total (eps 1e-3) | 1.52829e-4 | 0.000152828911949 |
| StageB5 rAU avgRatio | 2.50000000002 | 2.50000000002 |
| StageB5 GateB candB (eps 1e-5) | 1.05509085820e-09 | 1.0550908582e-09 |
| StageB4 dPhi_J L2 | 2.19206493915e-05 | 2.19206493915e-05 (L1207, L1597) |
| **Volume anchor projV_D1** | **0.1555186511144686** | **0.1555186511144686** (L1869) |

Volume anchor D2 = -0.08449667741658891, D3 = 0.06396237684432686 (production
replica; raw-design projection machinery regression green).

## 2. State-frozen / turbulence-frozen guarantee

Entry fingerprint (L1746) == exit fingerprint (L1880), byte-identical:
|U|2=16796.0670769, |p|2=7965118.75546, |phi|2=0.0041733366884,
|alpha|2=14587858210.1, |k|2=11566.7580246, |omega|2=13527763.4515,
|nutFrozen|2=0, |nuEffFrozen|2=0.00951359211546, |x|2=141.996485326,
|xp|2=141.993260856, |xh|2=142.223133214, eta5=0.748005161603, del=8.
No state variable re-solved; k/omega/nutFrozen/nuEffFrozen untouched; U/p/phi
rebuilt only in local copies inside rebuildResidual; filter solves act on
diagnostic local fields; no stageB6 field files written to the case (all probe
fields IOobject::NO_WRITE; verified `1/` contains no stageB6 files).

## 3. Gate results (executor + independent recompute; see verify_recompute_cycle3.out)

### RX-A — direct alpha derivative (fixed state), deltaAlpha-weighted analytic
| eps | R_U,alpha relL2 | R_U cos | R_P,alpha relL2 | R_P cos | FD R_P norm |
|---|---|---|---|---|---|
| 1e-3  | 9.238e-08 | 1 | 6.801e-07 | 1 | 6.22344978830e-10 |
| 3e-4  | 3.156e-07 | 1 | 2.333e-06 | 1 | 6.22344989835e-10 |
| 1e-4  | 9.081e-07 | 1 | 6.798e-06 | 1 | 6.22344999291e-10 |
| 3e-5  | 3.042e-06 | 1 | 2.269e-05 | 1 | 6.22344887681e-10 |
| 1e-5  | 9.280e-06 | 1 | 6.878e-05 | 1 | 6.22344766468e-10 |

R_P,alpha FD 5-eps plateau stable (6.22344978830e-10 … 6.22344766468e-10).
Weighted per-term (log L1751 + recompute):
|div(dphiHbyA_w)|=5.97626165832e-10, |-div(dflux_w)|=3.49379242410e-11,
|total_w|=6.22344977639e-10 == FD oracle; `stageB6_rpa_terms.mtx`
col1+col2==col3 to max|diff|=1.29e-26; col3 == `rxa_analytic` R_P block exactly.
(cycle-1 was unweighted-vs-weighted: relL2=1.906/cos=0.00066 for R_U,alpha —
B1 now fixed.)

### RX-B — xh-to-alpha closure (production alpha interpolation)
| eps | R_U,xh relL2 | R_U cos | R_P,xh relL2 | R_P cos | FD R_P norm |
|---|---|---|---|---|---|
| 1e-3  | 1.114e-06 | 1 | 5.614e-07 | 1 | 9.71259463259e-05 |
| 3e-4  | 1.003e-07 | 1 | 5.053e-08 | 1 | 9.71259108872e-05 |
| 1e-4  | 1.114e-08 | 1 | 5.615e-09 | 1 | 9.71259077717e-05 |
| 3e-5  | 1.003e-09 | 1 | 5.254e-10 | 1 | 9.71259074173e-05 |
| 1e-5  | 1.278e-10 | 1 | 4.501e-10 | 1 | 9.71259073862e-05 |

R_P,xh analytic |a|=9.71259073822e-05 vs FD 9.71259073862e-05 (eps=1e-5).
(cycle-1: R_P,xh relL2=1.012/cos=-0.176 — wrong drAU/dHbyA; B3 now fixed.)

### RX-C — raw-design D1/D2/D3 (filter→projection→alpha, U/p fixed)
xh-tangent / alpha-tangent / R_U,xd / R_P,xd relL2 (independent recompute):

| dir | quantity | eps=1e-3 | eps=3e-4 | eps=1e-4 | eps=3e-5 | eps=1e-5 |
|---|---|---|---|---|---|---|
| D1 | xh-tangent | 4.759e-4 | 8.100e-5 | 2.142e-5 | 1.209e-4 | 1.215e-4 |
| D1 | alpha-tangent | 3.899e-4 | 6.714e-5 | 1.954e-5 | 1.099e-4 | 1.101e-4 |
| D1 | R_U,xd | 4.543e-4 | 7.585e-5 | 2.298e-5 | 1.330e-4 | 1.339e-4 |
| D1 | R_P,xd | 2.031e-3 | 3.389e-4 | 4.141e-5 | 2.017e-4 | 2.014e-4 |
| D2 | xh-tangent | 7.231e-4 | 1.265e-4 | 2.865e-5 | 9.667e-5 | 9.657e-5 |
| D2 | alpha-tangent | 6.321e-4 | 1.109e-4 | 2.585e-5 | 9.209e-5 | 9.202e-5 |
| D2 | R_U,xd | 7.043e-4 | 1.217e-4 | 2.945e-5 | 1.002e-4 | 1.001e-4 |
| D2 | R_P,xd | 1.807e-3 | 3.419e-4 | 7.014e-5 | 1.561e-4 | 1.559e-4 |
| D3 | xh-tangent | 7.861e-4 | 1.412e-4 | 2.577e-5 | 9.702e-5 | 9.700e-5 |
| D3 | alpha-tangent | 6.742e-4 | 1.214e-4 | 2.264e-5 | 9.139e-5 | 9.137e-5 |
| D3 | R_U,xd | 8.659e-4 | 1.480e-4 | 2.768e-5 | 9.695e-5 | 9.692e-5 |
| D3 | R_P,xd | 3.601e-3 | 6.525e-4 | 1.100e-4 | 1.870e-4 | 1.867e-4 |

All cos ≥ 0.9999936. The eps=1e-3 R_P,xd points (1.8e-3…3.6e-3) and the
~1e-4…2e-4 floor at eps≤3e-5 are FD central-difference truncation + filter
solver tolerance; relL2 converges quadratically (1e-3→3e-4→1e-4) to the
eps-independent analytic side (the converged truth). R_P,xd at eps=1e-4:
D1=4.14e-5, D2=7.01e-5, D3=1.10e-4 — all below the task §13 diagnostic ceiling
relL2≤1e-3, preferred ≤1e-4 met for D1/D2 (D3 1.10e-4 marginally above).

### §9/§10 — λ-weighted importance and production completeness (analytic, eps-independent)
| dir | ‖R_U,x d‖ | ‖R_P,x d‖ | ratio ‖R_P‖/‖R_U‖ | λ_UᵀR_U,x d | λ_PᵀR_P,x d | D_momentum | D_pressure | D_total | prod gsensR·d | \|D_mom−prod\|/\|prod\| |
|---|---|---|---|---|---|---|---|---|---|---|
| D1 | 0.566318445188 | 0.000193062746532 | 3.409e-4 | 0.0486867554277 | 0.00754887795106 | −0.0486867554277 | −0.00754887795106 | −0.0562356333787 | −0.0486867553288 | 2.03e-9 |
| D2 | 0.861583739658 | 0.000447977917016 | 5.199e-4 | −0.229079662362 | −0.0685927693909 | +0.229079662362 | +0.0685927693909 | +0.297672431753 | +0.229079662485 | 5.40e-10 |
| D3 | 0.851652729879 | 0.000259542840587 | 3.048e-4 | 0.0301314016208 | 0.00621002821874 | −0.0301314016208 | −0.00621002821874 | −0.0363414298395 | −0.0301314015977 | 7.65e-10 |

- λ_P-weighted pressure-row fraction |D_pressure/D_momentum| = **15.5 % (D1),
  29.9 % (D2), 20.6 % (D3)** — NON-NEGLIGIBLE (raw norm ratio 3-5e-4 but
  λ_P-weighted effect 15-30 %).
- prod gsensR·d == D_momentum to ~1e-9 → the current production
  `gsenshPressureDrop` assembly is EXACTLY the momentum piece; the pressure-row
  piece −pcᵀR_P,x d is absent from production.

## 4. Scope / source integrity

- `git status --porcelain=v1` byte-identical to cycle-3/pre-state (only the two
  pre-existing BFINAL-003 production heads modified, hashes unchanged); probe
  header sha256 09832dd6... unchanged (audit+fix step made no further edit);
  no production math / sensitivity / residual / J / filter / projection /
  conductivity / MMA / thresholds touched; FLOW_JACOBIAN_FROZEN maintained
  (anchors byte-identical); no SIMPLE-transpose; no MMA; mmaUpdateEnabled=false
  (design not advanced).
- Artifacts: `cycle-3/artifacts/stageB6_*.mtx` (21 files, fresh 04:52,
  byte-identical to the case copies).

## 5. Routing

step4 acceptance MET: EXIT=0; RX-A/B/C + §9/§10 tables complete; volume anchor
0.1555186511144686 exact; no SIGABRT/calculatedFvPatchField; fingerprint
frozen; OFstream artifacts complete (all 21 .mtx present with fresh mtimes —
no tail truncation). → step5: FINAL_REPORT.
