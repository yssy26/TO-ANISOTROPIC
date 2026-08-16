# AG-SMOKE-001 — Consolidated smoke-test output (SMOKE-3)

Run: AG-SMOKE-001 — Read-only agent-group pipeline smoke test (plan -> pre_review -> execute -> post_review)
Step: SMOKE-3 — Consolidated record of command outputs + pre/post git status
Mode: VALIDATION_ONLY / READ-ONLY (no project artifact modified; evidence write only under `evidence/agent-group/AG-SMOKE-001/`)
Workspace: /home/ys/dsH/TO-ANISOTROPIC
Date: 2025-08-16 (repo-local clock)

Companion evidence: `smoke1-raw-output.md` (SMOKE-1 raw capture), `smoke2-verify.md` (SMOKE-2 diff/checksum proof).

---

## 1. Workspace identity (the four validation commands)

```
$ pwd
/home/ys/dsH/TO-ANISOTROPIC

$ git rev-parse --show-toplevel
/home/ys/dsH/TO-ANISOTROPIC

$ git rev-parse HEAD
df695c7373297425e2d76767d2c415a9f7f55b18

$ git status --short
 D CODEX_LOCAL_VALIDATION.md
 D LOCAL_FROZEN_GRADIENT_VALIDATION_REPORT.md
 D README.md
 M src/solveDiscreteFlowAdjoint.H
 M src/validateStageBRepeatability.H
?? build_minifd.sh
?? build_oracle.sh
?? check_dev2.sh
?? check_fix.sh
?? check_load.sh
?? check_oracle.sh
?? crosscheck_pd.py
?? evidence/
?? peek_log.sh
?? probe_oracle.sh
?? probe_python.sh
?? probe_windows_python.sh
?? project_pd.py
?? ps_check.sh
?? readme.md
?? run_b2diag.sh
?? run_b4_export.sh
?? run_b4_import.sh
?? show_cfg.sh
?? solve_oracle_pd.py
?? src/validation/
?? src/validation_stage_a_diagnostic_cleanup_build.log
?? src/validation_stage_a_processor_unit.log
?? src/validation_stage_a_rebuild.log
?? src/validation_stage_a_unit.log
?? src/validation_stage_b01_build.log
?? src/validation_stage_b02_build.log
?? src/validation_stage_b1_build.log
?? src/validation_stage_b2_build.log
?? src/validation_stage_b2diag_build.log
?? src/validation_stage_b2diag_build2.log
?? src/validation_stage_b2diag_build3.log
?? src/validation_stage_b4_build.log
?? src/validation_stage_b4b_build.log
?? src/validation_stage_b4c_build.log
?? src/validation_stage_b4d_build.log
?? src/validation_stage_b4e_build.log
?? src/validation_stage_b4f_build.log
?? src/validation_stage_b4g_build.log
?? src/validation_stage_b4h_build.log
?? src/validation_stage_b4i_build.log
?? src/validation_stage_b4j_build.log
?? src/validation_stage_b4k_build.log
?? src/validation_stage_b4l_build.log
?? src/validation_stage_b4m_build.log
?? src/validation_stage_b4n_build.log
?? write_C.sh
```

Interpretation:

- `pwd` == `git rev-parse --show-toplevel` == `/home/ys/dsH/TO-ANISOTROPIC` (exact match) — workspace identity confirmed (skill §12 requirement met).
- `git rev-parse HEAD` returns a non-empty 40-hex SHA `df695c7373297425e2d76767d2c415a9f7f55b18`, matching the approved plan's expected commit.
- `git status --short` shows the pre-existing dirty worktree, fully preserved: 3 deletions (README.md, CODEX_LOCAL_VALIDATION.md, LOCAL_FROZEN_GRADIENT_VALIDATION_REPORT.md), 2 modifications (src/solveDiscreteFlowAdjoint.H, src/validateStageBRepeatability.H), and 46 pre-existing untracked entries plus the permitted `evidence/` dir.

---

## 2. Pre-round git status snapshot (baseline)

Captured by the first PATCH-step pre-state baseline (before this round created anything):

- File: `/tmp/status_baseline.txt` (51 lines)
- sha256: `81b895c435e8f5fb8f54369fb0f946067f7114c961816b6f1d3888c4436eba1c`

Contents: identical to the listing in section 1 **minus** the single line `?? evidence/`
(51 pre-existing entries: 3 deletions + 2 modified src/ files + 46 pre-existing untracked files/dirs).

---

## 3. Post-round git status snapshots

| Snapshot | File | Lines | sha256 |
|---|---|---|---|
| SMOKE-1 post-capture | /tmp/status_post.txt | 52 | `294d1c458f4982524866ebdbaa3e089951b889aad1444730fc3dedd403ea06c8` |
| SMOKE-2 re-run | /tmp/status_smoke2.txt | 52 | `294d1c458f4982524866ebdbaa3e089951b889aad1444730fc3dedd403ea06c8` |
| SMOKE-3 pre-write | /tmp/status_smoke3_pre.txt | 52 | `294d1c458f4982524866ebdbaa3e089951b889aad1444730fc3dedd403ea06c8` |

All three post-round snapshots are byte-identical (`diff` exit 0).

---

## 4. Diff results

### vs pre-round baseline (`/tmp/status_baseline.txt`)

```
12a13
> ?? evidence/
[diff exit: 1]
```

The ONLY difference after the whole round is the single added line `?? evidence/` — the
permitted evidence directory created by SMOKE-1 and used by this round (plan
`files_allowed: evidence/agent-group/AG-SMOKE-001/**`). No tracked file changed; no other
untracked file added/removed.

### between post-round snapshots

```
[diff exit: 0]  (all three snapshots byte-identical)
```

---

## 5. Acceptance verification (SMOKE-1/SMOKE-2/SMOKE-3)

1. `pwd` == `/home/ys/dsH/TO-ANISOTROPIC` and `git rev-parse --show-toplevel` == `/home/ys/dsH/TO-ANISOTROPIC` — workspace identity confirmed.
2. `git rev-parse HEAD` returns a non-empty 40-hex commit `df695c7373297425e2d76767d2c415a9f7f55b18` (matches approved plan).
3. `git status --short` byte-identical across SMOKE-1 post-capture, SMOKE-2 re-run, and SMOKE-3 pre-write (sha256 `294d1c…` == `294d1c…` == `294d1c…`); vs the pre-round baseline the only delta is the permitted `?? evidence/` line → the four read-only commands modify no project file.
4. Evidence writes occurred only under `evidence/agent-group/AG-SMOKE-001/` (smoke1-raw-output.md, smoke2-verify.md, smoke-output.md); nothing outside that dir was written by this round.
5. Pre-existing user dirty worktree fully preserved; no `git reset --hard` / `git clean -fd` / destructive rollback used.

---

## 6. Conclusion

- VALIDATION_ONLY smoke round executed end-to-end: plan -> pre_review -> execute -> post_review ready for independent re-verification.
- All pass criteria for AG-SMOKE-001 satisfied from the Executor side; Post-Reviewer must independently re-run the same four read-only commands and confirm identical outputs before declaring PASS.
- Expected Post-Reviewer reproduction: pwd=`/home/ys/dsH/TO-ANISOTROPIC`, toplevel=`/home/ys/dsH/TO-ANISOTROPIC`, HEAD=`df695c7373297425e2d76767d2c415a9f7f55b18`, `git status --short` == 52-line snapshot (sha256 `294d1c458f4982524866ebdbaa3e089951b889aad1444730fc3dedd403ea06c8`).

Evidence file: `evidence/agent-group/AG-SMOKE-001/smoke-output.md`
