# AG-SMOKE-001 SMOKE-2 — Read-only re-verification evidence

Run: AG-SMOKE-001 (VALIDATION_ONLY smoke, plan -> pre_review -> execute -> post_review)
Step: SMOKE-2 — Re-run `git status --short` and diff against the pre-round snapshot
Mode: READ-ONLY (no project artifact modified; evidence write only under `evidence/agent-group/AG-SMOKE-001/`)
Date: 2025-08-16 (repo-local clock), workspace /home/ys/dsH/TO-ANISOTROPIC

## 1. Workspace identity (re-verified at start of this step)

```
pwd                        -> /home/ys/dsH/TO-ANISOTROPIC
git rev-parse --show-toplevel -> /home/ys/dsH/TO-ANISOTROPIC
git rev-parse HEAD          -> df695c7373297425e2d76767d2c415a9f7f55b18
```

- pwd == repo root == /home/ys/dsH/TO-ANISOTROPIC (exact match, skill §12 requirement met).
- HEAD is a 40-hex SHA and matches the approved plan's expected commit.

## 2. SMOKE-2 re-run of `git status --short`

Re-run output captured to `/tmp/status_smoke2.txt` (52 lines), identical to the SMOKE-1
post-capture snapshot `/tmp/status_post.txt`:

```
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

## 3. Diff results

### 3a. vs pre-round baseline `/tmp/status_baseline.txt` (51 lines)

```
12a13
> ?? evidence/
[diff exit: 1]
```

The ONLY difference is the single added line `?? evidence/` — the permitted evidence
directory created by SMOKE-1 and used by this round (plan `files_allowed`:
`evidence/agent-group/AG-SMOKE-001/**`). No tracked file changed; no other untracked
file added/removed.

### 3b. vs SMOKE-1 post-capture `/tmp/status_post.txt`

```
[diff exit: 0]
```

Byte-identical. The four read-only commands and this SMOKE-2 step modified nothing.

## 4. Checksums (byte-level proof)

```
sha256(/tmp/status_baseline.txt) = 81b895c435e8f5fb8f54369fb0f946067f7114c961816b6f1d3888c4436eba1c
sha256(/tmp/status_post.txt)     = 294d1c458f4982524866ebdbaa3e089951b889aad1444730fc3dedd403ea06c8
sha256(/tmp/status_smoke2.txt)   = 294d1c458f4982524866ebdbaa3e089951b889aad1444730fc3dedd403ea06c8
```

- pre-round baseline vs SMOKE-1 post-capture differ only by `?? evidence/` (checksum differs as expected);
- SMOKE-1 post-capture vs SMOKE-2 re-run: identical checksum -> no modification between SMOKE-1 and SMOKE-2.

## 5. Conclusion (SMOKE-2 acceptance)

- `git status --short` byte-identical before (SMOKE-1 post-capture) and after (SMOKE-2 re-run): PASS.
- No tracked or untracked project file changed: PASS (only permitted `evidence/` dir exists, unchanged since SMOKE-1).
- Pre-existing user dirty worktree (3 deletions, 2 modified src/ files, 46 pre-existing untracked entries)
  fully preserved; no git reset/clean used.
- The four read-only commands (pwd, git rev-parse --show-toplevel, git rev-parse HEAD,
  git status --short) are confirmed to modify nothing.
