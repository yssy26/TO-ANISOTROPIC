# TO-ANISOTROPIC — 接手与独立审阅任务书（转交给下一个 AI agent）

> 用法：把下面「任务指令」整段复制给接手的新 agent。本页其余内容是新 agent 需要了解的上下文背景，可一并附上或让新 agent 自己读仓库文档。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC（OpenFOAM-7 拓扑优化求解器，仓库 `yssy26/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）的**独立审阅者**。你的任务是两件事：

1. **先完整了解这个项目**——物理背景、代码架构、阶段历史、当前状态、尚未解决的问题；
2. **再客观审阅现有代码与未来路径**——不带着历史结论进来，独立判断代码质量、证据可信度、路径合理性，并指出风险。

这是一个**只读优先的审阅任务**：以阅读、比对、推理为主；只有在你认为一个小验证（例如运行某个静态测试、grep 某段代码、跑一次短编译）能显著改变结论时才可以做，并在报告中说明你做了什么、为什么。**不要修改任何源码，不要推送 GitHub。**

### 第一步：按顺序读这些文档（优先级从高到低）

1. **`AI_AGENT_HANDOFF.md`**（仓库根目录，权威交接文档）——项目物理、范围锁定、BFINAL-001..008 时间线、BFINAL-009 停止原因、BFINAL-010/011/012 任务规格、C1 并行化计划、P0–P4 问题排序、安全不变量。
2. **`OPENFOAM_USAGE_AND_FILES.md`**（仓库根目录）——编译/运行/环境操作手册、src 文件职责表、已知坑。
3. **`src/validation/stage_b/BFINAL_009_PRODUCTION_PRECONDITIONER.md`**——当前任务（生产 pressure-GAMG 预条件器）的正式规格与验收序列。
4. **`evidence/agent-group/BFINAL-008/`** 下的 `FINAL_REPORT.md`、`cycle-1/VERIFY_FRESH_6320943.md`、`cycle-1/s4_verify_fresh_unpinned.json`、`s3_gate_metrics_fresh.json`——最近一次已闭合的门（P1–P8）及其原始证据。
5. **代码**（重点）：`src/solveDiscreteFlowAdjointProduction.H`（生产伴随 + 新预条件器）、`src/solveDiscreteFlowAdjoint.H`（诊断/oracle）、`src/sensitivity.H`、`src/rxPressureRowTranspose.H`、`src/validateMmaUnlockGate.H`、`src/tests/test_stage_b_safety_gates.py`。
6. 历史 `src/validation/stage_b/STAGE_B_*` 报告——**只作历史背景**，不要因为某份报告写得很详细就认为它仍然有效。

### 第二步：审阅重点（逐项给出你的判断）

**A. 代码审阅**
- 生产路径与诊断路径是否清晰分离？`solveDiscreteFlowAdjointProduction.H` 是否保持可扩展（O(N)/O(nnz)），有没有把诊断机器混进生产？
- 已验证的算子（J_PU、J_PP、R_x、`p.needReference()` 参考语义）在本次提交（4c49903 引入 pressureGAMG）后是否被改动？有无回归风险？
- `PRODPRECGAMGSETUP diagRelL2=0.041` 的停止判定：审阅 `solveDiscreteFlowAdjointProduction.H` 中对角比较逻辑（裸 `prodPressurePrecMatrix.diag()` vs `prodJacobiDiag`），独立判断 handoff 中的「`addBoundaryDiag` 有效对角不一致」假设是否成立，以及 BFINAL-010 的 interior/boundary/effective/offDiag 四组分解方案是否完备。
- 静态测试（`test_stage_b_safety_gates.py`）与 CI（`.github/workflows/ci.yml`）覆盖是否充分？缺口在哪？

**B. 证据审阅**
- 区分「已验证事实」「强假设」「推测」三类，逐条列出。
- BFINAL-008 的 P1–P8 结论是否真的由原始证据支撑（读 json/日志，不只读摘要）？
- 生产求解器（diagonal 与 pressureGAMG）至今**都未收敛**到 1e-9——这对「冻结湍流梯度门」意味着什么？文档里对「可信路线是 oracle/直接求解」的表述是否被证据支持？

**C. 路径审阅**
- BFINAL-010（矩阵等价性诊断）→ BFINAL-011（生产外求解收敛）→ BFINAL-012（当前源 FD 幅值门）→ Stage C（MMA → 并行 → 完整优化）的顺序是否合理？有无遗漏的前置条件？
- P0（生产求解器扩展性/稳健性）、P1（FD 幅值门）、P2（并行伴随，需 processor-patch 转置交换）、P3（优化循环未证明）、P4（文档漂移）的风险排序是否恰当？
- 并行伴随：确认「先串行收敛、再做并行」的顺序是否必要，processor-patch 缺失项清单是否完整（U-U / U-P / P-U / P-P / frozen deviatoric 转置 + halo exchange）。

### 第三步：交付一份审阅报告（Markdown，中文）

报告结构建议：

1. **项目理解摘要**（200–400 字，证明你已掌握全局）；
2. **代码审阅发现**：按 严重（会误导/错误）→ 中等（风险/隐患）→ 轻微（风格/文档）分级，每条给出文件:行号与理由；
3. **证据审阅结论**：已验证事实 / 强假设 / 推测 三张表，并注明每个结论的证据来源；
4. **路径审阅结论**：对 BFINAL-010/011/012 → Stage C 顺序的独立判断，含风险与建议的调整；
5. **给下一轮执行 agent 的具体建议**（3–8 条，可执行、按优先级排序）；
6. **你本次审阅做了什么验证**（如有）。

要求：**客观、独立、可复核**。不要重复手记的结论而不验证；不要因为「文档写了」就当作「代码实现了」。凡是与 handoff 文档结论不一致的地方，明确标出并给出你的依据。

---

## 给新 agent 的背景速览（可一并附上）

- 仓库/分支：`yssy26/TO-ANISOTROPIC` @ `agent/dsH-stage-b-validation`，HEAD `6a0004b`（2026-08-19）。
- 当前宏观阶段：**B-final**——冻结湍流梯度/生产伴随闭合；**MMA 锁定**（`mmaUpdateEnabled=false`、`frozenGradientValidated=false`）。
- 最近状态：BFINAL-008 诊断路径 P1–P8 全部 PASS（J_PU/J_PP/R_x 闭合，物理矩阵无人工 identity 行，直接切向求解 ~1e-13）；生产迭代求解器停滞（cond~1e19），pressureGAMG 预条件器在首个检查点停止（`diagRelL2=0.041 > 1e-8`）。
- 环境注意：编译用干净 PATH + `source /opt/openfoam7/etc/bashrc` + `unset FOAM_SIGFPE`；运行用 `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`（勿用旧二进制）；详细见 `OPENFOAM_USAGE_AND_FILES.md`。
