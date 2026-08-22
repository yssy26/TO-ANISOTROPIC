# TO-ANISOTROPIC — 下一阶段任务书（BFINAL-016：双缺陷修正轮——先对齐量具，后修算子与源）

> 用法：把下面「任务指令」整段复制给执行 agent。**派发本任务书 = 用户授权触碰已验证算子层（J 的通量耦合语义）与伴随源装配层（b_TC 折叠）**；授权范围以文中明确列出的文件与条款为限。分支 `agent/dsH-stage-b-validation` @ `39e3219`（BFINAL-015 审阅见 `evidence/agent-group/BFINAL-015/cycle-1/REVIEW_BY_REVIEWER.md`，其中 §3 的两个量化缺口是本轮的先导步骤）。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC（OpenFOAM-7 拓扑优化求解器，`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）的**执行 agent**。当前任务：**BFINAL-016 — 双缺陷修正轮**。BFINAL-015 已用真值链（b^T·w_true ≡ FD，机器精度）把 BFINAL-012 失败拆为两个独立缺陷：(1) gDP ≈2.13× = J 算子（P 行/通量耦合语义）；(2) J 目标 = b_TC 源/消除链。本轮先补两个测量缺口、再修、最后完整重验。

### 第〇阶段：量具对齐（先于任何修复，两个必做实验）

**A1 — 偏应力对齐（裁决 C1）**：导出矩阵不含冻结偏应力转置（`validateMmaUnlockGate.H:50-68` 记载的已知缺口），而生产算子含。做法（零风险，离线）：在 `b15_t2t3.py` 中按 `applyProdFrozenDeviatoricJT`（solveDiscreteFlowAdjointProduction.H:105-231）的公式用 numpy 复现偏应力块的 matrix-free 作用，将其加入切线解（迭代格式：w ← splu(core)·(rhs − J_dev·w)，收敛到 1e-12）或直接组装进稀疏阵。判据：加入后离线 b_PD^T·w′ 是否移动到与生产 ADJ（−5.405/+1.106/+14.045）一致（≥3 位有效数字）——
- 一致 → 「生产 vs 导出差异」由偏应力解释完毕，产/导对齐完成，P 行集中度须**在完整算子上重测**后再作为修复靶点；
- 不一致 → 存在第三处差异，先定位（逐块差分导出 vs 产线装配），报告后再继续。

**B1 — 热介导项量化（清洁化缺陷 2 的量具）**：恒等式为 FD_J = 流介导（−b_TC^T·w_true）+ 热介导（经 T 方程）。用既有 T 伴随/或直接对 T 场做设计 FD，把 D2/D3 的热介导项算出来，更新缺陷 2 的真实幅度表（D1 的 2.63× 已干净）。判据：修正后的「流介导+热介导 vs FD_J」残差在三个方向上的一致解释 ≥90%，否则缺陷 2 还有未建模成分。

### 第一阶段：修复（每项独立提交，逐项验证）

**A2 — 算子修正（授权范围：`solveDiscreteFlowAdjointProduction.H` 的 J 装配、`solveDiscreteFlowAdjoint.H` 的对应诊断装配）**：以第〇阶段对齐后的完整算子为基准，针对 P 行失配推导正确语义——从 NS.H 实际不动点（rAtU=松弛后、pEqn 非正交循环、phi=phiHbyA−pEqn.flux）推导 dphi/dU、dphi/dp 的精确组合，与现行 BFINAL-003 relaxed-SIMPLE 切线逐项对照，修正差异项。**J_PU/J_PP/J_UU 各块作为残差-Jacobi 的已验证数学不得无故重写——只允许推导证明的语义修正，且每一处修正须引用推导条目**。
**B2 — 源修正（授权范围：`AdjNS_HT.H`/`AdjHeatTransfer.H` 的 b_TC 折叠装配）**：按 B1 清洁化后的量具审计热伴随消除链（dJ/dT、dJ/dphi → 通量转置折叠 → 流行 rhs），逐项对照推导，修正缺陷项。

### 第二阶段：完整重验（顺序执行，任何一步失败即停止定位）

1. 等价性门 PRODPRECGAMGCHECK 仍 1e-16 级（若 A2 改了 P 行装配，须先过该门再继续）；
2. 伴随双标签收敛 ≤1e-9（迭代数变化须可解释）；
3. T2/T3 复测：w′/w_true 三方向 cos>0.999、relL2<1e-2（或按 B1 清洁化口径达到与 FD 噪声底相称的水平）；
4. b_TC 量具：流+热介导 vs FD_J 三方向一致；
5. **完整 BFINAL-012 FD 门重跑**：gDP ADJ/FD ∈ 1±0.05（或明示 1±0.10 平台）、J 符号全对且幅值 <5%（或明示 5–10% 平台）、gV 不回归；
6. 静态测试全过（9+，相关断言同步更新）；P4/P5 点测、RX-A、P1 源点测试全部复跑通过。

### 纪律与不变量

- `frozenGradientValidated=false`、`mmaUpdateEnabled=false` 不动；解锁是用户在审阅后的决定；
- 禁止拟合因子（2.13/2.63 等不得出现在代码里）；每处修正引用推导；
- 先归档再清理；EXECUTOR_SUMMARY.md + EXECUTOR_DONE 协议；选择性提交不 push；
- 停止条件：任一验收检查点失败 → 停止并定位报告；A1 对齐出现第三处差异 → 先报告再继续。

### 证据

`evidence/agent-group/BFINAL-016/cycle-1/`：A1/B1 对齐与量化数据、推导对照文档、每项修复的独立 diff 与验证、最终全套重验数据。

---

## 审阅者备注（不需复制给执行 agent）

- 审阅重点：A1 是否真把偏应力块加进了离线解（而非又绕过）；P 行集中度是否在完整算子上复测后才作为 A2 靶点；B1 恒等式是否把热介导项放对了边；A2 的每一处修正是否有推导条目支撑（不是拟合）；最终 FD 门是否三目标全过；
- 风险提示：A2 触碰 BFINAL-003 锁定层——修正后 J_PU 的残差-FD 外部锚（stageB8 类）必须复跑；若推导显示 BFINAL-003 的 P←U 切线本身需要语义升级，这构成对「锁定里程碑」的正式修正，须在报告与 handoff 中显式声明；
- 全过后下一站：C0 串行 MMA smoke 准备轮（解锁决定留给用户）。
