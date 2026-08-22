# TO-ANISOTROPIC — 下一阶段任务书（BFINAL-017：装配审计 + 三缺陷修复轮）

> 用法：把下面「任务指令」整段复制给执行 agent。**派发本任务书 = 用户授权**：算子装配层（`solveDiscreteFlowAdjointProduction.H`/`solveDiscreteFlowAdjoint.H` 的 J）、伴随源层（`AdjNS_HT.H`/`AdjHeatTransfer.H` 的 b_TC 折叠）、收缩装配层（`sensitivity.H` 收缩项、`rxPressureRowTranspose.H` 语义对齐）。**`filter_chainrule.H`/过滤/投影/MMA 仍在禁改清单**——若审计把缺陷定位到那里，停下等用户单独授权。分支 `agent/dsH-stage-b-validation` @ `e67df00`（审阅见 `evidence/agent-group/BFINAL-016/cycle-1/REVIEW_BY_REVIEWER.md`）。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC（OpenFOAM-7 拓扑优化求解器，`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）的**执行 agent**。当前任务：**BFINAL-017 — 装配审计与三缺陷修复轮**。BFINAL-016 已把缺陷地图修正为三项（全部带干净外部量具）：①算子 P 行/通量耦合（1.769/0.856/2.274）；②装配链 vs λ-direct（D2 147%）；③b_TC 热消除折叠（86%/70%/符号翻转）。本轮顺序：审计装配矛盾 → 修复（装配→算子→源，各独立提交）→ 完整重验。

### 第〇阶段：装配矛盾审计（两个判别实验，先于一切修复）

**AD1 — 链式伴随点测试（隔离机理 b）**：对 `filter_chainrule.H` 做伴随性检验：取确定性场 g 与方向 d，验证 ⟨chain(g), d⟩ == ⟨g, chain^T(d)⟩（chain^T 按 filter/投影转置公式独立实现，勿复用生产代码）。若不闭合（>1e-10 相对）→ 链非线性（adaptive-eta 依赖被传播场）或转置实现错 → **缺陷在链上，属禁改层 → 停止并报告等授权**。若闭合 → 链无罪。
**AD2 — R_x 路径逐行对照（隔离机理 a）**：生产 R_x 收缩（rxPressureRowT·pc + gsenshMomentum 的 −(U&Uc)V·dAlphaDxh）与 stageB6 oracle 的 anRUa/assembleWeightedRPa（rxc 构造）逐单元对照（用已导出的量具数据，必要时补开关门控导出）。定位 147% 矛盾的具体项（动量行/压力行/边界项/V-因子/dAlphaDxh 裁剪）。
判别输出：装配缺陷的准确位置 + 修正后三项缺陷的最终形状（审计可能改变算子缺陷的表观幅度——若 AD2 揭示装配缺陷吞噬了部分算子信号，更新量具表）。

### 第一阶段：修复（每项独立提交，逐项引用推导/审计条目，禁止拟合因子）

1. **装配修正**（依 AD1/AD2 结论；授权范围：sensitivity.H 收缩项、rxPressureRowTranspose.H 语义对齐；**若结论指向 filter_chainrule.H → 停止等授权**）；
2. **A2 算子修正**：以干净口径复测 P 行/通量耦合块签名（未钉扎、含偏应力、完整算子），从 NS.H 实际不动点推导 dphi/dU、dphi/dp 精确语义，修正 J 装配差异项；修正后 J_PU 的残差-FD 外部锚（stageB8 类）必须复跑——**对 BFINAL-003 锁定里程碑的任何正式修正须在报告与 handoff 显式声明**；
3. **B2 源修正**：b_TC 热消除折叠逐项对照推导（dJ/dT、dJ/dphi → 通量转置折叠 → 流行 rhs），注意 B1 揭示的流/热近对消——修正后用「流+热 vs FD_J」恒等式验证三方向闭合。

### 第二阶段：完整重验（任何一步失败即停止定位）

1. PRODPRECGAMGCHECK 四值仍 1e-16 级；
2. 伴随双标签收敛 ≤1e-9（迭代数变化可解释）；
3. T2/T3 干净口径复测：−λ^T rxd / b^T w_true 三方向 ∈ 1±0.05（或明示 1±0.10 平台）；
4. 装配量具：sensitivity 装配 == λ-direct（三方向 ≤1e-8 相对）；
5. b_TC 量具：流+热介导 vs FD_J 三方向闭合（≤5%）；
6. **完整 BFINAL-012 FD 门重跑**：gDP、J 符号+幅值、gV 全部达标（口径同前）；
7. 静态测试（9+）+ RX-A、P1 源点测试、P4/P5 点测全部复跑通过。

### 纪律与不变量

- `frozenGradientValidated=false`、`mmaUpdateEnabled=false`；解锁是用户审阅后的决定；
- 先归档再清理；EXECUTOR_SUMMARY.md + EXECUTOR_DONE；选择性提交不 push；
- 每一阶段产物（审计表、推导条目、逐项验证）进 `evidence/agent-group/BFINAL-017/cycle-1/`。

---

## 审阅者备注（不需复制给执行 agent）

- 审阅重点：AD1 的 chain^T 是否独立实现（不得复用被审计对象）；AD2 是否逐单元而非只看投影和；修复是否逐项引用审计条目；A2 是否在干净口径复测签名后才动 J；最终 FD 门三目标是否全过；
- 风险提示：三处修复叠加后若 FD 门仍不过，须用三项量具分别归因（算子/装配/源各自的残差），不得笼统重试；
- 全过后下一站：C0 串行 MMA smoke 准备轮（解锁证据包 + 用户决定）。
