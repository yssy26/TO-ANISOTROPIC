# BFINAL-025 cycle-1 — 预注册（目标函数 Q 化 + FD 门复测轮）

先于本轮全部运行数据落盘（Phase 0 唯一产物）。授权凭证 = 派发文件
`NEXT_TASK_BFINAL025.md`（用户已授权，2026-08-21）。工作仓库
`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation` @ f0034c0（BFINAL-024）。

---

## 1. 目标函数 Q 定义与能量守恒一致性

### 1.1 病理回顾（本轮动机）
旧目标 `maximizeColdOutletTemperature`：
```
J = -N/(M·Tref),  N = Σ(phi_out·T_out),  M = Σ(phi_out),  Tmix = N/M
dJ/dphi_f = -(T_f - Tmix)/(M·Tref)
```
出口面温 T_f 是 Tmix 的面内插值，`(T_f − Tmix) ≈ 0` → dJ/dphi 近似相消，
辅以 1/M 缩小 → FD 门 15 轮反号/幅值病理（B24：D1 FD_J=+0.0437 反号，
ADJ_J=−0.035835，signJ=0）。

### 1.2 Q 定义
```
Q = N - T_in·M
T_in = 冷入口温度，读取自 case BC（不硬编码）：
      gAverage(T.boundaryField()[coldInletPatchID])
```
b24_gauge/b25_qgate `0/T`: `inlet { type fixedValue; value uniform 600; }`
→ T_in = 600 K。物理含义：以入口温度为焓参考的冷流焓增益（W）=
冷流实际吸收的热量 = 换热总量（total heat transfer）。

### 1.3 能量守恒一致性论证
`costfunction.H` L132-170 frozen-hot 能量诊断（既有代码，泛型）：
```
QcoldW = rhoc·(coldEnthalpyOut − coldEnthalpyIn)
coldEnthalpyOut = Σ( TEqn.flux() @ coldOutletPatchID )
coldEnthalpyIn  = −Σ( TEqn.flux() @ coldInletPatchID )
thermalBalanceError = |QhotW − QcoldW| / max(|QhotW|,|QcoldW|,SMALL)
```
`TEqn.flux()` 是装配好的温度方程通量（对流 + 扩散，匹配所选离散格式）。
- 出口 `zeroGradient` → 出口扩散通量 ≈ 0 → `coldEnthalpyOut ≈ Σ(phi_out·T_out) = N`；
- 入口 fixedValue 均匀 600 → `coldEnthalpyIn ≈ T_in·Σ(−phi_in) = T_in·M_in`；
  质量守恒（coldInletFlow/coldOutletFlow 诊断）→ `M_in ≈ M_out = M`；
- ⇒ `QcoldW ≈ rhoc·(N − T_in·M) = rhoc·Q`。

即 **Q = N − T_in·M 正是既有 frozen-hot 能量诊断所跟踪的冷流焓增益（对流归一化）**。
最大化 Q = 最大化实际换热总量。守恒锚：thermalBalanceError 由既有诊断逐状态监控
（mmaUpdateEnabled=false 下仅 Info，本案例不触发 FatalError）。

### 1.4 与旧目标关系（修复点）
```
J   = −N/(M·Tref) = −Tmix/Tref                 （旧）
J_Q = −(N − T_in·M)/(Tref·M_frozen)            （新）
快照下（M = M_frozen）：J_Q = −(Tmix − T_in)/Tref
dJ/dT   : −phi_f/(M·Tref)         →  −phi_f/(Tref·M_frozen)     （几乎同形，分母仅 M→M_frozen）
dJ/dphi : −(T_f − Tmix)/(M·Tref)  →  −(T_f − T_in)/(Tref·M_frozen)（相消 → O(温升)，非小）
```
dJ/dT 通道不变性保证热通道（F1、dJ/dT 路径）行为与 B24 基本一致（已知良好）；
dJ/dphi 通道消除相消是修复核心。

---

## 2. 归一化设计：J_Q 与 M_frozen 快照语义

### 2.1 定义
```
J_Q = −Q/Qref,   Qref = thermalReferenceTemperature × M_frozen
M_frozen = Σ( phi @ coldOutletPatchID )   每外层状态入口快照一次
Tref = thermalReferenceTemperature（案例 600.0，读 optProperties）
```
J_Q 无量纲。状态入口 M = M_frozen → J_Q = −(Tmix − T_in)/Tref。

### 2.2 快照放置与刷新时机（逐槽）
- **声明**：`createFrozenHotRegionFields.H` L16 附近新增
  `scalar thermalObjectiveMassFlowFrozen(0.0);`（热目标类型读取后即可使用）。
- **刷新**：`computeObjective.H` **入口顶部**（新类型分支）：
  `M_frozen = sum(phi@outlet); reduce; 若 mag(M_frozen)<SMALL → FatalError`。
  每次外层状态进入 computeObjective.H 刷新一次。
- **使用方（同一状态内全部 J 求值共用同一 M_frozen）**：
  - 主循环：`computeObjective.H`（MTO_HF.C L55）先于 `costfunction.H`（L103）
    → costfunction 的 MeanT/Q 与 normalizedThermalObjective 用新鲜快照；
  - 候选评估：`evaluateCandidate.H`（enableSSTAcceptanceCheck=false，分支仍实现）；
  - FD 门：`validateCommon.H evaluateObjective()` —— 全部 FD ±h 求值点共用状态级快照（见 2.3）；
  - CSV：`writeCSVLog.H`（enableCSVLog=true）；
  - F1/F2/F3 与 B2：见 2.3。
- **链内重包含 computeObjective.H 的唯一位置**（刷新点全集）：
  - `validateStageBRepeatability.H` runFullB1Chain L444（stageBEnabled=false，本案例不执行）；
  - `validateStageB2GradientAmplitude.H` runFullB2Chain L449（state-B 基线，L551 恰调用一次）。
- **FD 探针内禁止重包含 computeObjective.H**（无此代码路径；probeFD 用 evaluateObjective()）。

### 2.3 快照恒定性证明（审阅者注①：FD ±h 求值点不得偷偷用当前 M 重除）
逐路径代码级证据（本轮全部文件通读确认）：
- **B2**（`validateStageB2GradientAmplitude.H`）：`runFullB2Chain()` 在 L551
  恰调用一次（位于 ADJ 投影 + identity tsv + 全部 FD 扫描之前），其重包含
  computeObjective.H 仅刷新 state-B 基线快照。`probeFD`（L634-737）流程 =
  restoreStateInto(baselineState) → x±h·Dk → recomputeDesignChain() →
  `solveFrozenPrimalConverged(1000)` → `evaluateObjective()`/`evaluatePressureDrop()`/
  `evaluateVolumeConstraint()` —— **不重包含 computeObjective.H**。
  ⇒ 全部 FD ±h 求值点与全部 ADJ 投影共用**同一个** state-B M_frozen 快照。
  B2State 结构体不携带该量（探针从不修改它，无需保存/恢复）。
- **F1/F2/F3**（`validateFrozenGradient.H`）：solveThermalOnly /
  solveFrozenPrimalFixed 之后直接 `evaluateObjective()`（L688/707/876/895/1084/1107），
  无重包含。该文件（MTO_HF.C L125）位于主循环 computeObjective.H（L55）之后
  → 全程用 state-A 快照，±h 探针间不变。
- **B1**（`validateStageBRepeatability.H`）：runFullB1Chain L444 重包含
  （自洽；stageBEnabled=false，不执行）。
- **候选评估**（`evaluateCandidate.H`）：不重包含，用 state-A 快照。
⇒ 结论：**任何单个梯度/FD 求值内部，Qref = Tref·M_frozen 严格常数**；
M_frozen 仅在状态入口刷新一次。审阅者注①满足。

### 2.4 快照 vs 当前 M
J_Q = −(N − T_in·M)/(Tref·M_frozen)。M_frozen 视为常数：
dJ_Q/dphi 只剩分子 N − T_in·M 的显式依赖 → `−(T_f − T_in)/(Tref·M_frozen)`（非小）。

---

## 3. 导数与逐槽审计

### 3.1 解析导数（实现形式）
```
dJ_Q/dT[celli]  : 冷出口相邻胞 f：−phi_f/(Tref·M_frozen)；其余 0
                  （与旧型同形，分母当前 M → M_frozen）
dJ_Q/dphi[f]    : 冷出口面 f：−(T_f − T_in)/(Tref·M_frozen)；其余 0
                  T_f = T.boundaryField()[coldOutletPatchID].patchInternalField()
T_in            : gAverage(T.boundaryField()[coldInletPatchID])（从 case BC 读）
```
符号自检：T_f > T_in（出口被加热）→ dJ_Q/dphi_f < 0：增大出口通量降低 J_Q（改善）✓；
dJ_Q/dT_f < 0：升高出口温度降低 J_Q ✓。两通道同号、量级 O(温升/(Tref·M_frozen))，
无相消。

### 3.2 逐槽审计表（sensitivity.H / 伴随链中每个消费目标泛函的槽位）
标签：`需新分支` / `旧专用-新类型真空` / `泛型（自动正确，字段驱动）`。

| # | 槽位（文件:行） | 消费方式 | 标签 |
|---|---|---|---|
| 1 | computeObjective.H L28-57 dJ/dT | −phi_f/(M·Tref)，M=当前 → 新：M_frozen | 需新分支 |
| 2 | computeObjective.H L73-99 dJ/dphi | −(T_f−Tmix)/(M·Tref) → 新：−(T_f−T_in)/(Tref·M_frozen) | 需新分支 |
| 3 | computeObjective.H 入口顶部 | 新增 M_frozen 快照刷新 | 需新分支 |
| 4 | costfunction.H L9-41 MeanT | Tmix=N/M → 新：Q=N−T_in·M | 需新分支 |
| 5 | costfunction.H L71-74 normalizedThermalObjective | −MeanT/Tref → 新：−Q/(Tref·M_frozen) | 需新分支 |
| 6 | evaluateCandidate.H L72-92 candidateMeanT | Tmix → 新：Q_cand=N−T_in·M | 需新分支 |
| 7 | evaluateCandidate.H L140-143 candidateNormalizedThermal | → 新：−Q_cand/(Tref·M_frozen)，用 M_frozen 快照 | 需新分支 |
| 8 | validateCommon.H L16-34 evaluateObjective() | **硬编码旧型** −meanT/Tref（全部 FD 门使用） | 需新分支 |
| 9 | writeCSVLog.H L36-47 csvMeanT | 旧专用（新类型下算错值） | 需新分支（未列入派发清单，但 enableCSVLog=true，必须） |
| 10 | writeCSVLog.H L189-191 csvNormalizedThermal | 旧专用 | 需新分支 |
| 11 | createFrozenHotRegionFields.H L133-155 白名单 | frozen-hot 仅允许 maximizeColdOutletTemperature | 需新分支 |
| 12 | createFrozenHotRegionFields.H L16 附近 | 新增 thermalObjectiveMassFlowFrozen 声明 | 需新分支（新增） |
| 13 | createFrozenHotRegionFields.H L605-622 初始化区面积 dJ/dT | 旧专用；新类型下字段保持 0 → 空值检查平凡通过 | 旧专用-新类型真空 |
| 14 | validateDiscreteObjectiveDerivatives.H L21-42 面积 FD 检查 | 旧专用；新类型下 derivative=0 → error 恰为 0 → 平凡通过 | 旧专用-新类型真空 |
| 15 | AdjHeatTransfer.H L30 | thermalTransposeEqn.source() = thermalObjectiveDerivative | 泛型 |
| 16 | AdjHeatTransfer.H L150-155 | thermalObjectiveAdjointSource = dJ/dT / V | 泛型 |
| 17 | AdjHeatTransfer.H L224-310 | thermalDiffusionDerivativeDTCell（Tb、dDTDxh、T 面级） | 泛型 |
| 18 | AdjNS_HT.H L23-41 | discreteExternalFaceFluxAdjoint = 热残差项 + thermalObjectiveFluxDerivative@coldOutlet | 泛型（字段驱动） |
| 19 | solveDiscreteFlowAdjointProduction.H L340-419 | prodExternalTranspose 消费 discreteExternalFaceFluxAdjoint | 泛型（禁改清单） |
| 20 | sensitivity.H L65-66 | fsenshMeanT 动量项 | 泛型 |
| 21 | sensitivity.H L145-186 | thermalCouplingFaceFunctional（冷出口 outletFluxDeriv 字段驱动） | 泛型 |
| 22 | sensitivity.H L193-203 | rxPressureRowTranspose.H | 泛型（禁改清单） |
| 23 | sensitivity.H L239/255-259 | gsenshMeanTPressureRow / gsenshMeanTFluxDirect / fsenshMeanT | 泛型 |
| 24 | validateFrozenGradient.H F1/F2/F3 | 经 evaluateObjective()（validateCommon 新分支自动正确）；dfdxThermal 经 thermalDiffusionDerivativeDTCell | 泛型 |
| 25 | validateStageBRepeatability.H runFullB1Chain L444 | 重包含 computeObjective.H（自洽刷新） | 泛型；stageBEnabled=false |
| 26 | validateStageB2GradientAmplitude.H runFullB2Chain L551 + probeFD/scanDir | 单快照自洽（§2.3 证明） | 泛型 |

**审计结论**：
- `AdjHeatTransfer.H` / `AdjNS_HT.H` / `solveDiscreteFlowAdjointProduction.H` /
  `sensitivity.H` **无需新增分支** —— 目标泛函消费链完全字段驱动
  （经 thermalObjectiveDerivative / thermalObjectiveFluxDerivative 两个源字段），
  自动适配新类型；b_TC 折叠路径（#18+#19）同样字段驱动，无需改动。
- 旧类型两条路径（legacyMeanTemperature、maximizeColdOutletTemperature）**逐位不变**：
  全部新增代码位于新类型的 else-if 分支，旧分支文本不动。
- 未列入派发钩子清单但审计发现必须改的槽：**writeCSVLog.H**（#9/#10，enableCSVLog=true）。

---

## 4. 预注册预测（先于数据）

| # | 项 | 预测/门 |
|---|---|---|
| G1 | wmake WMAKE_EXIT=0 | 必过 |
| G2 | 静态测试 15/15 | 必过 |
| G3 | J_Q FD 门符号表：D1/D2/D3 signJ=1（3/3） | 必过（修复核心） |
| G4 | J_Q 幅值因子 ADJ/FD ∈ [0.9, 1.1]（正式步 relJ ≤ 0.10，目标 ≤ 0.05，ε 平台） | 必过 |
| G5 | baselineJ_Q = **−0.16745040316938**（由 B24 旧证据导出：MeanT_B=700.470241901628 K → (MeanT_B−600)/600；状态入口恒等式 J_Q=−(Tmix−T_in)/Tref 精确成立，M_frozen=state-B 入口 M） | 强数值交叉检验 |
| G6 | gV（体积约束方向导）FD_gV / ADJ_gV 与 B24 **逐位相同**（体积约束链与目标无关；mmaUpdateEnabled=false 下状态轨迹目标无关且确定性） | 必过（逐位） |
| G7 | 等价性门四值（identity 指纹：thermalCoupling gradProxy −510.7047296733922 / L2 2.451199924478948；pressureDrop gradProxy −315867.3295218884 / L2 97.47897160928267）与 B24 **逐位相同** | 必过（逐位） |
| G8 | gDP 同 B24（ADJ/FD 因子 ~2.0×：D1 2.061 / D2 2.007 / D3 2.100；已知状态，**非本轮失败判据**；passDP 判据 relDP≤0.10 不变） | 如实（与 B24 一致） |
| G9 | F1 热通道硬门：median ≤ 1%、max ≤ 5%、signMismatch=0 | 必过（dJ/dT 与旧型几乎同形，热通道已知良好） |
| G10 | 双标签伴随 trueRelRes（TC ~1e-13、PD ~1e-13）+ GRADSTABLE；FGMRES 迭代 ~1123/1323（算子未变，仅伴随源变化） | 必过 |
| G11 | 旧目标回归（b24_gauge 重跑，零输入改动）：FD 侧三列（J/gDP/gV）、ADJ_gV、gDP 因子、等价性门四值**逐位复现 B24 证据** | 必过（逐位） |
| G12 | 新增第 15 条静态测试（自读 5 个修改文件）：新类型存在 + 白名单 + M_frozen 冻结语义 | 必过 |

**归因预案**：若 G6/G7 非逐位 → 区分 (a) 目标渗漏 bug（J 相关量进入体积约束链）vs
(b) 良性浮点重排（轨迹无关）；定位后如实报告，不调参。若 G3 仍反号 → 停 + 定位 + 归因
（可能为 dJ/dphi 消费链某槽未按 §3.2 处理），不掩改。

---

## 5. 实现清单（Phase 1，最小侵入，只加新分支）

1. `createFrozenHotRegionFields.H`：L16 附近声明 `scalar thermalObjectiveMassFlowFrozen(0.0);`；
   L133-143 frozen-hot 白名单加 `&& thermalObjectiveType != "maximizeTotalHeatTransfer"`。
2. `computeObjective.H`：入口顶部新类型分支刷新 M_frozen（mag<SMALL → FatalError）；
   dJ/dT 新分支（−phi_f/(Tref·M_frozen)）；dJ/dphi 新分支
   （−(T_f − T_in)/(Tref·M_frozen)，T_in = gAverage(inlet)）。
3. `costfunction.H`：L9-41 新分支 MeanT = N − T_in·M；L71-74 新分支
   normalizedThermalObjective = −Q/(Tref·M_frozen)。
4. `evaluateCandidate.H`：L72-92 else 重构为 `else if (maximizeColdOutletTemperature)` +
   新分支（Q_cand=N−T_in·M）；L140-143 新分支 candidateNormalizedThermal
   = −Q_cand/(Tref·M_frozen)（用 M_frozen 快照）。
5. `validateCommon.H`：evaluateObjective() L16-34 新分支
   `−(mT − T_in·mFlow)/(Tref·M_frozen)`。
6. `writeCSVLog.H`：csvMeanT L36-47 新分支；csvNormalizedThermal L189-191 新分支。
7. `src/tests/test_stage_b_safety_gates.py`：第 15 条静态测试（自读 costfunction.H /
   computeObjective.H / evaluateCandidate.H / createFrozenHotRegionFields.H /
   validateCommon.H，断言新类型存在、白名单、M_frozen 冻结语义）。
8. 新案例 `/home/ys/dsH/b25_qgate`：b24_gauge 同构克隆，**仅**改
   `thermalObjectiveType maximizeTotalHeatTransfer;`（tolerance 1e-12、其余控制全同）。

---

## 6. 运行清单（Phase 2，严格顺序）

1. `wclean; wmake` → WMAKE_EXIT=0，存日志（~30-35 min 单 TU；后台 + 轮询；不用 tee）。
2. 静态测试 15/15。
3. 旧目标回归：重跑 b24_gauge（零输入改动）→ FD 三列 / ADJ_gV / gDP 因子 /
   等价性门四值逐位复现 B24（对照 §4 G5/G6/G7/G8/G11）。
4. 新目标 FD 门：跑 b25_qgate 全 BFINAL-012 协议（ε 阶梯 + 三方向），产出
   B24 同构 tsv（FD_J/ADJ_J/relJ/signJ、gDP、gV、双标签 adjoint trueRelRes、
   GRADSTABLE、P 行 identity 指纹）。
5. F1 热通道硬门通过（G9）。
6. 逐条裁决 §4 预注册预测（pass/fail/unexpected + 归因）；写 EXECUTOR_SUMMARY.md
   + 空 EXECUTOR_DONE。

---

## 7. 纪律声明

- 零拟合因子；不改任何验收阈值；frozenGradientValidated=false、mmaUpdateEnabled=false
  不动；不写 `FROZEN_GRADIENT_UNLOCK.txt`。
- 禁改清单遵守：`solveDiscreteFlowAdjointProduction.H` / `solveDiscreteFlowAdjoint.H`
  的 J/J^T 数学、`rxPressureRowTranspose.H`、filter_chainrule.H/过滤/投影、MMA、
  `validateMmaUnlockGate.H`、一切验收阈值。
- 旧类型两条路径逐位不变；旧目标回归必须逐位复现 B24 数字。
- 证据 append-only 进 `evidence/agent-group/BFINAL-025/cycle-1/`，不碰 B24 及更早。
- 任一门意外失败即停并定位报告，不调参掩盖。
- 结束：一次选择性 git 提交（风格仿 B24 f0034c0），绝不 push。

---

## 8. 修订 R1（协议修订，先于 FD 数据）—— Stage F / F1 结构性不触发的裁定

> 时间戳：2026-08-21。本修订先于全部 Phase 2 数据落盘（编译尚未完成、
> b24_gauge 回归与 b25_qgate 均未运行），属协议修订而非阈值调参。
> 授权来源：coordinator 对执行 agent 上报的「Phase 2.5 F1 门不可达」问题的裁决。

### 8.1 事实链（本修订依据，均为代码/日志证据）

1. `validateFrozenGradient.H` L17：整个 Stage F 块（F1/F2/F3）被
   `if (!frozenGradientValidated && opt >= 1 && !gradientValidated)` 守门。
2. 全部量具算例（b21_anchor、b22_gauge、b24_gauge、b24_diag、b25_qgate）
   optProperties 均 `gradientValidated true`（createFields.H L212-214 读取）。
3. 全谱系日志检索（`/home/ys/dsH/` 下全部 *.txt/*.log）：
   `Stage F: Frozen U-p-T` / `F1 (GATE4` / `medianRelErrorF1` **零命中**。
   唯一 FROZEN_GRADIENT_STATUS 行来自 Stage B2
   （validateStageB2GradientAmplitude.H L1064）。
   ⇒ **Stage F（含 F1）在本算例谱系结构性不触发，从未在任何轮次执行过**。

### 8.2 裁定

- b25_qgate **保持同构**：仅 `thermalObjectiveType` 一处不同；
  `gradientValidated` / `frozenGradientValidated` / `mmaUpdateEnabled` 一律不动。
- 任务书 Phase 2 第 5 步按**双判定线**执行（EXECUTOR_SUMMARY 与 tsv 后附表中
  两条独立呈现，不得合并）：
  - **主判定（战略判据，与 B24 同口径）**：signJ 3/3 + ADJ/FD 幅值因子 +
    ε 平台稳定性 —— 决定本轮过/不过。
  - **"F1-equivalent" 判定（第 5 步的忠实替代）**：对 b25_qgate 的 Stage B2
    全部扫描行（J_Q 的 relJ/signJ 列）施加 F1 三阈值——median relJ ≤ 1%、
    max relJ ≤ 5%、零符号失配（signJ=1 全行）——并显式注明
    「Stage F 在本谱系结构性不触发（证据：L17 守门 + 历史零命中），
    本判定为同阈值的 B2 数据替代」。
- 若 F1-equivalent 失败而主判定通过：禁止调参/改口径，两线如实并列上报，
  裁决权归审阅者。
- `FROZEN_GRADIENT_UNLOCK.txt`：B2 在 PASS 时写入算例目录属代码行为，
  不复制其内容进任何 optProperties，不提交该文件进仓库，EXECUTOR_SUMMARY
  中说明其存在即可。

### 8.3 对 §4 预注册预测的影响

- G9 口径更新为「F1-equivalent 判定」：对 B2 全部扫描行施加 F1 三阈值。
  其余 G1-G8、G10-G12 不变。
