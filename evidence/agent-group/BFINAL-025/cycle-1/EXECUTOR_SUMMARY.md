# BFINAL-025 cycle-1 — EXECUTOR_SUMMARY

任务：目标函数 Q 化（新增 `thermalObjectiveType maximizeTotalHeatTransfer`，
Q = N − T_in·M，J_Q = −Q/(Tref·M_frozen)）+ FD 门复测轮。分支
`agent/dsH-stage-b-validation`，起点 `f0034c0`。执行 2026-08-21。
授权：`NEXT_TASK_BFINAL025.md` + 协调人 2026-08-21 裁定（b25_qgate 同构、
双判定线、F1-equivalent 失败禁止调参、revision 先于 FD 数据、FROZEN_GRADIENT_UNLOCK.txt 处置）。

## 一行裁决

**J_Q FD 门双线均 FAIL。** 主判定（与 B24 同口径）：D1 signJ=0（FD 仍正号
+0.0424，ADJ 负号 −0.0358），2/3 方向符号通过；幅值 relJ 2.185/0.852/0.963
≫ 0.10。F1-equivalent 判定（R1 §8.2，对 B2 全部 13 扫描行施加 F1 三阈值）：
median 215.3% > 1%，max 220.8% > 5%，signMismatch=7 ≠ 0，三阈值全败。
修复目标的近相消病理未复现为通过——D1 的符号反号 + ~2.2× 幅值病理位于
dJ/dT（热伴随）通道，而 Q 化只改写 dJ/dphi（流伴随通道），后者沿
D1/D2/D3 对 projJ 的贡献 ~1e-10 相对量级，故 J_Q 门与旧型同态。实现、
传播、回归全链路已验证（G1/G2/G5/G6/G10/G11/G12 全过），失败为诚实的
目标语义级失败，非装配/接线 bug。按任务书判定树「J_Q 门失败 → 回到战略
菜单」，本轮不再进入 BFINAL-026 冒烟；战略裁决权归协调人/审阅者。

## 双判定线（R1 §8.2 两条独立呈现，不合并）

| 判定线 | 判据 | b25_qgate 实测 | 结论 |
|---|---|---|---|
| **主判定**（战略，与 B24 同口径） | signJ 3/3 + ADJ/FD 幅值因子 + ε 平台 | signJ D1=0/D2=1/D3=1（2/3）；bestRelJ 2.185/0.852/0.963；plateauJ 3/3 但全行 passJ=0 | **FAIL** |
| **F1-equivalent**（Stage F 结构性不触发的同阈值替代） | median relJ ≤ 1%、max relJ ≤ 5%、signMismatch=0 | median **215.29%**、max **220.80%**、signMismatch **7**（D1 全部 7 行 signJ=0） | **FAIL** |

Stage F 在本谱系结构性不触发（证据：`validateFrozenGradient.H` L17 守门 +
全部量具 `gradientValidated true` + 全谱系日志 F1 零命中；R1 §8.1）。
F1-equivalent 线为对 B2 数据的同阈值替代，已显式注明。两线同败，无
「F1 败主判定过」的裁定冲突情形；按纪律未做任何调参/口径修改。

## G1-G12 预注册预测裁决

| 门 | 预注册 | 实测 | 裁决 |
|---|---|---|---|
| G1 | wmake WMAKE_EXIT=0 | WMAKE_EXIT=0 | **PASS** |
| G2 | 静态测试 15/15 | 15/15（fresh rerun OK） | **PASS** |
| G3 | J_Q FD 符号表 3/3 signJ=1 | D1 signJ=0（FD +0.0424 vs ADJ −0.0358）；D2/D3 signJ=1 | **FAIL** |
| G4 | ADJ/FD ∈ [0.9,1.1]，relJ ≤ 0.10 | bestRelJ 2.185/0.852/0.963 | **FAIL** |
| G5 | baselineJ_Q = −0.16745040316938 | −0.1674504032463462（rel ~4.6e-10，印截断精度内：Tmix 差 4.6e-8 K） | **PASS**（精度注记） |
| G6 | gV 与 B24 逐位相同 | FD_gV / ADJ_gV 逐位相同 | **PASS**（逐位） |
| G7 | 等价性门四值逐位相同 | PD 两值逐位同；TC gradProxy 差 ~6e-11（b_TC 折叠携带新 dJ/dphi 源 68.25→108.82，预期良性） | **PARTIAL**（逐位同 + 良性解释） |
| G8 | gDP 同 B24（因子 ~2.0× 已知状态） | gDP 全部逐位同 B24；因子 D1 2.061/D2 2.007/D3 2.100 | **PASS**（如实，非失败判据） |
| G9 | F1 热通道硬门：median ≤1%、max ≤5%、signMismatch=0 | median 215.3%、max 220.8%、signMismatch 7 | **FAIL** |
| G10 | 双标签 trueRelRes ~1e-13 + GRADSTABLE | TC 9.81e-13（r1）/9.79e-13（B）；PD 9.48e-13（r1）/9.69e-13（B）；GRADSTABLE 4/4 relChange ≤ 1.85e-9 | **PASS** |
| G11 | 旧目标回归逐位复现 B24 | b24_gauge 重跑全 log diff 仅 Time/PID；stageB2 tsvs 逐位同 B24 证据 | **PASS**（逐位） |
| G12 | 新增第 15 条静态测试 | 15 项含新类型存在 + 白名单 + M_frozen 冻结语义 | **PASS** |

## Phase 2.4 定位与归因（G3/G4/G9 失败根因，不掩改）

数据链（全部来自 `/home/ys/dsH/b25_qgate`，对照 B24 证据）：

1. **projJ 是 FULL 设计梯度 dfdx 的方向投影**（validateStageB2GradientAmplitude.H
   L464-471/L591-593），dfdx = fsensMeanT（sensitivity.H L465）。
2. **BFINAL-019 装配分解**（sensitivity.H L110-280）：fsenshMeanT = momentum +
   pressureRow + fluxDirect + thermalC。state-B 实测（b25 L10409）：
   |momentum|L2=0.002267 |pressureRow|L2=0.000409 |fluxDirect|L2=0.000323
   |thermalC|L2=0.001113。
3. **dJ_Q/dphi 修复确实传播**：B0.2 state-B maxObjectiveFluxDerivative
   68.252（B24）→ 108.816（b25），+59%，且 b25 maxCombined=108.816
   （目标通量 > 内部残差 98.965），B24 maxCombined=98.965（内部主导）。
   修正仪器（AdjNS_HT.H B0.2 打印）确认槽 #18（离散外部面泛函折叠）正确消费。
4. **但 projJ 几乎未动**：b25 vs B24 ADJ_J 差 Δ = +4.19e-12（D1）、+1.49e-12
   （D2）、+1.11e-11（D3）——相对 ~1e-10。
5. **根因**：state-B 入口 M_frozen ≡ M，故 dJ_Q/dT ≡ dJ/dT 逐位（M→M_frozen
   仅替换同一快照值），热伴随 λ_T 及其下游 thermalC 不变；dJ/dphi 的修复
   只经流伴随通道（b_TC 折叠 → Ub → momentum/fluxDirect）进入 projJ，沿
   D1/D2/D3 贡献 ~1e-10 相对量级（数据 (2)(4) 一致）。
6. **结论**：D1 符号反号 + ~2.2× 幅值病理位于 dJ/dT（热伴随）通道，Q 化
   按设计不改该通道。预注册假设（近相消为根因 → G3/G4 过）被数据**证伪**。
   无接线 bug（槽 #18/#19 消费验证正确）。零调参。

## 关键数字

| 量 | 值 | 对照 |
|---|---|---|
| thermalObjectiveType | maximizeTotalHeatTransfer | B24: maximizeColdOutletTemperature |
| T_in / Tref | 600 / 600 K | 两轮相同 |
| Stage B2 baseline J_Q | −0.167450403246 | 预注册 −0.16745040316938（~4.6e-10） |
| B0.2 maxObjectiveFluxDerivative (state-B) | 108.816 | B24 68.252（+59%） |
| B0.2 maxCombinedFluxAdjoint (state-B) | 108.816 | B24 98.965（内部主导） |
| ADJ_J (D1/D2/D3) | −0.0358350400486 / −0.0249259903898 / −0.0120030122699 | B24 Δ ~1e-12 |
| FD_J D1 (ε=1e-5) | +0.042446540875 | B24 +0.043667921212（−2.8%，仍正号） |
| bestRelJ D1/D2/D3 | 2.185/0.852/0.963 | passJ=0 all |
| F1-equivalent median/max/signMismatch | 2.15292 / 2.20804 / 7 | 三阈值全败 |
| gDP / gV 列 | 逐位同 B24 | 体积链目标无关 |
| 等价性门四值 | PD 逐位同；TC gradProxy −510.7047296429894 vs B24 −510.7047296733922 | ~6e-11 |
| 双标签 trueRelRes + GRADSTABLE | TC 9.79-9.81e-13；PD 9.48-9.69e-13；GRADSTABLE 4/4 | 全 ≤ 1e-6 |
| FGMRES 迭代 | TC 1112/1108；PD 1323/1102 | 算子未变，仅伴随源变化 |
| Stage B2 状态 | FROZEN_GRADIENT_STATUS=FAIL；formalConvergenceOK=1；plateauOK=1 | 与 B24 同为 FAIL |
| ExecutionTime | 1576.43 s | B24 1396.48 s（FD 迭代更多） |

## FROZEN_GRADIENT_UNLOCK.txt

B2 判定为 FAIL，代码不写 `FROZEN_GRADIENT_UNLOCK.txt`（已确认 b25_qgate 与
b24_gauge 两算例目录均无该文件）；未将其内容复制进任何 optProperties，
未提交该文件进仓库。`frozenGradientValidated=false` / `mmaUpdateEnabled=false`
维持未动。本文件仅按 R1 §8.2 要求说明其存在/不存在的语义。

## 改动清单（本轮唯一意图变更）

- `src/computeObjective.H`：`maximizeTotalHeatTransfer` 分支——M_frozen
  快照（L36-51，冷出口 phi 求和 + FatalError 守卫）+ dJ_Q/dT
  （L86-103，−phi_f/(Tref·M_frozen)）+ dJ_Q/dphi（L152-172，
  −(T_out−T_in)/(Tref·M_frozen)）。
- `src/costfunction.H`：J_Q = −(N − T_in·M)/(Tref·M_frozen) 目标值分支。
- `src/createFrozenHotRegionFields.H` / `src/evaluateCandidate.H` /
  `src/validateCommon.H` / `src/writeCSVLog.H`：类型传播/日志。
- `src/tests/test_stage_b_safety_gates.py`：第 15 项静态测试（新类型存在 +
  白名单 + M_frozen 冻结语义）。
- 禁改清单（J/J^T 数学、rxPressureRowTranspose.H、filter_chainrule.H、
  MMA、validateMmaUnlockGate.H、验收阈值）零触碰；旧目标两路径逐位不变
  （b24_gauge 全 log 逐位复现）。

## 证据产物（cycle-1/）

- PREREGISTRATION.md（含 R1 §8.2 修订）、wmake_b25.log
- EXECUTOR_SUMMARY.md（本文件）、EXECUTOR_DONE（空标记）
- stage_b2_fd_scan.tsv、stage_b2_summary.tsv、stage_b2_adjoint_identity.tsv
  （从 b25_qgate 算例复制）
- f1_equivalent_metrics.tsv（B2 全行 × F1 三阈值）
- 原始日志：`/home/ys/dsH/b25_qgate/Log.verify_b25_qgate.txt`（81028 行，
  1576.43 s）

## 结论

J_Q FD 门双线 FAIL，诚实上报。实现与传播全链路正确（G1/G2/G5/G6/G10/
G11/G12 过，G7 良性注记），失败根因已定位：病理在 dJ/dT 热伴随通道，
Q 化（dJ/dphi 通道）不改该通道，故门结果与旧型同态。按判定树回战略菜单，
战略决策（是否改 dJ/dT 通道 / 接受旧型同态 FAIL / 换判据）归协调人与审阅者。
