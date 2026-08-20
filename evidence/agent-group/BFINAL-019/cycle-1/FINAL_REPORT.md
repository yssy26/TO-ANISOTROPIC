# BFINAL-019 cycle-1 — FINAL REPORT（J 装配完成轮）

任务：BFINAL-018 审阅派发的完成轮——sensitivity.H 装配补 `−pb^T·R_P,x` 与
`Gx = g^T·φ_x`（DERIVATION_FIX3 §1/§5 的两个不可被 (U,p) 源吸收项），一次编译，
验证序列 1–5；J(D1) 仍反则停止并按三分量框架定位。分支
`agent/dsH-stage-b-validation`，起点 `caba940`。执行 2026-08-20。

---

## 0. 结果速览

| 项 | 结果 |
|---|---|
| 遗留改动审计 | **结论 A（方向正确，直接沿用；补静态测试第 12 项）**，见 §1 |
| 两项装配 | 完成：pb 压力行收缩 + Gx 直接项，均经门控/逐位不变性验证，见 §2 |
| 不变性门 | 全部逐位通过：四值、gDP/gV 投影与 FD 列、双标签 1e-12、GRADSTABLE 4/4、静态 12/12、WMAKE_EXIT=0、MTO_RC=0 |
| FD 门（功能口径） | **J 符号 2/3（D1 仍反：FD +0.0429 vs ADJ −0.03475）**；gDP 符号 3/3（比值 2.130/2.113/2.169 逐位复现）；gV 符号 3/3（逐位复现） |
| 停止条款 | 已触发并执行：三分量定位见 §4，无重试/无调参/无拟合 |

---

## 1. 遗留改动审计（结论 A：方向正确）

工作区有前任实例的未验证改动（`src/rxPressureRowTranspose.H` +173、
`src/sensitivity.H` +220）。逐块审计如下（判定依据 = DERIVATION_FIX3 +
B18 离线工件 + 本轮可执行验证）：

### 1.1 rxPressureRowTranspose.H 宏参数化

| 审计点 | 判定 | 依据 |
|---|---|---|
| 默认展开 token 级兼容（gDP 逐位不变的前提） | **成立** | `gcc -E -P` 新旧对照（`audit_token_compat.txt`）：仅 (i) `rxT1[celli] += ...` 一处换行（token 序列相同）、(ii) 门控诊断 Info 标签字符串 `"RxPressureRowTranspose[" "pc" "]"`（编译期拼接，仅日志文本）两处差异；数值路径零差异。运行级证据：b19_main 全日志 diff 中 dV/dDP 范围、PRODPRECGAMGCHECK、GRADSTABLE、伴随残差逐位相同 |
| pb 收缩与 pc 路径同构 | **成立** | 同一文件二次包含，同一 `{...}` 作用域、同一 LOCKED 重建（UEqn/rAU/HbyA/drAU/dHbyA/g0），仅 `RX_ADJ_PRESSURE=pb`；宏尾部 `#undef` 自清理。辅助头未定义任何全局状态，二次包含无副作用（两次包含之间无任何场变异——g 复刻只读） |
| 用了正确的伴随压力场 | **成立** | `pb` 于 `readTransportProperties.H:125` 声明、`AdjNS_HT.H:89` 作为热耦合离散伴随压力别名进入 `solveDiscreteFlowAdjointProduction.H` 求解；后续 `AdjNS_PD.H:6-7` 只用 `Uc/pc`，`AdjNS_FF.H`/`stageB7`/`costfunction` 不写 pb；MTO_HF.C 包含顺序（AdjNS_HT→…→sensitivity）保证 sensitivity 时 pb 为收敛热伴随压力 |
| Gx 的 φ_x 公式与 DERIVATION_FIX3 一致 | **成立** | φ_x[f] = Sf&(interp(dHbyA·wc)) − interp(drAU·wc)·g0[f]（内部）+ Sf_b&dHbyA[cell]·wc[cell]（U-assignable 边界）——与文件头 BFINAL-004 前向算子的 `dphiHbyA_w − dflux_w` 通道逐槽位对应、与 `b18_folding_lab.py Gx_of`（含 `am=U_assign` 边界处理）完全一致。逐面代数复核：g^T·φ_x = Σ_f g[f]·[wf((Sf&dHbyA[own])−g0·drAU[own])·w_own + (1−wf)(...nei...)·w_nei] + Σ_b g[b]·(pSf&dHbyA[cell])·w_cell，与代码 `rxFx` 累加逐项相等（精确收缩，非近似） |
| 符号 | **成立** | `dJ/dx = −λ^T R_x + C + g^T φ_x` → `gsenshMeanTPressureRow = −rxPressureRowTb·dAlphaDxh`（无 V，同 gDP 侧压力行约定）、`gsenshMeanTFluxDirect = +rxFluxDirectT·dAlphaDxh`（正号） |

### 1.2 sensitivity.H 侧

| 审计点 | 判定 | 依据 |
|---|---|---|
| g 复刻逐值正确 | **成立** | 内场面 `−mask·Tb_down·(T_n−T_o)`（mask≤0 跳过）+ 冷出口 `+dJ/dφ`，与 `AdjNS_HT.H:5-41` 的 `discreteExternalFaceFluxAdjoint` 逐行同构；依赖场（phiThermal/coldFaceMask/thermalObjectiveFluxDerivative/coldOutletPatchID/T/Tb）均在 main 作用域（createFrozenHotRegionFields.H/updateThermalFlux.H/computeObjective.H），且 AdjNS_HT→sensitivity 之间无写者 |
| 接入 fsenshMeanT 链 | **成立** | 两项在 momentum（line 65）之后、freeze 零化之前、C 项（行末 forAll）之前加入；freezeColdFlowForValidation 分支同步零化新分解场；`#include "rxPressureRowTranspose.H"` 共两次、宏定义块恰一次（静态测试断言） |
| 分解打印齐全 | **成立** | 无条件 L2 分解打印（momentum/pressureRow(pb)/fluxDirect(Gx)/thermalC）+ stageB6RxDesignOracle 门控 4 文件 mtx 导出（rxj_prod_*.mtx，与 gDP 侧同型） |
| 缺陷 | **仅一项：静态测试第 12 项未写** | 本轮补 `test_bfinal019_j_assembly_pressure_row_and_flux_direct`（断言双包含结构、宏块唯一性、两项生产公式与接线、g 复刻、freeze 零化、helper 契约与 pc 项原样保留），12/12 通过 |

**审计结论：A。** 遗留改动方向正确、与推导逐条对应，除缺静态测试外无需修正；
未发现原则性错误，故不回退。

### 1.3 与 B18 离线预测的一致性（装配正确性旁证）

`b18_decomp.json`（冻结态 λ，另一线性化状态）预测两项净贡献
D1 +0.004714 / D2 +0.005444 / D3 −0.003566；本轮 in-pass（同 λ 对照）
实测 **+0.002735 / +0.000905 / −0.000929**——三个方向净符号全部一致，
幅度差 2–4×与 B18 已实证的状态敏感性（1.8% 状态差使 b_TC^T·w(D1) 塌缩 84%）
相容。若装配有符号/槽位错误，难以在两个独立状态上重现同号净效应。

---

## 2. 实现内容（最终代码）

1. `src/rxPressureRowTranspose.H`：
   - 宏参数化（RX_ADJ_PRESSURE/RX_OUTPUT_FIELD/RX_OUTPUT_NAME/RX_LABEL，
     默认值展开 = 历史文本；`#ifndef` 守卫 + 文件尾 `#undef` 自清理）；
   - `RX_FLUX_DIRECT_SOURCE`（可选）：Gx 基场 `rxFluxDirectT`，逐槽位镜像
     前向算子设计通道（BFINAL-004-LOCKED per-cell 导数 drAU/dHbyA 与 g0 复用，
     边界通道含 U-assignable 判定与 g0/drAU 通道零边界贡献的既有论证）。
2. `src/sensitivity.H`（discrete 分支）：
   - `thermalCouplingFaceFunctional`（g 逐值复刻，nFaces 尺寸）；
   - pb 参数化二次收缩（产出 `rxPressureRowTb` 与 `rxFluxDirectT`）；
   - `gsenshMeanTMomentum`（重复 line-65 项，仅分解用，不重复相加）、
     `gsenshMeanTPressureRow = −rxPressureRowTb·dAlphaDxh`、
     `gsenshMeanTFluxDirect = +rxFluxDirectT·dAlphaDxh`；
   - `fsenshMeanT += 两项`（C 项仍在其后原样相加）；
   - freeze 零化扩展；L2 分解打印；门控 `rxj_prod_*.mtx` 导出。
3. `src/tests/test_stage_b_safety_gates.py`：+1 项（共 12）。

未触碰（授权外）：J 算子语义、两伴随求解器 rhs 折叠、
filter_chainrule/过滤/投影/MMA/目标定义、等价性门阈值、
`frozenGradientValidated=false`、`mmaUpdateEnabled=false`。

---

## 3. 验证序列（一次编译，按序执行）

| # | 步骤 | 结果 |
|---|---|---|
| 1a | 静态测试 | **12/12 OK**（unittest 实跑） |
| 1b | wmake（干净 PATH→bashrc→FOAM_USER_APPBIN→unset FOAM_SIGFPE→wclean→wmake） | `wmake_b19.log`，**WMAKE_EXIT=0**，0 error；二进制 21:51 新于全部源，含 `BFINAL-019 J assembly decomposition` 标记 |
| 2a | 等价性门四值（b19_main 单 pass，153 s，MTO_RC=0） | 两标签均 `1.13142348768e-16 / 8.25782114104e-17 / 1.42213008543e-16 / 1.34302013134e-16`——与 B18 逐位相同；全日志 diff 仅：时间戳/路径、新增 J 分解行、fsensMeanT 过滤残差（2.6403e-10→2.6341e-10）、dJ 活动范围（J 侧，允许）；**dV=[0.000198412698111, 0.000198412698776] 与 dDP=[-0.0152865454664, 0.0106732129219] 逐位相同** |
| 2b | gDP 侧投影与 B18 逐位相同（b19_fix FD 门） | projDP_D1/D2/D3 = −5.405058716169706 / +1.10592702024591 / +14.04471433326105（逐位）；ADJ_gDP 13 行逐位；FD_gDP 13 行逐位 |
| 2c | gV 逐位相同 | projV = 0.1555186511144686 / −0.08449667741658891 / 0.06396237684432686（逐位）；bestRelV 1.58353754597e-06 / 3.60173702753e-06 / 6.82658565947e-06（逐位） |
| 3 | 双标签 ≤1e-12 + GRADSTABLE | 主循环轮 TC 898 iter 9.80383374526e-13 / PD 1130 iter 9.68035507374e-13；B2 基线轮 TC 894 iter 9.82115216977e-13 / PD 855 iter 9.36148075111e-13；GRADSTABLE 4/4（relChange 5.14e-11 / 3.53e-13 / 3.65e-10 / 2.82e-12）；0 次 NOT stable；每运行仅 1 处启动 FOAM Warning（B18 同款）；MTO_RC=0 ×2 |
| 4 | FD 门（b19_fix，cp -a 自 b18_fix，1299 s，MTO_RC=0） | **J 符号 2/3**（下表）；FD 平台稳定（D1 阶梯 1e-5…1e-2 FD_J=0.0416–0.0437 恒正，与 B18 逐位相同）；formalConvergenceOK=1、plateauOK=1；FROZEN_GRADIENT_STATUS=FAIL（幅度判据，与 B18 相同，属缺陷①） |

### FD 门 J 符号表（h=0.001；FD 列与 B18 逐位相同）

| 方向 | FD_J | ADJ_J(B18) | ADJ_J(B19) | 新两项净效 | signJ(B19) |
|---|---|---|---|---|---|
| D1 | +0.042856 | −0.037489 | **−0.034754** | +0.002735 | **0 ✗** |
| D2 | −0.003469 | −0.024719 | −0.023814 | +0.000905 | 1 ✓ |
| D3 | −0.000283 | −0.012713 | −0.013641 | −0.000929 | 1 ✓ |

运行内两项分解（B2 基线态）：|momentum|L2=2.2712e-3、|pressureRow(pb)|L2=4.1494e-4、
|fluxDirect(Gx)|L2=3.2261e-4、|thermalC|L2=1.1131e-3（主循环态：3.5092e-3/4.3447e-4/2.8471e-4/1.1131e-3）。

---

## 4. J(D1) 仍反的三分量定位（按条款 5 停止，不重试不调参）

**前提**：装配层已闭合——恒等式 `dJ/dx = −λ_U^T R_U,x − λ_P^T R_P,x + C + g^T φ_x`
的四个生产项全部就位，其中三项（momentum/pressureRow/Gx）是已验证算子片段的
精确转置、C 经 B18 双重裁定无罪；λ 干净（1e-12 + GRADSTABLE）；全部不变性门
逐位通过。因此 D1 的残余差（ADJ −0.0348 vs FD +0.0429，全距 0.0776）只能由
以下三个分量解释：

1. **新两项量级（已定量，不足翻号）**：D1 净贡献 +0.00274（in-pass 同 λ 对照），
   为翻号所需（+0.0348）的 **7.9%**、全距的 3.5%。两项各自的 L2（4.15e-4/3.23e-4）
   对 momentum 项（2.27e-3）是 18%/14% 的修正。B18 §4-1「不足以独解释 0.08
   缺口」的预判在 in-pass 得到定量证实。
2. **缺陷①（算子 P 行/通量耦合语义）在 TC 源上的作用（主导分量）**：本轮 gDP
   逐位复现端到端单因子 2.130/2.113/2.169——λ 干净 + 装配完整后，该因子即算子
   切线与真实 SIMPLE 响应失配的纯表现。同一失配通过 λ_TC = J_flow,op^{-T}·b_TC
   作用于 J 链：post-B19 的 J(D1) 流部分（projJ − C ≈ −0.0348 + 0.0303 =
   −0.0044）对真实流部分（FD − C − Gx = +0.0754，B18 量具）是**符号级失配**
   （≈6% 幅度且反号），而 PD 链同缺陷仅表现为 47%（1/2.13）幅度失配且同号——
   差异来自 b_TC 泛函的近对消结构（B0.2：内面项 max 107 vs 出口项 max 114，
   同量级对消），使算子失配在 D1 上被放大而非折半。这同时解释 B18 同-pass 量具
   b_TC^T·w = +0.018 vs 目标 +0.075 的 76% 缺口（PD 侧对应缺口仅 53%）。
3. **状态敏感性（放大器，B18 已实证，本轮复现旁证）**：同一代码、同一设计，
   线性化态差 1.8%（主循环 vs B2 full-SST 基线）使 b_TC^T·w(D1) 从 +0.1146 塌缩
   到 +0.0180（−84%）；本轮两项净效在两个状态间幅度差 2–4×（符号稳定）。
   D1 是 C(−0.030) 与流部分(+0.075)的近对消方向，任何 % 级状态/算子误差都可
   翻符号——FD 平台恒正（+0.0429）而 ADJ 对状态敏感，是符号分歧的结构性土壤。

**结论**：B19 把 J 装配推到恒等式完备（其余三项+两项全部为已验证片段的精确
转置），D1 未翻号不再有装配层解释；残余符号分歧的量化归属为缺陷①作用于近对消
TC 泛函（主导，gDP 2.13× 单因子的同源表现）+ 状态敏感性（放大器）。下一轮
若授权缺陷①（算子通量图语义），J(D1) 与 gDP 幅度、b_TC 量具 76% 缺口应同时
收敛——三者共享同一通量图。

---

## 5. 工件清单

- 编译：`wmake_b19.log`（WMAKE_EXIT=0）。
- 运行：`run_b19main/Log.verify_b19main.txt`（单 pass，153 s，MTO_RC=0）、
  `run_b19fix/`（FD 门 1299 s，MTO_RC=0；Log、optProperties、
  adjointCheckpoint_*.tsv、stage_b2_{summary,fd_scan,adjoint_identity}.tsv）。
- 对比/审计：`fdgate_b19_vs_b18.md`（13 行逐位对照 + 新项净效 vs B18 离线估计）、
  `audit_token_compat.txt`（gcc -E token 级兼容证明）、
  `workspace_snapshot.txt`。
- 上游依据：`../BFINAL-018/cycle-1/`（DERIVATION_FIX3.md、b18_gauge.json、
  b18_decomp.json、REVIEW_BY_REVIEWER.md）。

## 6. 纪律与不变量

- 无拟合因子：所有实现逐条引用 DERIVATION_FIX3 §1/§5 与 b18_folding_lab.py
  的既有定量；两个新项的幅度未做任何调整。
- 未 push；git 选择性提交（代码 3 文件 + 本证据目录）。
