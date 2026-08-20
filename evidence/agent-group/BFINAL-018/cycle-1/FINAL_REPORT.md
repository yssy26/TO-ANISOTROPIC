# BFINAL-018 cycle-1 — FINAL REPORT（合并修复轮：缺陷②+③）

任务：NEXT_TASK 派发书（BFINAL-017 审阅通过后的合并修复轮）。分支
`agent/dsH-stage-b-validation`，起点 `61abbb0`。执行 2026-08-20。

---

## 0. 结果速览

| 项 | 结果 |
|---|---|
| 修复②（工程） | **完成**：容差 1e-12、GRADSTABLE 泛函稳定判据（默认开）、B2 λ 身份协议、门控 rhs 导出 |
| 修复③（源） | **完成（按数学正确口径）**：b_TC 折叠 = 算子一致通量图的精确转置；实现经独立块回归验证 |
| J 算子/预条件/过滤/MMA | **零改动**（等价性门四值 1e-16 逐位不变；`frozenGradientValidated=false`、`mmaUpdateEnabled=false` 未动；无任何因子拟合） |
| FD 门（功能口径） | gDP 符号 3/3 ✓（单因子 2.130/2.113/2.169）；J 符号 2/3（**D1 仍反**）；gV 与 B12 逐位一致 |
| b_TC 量具 | 未闭合（76%/758%/150%），定位：缺陷① + 装配缺口 + 泛函病态（详见 §4） |

---

## 1. 修复② —— 工程层

### 1.1 实现
1. **容差**：`discreteFlowAdjointTolerance` 代码默认 1e-9→1e-12
   （`solveDiscreteFlowAdjointProduction.H:1499`、`solveDiscreteFlowAdjoint.H`）；
   算例侧同步；maxIter 4000 保持（实测 855–1130 充裕）。
2. **泛函稳定判据**（防 λ_679 类假收敛）：FGMRES 收敛后须最后两个重启
   cycle 的 gradProxy 相对变化 ≤ `discreteProdGradientStabilityTolerance`
   （默认 1e-6）方准出；开关 `discreteProdGradientStabilityCheck` 默认开；
   不满足则 Warning 并继续迭代至 maxIter；循环出口加二次 Warning。
   BCGS 路径不受影响（文档注明判据按任务定义只作用于 FGMRES）。
3. **B2 λ 一致性协议**：`validateStageB2GradientAmplitude.H` 在 B2 基线
   伴随解后写 `stageB2/stage_b2_adjoint_identity.tsv`（两标签 gradProxy
   指纹），并加协议注释：该模块重解两伴随并重跑 sensitivity → 盘上
   灵敏度族属**当轮** λ；离线量具只能同-pass 收缩，禁止与主循环写盘
   Uc/pc/Ub/pb 混用（λ_679/λ_902 混用路径就此关闭）。
4. **门控导出** `stageB18RhsExport`（默认关）：生产路径 rhs 物理单位写盘
   （`b18rhs_<label>.mtx`，4 列/单元），供量具复算。

### 1.2 验证
- 静态：**11/11 OK**（原 9 项 + 新增 2 项：`test_bfinal018_tolerance_and_
  functional_stability_gate`、`test_bfinal018_source_folding_matches_
  operator_flux_map`）。
- 编译：wclean+wmake 全量，`WMAKE_EXIT=0`（wmake_b18.log）；二进制含
  `GRADSTABLE`/`B18RHSEXPORT` 标记，新于全部改动源。
- **等价性门逐位不变**：两标签×两轮 PRODPRECGAMGCHECK
  interior/boundary/effective/offDiag 全部 1.1–1.5e-16。
- **双标签双轮 1e-12**：thermalCoupling 898/894 迭代（9.80e-13/9.82e-13）、
  pressureDrop 1130/855 迭代（9.68e-13/9.36e-13）；**GRADSTABLE 4/4**
  （relChange 5.1e-11/3.5e-13/3.6e-10/2.8e-12）；全程 0 次「NOT stable」
  Warning；`discreteProdConvergeFatal false` 保持，MTO_RC=0。

### 1.3 λ 身份的修正认识（重要）
两轮 gradProxy **不**一致且**持久**：PD 主循环 −205194.141009（=B17
SuperLU 真值 λ* 的 −205194.1410）vs B2 第二轮 −219717.066521（与 b16
时代 −219717.066609 同源）。1e-12+稳定判据下仍相差 7.1% ⟹ **两轮解的
是两个不同线性化系统**：B2 模块 phase-1 的 full-SST 再基线改变了流态
（主循环 J=−1.14640 / B2 基线 J=−1.16745，b16/b18 两轮逐位复现）。
B17 的「λ_679 收敛不足」叙述据此精化：**当时 20–147% 泛函差 = 状态差
+ 收敛差的混合**；容差收紧消除了后者，状态差按 B2 门设计保留（FD 探针
与 ADJ 同在 B2 基线态，门内部自洽）。该认识已写入协议注释与 §4。

---

## 2. 修复③ —— b_TC 折叠（推导全文见 `DERIVATION_FIX3.md`）

### 2.1 定位：缺陷在路由，不在 g 公式
- g 公式（内部面 `−mask·Tb_down·(T_n−T_o)`、出口 `+dJ/dφ`、入口跳过）
  是 B0.2 以 1e-10 逐面 FD 验证过的精确转置——**未改**。
- 旧路由是「`dφ = interp(dU)&Sf + kf(p_n−p_o)`」的转置：(i) 缺
  `αrel·rAU·dH` 的 H^T 非局部块；(ii) 缺 (1−αrel) 因子；(iii) **内部 kf
  符号与算子自身 J_PP 块相反**。三者合计即 B16 量具的 86%/70%/符号翻转。
- 新路由 = 算子 P 行通量切线（BFINAL-003）的逐槽位转置：
  `dφ_f = Sf&(αrel(w·rAU_o·dH_o+(1−w)rAU_n·dH_n)) + (1−αrel)Sf&(w dU_o+
  (1−w)dU_n) − kf_f(dp_n−dp_o)` → hA(αrel·rAU_u)+deltaH^T（含边界 H 对角
  块）+ direct((1−αrel)) + kf(+own/−nei) + 出口 kf_b。生产/诊断两文件
  同构；删除旧 pressure-flux-correction-transpose 辅助。

### 2.2 离线定量（`b18_folding_lab.py`，15 部分）
- **重建验证**：几何 1.4e-19；旧 b_TC 复现导出 rhs 2.4e-5；热算子
  A_T·T−src−Q = 2.7e-7；A_T^T·Tb−dJdT = 3.4e-3。
- **C 真值裁定**（不经生产公式）：直接解 `A_T δT = −(A_x dx)T` 得
  C_true = −0.0303/−0.0083/−0.0190 ⟹ 生产 `thermalDiffusionDerivativeDTCell`
  **无罪**（与我的公式转写一致；B13 in-pass 探针 D 配对 0.4% 互证）。
  期间澄清：b15 写盘 fsenshMeanT 含陈旧导入 λ_TC 的动量项（b15 算例 1/
  为 cp 残留、Ub=uniform(0) 与内存不符），不可用于反演 C。
- **恒等式口径**：FD_J = A+B+C（A=出口通量显式、B=T 对通量响应、C=直接
  热项）；T-消除源的正确量具目标为 **FD−C−Gx**（Gx=通量直接 α 项，
  g^T·φ_x，离线精确定量 −0.0021/+0.0024/+0.0076）。B16 预注册比较子
  FD−thermal_gauge(=A) 只能被「仅出口项折叠」满足——生产不可实现。
- **冻结态（主循环态）路由对照**：旧 b_TC^T w = (+0.1125, −0.0256,
  −0.0448)；numpy V1 = (+0.0531, −0.0011, +0.0389)（D3 符号修复、D1
  误差减半）；目标 FD−C−Gx = (+0.0753, +0.0023, +0.0119)。
- **实现保真**：单 pass 主循环运行导出的 b18rhs vs numpy V1 块回归：
  U 块系数 hA=0.9947 / direct=1.0032（残差 6.2% = 我方重建输入精度，
  非正交修正等）；P 块 kf=1.0010 ⟹ **C++ 折叠即 V1 公式**。

### 2.3 同-pass 量具（b18 运行，第 2 轮数据）
| 方向 | b_TC^T w_true | FD−C−Gx | 误差 |
|---|---|---|---|
| D1 | +0.01797 | +0.07536 | 76% |
| D2 | −0.01487 | +0.00226 | 758%（近对消区） |
| D3 | −0.00596 | +0.01190 | 150% |

未闭合。三个已定量的贡献因子：
1. **泛函病态**：同一代码、同一设计、线性化态差 1.8%（主循环 vs B2
   full-SST 基线）使 b_TC^T w(D1) 从 +0.1146 → +0.0180（−84%）——
   b_TC 是流/热近对消结构（B16 已示），对状态极度敏感；
2. **缺陷①**：算子 P 行/通量耦合语义（本轮禁改）——gDP 端到端 2.13×
   单因子即其同源表现；
3. **装配缺口**（sensitivity.H，授权外）：J 装配缺 `−pb^T·R_P,x` 与
   `Gx` 直接项（拉格朗日推导证明不可被源向量吸收）。

---

## 3. FD 门（功能口径，完整三方向三目标）

FD 平台稳定（D1 阶梯 1e-5…1e-2 FD_J = 0.0416–0.0437 恒正）。
所有 FD 点有效（formalConvergenceOK=1、plateauOK=1）。

| 方向 | FD_J | ADJ_J | signJ | FD_gDP | ADJ_gDP | signDP | 比值 | gV |
|---|---|---|---|---|---|---|---|---|
| D1 | +0.042856 | −0.037489 | **0 ✗** | −2.53746 | −5.40506 | 1 ✓ | 2.130 | sign ✓, rel 1.6e-6 |
| D2 | −0.003469 | −0.024719 | 1 ✓ | +0.52341 | +1.10593 | 1 ✓ | 2.113 | sign ✓, rel 3.6e-6 |
| D3 | −0.000283 | −0.012713 | 1 ✓（b16 为 0 ✗，**修复**） | +6.47722 | +14.04471 | 1 ✓ | 2.169 | sign ✓, rel 6.8e-6 |

- **gDP：符号 3/3，幅度= 算子单因子 2.130/2.113/2.169**（±1.3% 内
  近常数）——λ 干净后 20–147% 的方向依赖散布消失，剩余为缺陷①纯算子
  因子（与任务预期「符号全对即本轮合格；幅度残留属①」一致；注意
  B16 的 1.77/0.86/2.27 是 λ-direct 口径，与端到端口径不同）。
- **J：D3 符号修复（+0.0032→−0.0127 对 FD −0.000283），D2 保持；D1
  仍反**（FD 平台恒正 +0.043，ADJ −0.0375）。定位见 §4。
- **gV：与 B12 逐位相同**（bestRelV 1.58e-6/3.60e-6/6.83e-6，符号 3/3）
  ——零回归。
- 正式门 FROZEN_GRADIENT_STATUS=FAIL（幅度判据未过——预期内，幅度属①）。

## 4. J(D1) 未恢复的定位（按任务条款停在此，不猜）

1. **装配缺口（下一轮授权候选）**：J 的 xh 装配
   `fsenshMeanT = −dAlphaDxh(U&Ub)V + C` 相对完整恒等式
   `dJ/dx = −λ_U^T R_U,x − λ_P^T R_P,x + C + Gx` 缺 `−pb^T R_P,x`
   （rxPressureRowTranspose.H 目前只对 pc 常驻）与 `Gx`（两项均
   sensitivity.H 层，本轮授权范围 AdjNS_HT/AdjHeatTransfer 之外）。
   冻结态 λ 粗估 pb 项量级 0.003–0.011、Gx −0.002…+0.008（精确值
   Gx 已定量）——不足以独解释 D1 的 0.08 缺口。
2. **缺陷①的 TC 源表现**：算子切线 w_op 与 w_true 的失配（gDP 2.13×
   即其表现）同样作用于 TC 源泛函；b_TC^T w_op 与 b_TC^T w_true 可差
   数倍。
3. **泛函病态**（§2.3-1）：D1 的 J 是 C(−0.030) 与流部分(+0.075 目标)
   的近对消，任何 % 级状态/算子误差都可翻符号。

## 5. 纪律与不变量

- 未触碰：J 算子（applyProdFlowJT/预条件/两诊断算子段）、
  filter_chainrule/过滤/投影/MMA/目标与约束定义、等价性门 1e-8 阈值、
  `frozenGradientValidated=false`、`mmaUpdateEnabled=false`、
  validateMmaUnlockGate；无任何因子拟合（所有数字来自推导或既有量具）。
- 静态 11/11；WMAKE_EXIT=0；两运行 MTO_RC=0；未 push。
- 工件清单见 `EXECUTOR_SUMMARY.md`；推导全文 `DERIVATION_FIX3.md`；
  实验室日志 `b18_folding_lab.log`（parts 1–15）、量具 `b18_gauge.json`。

## 6. 给下一轮（缺陷①轮）的输入

1. 算子通量图语义是 b_TC 量具与 gDP 2.13× 的共同嫌疑；建议先做
   「算子 P 行 vs 真值通量响应」的直接量具（需导出 FD φ 场——本轮
   证实 φ 不在写盘集合中，需加开关）。
2. sensitivity.H 装配补 `−pb^T R_P,x` 与 `Gx`（或等价的 g 泛函直接项）。
3. B2 门的两态设计（主循环态 vs full-SST 基线态）建议显式文档化或统一，
   消除 λ 身份类混淆的土壤（本轮已加 identity 文件缓解）。
