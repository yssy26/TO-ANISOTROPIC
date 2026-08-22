# BFINAL-026 — 热伴随链（dJ/dT 通道）定向诊断轮 — 预注册

> 分支 `agent/dsH-stage-b-validation` @ `055db8c`。本文件在**任何数据计算之前**落盘（Phase 0 纪律）。
> 判别矩阵、预期分支数字、仪器协议全部先于此文件写入；后续测量只追加证据，不改本矩阵。
> 生产数学零改动；所有新增导出开关门控默认关闭；零拟合因子。

---

## 0. 任务与背景（B25 已确立事实）

- J 的设计梯度**全部**由 dJ/dT 通道承载（dJ/dphi 通道贡献 ~1e-10 相对）。
- J 门失败：D1 反号（signJ=0）、D2/D3 幅值超配 ~6.7×/28×（signJ=1）。
- 病灶候选三段：
  - **(a) 热伴随求解**：λ_T（`AdjHeatTransfer.H` + adjointOutlet*Heat BC）——从未被外部验证过；
  - **(b) thermalC 收缩**（`sensitivity.H` BFINAL-019 分解第四项，L2=0.00111307150907）；
  - **(c) T-消除折叠**：dJ/dT → b_TC → (U,p) 行 → λ_TC → momentum/pressureRow/fluxDirect（L2 合计 ~0.0023，momentum 主导）。
- 核心机会：**源层恒等式 b_TC^T w_true vs FD_J** 自 B15 后从未用清洁仪器复测。
  B15 旧读数 2.6×/7.4×/158× 因三项仪器缺陷（A-B 线性化点错配、stale kf、离线基符号翻转）不可信。

## 1. 仪器协议（清洁量具，state-B 单点）

### 1.1 数据源（全部 B25 已落盘，不重跑）

| 项目 | 源 | 用途 |
|---|---|---|
| FD_J(Q) 真值 | `/home/ys/dsH/b25_qgate/stageB2/stage_b2_fd_scan.tsv` | M1/M4 对照（**直接取 tsv，不得重跑重定义口径**） |
| 生产 λ_T（Tb/pb/Ub） | `/home/ys/dsH/b25_qgate/1/` (Tb, pb, Ub …) | M2 生产侧对照 |
| 生产 b_TC 导出 | `/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx`（4 列/cell，按 prodMomentumScale/prodAreaScale 缩放，需反缩放） | M1 源层 |
| 全部 state-B 场 | `/home/ys/dsH/b25_qgate/1/`（52 文件：T, U, p, phi, phiThermal, Tb, alpha, dAlphaDxh, dDTDxh, x, xh, xp, designMask, coldFaceMask, …） | 离线重算基、方向、Gx、C |
| wstate 基线 | `/home/ys/dsH/b25_qgate/stageB2/wstate_baseline_*.mtx` | b_TC 反缩放辅助、状态指纹 |
| B25 日志 | `/home/ys/dsH/b25_qgate/Log.verify_b25_qgate.txt` | λ_T 收敛审计（L936-937）、B19 分解（L10409） |

### 1.2 离线基（B23 修正符号版本，全案唯一可用基）

- 使用 `/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-023/cycle-1/b23_t9_abform_corrected.py` 的 `build_basis` 逻辑。
- **SB 必须重指向 `/home/ys/dsH/b25_qgate/stageB2`**（B23 原脚本指向 b22_gauge，对 B26 状态错配，须改）。
- 值门限：`norm(mob − mob_ex)/norm(mob_ex) < 1e-12`（B23 惯例）。
- 装配数据来自 b25_qgate 的 wstate_baseline_phi/phiB/alpha/nutFrozen/nuEffFrozen/U/p + b16_mesh。
- **若基门限失败 → BLOCKED**（B23 教训：基符号翻转是全案历史读数的污染源）。

### 1.3 D1/D2/D3 方向（离线精确重建）

- 公式：`xx=Cx/Lx, yy=Cy/Ly, zz=Cz/Lz`（Lx/Ly/Lz = mesh C 分量极差）；
  - D1 = sin(2πxx)cos(πyy) + 0.5·sin(3πzz)
  - D2 = sin(4πxx) + 0.3·cos(2πxx)
  - D3 = cos(2πyy)sin(2πzz)
  - mask：designMask>0.5 之外置 0；x<0.02 或 x>0.98 置 0（物理坐标 x，用场 `x`，非 xx）；max-范数归一化。
- 重建校验：与 B25 `Stage B2 adjoint projections` 的自洽性通过 `projV`（gV 投影）交叉验证（D1 0.1555186511、D2 −0.0844966774、D3 0.0639623768，gV 场 1 或 0 结构，投影对设计 mask 敏感度低，是方向重建的独立锚）。

### 1.4 M1 恒等式与链式规则（关键口径）

- probeFD（`validateStageB2GradientAmplitude.H` L634）扰动的是**物理密度 x**（`x.primitiveFieldRef()`），所以 FD_J 是 dJ/dx（物理密度方向导）。
- 离线侧各项必须把方向 D 经 x→xh 链式规则转成 Dxh（`filter_chainrule.H`：Heaviside 投影 drho + 滤波 + designMask）：
  - `Dxh = d(xh)/dx · D`（在 state-B 处逐 cell 求 Jacobian，含体积保持的 eta5 隐式项）。
- **M1 源层恒等式**：
  ```
  FD_J(Q)·D = b_TC^T w_true + C·Dxh + Gx·Dxh
  其中 b_TC^T w_true = (J^{-1} b_TC)^T R_x·Dxh = λ_TC^T R_x·Dxh
      = momentum 段 + pressureRow 段（流介导全部）
  C = thermalDiffusionDerivativeDTCell（热直接项，离线可由 AdjHeatTransfer.H
      面公式从 Tb/T/dDTDxh/网格重算）
  Gx = rxFluxDirectT·dAlphaDxh（通量直接-alpha 项，正号）
  ```
- **w_true 获取方式（B21 惯例）**：`J w_true = −R_x·Dxh`；`M = J^T`（B23 修正基 + applyProdFlowJT 公式离线装配的稀疏矩阵）；`w = splu(M).solve(−rxd, trans='T')`；**清洁未钉扎**；闭合验证 `r_P = (M.T@w)[P0:] + rxd[P0:]`，rel-norm ~1e-11 级。若闭合 > 1e-9 → BLOCKED。

### 1.5 λ_T 求解收敛审计（M2 前置，风险项）

B25 日志证据（L936-937，已验证）：
- `Tb` DILUPBiCGStab：**final residual 6.57525213585e-13**，66 iters；fvSolution T 求解器 tolerance 1e-12、relTol 0 → 绝对收敛，**PASS**。
- 转置点积测试：**1.09689091013e-13 < 1e-10**，**PASS**。
- 结论：λ_T 收敛审计 PASS，无 BLOCKED 触发。M2/M3 对照可用 B25 生产 Tb 作为生产侧场。

## 2. 判别矩阵（先于数据）

| 分支 | 观测 | 裁决 |
|---|---|---|
| **B1** | 源层恒等式闭合：`b_TC^T w_true` ≈ `FD_J(Q) − C − Gx`，**三方向 ≤5% 量级** | 折叠 (c) 无罪；病灶在 (a) 或 (b)，进入 B3 |
| **B2** | 源层失败：至少一个方向偏差 > 5% | 病灶在 (c) 的折叠语义；按方向特征归档（哪些方向、超配还是反号） |
| **B3a** | M2 λ_T 外部锚：生产 λ_T 与直接解转置热算子 relL2 大偏差；或切线恒等式 ⟨λ_T, R_T,x·d⟩ vs 热 FD 失配 | 病灶在 (a) 热伴随求解 |
| **B3b** | M2 锚 PASS 且 M3 热伴随重算 thermalC 与生产 thermalC 偏差显著（符号/尺度/支承异常） | 病灶在 (b) thermalC 收缩 |
| **B3c** | M2 锚 PASS 且 M3 thermalC 重算与生产一致（relL2 小） | 三段全清洁 → 缺陷在更高链层/共同模态（B24 假设升级） |

**先验裁决偏好**（M4 分段投影补充）：B25 分解 momentum L2=0.00227 主导、(c) 段的动量行是 ADJ_J 主体，若 B1 失败最可能 D1（反号）由 (c) 某行符号语义或 A-B 链式规则错配引起；D2/D3 幅值超配可能由 (c) 中某行幅值通道（alphaRel/prodRAU/weight 通道）错误引起。

## 3. 预期数字（幅值量级预注册）

### 3.1 FD_J(Q) 真值（B25 tsv 直接取，主参考行）

| 方向 | FD_J(Q)（主参考） | FD_J 范围（h 扫描） | ADJ_J（生产投影） | signJ | relJ |
|---|---|---|---|---|---|
| D1 | **+0.04244654087520727** (h=1e-5) | +0.0413…+0.0433 | −0.0358350400485765 | 0（反号） | ~2.18 |
| D2 | **−0.003699503574594587** (h=1e-4) | −0.00370…−0.00353 | −0.02492599038981423 | 1 | ~0.85 |
| D3 | **−0.0004211702528400529** (h=1e-4) | −0.000421…−0.000445 | −0.01200301226991504 | 1 | ~0.96 |

主参考行选择：D1 取 h=1e-5（FD 收敛区，h=1e-5 与 1e-3 间漂移 <2%）；D2/D3 取 h=1e-4（可解析区）。所有 h 行的结果在报告中等价列出，主判据用主参考行。

### 3.2 各段量级预期（B25 分解，relJ 对照）

B19 分解（L10409）：|momentum|L2=0.002266994、|pressureRow(pb)|L2=0.000409452、|fluxDirect(Gx)|L2=0.000322614、|thermalC|L2=0.001113071。

按典型支撑（设计区 active cells ~5040，D 归一化）外推各段对 D 投影的量级：
- **momentum 段**（流介导主体，方向投影后）预期主导，D1 反号或 D2/D3 超配的主嫌疑；
- **pressureRow 段**量级约为 momentum 的 ~1/5.5（L2 比）；
- **fluxDirect(Gx)** 量级约为 momentum 的 ~1/7（L2 比）；
- **C（热直接项）**量级约为 momentum 的 ~1/2（L2 比）。
- ADJ_J 与各段和的自洽：`ADJ_J = momentum + pressureRow + fluxDirect + C`（M4 核心自洽校验）。

### 3.3 M2 预期

- 生产 Tb 与直接解转置热算子 relL2：**预期 ≤ 1e-6**（若两个算子的离散一致且求解收敛）。
  若 > 1e-4 → (a) 内部不一致强信号；若 > 1e-2 → (a) 定案。
- 切线恒等式 ⟨λ_T, R_T,x·d⟩ vs 热通道 FD（冻结 U,p,phi 只解 T）：**预期 ≤ 5%**。
  > 注：F1 型热 FD 在 B25 无独立导出（冻结 U,p,phi 的 T-only FD 未跑）；若 b26_diag 门控运行不可行（时间/崩溃风险），M2 以「λ_T 外部锚 relL2」为主判据，切线恒等式为强假设级证据（用 B25 wstate 的 T+/T− 差与冻结口径核对）。

### 3.4 M3 预期

- 离线重算 thermalC 与生产 thermalC（B19 分解打印 + sensitivity.H 同公式）：
  - 若 relL2 ≤ 1e-6 → (b) 收缩一致，病灶不在 (b)；
  - 若 relL2 ∈ (1e-6, 1e-2] → 微差，按符号/尺度/支承归档，待 B1/B3 联合判定；
  - 若 relL2 > 1e-2 → (b) 收缩缺陷，符号/尺度/支承分类。

### 3.5 自洽校验（审阅点 ⑤）

- M4 分段投影：`ADJ_J = momentum + pressureRow + fluxDirect + C`（各段之和 == 总 ADJ_J，逐方向）必须成立到 ~1e-10 相对（线性恒等式，若破则仪器不自洽，报告 BLOCKED）。
- 方向重建交叉校验：`projV`（gV 投影）重建值与 B25 tsv 一致到 ~1e-6 相对。

## 4. 风险与执行纪律

- **legacy-ILU NaN 崩溃路径**（b14/b24_diag 模板）：仅 `stageB4JacobianProbe=true → solveDiscreteFlowAdjoint.H` 触发。本轮**不设 stageB4JacobianProbe=true**；所有交付物在崩溃点前完成。如需门控导出运行，走生产 FGMRES 路径（b25_qgate 已验证无崩溃）。
- 禁改清单：两个伴随算子文件的 J/J^T 数学、AdjHeatTransfer.H 生产数学、rxPressureRowTranspose.H、filter_chainrule.H、MMA、validateMmaUnlockGate.H、一切验收阈值、两种既有目标类型路径。
- frozenGradientValidated=false、mmaUpdateEnabled=false 不动；不写解锁文件（FROZEN_GRADIENT_UNLOCK.txt 绝不创建/复制/提交）。
- 静态测试 15/15 不得回归（本轮纯诊断，若无编译运行则基线验证仍须引用 B25 的 15/15 记录并声明未触碰生产源码）。
- 证据追加式写入本目录；结束：EXECUTOR_SUMMARY.md + EXECUTOR_DONE；一次选择性 git 提交（B25 08ff950 风格），不 push。
- 编译环境按 `OPENFOAM_USAGE_AND_FILES.md`（OpenFOAM-7）；wmake 与运行日志完整保留；长跑后台/日志轮询。

## 5. 交付物清单

1. PREREGISTRATION.md（本文件，先于一切数据）
2. 离线仪器脚本（基装配、方向重建、J 稀疏矩阵、w_true 求解、C/Gx 重算、M1/M3/M4 计算）
3. M1 源层恒等式结果（三方向，数字）
4. M2 λ_T 外部锚（relL2 + 切线恒等式）
5. M3 thermalC 重算对照（符号/尺度/支承分类）
6. M4 分段投影表（自洽校验）
7. FINAL_REPORT（判别矩阵逐项裁决 + 分级：已验证/强假设/推测）
8. EXECUTOR_SUMMARY.md + EXECUTOR_DONE + 选择性 git 提交

*预注册落盘时间：2026-08-21（Phase 0 完成；此后再计算任何数据）。*
