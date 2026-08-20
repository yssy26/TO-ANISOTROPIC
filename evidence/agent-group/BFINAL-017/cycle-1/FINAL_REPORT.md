# BFINAL-017 cycle-1 — 第〇阶段（装配矛盾审计）FINAL REPORT

**任务：NEXT_TASK_BFINAL017.md 第〇阶段（AD1 + AD2 + 判别输出），仅审计、零修复。**
**结论摘要见 EXECUTOR_SUMMARY.md；本报告给出完整证据链。**

数据基线（三次运行共享的冻结验证态，已逐位验证）：
`b15_export/1`、`b16_states/1` 的 x/xp/xh/alpha/U/p/phi/dAlphaDxh/nuEffFrozen
全部 maxabs=0；网格与算子输入（fvSolution/fvSchemes）相同；λ_902（b12/b13/b16
写盘 Uc/pc）三运行逐位相同。

---

## 1. AD1 — 链式伴随点测试（`b17_phase0_audit.py`，日志 `b17_console.log`）

### 1.1 离线复现基础（两条路互证）

- 过滤矩阵 A = laplacian(designFilterFaceMask) − diag(b·V)（零梯度 BC → 无边界项；
  b = 1/(filterR·len/3.464)² = 5333020.444，源 = −b·V·x，V 因子由生产数据对裁定）：
  - 前向对 |A·xp + b·V·x|/|b·V·x| = **9.8e-10**（= 生产 PCG 容差）；
  - 伴随对 |A·gsensPD + b·V·gsensh|/|·| = **6.5e-10**（写盘场对）。
- 投影数据：η5 二分复现 0.748005161574（b15 日志 0.748005161603）；
  projectionEtaDenominator = −3.138599681023e-06（日志 −3.1385996818e-06）；
  drho 与 `stageB6_proj.mtx` maxabs 7.5e-9。
- 全链离线复现：L(g15_tot) vs `rxpr_prod_gsens.mtx`（生产链后导出）rel **7.1e-10**。

### 1.2 独立转置 L^T 与伴随点测试（≥3 组 (g,d)）

L = Mask ∘ FilterSolve ∘ E，其中 E(g) = drho·(mask·g) + mask·V·(1−drho)·C(g)，
C(g) = ⟨mask·dPdeta, g⟩/⟨mask·dPdeta·V⟩（filter_chainrule.H:74-161 的 numpy 重实现）。
**独立转置**：L^T(d) = E^T(FilterSolve(mask·d))，E^T(v) = drho·mask·v +
mask·dPdeta·⟨mask·V·(1−drho), v⟩/D（rank-one 项转置，未复用任何生产代码）。

| (g, d) | ⟨L(g), d⟩ | ⟨g, L^T(d)⟩ | 相对闭合 |
|---|---|---|---|
| (g15_tot, D1) | −5.623563345872e-02 | −5.623563345872e-02 | **4.9e-16** |
| (g15_tot, D2) | +2.976724315086e-01 | +2.976724315086e-01 | **1.7e-15** |
| (g15_tot, D3) | −3.634142961846e-02 | −3.634142961846e-02 | **3.1e-15** |
| (g_rand(seed 20260820), d_rand) | −7.126010077501e+01 | −7.126010077501e+01 | **4.2e-15** |
| (g_rand, D2) | +2.267023778463e+02 | +2.267023778463e+02 | **7.5e-16** |

全部 ≪ 1e-10 阈值 → **生产链是其自身公式的精确转置**。

### 1.3 线性（adaptive-eta / rank-one 检查）

- |L(2·g15_tot) − 2·L(g15_tot)| = **0.0**（严格零）；|L(2·g_rand) − 2·L(g_rand)| = **0.0**。
- 机理：η5 在 filter_x.H 由 xp（几何/设计）二分确定，**不依赖被传播场**；
  eta 修正是固定 rank-one 线性算子 → L(2g)=2L(g) 严格成立。**机理 (b) 否决**。

### 1.4 真实前向切线 z 的独立推导与验证

由 `diff.c` 的保体积约束 F(xp,η)=Σ_mask (xp−xh(xp,η))·V = 0 隐函数微分：
∂F/∂xp_j = mask_j·V_j·(1−drho_j)、∂F/∂η = −D ⟹ dη = +⟨mask·V·(1−drho)·y⟩/D，
z = drho·y + dPdeta·dη（y = filterTangent(d)）——与生产 (1−ρ) 公式**代数同一**。

数值验证：z_mine vs `stageB6_rxc_z_analytic.mtx` rel 5.6e-10/6.2e-10/7.8e-10；
z_analytic vs 设计 FD（`stageB6_rxc_xh_FD_all_eps.mtx`，eps=1e-3→1e-4）rel
4.8e-4→2.1e-5（D1）、7.2e-4→2.9e-5（D2）、7.9e-4→2.6e-5（D3），随 eps 收敛
（FD 噪声底 ~1e-4）→ **z 公式 = 真实前向切线；生产伴随链 = 生产前向链的精确转置**。

### 1.5 链环节的矛盾份额

δ_chain = 装配 − ⟨g_tot, z⟩ = **2.8e-10 / 4.0e-9 / 2.9e-9**（D1/D2/D3）——
机器/求解器精度。**B16 观察到的 20.5%/147%/4.9% 不在链上。**

---

## 2. AD2 — R_x 收缩路径逐单元对照

### 2.1 生产装配数字的坐实

⟨写盘 gsensPressureDrop, D_k⟩ = −5.4050587339 / +1.1059269805 / +14.0447142088
——复现 B12 ADJ 表 10 位（`sanity_v0.log`）；⟨写盘 gsensh, D⟩ = −5.9951/+1.0992/
+16.3348 复现 B13 P2 "chained-xh" 8 位（`b17_followup3.log`）。

### 2.2 逐项恒等式（任意 λ 均成立的代数恒等式）

- 动量行：⟨−dAlphaDxh·(U&Uc)·V, z⟩ vs −Uc·anRUd：差 **1e-13**（三方向），
  逐单元 L2 1.4–2.4e-12（`b17_phase0_results.json` AD2_localisation）。
  公式复现：−dAlphaDxh·(U&Uc15)·V vs `rxpr_prod_gsensh_momentum.mtx` maxabs 3.5e-14。
- 压力行（b15 同-pass 单-λ 恒等式，`b17_followup2.log`）：
  ⟨g15_pr, z_k⟩ vs −pc15·anRPd_k 三方向差 **1.7e-17 / 5.6e-17 / 1.0e-16**。
- 全等式（XID，同上日志）：⟨g15_tot, z_k⟩ vs −λ15ᵀrxd_k 三方向 1.2e-16/3.7e-16/5.0e-15。

### 2.3 XID 与 AD2 的配对定义——为何二者不冲突（审阅要点 1）

**XID 配对（单-λ，闭合）**：
- 左 = ⟨g15_tot, z_k⟩，其中 g15_tot = `rxpr_prod_gsensh_total.mtx` =
  b15 运行**同一 pass** 内生产的链前总收缩场 = 动量(Uc15) + 压力行(pc15 经
  rxPressureRowT)，Uc15/pc15 即该 pass 的伴随（后被导出为 `stageB6_lambda.mtx`）；
- 右 = −λ15ᵀ rxd_k，rxd_k = oracle 的 [anRUa·w; J_P·w]（**λ 无关**，
  只依赖状态），w = dAlphaDxh·z_k，z_k = oracle 解析切线（λ 无关）。
- 两侧用**同一个 λ15** 与同一状态 → 逐项恒等式以机器精度闭合。
  数字小（−0.0562/+0.2977/−0.0363）是因为 λ15 是陈旧显式解（见 §4），其收缩
  值恰小——与量具真伪无关（恒等式对任意 λ 成立）。

**AD2 配对（混-λ，矛盾）**：
- 装配侧 = ⟨写盘 gsensPressureDrop, D_k⟩ = L(g_total(λ_679))——B2 模块
  第 2 轮 pass（λ_679，warm-start FGMRES 679 迭代）的输出；
- λ-direct 侧 = −λ_902ᵀ rxd_k，λ_902 = 写盘 Uc/pc（第 1 轮，902 迭代）。
- **两侧 λ 不同** → −0.919/+0.658/−0.686 的差 = λ_679 与 λ_902 的泛函差。

即：XID 证明「同一 λ 下装配 == λ-direct」（恒等式成立）；AD2 的矛盾只是把
两个不同伴随迭代的数字放在了一起。二者不冲突，共同指向 §3 的 λ 身份机制。

### 2.4 λ 身份机制的直接证据（`b17_followup3/4.log`）

1. gradProxy(写盘 Uc) = **−205194.1411241132** == b16 日志第 1 轮
   `[GRADPROXY-FINAL pressureDrop]`（902 迭代）逐位；≠ 第 2 轮 −219717.066609。
   （gradProxy ≡ Σ(U·Uc)·D1，`solveDiscreteFlowAdjointProduction.H:1530-1546`。）
2. B13 P2 pass 内动量 ⟨−dAlphaDxh(U&Uc)·V, D⟩ = −3.5036/+1.5889/+2.2240
   ≠ 写盘 λ 公式 −3.2626/+1.4745/+2.2924 → **装配 pass 的 λ ≠ 写盘 λ**。
3. 写盘顺序：`MTO_HF.C` 主循环 costfunction → **writeOptimizationState.H:9
   `runTime.write()`** → sensitivity（第 1 轮，λ_902 已在 Uc/pc → 写盘）；
   随后 `validateStageB2GradientAmplitude.H` 的第 2 轮链（711/679）重解伴随并
   重跑 sensitivity（`finalSolvedIteration` 仍真）→ **只覆写灵敏度族文件**
   （gsensh/gsens 等），Uc/pc 文件保留 λ_902。B2State（该模块的保存/恢复结构）
   不含 Uc/pc，扫描结束的 `restoreStateInto(entryState)` 也不触碰它们。
4. b16 日志结构（行号）：3957 PD-902 → 3970 sensitivity#1 → 3991 Stage B2
   phase1 → 7380/8734 TC-711/PD-679 → 8735 sensitivity#2 → 8751 基线投影表
   → 8755 FD 扫描 → 68129 模块结束。

### 2.5 机器精度终裁（`b17_lambda_true.{py,log,json}`）

显式矩阵 M = `explicitJT.mtx`（未钉扎，467MB）SuperLU 直解：
λ*：|M·λ*−bPD|/|bPD| = **2.1e-11**；w*（J·w=bTC）4.8e-11。

```
                     D1            D2            D3
λ-direct(λ*)     −4.4859028311  +0.4478786970  +14.7302689231
λ-direct(λ_902)  −4.4859028556  +0.4478786666  +14.7302687819   (≈λ*, ~1e-8)
装配表(λ_679)     −5.4050587339  +1.1059269805  +14.0447142088   (偏离 λ*)
|M⁻¹rxd|/|rxd|    6.49e7         4.65e7         1.14e8
gradProxy(λ*) = −205194.1410 ≈ 第1轮(λ_902)；第2轮 −219717.0666 偏离 7%
```

- **λ_902（写盘）已是该泛函的机器精度解**——B16 的 λ-direct 数字**正确**；
- λ_679（装配 pass 的 warm-start 迭代）偏离真值 → B12 系装配数字携带
  λ 收敛误差；1e-9 的 FGMRES 容差 × ~1e8 的 |M⁻¹rxd| 放大 = O(1) 泛函误差，
  与 20.5%/147%/4.9% 的观察一致。

---

## 3. 生产 rxPressureRowT 与 oracle assembleWeightedRPa 的子项对照（审阅要点 2）

**裁决：两侧构造逐子项一致；不存在边界 dflux_b、drAU 边界值、uAssignable
集合或 dAlphaDxh 裁剪上的差异。** 证据：

1. **b15 同-pass 机器精度闭合**（§2.2）：生产压力行场 g15_pr（由
   rxPressureRowTranspose.H 的 T1/T2 面环构造）与 oracle 的 anRPd（由
   assembleWeightedRPa 的前向构造）在 ⟨·, z_k⟩ 收缩下差 ≤1e-16——同一 λ、
   同一状态下两条独立实现的转置恒等式严格成立。
2. **第一性原理重建**（`b17_rebuildT.py`，`b17_rebuildT.log`）：从 polyMesh
   重建 UEqn（bounded Gauss upwind 对流 + corrected 拉普拉斯 + Sp(α) +
   relax(0.4) + (−grad p − div dev2) 源）→ drAU/dHbyA/g0 → T1/T2：
   - rAU avg 2.84e-7；g0 边界 max **6.506204e-02** vs b15 日志
     **0.0650620425279**（6 位吻合）；
   - V1：T(pc15) vs b15 原始压力行导出 rel **6.3e-3**（含我的格式近似）；
   - V2：J_P·w_k vs oracle anRPd relL2 **7.6–8.4e-3**、cos≈0.9999965；
     sin 权重 7.1e-3；
   - 用 pc16 收缩：⟨−T(pc16)·dA, z⟩ = −0.0919/+0.0390/−0.1511 vs
     oracle −pc16·anRPd = −0.0917/+0.0386/−0.1468（**与 oracle 一致**，
     差 0.3–3% = 重建格式近似）。
3. **曾经观察到的"97% 压力行失配"是提取伪影**：我们从写盘 gsensh 反演
   g_total(λ_679)，再减去用**写盘 Uc(λ_902)** 算的动量公式 → 所得"压力行"
   = 真压力行 + [动量(λ_679) − 动量(λ_902)] 的混合场。变体扫描
   （`b17_variant_sweep.log`：αRel ∈ {0.3,0.4,0.5,0.7,1.0}、去 (1−α)U 项、
   去 bounded 修正，全部 rel≈0.97 不匹配；比值中位数 1.62、p10 −4.4——
   非任何标度/符号变体）排除了"重建参数差异"解释；§2.4 的 λ 身份证据
   给出正解。**生产压力行收缩（sensitivity.H:115 与 rxPressureRowTranspose.H）
   无缺陷。**
4. oracle rxc 的外部验证链（引证）：BFINAL-014 RX-A ∂R_P/∂α vs 设计 FD
   ~1e-6、cos=1（历史结论）；本轮 λ* 直解独立证实其收缩口径正确。
5. **差异峰值单元的空间特征**（`b17_spatial.log`）：5610 = (0.00525, 0.00025,
   0.00265)、6650 = (0.00525, 0.00675, 0.00265)——设计区（x∈[0.0052,0.0348]）
   **x 最小（入口侧）边缘**与**底壁/顶壁**相交的角点单元，边界相邻；
   |dAlphaDxh| = 2.769e5（全场最大），|z| 1.8–5.8。即 λ_679−λ_902 差被
   最高敏感度单元放大的位置——与任何边界项语义无关（边界项在 T1 中仅
   outlet 的 84 个可赋值面，量级由 g0_bnd 与点测约束）。

---

## 4. 附带发现：b15 导出运行的伴随是陈旧显式解（数据质量警示）

`Log.verify_b015_export.txt:1757-1762`：
`Loaded explicit solution "/home/ys/b2_case_smoke/explicitSol_pressureDrop.mtx"
relRes=9079.96` / `GMRES iterations=0`（thermalCoupling 同样 relRes=1）。
b15 的 `discreteUseExplicitSolution true` + 陈旧文件 → 该运行所有 λ 相依
产物无效：`stageB6_lambda.mtx`（λ15 与 λ16 relL2=1.0）、§9/§10 in-run
D_momentum/D_pressure/D_total（其中 prodGsen==D_momentum 的 1e-9 闭合仍是
有效的**恒等式**检验——恒等式对任意 λ 成立）。**λ 无关产物（rxc/z/dirs/
celltype/proj/dalphadxh/rpa_terms/rxa/rxb）有效**（状态逐位一致已验证）。
后续任何量具不得再用 `stageB6_lambda.mtx` 或 b15 in-run D_* 作伴随侧。

## 5. 修正后三项缺陷量具表（判别输出）

| 缺陷 | 审计后定位 | 干净量具 | 相对 B16 的变化 |
|---|---|---|---|
| ① 算子 | 不变：单一共享算子（产=导）P 行/通量耦合语义 | −λ*ᵀrxd / bPDᵀw_true = **1.768 / 0.856 / 2.274** | 无实质变化（λ_902≈λ*，1.769/0.856/2.274 → 1.7679/0.8557/2.2742） |
| ② 装配 | **重定性**：装配链与收缩逐项精确（本轮 §1/§2/§3）；B12 系装配数字的 20.5%/147%/4.9% 偏差 = λ_679 迭代误差（1e-9 容差 × ~1e8 泛函放大）+ 写盘 λ(λ_902) 与 pass λ(λ_679) 不一致的量具伪影 | 装配/λ-direct 对比必须同-λ；修复轮动作：以 λ*（或收紧 `discreteFlowAdjointTolerance`）重出装配表，并修正写盘顺序或量具协议 | 「独立的装配链缺陷」**撤销**；替换为「伴随容差不足 + 量具 λ 一致性协议」——处置位于授权层（sensitivity/量具协议），**不触发停止条款** |
| ③ 源 | 不变：b_TC 热消除折叠 | b_TCᵀw_true vs (FD_J−热介导) = 86%/70%/符号翻转 | 本轮未触及（B1 结论维持） |

**λ-direct 口径算子量具的稳定性说明**：λ-direct 泛函的伴随残差放大系数
|M⁻¹rxd|/|rxd| ≈ 4.7e7–1.1e8；λ_902（残差 ~1e-10 等效）给出与 λ* 一致到
~1e-8 的值，而 λ_679（残差 9.8e-10）偏离 20%/147%/5%。算子缺陷比值本身
（对 λ*）稳定为 1.768/0.856/2.274——**算子缺陷为真且量具有效**；但任何
未来 λ-direct 量具都应报告 |M⁻¹rxd| 加权的不确定度或直接用 λ*。

## 6. 纪律与不变量

- 零修复：`git status` 仅未跟踪的证据文件；无任何算子/源/装配/链/MMA 语义
  改动；未新增导出开关（全部用既有数据完成）；未重编译、未跑求解器。
- `frozenGradientValidated=false`、`mmaUpdateEnabled=false` 未动。
- 第一/第二阶段未启动。修复轮输入 = 本报告 §5 的修正后量具表。

## 7. 工件清单

| file | 内容 |
|---|---|
| `b17_phase0_audit.py` + `b17_console.log` + `b17_phase0_results.json` | AD1/AD2 主审计（链复现、伴随点测试、线性、z 验证、装配/λ-direct 分解） |
| `b17_followup1.log` | 支撑集/相关性/sin 权重恒等式 |
| `b17_followup2.log` | E 步对 b15 原始导出验证；b15 单-λ 压力行恒等式（1e-16） |
| `b17_followup3.log` | 写盘 gsensh == B13 pass（8 位）；写盘 Uc 公式 ≠ B13 pass 内动量 → λ 身份分裂 |
| `b17_followup4.log` | gradProxy(写盘 Uc) == 第 1 轮；b12/b13/b16 写盘 Uc/pc 逐位相同 |
| `b17_rebuildT.py` + `b17_rebuildT.log` | 第一性原理 UEqn 重建 + V1/V2 验证 + pc16 收缩 == oracle |
| `b17_variant_sweep.py` + `b17_variant_sweep.log` | 重建变体扫描（排除参数差异解释） |
| `b17_lambda_true.py` + `b17_lambda_true.log` + `b17_lambda_true.json` | λ*/w* 机器精度直解与终裁、放大系数、量具表 |
| `b17_spatial.log` | 5610/6650 空间特征 |
| `sanity_v0.log` | 写盘 gsensPD 复现 B12 表（10 位） |
