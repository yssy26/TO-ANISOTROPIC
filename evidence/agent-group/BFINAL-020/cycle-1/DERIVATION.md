# BFINAL-020 cycle-1 — 第〇阶段推导：真实不动点 Jacobi 的 U 行通量耦合语义

日期：2026-08-20。基线 `3d39aa9`。授权：用户显式授权触碰 BFINAL-003 锁定算子层（本派发即凭证）。
范围：缺陷①（P 行/通量耦合语义）——同时解释 gDP 2.13× 幅度与 J(D1) 反号的最后一块数学拼图。

## 0. 方法与证据来源

从 NS.H 实际前向不动点出发，逐行核对 OpenFOAM-7 源码的精确装配语义
（`boundedConvectionScheme.C`、`gaussConvectionScheme.C`、`fvmSup.C`、
`lduMatrixOperations.C::negSumDiag`、`fvMatrix::flux`、
`zeroGradientFvPatchField::valueInternalCoeffs`），推导消去 phi 后的完整
固定点残差及其 Jacobi，逐块对照现行 J（`solveDiscreteFlowAdjointProduction.H`
的 `applyProdFlowJT`、`solveDiscreteFlowAdjoint.H` 的
`applyDiscreteFlowJ`/`applyDiscreteFlowJT`/CSR `addCoo`）。

## 1. 前向不动点（NS.H 逐行事实，含算例配置）

算例（b19_main 系）配置事实：
- `divSchemes: div(phi,U) bounded Gauss upwind`（fvSchemes:32）——**bounded 关键字实际生效**；
- `SIMPLE` 无 `consistent` → `rAtU = rAU`（松弛后 1/A）；`nNonOrthogonalCorrectors 1`；
- U：inlet/hotInlet/hotOutlet/墙 fixedValue（noSlip），outlet zeroGradient（唯一 assignable）；
- p：outlet fixedValue 0，其余 zeroGradient（→ `p.needReference()=false`，无规范行）；
- 松弛：U 0.4（α_rel=0.4）、p 0.2；fvOptions 空；网格 33600 hex（B8 P1 验证 J_PP kf-only
  与实际 pEqn.flux() FD 1.44e-8 闭合 → 网格正交，非正交修正为零）。

每 corrector（收敛不动点处）：
1. `UEqn = fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U)`，RHS `−fvc::grad(p)`；
2. `UEqn.relax()` → `solve(UEqn == −fvc::grad(p))`；
3. `rAU = 1/UEqn.A()`（松弛后），`HbyA = rAU·UEqn.H()`（constrainHbyA 钉扎固定 U 边界），
   `phiHbyA = fvc::flux(HbyA)` + adjustPhi；
4. 非正交循环 `pEqn: fvm::laplacian(rAtU,p) == div(phiHbyA)`（2 遍）；最终
   `phi = phiHbyA − pEqn.flux()`；
5. `p.relax(); U = HbyA − rAtU·grad(p)`。

不动点消去 phi：装配时的 phi（上一 corrector 末端）在不动点等于闭环通量
**phi\*(U,p) = phiHbyA(U) − pEqnFlux(p)**（phiHbyA 经松弛映射
`HbyA = α_rel·rAU_u·H_u`，含 relax-source 通道）。

## 2. 关键源码语义（逐字核实）

**(a) bounded Gauss upwind**（boundedConvectionScheme.C:76-81）：
```
fvmDiv(bounded) = gaussUpwind.fvmDiv(phi,U) − fvm::Sp(fvc::surfaceIntegrate(phi), U)
```

**(b) gauss upwind 矩阵**（gaussConvectionScheme.C:92-97）：`lower = −w·phi`，
`upper = lower + phi`，`w = pos0(phi)`；`negSumDiag`（lduMatrixOperations.C:71-72）：
`Diag[lowerAddr] -= Lower; Diag[upperAddr] -= Upper`。行约定
`(Ax)_own = diag·x_own + upper·x_nei`。⇒ 矩阵部分**恰为守恒上风**：
phi≥0 时 own 行 `+phi·U_own`、nei 行 `−phi·U_own`；phi<0 时对称。

**(c) fvm::Sp(sp,U)**（fvmSup.C:113-118）：`diag += V·sp`（无源项）。
⇒ bounded 修正是**纯对角**项 `−V·div(phi)·U_P`。

**(d) zeroGradient 边界**（zeroGradientFvPatchField.C:107-131）：
`valueInternalCoeffs = 1`，`valueBoundaryCoeffs = 0`（与 pw 无关）；
fixedValue：0 / 1。⇒ zeroGradient U 边界的对流向
`internalCoeffs_b = phi_b`（系数 1）。

## 3. 固定点残差与逐块 Jacobi

残差定义（与 rxd/rxc 量具一致）：
```
R_U(P) = [A_upw(phi*)·U]_P − V_P·div(phi*)_P·U_P − [L·U]_P + α_P·V_P·U_P + V_P·grad(p)_P [+RANS dev]
R_P(P) = Σ_{f∋P} s_{P,f}·phi*_f
```
切线通量（**不变**，BFINAL-003/008 已定）：
```
dphi*_f = Sf&(α·rAU_u·dH + (1−α)·dUbar) − kf·(dp_n − dp_o)      ≡ deltaPhiFacePU
```

### 3.1 U 行对 phi 的反馈（缺陷①本体）

对每个内部面 f（own=o, nei=n, phi_f≥0 记 down=n；否则 down=o；jump = U_n − U_o）：

| 来源 | ∂R_U(o)/∂phi_f | ∂R_U(n)/∂phi_f |
|---|---|---|
| 守恒上风矩阵 | +U_up | −U_up |
| bounded Sp（−V·div(phi)·U） | −U_o·(+1) | +U_n·(−1)·(−1)=+U_n |
| **合计（phi≥0）** | U_o−U_o = **0** | −U_o+U_n = **jump** |
| **合计（phi<0）** | U_n−U_o = **jump** | −U_n+U_n = **0** |

即：**净反馈 = (U_n − U_o)·dphi\*_f，只落在下风单元行**（等价紧凑式：贡献
= phi_f·(U_P − U_up)，对下风行 = phi_f·jump，对上风行 = 0）。

**边界行：精确为零。** zeroGradient-U（outlet，唯一 assignable）：矩阵 +U_c·dphi_b
（vic=1）与 Sp −U_c·dphi_b 抵消；fixedValue-U 补丁：vic=0 且
dphi_b=0（fixedFluxPressure/zeroGradient-p 边界 phi_b 状态无关；inlet 通量钉扎）。
⇒ 边界无任何新 U 行项。

### 3.2 U 行对 U 的对角（Sp 通道二阶项）

`∂[−V·div(phi*)·U_P]/∂U_P = −V·div(phi*)·dU_P`：对角系数 = −contRes_P
（= Σ±phi_f + Σphi_b，收敛处 ~ 连续性残差，量级可忽略但可精确包含）。

### 3.3 P 行

`dR_P(P)/dx = Σ s·dphi*_f`——现行算子的 relaxed-SIMPLE 闭合
（J_PU/J_PP，BFINAL-003/008 外部验证：B8 P1 relL2 1.44e-8）**完整，无需改动**。
J_PP（diag +kf、off −kf、outlet 边界 +kf_b）不动。

## 4. 对照现行 J：缺陷清单（唯一缺陷的三层表现）

现行（内部面，forward J）：
```
output[down,c] += jump_c · deltaPhiFace_old(w),
deltaPhiFace_old = Sf&rAUdH + kf·(p_n − p_o)        // 旧、未松弛、kf 符号相反
```
正确：
```
output[down,c] += jump_c · deltaPhiFacePU(w)          // 与 P 行同一闭环切线
```

| # | 差异 | 推导依据 | 量级 |
|---|---|---|---|
| D1 | 通量图：hA 部分缺 α（0.4）缩放 | §3.1（dphi* 的 hA 部分带 α） | 同项量级 ×0.6 |
| D2 | 通量图：缺 (1−α) relax-source 直接部分 | §3.1 | 同项量级 ×0.6 |
| D3 | 通量图：kf 部分符号相反（+kf vs −kf） | §3.1（phi* = phiHbyA − pEqnFlux） | 该子项反号 |
| D4 | Sp 对角 −contRes_P·dU_P 缺失 | §3.2 | ~连续性残差（≤1e-8 相对，含入但预期不可见） |

**转置语义**（J^T，面伴随标量）：
```
φ_f = lambdaPdiff + convTerm        // convTerm = λ_down & jump（现行定义不变）
路由：hA += α·φ_f·w·rAU·Sf（现行 convTerm 无 α）；direct += (1−α)·φ_f·w·Sf（现行仅 lambdaPdiff）；
      kf：P(o) += kf·φ_f；P(n) −= kf·φ_f（现行 kf·(lambdaPdiff − convTerm)）
```
注意 kf 路由中 convTerm 的符号由 −（现行）翻为 +（正确），且 hA/direct 路由对
convTerm 补 α/(1−α)。

**P 行转置（含边界 hB/direct）不变**；源折叠（B18 fix 3）不变（外部面泛函只乘
dphi\*/dx 本身 = P 行通量图，不经过 λ_down）；J_PP 不变 ⇒ 等价性门四值
（PRODPRECGAMGCHECK）逐位保持。

## 5. 与既有实测的一致性（签名预测）

1. **B8 stageB8 dir1（dp-only）U 行 relL2=0.123**：冻结-phi 残差 FD 无任何 phi 反馈，
   而 Jv 含旧 velocityJump·kf·dpdiff → 12% 伪内容。此锚在冻结语义下本来就测不了
   耦合项——与新推导一致（该 FD 不是耦合项的真值）。
2. **B15 P 行集中失配（relP 0.72–1.13 > relU 0.18–0.53）**：w′=J⁻¹(−rxd) 是耦合解，
   U 行缺陷经 J_PU^T（hA/kf 路由）污染 λ_P；B15 的「P 行集中」不排斥 U 行成因
   （B15 证伪的是「U 行集中」的预测，不是 U 行成因本身）。
3. **方向依赖比值 1.768/0.856/2.274**（λ\* 口径）：三项修正（α 缺失、relax-source
   缺失、kf 反号）都是结构性（非标量）修正，方向依赖的残余谱与之一致；
   λ-direct 比值非纯标量（1.3% 跨方向离散，B14 已证）。
4. **J(D1) 反号**：D1 是近对消方向，旧项 kf 子项反号贡献 O(主信号) 量级
   （B19 测得新装配两项净贡献 +0.00274 仅 7.9%——那是装配层项，本轮改的是算子层
   不同的、大得多的通道）。

## 6. 验证前置：离线 ΔJ 预验证（先于任何代码改动）

用 b17 第一性原理重建（rAU/HbyA，已对 oracle 0.6–3% 验证）+ b16 mesh 导出 +
b15 w_true 状态差分，构造 ΔJ（=J_new − J_old）稀疏矩阵，SuperLU 直解：

- **V1（重建自检）**：kf 重建 vs explicitJT P-P offdiag 槽位逐面机器精度；
  P 行槽位 M[P(o),U(down,c)] = −kf·jump_c ∓ w·sf_c 逐面对照；
- **V2（决定性）**：w′_new vs w_true 三方向 cos/relL2/relU/relP
  （门：cos>0.999、relL2<5e-2；旧值 0.72–1.13）；
- **V3（gDP 因子预测）**：bPD^T w′_new / bPD^T w_true → 1±0.10（目标 1±0.05）；
- **V4（ADJ 侧预测）**：λ\*_new = M_new⁻¹ bPD 的 λ-direct vs 真值。

若 V2/V3 不收敛 → 实测与推导矛盾 → 停止报告（任务条款）。

## 7. 结论（第〇阶段推导稿，已被实验部分修正——见 §8）

缺陷① = **U 行对流反馈的通量图语义**：bounded Gauss upwind 的 U 行 phi 反馈净结构
（jump·dphi\* 仅落下行）与现行 velocityJump 结构一致，但现行乘的是
**旧未松弛、kf 反号的 deltaPhiFace**，而非 P 行所用的闭环切线 deltaPhiFacePU。
修正 = 三个实现点（生产 JT、诊断 J/JT、CSR 导出）统一把该反馈乘子换成
deltaPhiFacePU + 补 −contRes 对角（D4）。J_PP/J_PU/边界/源折叠全部不动。

## 8. 第〇阶段实验记录（2026-08-20 夜，全部离线、零代码改动）

### 8.1 ΔJ(D1–D4) 离线预验证（b20_derivation_check.py/.log/.json）

构造 ΔJ（velocityJump 换闭环切线 + kf 反号修复 + Sp 对角），SuperLU 直解
（基线 1352 s 复现 b17 的 λ-direct 1.7679/0.8557/2.2742 ✓）：

| 方向 | w′ relL2 旧→新 | gDP 因子 旧→新 |
|---|---|---|
| D1 | 0.7247 → 0.7198 | 1.7049 → 1.6923 |
| D2 | 1.0223 → 1.0393 | 0.8142 → 0.8966 |
| D3 | 1.1292 → 1.1284 | 2.2003 → 2.1802 |

**D1–D4 不是主导缺陷**——修正后失配基本不动。推导的 U 行反馈通道本身正确但量级不足。

### 8.2 既有量具的重新解读（关键转折）

`b15_export` 运行日志中沉睡着 GateOracle/GatePR 量具（stageB2 gated）：
- **A：dPhiHfd vs dPhiFD = 1.65e-11**（未松弛 H-FD 与 phi-FD 机器精度互证）；
- **B：dPhiJ vs dPhiHfd = relL2 1.3638, cos 1.0587**（h 无关、系统性）；
- GatePR P 行：|J_P|=2.80e-6 vs |FD_P|=3.69e-6，relL2 0.858，cos 0.580。

即：**现行 P 行通量切线（BFINAL-003 松弛基 Candidate B）与未松弛重建 FD 切线
相差 relL2 1.36 / cos 0.61（|J|/|FD|≈1.72）**——这才是与 λ-direct 1.77/0.86/2.27
及 w′ 失配 P 行集中（relP 0.72–1.13 ≫ relU）同量级、同形状的缺陷。
BFINAL-003 当时的 G1/G2 裁决中 G2（1.06e-11）是自比（重构与 FD 同基），
G1（2.99e-5 平台）只验证了松弛基在自身基内的闭合。

### 8.3 通道分解（b20_flux_tangent_decomp.py/.json）

Python 复算（算子模型 vs 导出 J：relL2 1.47e-3 ✓；重建含 dev/bsrc/relax-src），
对「算子 − 重建真值」的通道量化（方向 dUdir）：

| 通道 | \|v\| | 相对失配占比 |
|---|---|---|
| basis（α 缩放 hA + (1−α) direct vs 未松弛） | 1.732e-5 | ~100%（即 GateOracle-B 的全部） |
| bsrc（outlet 零梯度 U 边界源 ∂(bc·U_b)/∂U） | 3.207e-7 | ~1.8%（relL2 ~0.025） |
| dev（RANS 偏应力源导数） | 3.801e-9 | 可忽略 |

### 8.4 变体扫描终局（b20_variant_sweep.py/.log/.json，w_true 仲裁，全部完成）

| 变体 | relL2 D1/D2/D3（旧 0.72/1.02/1.13） | gDP 因子 D1/D2/D3（旧 1.70/0.81/2.20） | 判定 |
|---|---|---|---|
| D1–D4（U 行闭环切线 + kf 反号 + Sp 对角） | 0.720/1.039/1.128 | 1.692/0.897/2.180 | 无效（~1%） |
| R3（未松弛基 + dev + bsrc） | 0.855/0.992/1.016 | 2.029/0.547/1.973 | **更差** |
| R2（纯未松弛基 = GatePR FD 基） | 0.966/1.010/1.029 | 2.138/0.745/1.985 | **更差** |
| S1（松弛基 + dev + bsrc） | 0.502/1.009/1.168 | 1.418/0.241/2.079 | 部分改善、极端不均匀，未过门 |

门（cos>0.999、relL2<5e-2、gDP∈1±0.10）：**全部变体失败**；两个「换基」变体使
量具恶化——测量与推导候选矛盾成立。

### 8.5 深层诊断（b20_augmented_check.py/.json/.log + 增补计算）

1. **r = J·w_true + rxd**（真 ΔJ 必须精确满足 ΔJ·w_true = −r）：
   - U 行失衡 |r_U|/|rxd_U| = 32/68/35%（D1/D2/D3），但相对主导块
     |frozenA·dU|（≈1233，D1）仅 **1.5e-4**——精巧近对消的残差；
   - P 行失衡 |r_P|/|rxd_P| = 166/214/151%，峰值集中于设计区角点单元
     （5610/6650 等，与 B17 峰值单元一致，成对对称）。
2. **φ-闭环量化**：dφ 的 (I−Ψ_φ)⁻¹ 修正 relL2 仅 **4.2–7.1%**——非 O(1) 主因。
3. **增广 P 行恒等式失败**：div(dφ_state) = −rxd_P 应成立，实测
   cos(div, −rxd_P) = **−0.59/−0.78/−0.55**（符号级反向），|div|≈1.2×|rxd_P|；
   含闭环修正的 rxd_P_c（14% 幅度，cos 0.93）也不能闭合（差 9×）。
4. **rxd 自身可信**：rxd_U = V·U·(dAlphaDxh·z) 机器精度复现（2.4e-12）；
   rxd_P 单位方向归一化后与 NS.H Rx-B FD 探针一致（1.087e-9 vs 1.152e-9/unit，差 5.6%）。
5. **无通道解释 r_U**：全部候选（jump×基差 |v|=0.106 cos 0.18、jump×φ_α 0.4%、
   −U·div(φ_α) 3–4%、dev 1%）投影 cos ≤ 0.19。

### 8.6 最终结论（停止裁决）

**缺陷①不（只）位于本轮授权的算子耦合层内。** 七个第一性原理候选（含 OF7 源码
级语义支撑的 U 行 bounded-upwind 通道、松弛/未松弛基、dev/bsrc、φ-闭环）经 w_true
外部仲裁全部失败，其中两个使量具恶化。深层诊断显示：

- 算子 U 行与真值在相对意义下已经接近（失衡 1.5e-4 相对主导块）；
- P 行切线恒等式对 (J_P, rxd_P) 对在 w_true 上**符号级不闭合**，且 rxd_P 本身
  经独立 FD 验证可信——矛盾指向 **J 的 P 行切线与 w_true 所响应的真实系统之间
  存在表示级（representation-level）差异**，其可能根源包括：w_true 状态差分对
  P 行恒等式的精度不足（近对消放大）、J 与 R_x 残差表示的不完全一致性、或
  更深的未建模系统语义（如 SIMPLE 外环收敛判据引入的映射差异）。

按任务条款（「若实测与推导矛盾，停下报告」「宁可报告不要猜」）：
**第〇阶段后停止；零代码改动；无 BFINAL-003 升级声明（无语义变更落地）。**

### 8.7 对下一轮的建议（按信息价值排序）

1. **提高 w_true 精度再仲裁**（最高优先）：用更大的状态导出 eps 阶梯（1e-2…1e-4）
   直接测 P 行恒等式 div(dφ) = −rxd_P 的收敛性——区分「恒等式真失败」vs
   「w_true 噪声放大」。零代码（复用 stageB2 开关），一次运行。
2. **T1 级外部锚定 J_P 行**：对 P 行做 x±h·d 的**全链重收敛**（不是冻结态重建）
   的 R_P FD，与 J_P·w_true 直接对照（现无此量具）。
3. **表示一致性审计**：J（算子）与 stageB6 rxc（R_x oracle）是否为同一残差
   表示的导数——逐项对照 assembleWeightedRPa 与算子 P 行的前向作用
   （两者都声称是「flux-map 切线」但结构不同源）。
4. 若 1–3 闭合后缺陷仍在算子层，再评估增强 3-块伴随（(U,p,φ) 显式化，
   (I−Ψ_φ)⁻¹ 隐式化）的可行性（本轮已量化其修正仅 5–7%，优先级最低）。
