# BFINAL-020 cycle-1 — FINAL REPORT（算子修复轮：第〇阶段停止）

任务：BFINAL-020 —— 缺陷①（P 行/通量耦合语义）修复轮。用户授权触碰
BFINAL-003 锁定算子层。分支 `agent/dsH-stage-b-validation`，起点 `3d39aa9`。
执行 2026-08-20/21。**结果：第〇阶段（推导+复测）完成后触发停止条款——
七个修复候选全部被外部仲裁否决，缺陷定位超出授权算子耦合层；零代码改动。**

---

## 0. 结果速览

| 项 | 结果 |
|---|---|
| 第〇阶段推导 | 完成（DERIVATION.md：OF7 bounded Gauss upwind 精确语义 + 不动点 Jacobi 全块） |
| 干净签名复测 | 完成：基线复现 b17 λ-direct 1.768/0.856/2.274；V1 槽位自检（kf/P 行 U 槽 1.5e-10 中位） |
| 候选修复仲裁 | **7/7 失败**（§2 表）；2 个使量具恶化 → 「实测与推导矛盾」成立 → 停止 |
| 深层诊断 | P 行切线恒等式符号级失败；φ-闭环 5–7%；U 行失衡 1.5e-4（相对主导块）；rxd FD 可信 |
| 代码/编译/门 | 零改动、零编译、无门运行（无变更可验）；无 BFINAL-003 升级声明 |

## 1. 第〇阶段方法与自检

1. **推导**（DERIVATION.md §1–4）：从 NS.H 实际不动点（bounded Gauss upwind 的
   bounded 路由、relax(0.4)、rAtU=松弛后、非正交 1 次网格正交、adjustPhi/
   constrainPressure、phi=phiHbyA−pEqn.flux）推导消去 φ 的系统 Jacobi 全块，
   逐字核实 OF7 源码：boundedConvectionScheme.C（bounded = 上风矩阵 −
   fvm::Sp(surfaceIntegrate(phi),U)）、fvmSup.C（Sp 对角 += V·sp）、
   gaussConvectionScheme.C+negSumDiag（守恒上风）、fvmDiv 边界
   （zeroGradient vic=1）、fvMatrix::flux（含 faceFluxCorrection；本网格正交为零，
   B8 P1 旁证 1.44e-8）。
2. **离线量具自检**：python 重建算子 vs 导出 J relL2 1.5e-3；P 行 U 槽位中位
   1.5e-10；rxd_U = V·U·(dAlphaDxh·z) 机器精度（2.4e-12）；基线 SuperLU 复现
   b17 全部三位比值。

## 2. 候选修复仲裁（w_true = b15 状态差分真值；门 cos>0.999/relL2<5e-2/gDP 1±0.10）

| # | 候选（推导依据） | gDP 因子 D1/D2/D3 | relL2 D1/D2/D3 | 判定 |
|---|---|---|---|---|
| 0 | 现行算子（基线） | 1.705/0.814/2.200 | 0.725/1.022/1.129 | 缺陷在 |
| 1 | D1–D4：velocityJump×闭环切线 deltaPhiFacePU + kf 反号 + Sp 对角 −contRes | 1.692/0.897/2.180 | 0.720/1.039/1.128 | 无效 |
| 2 | R3：通量图换未松弛基（去 α/去 direct）+ dev + bsrc 通道 | 2.029/0.547/1.973 | 0.855/0.992/1.016 | 恶化 |
| 3 | R2：纯未松弛基（= GatePR FD 基） | 2.138/0.745/1.985 | 0.966/1.010/1.029 | 恶化 |
| 4 | S1：松弛基 + dev + bsrc | 1.418/0.241/2.079 | 0.502/1.009/1.168 | 不均匀、未过门 |

（每个候选 = 完整 ΔJ 稀疏装配 + SuperLU 直解 + 三方向全指标；见
b20_variant_sweep.log / b20_variant_*.json / b20_derivation_check.log。）

## 3. 深层诊断（缺陷的精确定位尝试）

1. **r = J·w_true + rxd**（真 ΔJ 须满足 ΔJ·w_true = −r）：
   |r_U|/|rxd_U| = 32/68/35%，但 |r_U|/|frozenA·dU| ≈ **1.5e-4**（近对消残差，
   D1 主导块范数 1233）；|r_P|/|rxd_P| = 166/214/151%（O(1)），峰值在设计区
   角点单元 5610/6650…（与 B17 峰值一致，成对对称）。
2. **φ-闭环**：(I−Ψ_φ)⁻¹ 对 dφ 修正 relL2 4.2–7.1%（b20_augmented_check）。
3. **P 行恒等式**：div(dφ_state) = −rxd_P 应成立；实测 cos = −0.59/−0.78/−0.55
   （符号级反向），含闭环修正亦不闭合（差 9×）。
4. **rxd 可信**：rxd_P 单位方向归一化 1.087e-9 vs NS.H Rx-B FD 探针 1.152e-9
   （差 5.6%）——R_x 层无 5.6% 以上缺陷。
5. **无通道解释 r_U**：全部候选投影 cos ≤ 0.19（jump×基差 |v|=0.106 cos 0.18；
   jump×φ_α 0.4%；−U·div(φ_α) 3–4%；dev ~1%）。
6. **既有量具重解读**：b15 日志 GateOracle-B（dPhiJ vs dPhiHfd relL2 1.364，
   h 无关）+ GatePR（|J_P|/|FD_P|=0.76，cos 0.58）——P 行通量切线 O(1) 失配
   的直接历史证据；本轮证明其不可由任何单一候选解释。

**结论**：缺陷为 **J 的 P 行切线与真实系统响应之间的表示级差异**
（P 行恒等式符号级不闭合 + 换基恶化 + 通道全解释失败）。授权范围（算子耦合
补项）内无可实现的修复；继续即猜测——停止。

## 4. 停止后的建议（DERIVATION.md §8.7）

1. **w_true 精度阶梯**（零代码）：复用 stageB2 状态导出，eps 1e-2…1e-4 重测
   P 行恒等式收敛性——区分「恒等式真失败 vs 近对消噪声放大」。
2. **T1 级外部锚**：全链重收敛的 R_P 设计-FD vs J_P·w_true（现无此量具）。
3. **表示一致性审计**：算子 P 行 vs stageB6 assembleWeightedRPa 的残差表示
   是否同源一致。
4. 增广 3-块伴随（(U,p,φ) 显式化）仅当 1–3 闭合后仍指向算子层再评估
   （本轮量化其修正幅度 5–7%）。

## 5. 工件清单

- 推导与实验记录：`DERIVATION.md`（§1–7 推导稿 + §8 实验全记录与停止裁决）。
- 脚本（python3.8 + numpy1.24.4/scipy1.10.1，pip --user 引导）：
  `b20_derivation_check.py`（ΔJ(D1–D4) 仲裁 + V1 自检）、
  `b20_flux_tangent_decomp.py`（通道分解）、`b20_variant_sweep.py`（R2/R3/S1）、
  `b20_augmented_check.py`（φ-闭环 + 恒等式）。
- 数据：`b20_derivation_check.json/.log`、`b20_flux_decomp.json`、
  `b20_variant_{R2,R3,S1}.json`、`b20_augmented_check.json`、
  `b20_variant_sweep.log`。
- 上游数据（只读）：/home/ys/dsH/b8_verify_diag/explicitJT.mtx、
  b15_export/*、b16_mesh/*、b16_states/1/*、b13_probe3/stageB2/*。

## 6. 纪律与不变量

- `git status`：无源码改动（本证据目录为唯一新增）；
- 未重编译、未跑求解器；锁定开关/阈值未动；无拟合因子；
- 未 push；git 选择性提交（仅证据目录）。
