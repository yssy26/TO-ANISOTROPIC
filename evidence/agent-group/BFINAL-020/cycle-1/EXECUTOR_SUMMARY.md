# BFINAL-020 cycle-1 — EXECUTOR SUMMARY（结论先行）

**结论：按停止条款停止。第〇阶段（推导 + 干净签名复测 + 七候选离线仲裁）完成并
带来决定性新证据：缺陷①不在本轮授权的算子耦合层内。零代码改动、零编译、
无 BFINAL-003 升级声明（无语义变更落地）。**

## 三个决定性事实

1. **七个第一性原理修复候选全部被 w_true 外部仲裁否决**（SuperLU 直解 +
   b15 真值链，基线复现 b17 的 λ-direct 1.768/0.856/2.274 ✓）：

   | 候选 | gDP 因子 D1/D2/D3（门 1±0.10；旧 1.70/0.81/2.20） | relL2 D1/D2/D3（门 <5e-2） |
   |---|---|---|
   | D1–D4：U 行 bounded-upwind 闭环切线 + kf 反号 + Sp 对角 | 1.69/0.90/2.18 | 0.72/1.04/1.13（不变） |
   | R3：未松弛基 + dev + bsrc | 2.03/0.55/1.97（恶化） | 0.86/0.99/1.02（恶化） |
   | R2：纯未松弛基（GatePR FD 基） | 2.14/0.74/1.99（恶化） | 0.97/1.01/1.03（恶化） |
   | S1：松弛基 + dev + bsrc | 1.42/0.24/2.08（不均匀） | 0.50/1.01/1.17 |

2. **深层诊断把缺陷逼到「表示级」**：r = J·w_true + rxd 显示 U 行失衡仅为主导块
   （|frozenA·dU|≈1233）的 1.5e-4，而 **P 行切线恒等式 div(dφ_true)=−rxd_P
   符号级失败**（cos −0.59/−0.78/−0.55）；φ-闭环 (I−Ψ_φ)⁻¹ 仅 5–7%；
   rxd 自身经 Rx-B FD 探针独立验证可信（单位方向归一化差 5.6%）。
   即：算子的 P 行切线与 w_true 所响应的真实系统之间，差的不只是可枚举的耦合项。

3. **既有量具重解读**：b15 运行日志中沉睡的 GateOracle/GatePR 早已记录了
   P 行通量切线的 O(1) 失配（dPhiJ vs dPhiHfd relL2=1.364，h 无关系统性）——
   本轮首次把它与 λ-direct 1.77/0.86/2.27、w′ relP 集中连接成完整证据链，
   并证明它**不能**由松弛基/未松弛基/附加通道/闭环中任何一个解释。

## 推导贡献（保留资产）

- OF7 源码级语义确认：`bounded Gauss upwind` = 守恒上风矩阵 − Sp(div φ)·U
  （boundedConvectionScheme.C/fvmSup.C）；U 行 φ 反馈净结构 = jump·dphi 仅落
  下风行（矩阵与 Sp 在上风行精确抵消；边界行全零）；negSumDiag/flux()/
  zeroGradient valueInternalCoeffs=1 等逐字核实（DERIVATION.md §2–3）。
- 完整可复用的离线仲裁工具链（b20_derivation_check / b20_flux_tangent_decomp /
  b20_variant_sweep / b20_augmented_check，全部 python3.8+numpy/scipy，
  pip --user 引导记录在案）。

## 停止依据与下一步

任务条款「若实测与推导矛盾，停下报告」「宁可报告不要猜」：候选使量具恶化 =
矛盾成立。建议下一轮（DERIVATION.md §8.7）：(1) w_true 精度阶梯重测 P 行恒等式
（区分「真失败 vs 噪声放大」，零代码）；(2) T1 级全链重收敛 R_P FD 外部锚；
(3) J 与 stageB6 rxc 的残差表示一致性审计。本轮不触碰源码。

## 纪律

- `git status` 无源码改动；全部产物 = 离线脚本 + 既有导出数据；
- 未重编译、未跑求解器；`frozenGradientValidated=false`、`mmaUpdateEnabled=false` 未动；
- 证据（脚本/日志/json/推导）全部在 `evidence/agent-group/BFINAL-020/cycle-1/`；
- 无拟合因子：所有候选均为推导量的直接组装，LU 求解为 SuperLU 机器精度。
