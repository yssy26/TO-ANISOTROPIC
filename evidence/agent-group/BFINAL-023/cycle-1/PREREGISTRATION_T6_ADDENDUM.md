# BFINAL-023 cycle-1 — PREREGISTRATION ADDENDUM T6（第二扩展，先于 T6 数据）

## 已落盘事实（T5，b23_t5_matrix_motion.log）

- H5a **未过门**：修正基值级重建 internal relL2 = 3.0%（outlet 29%，
  p99 相对 3.8）——mobility 机器精确但 phi 映射模型仍缺件。
- H5b **否证**：B-form（含矩阵运动）D2 relL2 0.478 > A-form 0.312；
  |motion|/|truth| = 0.204 且 corr(motion, residA) = +0.70（同向，
  加上更糟）。矩阵运动由 α 驱动（0.203/0.204，即 da 通道的值形式），
  φ 驱动仅 0.030。D1/D3 的 B-form 略好（corr −0.69/−0.75，量级 5-6%）。
- H5c：D2 残差与 du 通道投影 1.016（几乎平行）；leave-one-out：
  去掉 da 通道 0.454→**0.312**（da 通道在 D2 上净有害）；du 通道
  必需（去掉 0.95）；dp 通道小有益。
- H5d：上风施主翻转不在载体上（carrier flip 0.000）——施主翻转
  不是 D2 机制。

## 新事实（源码层，T6 基础）

生产算子的动量-H 转置（hA 路径，`solveDiscreteFlowAdjointProduction.H`
applyProdFlowJT 内，line ~671）以 `applyProdFrozenDeviatoricJT` 收尾
——即生产 P 行通量切线的 dH 通道**包含**显式冻结偏应力响应
（NS.H:39-43 `UEqn -= fvc::div(nuEffFrozen*dev2(T(fvc::grad(U))))`，
源项 S_phys 线性于 U）。**b16 谱系离线仪器（B15-B22 与本轮 T3/T5）
从未建模该项**——H5c 的「残差平行于 du 通道」与之直接相容。

## H6（本轮第二扩展假设）

**H6：离线（及历轮）通量切线与值映射缺失显式偏应力项
S_phys = −fvc::div(nuEff·dev2(T(grad U)))，该项携带 D2 面锚的
大部分失配（以及 H5a 的 3% 值差）。**

离线模型（Gauss linear，与 fvSchemes 一致）：
- G_c = (1/V)Σ_f Sf·U_f（内部 w 插值；壁面/inlet/hotInlet/hotOutlet
  边界 U_b 取 BC 值：墙 0、inlet (193.97,0,0)；outlet zeroGradient
  → U_c）；
- X_c = nuEff_c·dev2(T(G_c))；X_f = w·X_o+(1−w)·X_n（内部）；
- divDev_c = (1/V)[Σ_int Sf·X_f + Σ_bnd Sf_b·X_b]，边界 X_b 取
  相邻胞值（calculated）；墙面单侧梯度变体作括号报告；
- 值映射：HbyA += −mob·divDev；切线：dHbyA += alphaRel·rAU_u·
  (−divDev(dU))（与生产 hA 路由同因子）。

### 判别实验
- **H6a（值门）**：加偏应力后 H5a 值级 internal relL2 从 3.0% 降到
  <1%（判「过门」）；1–2% 记「部分」；>2% 记「未过」。
- **H6b（D2 主判）**：加偏应力的切线（生产语义 VP+dev）在 D2 面锚
  relL2 从 45.3% 降到 **<15%** 判「携带大部分」（≥1.5× 以上改善且
  过半能量被解释）；15–30% 记「部分携带」；>30% 记「否证」。
  D1/D3 同报（预计从 3.3% 降）。
- **H6c（结构自证）**：偏应力修正量（VP+dev − VP）与 T5 残差
  （VP − truth）的相关应显著为负（修正方向对）；|修正|/|残差|
  报告。若相关为正而 relL2 改善——记「部分携带（方向对但幅度/
  结构不全）」并如实报告。

### 判读
- H6 证实（H6b <15%）→ SLOT-5 结论改写：**生产算子本身含该槽位，
  历轮 D2「45% 真信号」的相当部分是离线仪器缺槽伪影**；生产算子
  面级误差需用修正仪器重测给出（本轮给出修正读数），系统级 O(1)
  失配（gDP 2.15×/D1 符号）的归因在此新仪器下重新表述；实现轮
  候选 = 无（生产无变更）或按新读数另立。
- H6 部分携带 → 报告剩余量与通道归因，评估是否还有单一缺槽
  （无新数据不扩展）。
- H6 否证 → 退出坡道：如实报告无明确可修机制。

无拟合因子；边界 BC 值取自 case 文件与 bnd 表。
