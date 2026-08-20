# BFINAL-023 cycle-1 — PREREGISTRATION ADDENDUM T7（第三扩展，先于 T7 数据）

## 已落盘事实（T6 与 H5a 取证）

- H6 否证：偏应力通道 |devchan|/|tot| ≈ 0.004，加与不加 relL2 不动。
- H5a 取证（b23_h5a_forensics.log）：网格完全正交（1−n·d̂ = 0）、
  weights≈0.5；3.0% 值差能量 99.7% 在 x-法向面、80.5% 边界邻胞、
  73% 流侧 owner；|div(重建)|/|div(真值)| = 1.008（非整体失衡，
  是重分布型误差）；top-2000 误差面 100% x-法向、40% 边界邻、
  中位 |phi| 比全场大 3 个量级（大通量面）。
- 边界 ic 全部均匀分量（fixedValue: 梯度 ic = −δ|Sf|γ·(1,1,1)、
  对流 ic = 0；outlet: (phi_b,phi_b,phi_b)）→ H() 边界对角项严格
  相消，mobility 机器精度验证与此一致。H() 边界项不是缺口。
- 源码事实：NS.H:28-33 构造 `... == −fvc::grad(p)`（operator==
  语义 source += V·su，fvMatrix.C:1401-1411）→ UEqn 自身源含
  −∇p·V 一次；solve(UEqn == −fvc::grad(p)) 的第二个 −∇p 只进
  临时拷贝（求解用），不留在 UEqn.H()。因此 NS.H:114 的
  HbyA = rAU·UEqn.H() 含 −rAU·∇p，phiHbyA = fvc::flux(HbyA) 含
  **−interp(rAU·∇p)·Sf**；NS.H:1032 phi = phiHbyA − pEqn.flux()。
  **离线仪器与生产算子的通量切线都只建 −kf(p_n−p_o)（pEqn.flux
  通道），均缺 −interp(rAtU·∇p)·Sf 通道。**

## H7（本轮第三扩展假设）

**H7：通量映射缺 phiHbyA 压力通道 −interp(rAtU·∇p)·Sf（值级与
切线级，生产与离线同缺）。该项携带 H5a 的 3.0% 值差的大部分与
D2 面锚失配的相当部分。**

∇p 取 Gauss linear（正交网格，内部面 ∇p_c 由 Green-Gauss；
边界 p：outlet fixedValue 0（p 出口 BC——从 case 0/p 确认）、
壁面/inlet zeroGradient → 胞值）。

### 判别实验
- **H7a（值门）**：加该通道后值级 internal relL2 3.0% → **<1.5%**
  判「携带大部分」；1.5–2.5% 记「部分」；≥2.5% 否证。
- **H7b（D2 主判）**：切线加 dp-通道扩展 −interp(rAtU·∇dp)·Sf
  （及 du/da 不变；kf 通道不变）后 D2 relL2 45.3% → **<20%** 判
  「携带大部分」（≥2× 改善）；20–35% 记「部分」；≥35% 否证。
  D1/D3 同报。方向自证：corr(修正量, 旧残差) 显著为负。
- **H7c（符号判别）**：同时跑 +/− 两种符号的值级测试——只有正确
  符号能降值差（防符号误判）。

### 判读
- H7 证实 → **SLOT-5 候选机制成立且为生产算子真缺陷**（J 的 P 行
  通量切线缺 phiHbyA 压力通道）；实现轮规格 = 生产 P 行（与镜像
  转置）补 −interp(rAtU·∇(·))·Sf 通道（值级 HbyA 源 −∇p 与其
  p-切线），需用户授权，本轮不实现。
- 部分携带 → 报告剩余量，评估下一缺槽（无新数据不扩展）。
- 否证 → 退出坡道。

无拟合因子；边界 p 值取自 case BC。

## T7 落盘结果与 T8（附注，先于 T8 数据）

T7 数据（b23_t7_phihbya_pressure.log / b23_t7b_full_channel.log）：
- H7a：3.01% → 2.60%（正确符号；+号变体 3.80% 恶化——T7 脚本输出
  标签 minus/plus 恰好写反，代码 extra 语义：extra=+1 为推导的
  HbyA −= mob∇p 通道）。按预注册阈值 2.60% ≥ 2.5% 记「否证边缘/
  部分携带」，但 corr(修正,残差) = −0.996/−0.879/−0.948（三方向
  全改善、方向正确）→ 通道为真、必要、非主导（|通道|/|残差| ≈ 0.21）。
- H7b：D2 45.3% → 35.7%（仅 dp 扩展）/ 38.2%（含 da 扩展全线性化）；
  LOO 去掉 da 通道 → 22.2%。da 通道（map 的设计偏导）在 D2 上仍
  净有害。D1/D3 → 2.68%/2.45%。
- 值残差定位：加通道后 2.60%，**98.1% 边界邻胞、98.1% 流侧 owner、
  99.8% x-法向**——剩余值差集中在边界邻流胞。

### T8（下一缺槽评估，由 T7 残差定位数据授权）

**H8：relax 源系数的边界后段缺失。** OF7 relax() 中源补偿用
S += (D − D0)·psi，其中 D 为 relax 内部对角 = clamp/alphaRel −
cmptMin(ic)（Ufixed 墙：+ab；outlet：−phi_b；ab = nuEff·|Sf|·δ），
而 D0 为纯内部对角。离线/生产的松弛源通道用 D_rel = clamp/alphaRel
（无边界后段）。边界邻胞上差 = bnd_post·U（值）与 bnd_post·dU（切线
du 通道），恰在 T7 残差定位处。
判别：H8a 值门 2.60% → **<1.5%**（否则部分/否证）；H8b D2（du+dp+
press 配置）22.2% → **<12%** 记携带过半；D1/D3 同报；方向自证
corr 为负。若 H8a/H8b 达标 → 与 H7 通道合并为「生产算子 P 行通量
切线缺失两槽：phiHbyA 压力通道 + 松弛源边界后段」，实现轮规格
据此立；仍未达标 → 按退出坡道收束（本轮不新开假设）。
