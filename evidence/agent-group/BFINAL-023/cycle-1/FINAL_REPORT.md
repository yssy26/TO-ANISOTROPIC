# BFINAL-023 cycle-1 — FINAL_REPORT（SLOT-5 专项轮）

## 0. 执行概况

预注册主假设 H-CLAMP + 退出坡道；零生产改动，全程离线（B22 导出
+ b16 网格）。四份预注册文件先于各自数据；判别表见 VERDICT_TABLE.md。
检验链条：T0 复现控制 → T1/T2/T3 主假设 → H5a 值门驱动 → T5/T6/T7/
T8/T9 扩展假设（每步先预注册判别实验再跑数据）。

## 1. H-CLAMP 裁决：否证

结构性前提在源码层成立：生产算子重组动量矩阵（`fvm::div(phi,U) −
fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U)`）后**不调用
.relax()**，`prodRAU = 1/A_unrelaxed`（solveDiscreteFlowAdjoint
Production.H:144-173；模块版 solveDiscreteFlowAdjoint.H:22-81 同），
P 行通量切线 U 通道 = `alphaRel*prodRAU`、kf 通道 = 精确钳制导出
`primalPressureMobility`（B22 W1 后）。

但语义推演+数据否证了「钳制携带 D2」：

1. **真矩阵 0/33600 胞钳制绑定**。用生产符号+对称插值 γ_f 重建的
   基把精确导出 mobility 复现到 relL2 = 2.19e-15（逐胞，全区域
   frac>1% = 0）——该基下 sumOff 处处 ≤ |D|（占优余量见 json），
   fvMatrix::relax 的钳制分支永不触发。
2. 因此真 rAtU = alphaRel·V/|D0| 处处成立，生产 `alphaRel·rAU_u`
   **恰为真 mobility 因子**；(1−alphaRel) 直通系数精确；da 通道
   drAU = −mob²·da/alphaRel、dHsrc = (1/alphaRel−1)·da·U 精确。
   六变体（含全钳制语义 VC 与只修 (b)/(c) 的 VC1/VC2、生产镜像
   VP）数字完全相同（D2 relL2 0.4534, cos 0.9134）。
3. 审阅者引用的三组指纹重合（12.9% 差集中 nut 活跃钳制胞 / D2
   流向集中 / 固侧更差）由**错符号离线基的伪钳制结构**统一解释：
   b16 谱系仪器把扩散装配成 diag −= kfv、offdiag += kfv（生产为
   diag += kfv_f、offdiag −= kfv_f、γ_f 对称插值），在 nut 活跃胞
   人为制造非占优（14.1% 伪绑定）。**H-CLAMP 的经验基础（钳制绑定
   分布）本身是量具伪影。**

裁决依据：T1 钳制集空 + T3 六变体恒等 + T9 复核（修正基下 D1/D3
闭合见 §4）。预注册阈值「VC 改善 ≥1.5×」远未达到（改善 = 0）。

## 2. W-B（12.9% 根因）：证实 = 离线仪器扩散符号翻转

- 修正基（生产符号 + 对称 γ_f）→ mob_exact 重建 2.19e-15。
- 旧基 12.93%（T0 逐位复现 B22）。固体区两符号不可见（α·V 支配），
  偏差集中于 nut 活跃胞——与 B22 取证分布一致。
- 连带修正：B15-B22 全部以旧基给出的**绝对**读数继承该缺陷；
  相对结论（同基比较）大多不变（如 D2 44.9→45.3%）。

## 3. W-C 通道分解与扩展假设链（每步预注册）

| 假设 | 裁决 | 关键数据 |
|---|---|---|
| H5a 值门 | 未过（3.0%→2.0% 渐近） | 误差 99.7% x-法向、80% 边界邻、73% 流侧；div 比对 1.008（重分布型） |
| H5b 矩阵运动携带 D2 | 否证 | B-form 0.478 > A-form 0.312；corr(motion,residA)=+0.70；α 驱动 0.203，φ 驱动 0.030 |
| H5c 固定矩阵通道归因 | 完成 | 残差∥du（proj 1.016）；LOO 去 da：0.454→0.312（da 净有害） |
| H5d 上风施主翻转 | 否证 | 载体面翻转 0.000 |
| H6 显式偏应力槽 | 否证 | 通道 |dev|/|tot|≈0.004；锚不动；生产本有该槽（applyProdFrozenDeviatoricJT） |
| **H7 phiHbyA 压力通道** | **证实真缺失槽，部分携带** | 值 3.01→2.60%（反号 3.80% 恶化）；D2 45.34→35.68；corr(修正,残差) −0.996/−0.879/−0.948；\|通道\|/\|残差\|≈0.21 |
| H8 松弛源边界后段 | 值级部分；切线正交 | 值 2.60→2.00%；D2 0.2216→0.2224；corr 0.001 |
| T9 终局 A/B-form | D1/D3 闭合；D2 不闭合 | D1 2.77%、D3 2.40%（B-form）；D2 22.12%（A-form）→39.48%（+运动，α 驱动 39.6） |

H7 机理（生产+离线同缺）：构造 `== −fvc::grad(p)` 使 UEqn 源含
−∇p·V 一次（solve 的第二个只进临时拷贝）→ HbyA = U + rAtU·∇p
（NS.H:1039 互证）→ 真 phi 含 **+interp(rAtU·∇p)·Sf**，即相对
「offline 基映射」缺 **−interp(mob·∇p)·Sf** 通道；生产 J_P 的
kf 通道（pEqn.flux 转置）不含它。

## 4. 面锚终局与系统级含意

| 配置 | D1 | D2 | D3 | 值门 |
|---|---|---|---|---|
| B22 旧仪器 | 8.68% | 44.90% | 8.99% | — |
| 修正基（=生产语义） | 3.32% | 45.34% | 3.35% | 3.01% |
| +H7 | 2.74% | 35.68% | 2.79% | 2.60% |
| +H7+H8（du+dp 值 FD / +运动） | 2.77% | 22.12 / 39.48% | 2.40% | 2.00% |

- **D1/D3 面级闭合**（2.4-3.3% ≈ 值门 2.0% 水平）。E3 的「仪器
  地板 3.4%」与新值门一致；旧 8.7/9.0% 大头为仪器伪影。
- **D2 为唯一大方向特异缺陷**：残余 22.1%（du+dp 型），且「映射
  设计偏导（da 通道 / α 驱动矩阵响应）在 D2 净有害、在 D1/D3 必需」
  ——机制未明，是本轮留下的核心未解事实。
- 系统级 O(1)（gDP 2.15×、J 符号 D1 反）与面级的新关系：D1/D3
  面级已 2-3%，O(1) 放大不能来自 D1/D3 面级缺陷；候选 = D2 集中
  缺陷经 div/压力反演的放大（B21 已有机制），或 J 链更高层。本轮
  未重跑 w* 闭包链（预算原因），列为下一轮首选。

## 5. 下一轮建议（供用户战略决策）

A. **实现轮（需授权）**：生产 P 行通量切线补 H7 压力通道（及
   J^T 折叠、边界通道），量级面级 D2 −9.7pp、D1/D3 −0.6pp；
   同时以修正仪器重跑 w* 闭包链与 FD 门，重估 gDP 因子与 J 符号
   的归因。风险：H7 单独不闭合 D2（22.1% 残余）。
B. **退出坡道**：D2 残余无已识别可修机制；记录「D2 方向特异的
   通量切线设计响应缺陷，机制未明」并转入用户战略决策。

## 6. 产物清单

PREREGISTRATION.md、PREREGISTRATION_T5_ADDENDUM.md、
PREREGISTRATION_T6_ADDENDUM.md、PREREGISTRATION_T7_ADDENDUM.md
（含 T8 节）、VERDICT_TABLE.md、b23_hclamp_gauge.{py,log,json}、
b23_t5_matrix_motion.{py,log,json}、b23_t6_deviatoric.{py,log,json}、
b23_h5a_forensics.log、b23_t7_phihbya_pressure.{py,log,json}、
b23_t7b_full_channel.{py,log,json}、b23_t8_relaxsrc_bnd.{py,log,json}、
b23_t9_abform_corrected.{py,log,json}、sha256_artifacts.txt、
EXECUTOR_SUMMARY.md、本报告、EXECUTOR_DONE（空）。
