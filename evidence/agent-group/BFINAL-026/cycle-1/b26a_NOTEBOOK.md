# BFINAL-026 尝试 2 — b26a NOTEBOOK（执行磁盘笔记本，每小步追加）

> 任务：M1 源层恒等式测量（纯离线 Python，零编译、零求解器运行、零生产代码改动）。
> 仓库 `/home/ys/dsH/TO-ANISOTROPIC`（即 `/home/ys/dsH` 下真实目录；`/home/ys/TO-ANISOTROPIC` 为无关空壳，勿用）。
> 分支 `agent/dsH-stage-b-validation`，HEAD `055db8c03c25b439b3ed5dcad864afcdf3d1a7cf`。
> 权威任务书：`NEXT_TASK_BFINAL026_ATTEMPT2.md`。判别矩阵原样引用见 `b26a_PREREG_CONFIRM.md`。

## 0. 时间轴

- 2026-08-21（继续会话）：核对全部输入；写本 NOTEBOOK + PREREG_CONFIRM；建仪器；跑 M1。

## 1. 输入核验表（全部通过）

| 输入 | 路径 | 核验结果 |
|---|---|---|
| w_true(H7) | `evidence/agent-group/BFINAL-026/cycle-1/b26_wstar_h7.npz` | D1/D2/D3 各 (134400,) float64 全有限；maxabs 295844/573932/714579 |
| 验收记录 | `b26_wstar_h7.json` | \|r_P\|/\|rxdP\| = 1.037e-11 / 9.649e-12 / 2.428e-11（PASS，闭合 <1e-9 无 BLOCKED）；relL2 vs w_state = 0.6849/0.9023/1.0893；cos ≈ 0.99 |
| 求解日志 | `b26_wstar_h7.log` | M nnz 6395416（H7）；MMD 失败 → COLAMD 回退；耗时 1006.2s；协调人执行，勿重跑 |
| w_true(pre-H7) 对照 | `evidence/agent-group/BFINAL-021/cycle-1/b21_wstar_cache.npz` | (134400,) 全有限；b21_wstar_clean.json：闭合 1.76e-11/1.25e-11/2.99e-11；relL2 0.7593/1.0244/1.1722 |
| b_TC（Q 型，B 态） | `/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx` | 33600 行 × 4 列，纯空白分隔（Ux Uy Uz P 每 cell）；**PHYSICAL/未缩放**，NO descaling（stageB18RhsExport L870-898） |
| FD_J(Q) 真值 | `/home/ys/dsH/b25_qgate/stageB2/stage_b2_fd_scan.tsv` | 14 行（7 D1 + 4 D2 + 3 D3），勿重跑 |
| FD 判据行（任务书 h=1e-3） | 同上 tsv | D1 +0.04252469623637622 (signJ=0)；D2 −0.003528972517971574 (signJ=1)；D3 −0.0004449511895182612 (signJ=1) |
| FD 主参考行（PREREG §3.1） | 同上 tsv | D1 +0.04244654087520727 (h=1e-5)；D2 −0.003699503574594587 (h=1e-4)；D3 −0.0004211702528400529 (h=1e-4) |
| FD 平台稳定性 | 同上 tsv | D1 h 扫描 +0.0413…+0.0433（1e-5..1e-2）；D2 −0.00353…−0.00370；D3 −0.0004148…−0.0004450。h=1e-3 vs 主参考漂移：D1 +0.18%、D2 −4.6%、D3 +5.6% |
| ADJ_J（生产投影） | 同上 tsv | D1 −0.0358350400485765；D2 −0.02492599038981423；D3 −0.01200301226991504 |
| M4 分段导出 | `/home/ys/dsH/b8_verify_diag/` | rxpr_prod_gsensh_{momentum,pressurerow,total}.mtx（各 33600，nnz 5040）；stageB6_rxc_z_analytic.mtx（100800=3×33600，maxabs 8.17）；stageB6_dirs.mtx；stageB6_dalphadxh.mtx |
| state-B 场 | `/home/ys/dsH/b25_qgate/1/` | T,U,p,phi,phiThermal,Tb,alpha,dAlphaDxh,dDTDxh,x,xh,xp,designMask,coldFaceMask 等（52 文件） |
| b8_verify_diag == b25_qgate | — | 全部 List 场 bit-identical（maxdiff=0，relL2=0），见尝试 1 NOTEBOOK |

## 2. 符号推导与链式规则（M1 核心口径，含符号约定）

约定：`J` 为雅可比（（U,p）残差对状态），`M = J^T`（显式装配的 explicitJT 稀疏阵）。
`w_true` 求解（b21/b26_wstar_h7 惯例）：

```
J w_true = −R_x·Dxh
w = splu(M).solve(−rxd, trans='T')      # M = J^T
```

闭合验证 `r_P = (M.T@w)[P0:] + rxd[P0:]`，rel-norm ~1e-11（协调人验收 PASS）。

源层恒等式（PREREG §1.4）：

```
FD_J(Q)·D = b_TC^T w_true + C·Dxh + Gx·Dxh
b_TC^T w_true = (J^{-1} b_TC)^T R_x·Dxh = λ_TC^T R_x·Dxh
             = momentum 段 + pressureRow 段（流介导全部）
```

其中 `λ_TC = J^{-1} b_TC`，且 `w_true = −J^{-1} R_x·Dxh` ⇒ `b_TC^T w_true = −λ^T R_x·Dxh`。
符号组合（本次会话再确认）：w_true 解 `J w_true = −R_x·Dxh`，故 `b_TC^T w_true = −(J^{-T} b_TC)^T R_x·Dxh`，即
`b_TC^T w_true = momentum段 + pressureRow段` 与 PREREG §1.4 一致（正号在 rhs 侧）。

`C = thermalDiffusionDerivativeDTCell`（热直接项，AdjHeatTransfer.H 面公式，从 Tb/T/dDTDxh/网格离线重算，**不乘 V**，sensitivity.H L396-399）。
`Gx = rxFluxDirectT·dAlphaDxh`（通量直接-alpha 项，**正号**）。

D → Dxh 链式规则：`Dxh = d(xh)/dx · D`（Heaviside 投影 drho + 滤波 + designMask，在 state-B 处逐 cell 的 Jacobian，含 eta5 隐式项）。离线重建方向 D1/D2/D3 见 PREREG §1.3。

## 3. H7 vs pre-H7（控制组自洽锚）

- explicitJT_H7.mtx nnz 6395416；explicitJT.mtx（pre-H7）nnz 6209576；834930 个 diff 条目（H7 nnz 的 13%），max abs diff 1.49e-09。
- H7 份额锚：`PRODH7SHARE = 0.701947825646` = |H7 dp-channel|L2 / |kf dp-channel|L2（b25 Log）。
- 预判：`b_TC^T w_true(pre-H7)` 与 `b_TC^T w_true(H7)` 之差应自洽于 ~0.7 份额（H7 与 kf 算子差异集中在 dp 通道）。

## 4. b18_folding_lab.py 模板 DEFECTS（引用时必须修正的清单）

来源：尝试 1 NOTEBOOK + 本次复核 b18_folding_lab.py 源码。修正后使用并在报告列出。

1. **L164-165 顶点均值 centroid**（`Ccen[c] = pts[sorted(vs)].mean`）→ 改用 EXACT pyramid centroid（primitiveMeshCellCentresAndVols.C L95-140 公式：cEst=面心均值；内部面 pyr3Vol=fAreas&(fCtrs−cEst[own])，pc=(3/4)fCtrs+(1/4)cEst[own]，FINAL /= cellVols）。
2. **L38 B16 FD 引脚** → 全部用 B25 tsv（`stage_b2_fd_scan.tsv`）FD_J 行；勿引用 B16 数字。
3. **L373-374 用 E15=b15_export z3** → 改用 `/home/ys/dsH/b8_verify_diag/stageB6_rxc_z_analytic.mtx`（= 链式规则后的 Dxh 方向，reshape(3,N)，maxabs 8.17）；b8==b25 状态身份已验证。
4. **L235 ST=b16_states/stageB2** → 改指 `/home/ys/dsH/b25_qgate/stageB2`（wstate_outletMeta.mtx 等）。
5. **L181-189 硬编码 B16 patch 面数** → 用 polyMesh/boundary 正则（b23 逻辑）读 patch 尺寸。
6. **b_TC 读取**：按 route_V0 布局 `m=matrix.reshape(N,4)`；`bexp[0:NV:3]=m[:,0]`（Ux）、`bexp[1:NV:3]=m[:,1]`（Uy）、`bexp[2:NV:3]=m[:,2]`（Uz）、`bexp[NV:]=m[:,3]`（P）。NV=100800, NUNK=134400, N=33600。

可复用（无缺陷）：read_field（L44-127，括号剥离）、mesh build 骨架（Sf via 0.5*cross(fv,roll).sum）、C_field 面公式（L376-408）、build_dHbyA_drAU（L462-513）、Gx_of（L522-533）、g 面泛函（L229-240）。

## 5. M4 口径声明（弹性目标，本会话确认）

- **b24_diag 的 J-functional rxj 系列不可用**：/1/Ub、/1/pb 为 NaN，rxj_prod_gsensh_*.mtx 全 NaN（NaN-crash 运行残留）。M4 只能用 b8_verify_diag 的 **rxpr（压力降泛函）系列**（sensitivity.H L334-360 导出 gsenshPressureDrop*，有限）。
- **rxpr vs ADJ_gDP 不匹配（开放差异，声明不掩盖）**：dot(rxpr_{momentum,pressurerow,total}, z) = (−0.04869/0.22908/−0.03013)、( −0.00755/0.06859/−0.00621)、( −0.05624/0.29767/−0.03634)；而 ADJ_gDP = −5.2303/1.0505/13.6034。非均匀比值 → 不是尺度因子。projV 锚验证（dot(gsensVol,D)=dot(gsensVol,z)=projV anchors，~1e-8 级吻合）确认 z 是正确链式方向。假设（不掩盖）：生产 ADJ_gDP 投影路径可能用了 filter 链后的 z，或 b8 rxpr 导出早于 Uc 完全一致；M4 表中声明为开放差异。
- thermalC 段 = total − momentum − pressureRow 的框架适用于 J 泛函分解；本表声明实际用的是压力降泛函 rxpr 系列。

## 6. 执行日志

- 2026-08-22 03:43 M1 运行（仪器未打补丁，PRODH7SHARE=0.701947825646 为 b25 Log 锚；偏差见 EXECUTOR_SUMMARY §4）：
  - D1：lhs=−4.884e-02，FD=+4.252e-02，C=−3.034e-02，Gx=−1.570e-03，rhs=+7.444e-02，ratio=−0.656，sign −/+，reldev=165.6% → O(1)-MISMATCH（c 段）
  - D2：lhs=−6.441e-03，FD=−3.529e-03，C=−8.211e-03，Gx=+4.043e-03，rhs=+6.388e-04，ratio=−10.08，sign −/+，reldev=1108% → O(1)-MISMATCH（c 段）
  - D3：lhs=+1.157e-02，FD=−4.450e-04，C=−1.979e-02，Gx=+8.407e-03，rhs=+1.094e-02，ratio=+1.057，sign +/+，reldev=5.7% → SOURCE-CLEAN（a/b 段）
  - pre-H7 控制组：lhs_pre/lhs = +1.127 / −0.418 / +1.397（D2 反号，节线敏感，见总结 §3.4）
  - 门控：mobility 2.19e-15、V 2.72e-14（机器精度 PASS）；route_V1 vs b_TC relL2=0.397（cos 0.9583）、route_V0 relL2=0.323（cos 0.9584）
  - 裁决：**混合裁决**——D1/D2 指向 (c) T-消除折叠，D3 指向 (a)/(b)；M2/M3（recon）未执行，B3a/b/c 开放
- 产物：b26a_m1_identity.tsv、b26a_m1_gates.tsv、b26a_m1_results.json、b26a_EXECUTOR_SUMMARY.md、b26a_EXECUTOR_DONE
