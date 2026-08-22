# BFINAL-029 — R1 修复轮：c 段逐槽隔离与折叠修复（执行磁盘笔记本，逐小步追加）

> 任务书：`/home/ys/dsH/TO-ANISOTROPIC/NEXT_TASK_BFINAL029.md`。
> 仓库 `/home/ys/dsH/TO-ANISOTROPIC`（真实仓库；`/home/ys/TO-ANISOTROPIC` 为无关空壳，勿用）。
> 分支 `agent/dsH-stage-b-validation`。预注册（含第二阶段定罪判据与预期，**数据之前**）：`b29_PREREGISTRATION.md`。

## 0. 时间轴

- 2026-08-21（继续会话）：完成全部侦察（g 层、生产/诊断折叠槽位、目标类型、数据可用性）；**关键发现**：
  b26a route_V1 的出口 dJ/dphi 用了 `-(T-Tmix)`（maximizeColdOutletTemperature 公式），而 b25_qgate 案例实际运行
  `maximizeTotalHeatTransfer`（生产 computeObjective.H L152-172 用 `-(T-T_in)`，T_in=gAverage(coldInlet)=600）。
  离线量化：Tmix=687.840，TIn=600；dJdphi(Tmix) L2=4.949e2 vs dJdphi(TIn) L2=6.668e2，relL2=0.494、L2 比 0.742。
  即 B26 的 0.397 中有一部分是 g1 槽位离线公式偏差，需在 Phase 1 中分解。
  另确认：b26a 跳过的 T8 必须作为可切换槽位实测（本案例内生产 T8 循环无 U-fixed 跳过；
  但 noSlip/fixedValue-U 补丁 internalCoeffs 各分量相等 → (−ic+cav)=0；outlet zeroGradient-U → ic=0；
  故本案例预期 T8≡0，**必须实测而非假设**）。

## 1. 输入核验表（全部通过）

| 输入 | 路径 | 核验结果 |
|---|---|---|
| w_true(H7) | `evidence/agent-group/BFINAL-026/cycle-1/b26_wstar_h7.npz` | D1/D2/D3 各 (134400,) float64 全有限；maxabs 295844/573932/714579 |
| 求解日志/闭合 | `b26_wstar_h7.json/.log` | 闭合 1.037e-11/9.649e-12/2.428e-11（<1e-9 PASS）；耗时 1006.2s，**勿重跑** |
| b_TC 导出 | `/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx` | 33600×4（Ux Uy Uz P），PHYSICAL/未缩放 |
| FD_J(Q) | `stageB2/stage_b2_fd_scan.tsv` | 判据行（h=1e-3）：D1 +0.04252469623637622，D2 −0.003528972517971574，D3 −0.0004449511895182612 |
| M1 rhs 三件套 | `b26a_m1_results.json` | rhs=FD−C−Gx：D1 0.07443678641907196，D2 0.0006387577668902512，D3 0.01094029857921596 |
| 网格 | `/home/ys/dsH/b16_mesh/` | owner/neighbour/boundary/sf/weights/b16mesh_V_nueff.mtx 齐备 |
| state-B 场 | `/home/ys/dsH/b25_qgate/1/` | T,U,p,phi,phiThermal,Tb,alpha,dAlphaDxh,dDTDxh,coldFaceMask 等 |
| SB 基准 | `stageB2/wstate_baseline_*.mtx`, `wstate_outletMeta.mtx` | Mflow=4.07337e-3, Tref=600 |

## 2. 目标类型与 g1 公式（本会话新确认）

- 案例（b25_qgate Log L115）：`thermalObjectiveType=maximizeTotalHeatTransfer`。
- 生产 `computeObjective.H` L152-172：`dJ/dphi_f = -(T_out - T_in)/(M_frozen*Tref)`，`T_in = gAverage(T[coldInlet])`。
- 案例 `constant/thermalProperties` L30-31：`coldInletPatch inlet; coldOutletPatch outlet`。
- 案例 `1/T` L33628-33632：inlet fixedValue 600 → **T_in = 600.0**。
- b26a instrument L289-291 用 `Tmix=sum(ophi*Toc)/Mflow=687.840` 和 `-(Toc-Tmix)/(Mflow*Tref)`：
  **与生产公式不符**（那是 maximizeColdOutletTemperature 的公式）。
- b29 默认 route 采用生产公式 g1_TIn：`dJdphi_out = -(Toc - 600.0)/(Mflow*Tref)`；
  g1_Tmix 作为对照变体保留（复现 B26 0.397）。

## 3. 槽位表（17+2，含 b29 默认公式与变体）

| 槽 | 生产文件:行 | 内容 | 符号 | b29 状态 |
|---|---|---|---|---|
| g0 | AdjNS_HT.H L23-25 | 内部面 g_f = −coldmask·adjDown·jumpT | − | 实现 |
| g1_TIn | computeObjective.H L152-172 | 出口 dJ/dphi = −(T−T_in)/(M·Tref) | − | **默认** |
| g1_Tmix | b26a L289-291（对照） | 出口 dJ/dphi = −(T−Tmix)/(M·Tref) | − | 变体 |
| T1 | Production:359-362 | hA += αRel·rAU·w·(Sf·g) | + | 实现 |
| T3 | Production:363-369 | U += (1−αRel)·w·(Sf·g) | + | 实现 |
| T4 | Production:370-378 | P(own) += kf·g；P(nei) −= kf·g | +/− | 实现 |
| H7s1 | Production:379-385 | h7G += mob·w·g·Sf | + | 实现 |
| T5 | Production:403-404 | 边界 hA += αRel·rAU·(Sf_b·g_b) | + | 实现 |
| T6 | Production:405-409 | 边界 U += (1−αRel)·(Sf_b·g_b) | + | 实现 |
| H7s1b | Production:410-412 | 边界 h7G += mob·g_b·Sf_b | + | 实现 |
| T7 | Production:413-420 | 边界 P += mob·δ_b·\|Sf_b\|·g_b（pressureFixed） | + | 实现 |
| T2 | Production:428-441 | δH^T：U(nei) −= upper/V_o·hA_o；U(own) −= lower/V_n·hA_n | − 两侧 | 实现 |
| T8 | Production:442-465 | 边界 H 对角：(−ic+cav)·Vinv·hA(c) | + | **实测（预期≡0）** |
| H7s2 | Production:469-481 | q=Sf&(h7G_o/V_o−h7G_n/V_n)；P(own) −= w·q；P(nei) −= (1−w)·q | − 两侧 | 实现 |
| H7s2b | Production:482-497 | P(c) −= (h7G&Sf_b)/V_c（仅 zeroGradient-p） | − | 实现 |
| A1 | Production:499-508 | U/P 汇总 | | 实现 |
| A2 | stageB18RhsExport L870-898 | 导出布局（bexp 排列） | | 实现 |
| C（外挂） | sensitivity.H L396-399 | thermalDiffusionDerivativeDTCell（不乘 V） | | M1 rhs 已有 |
| Gx（外挂） | — | g^T·φ_x（真缺失，超范围） | | M1 rhs 已有 |

## 4. 执行计划（与预注册一致，数据前已写）

- Phase 1（纯离线）：b29_slot_instrument.py；逐槽独立贡献 + 线性叠加门；leave-one-out relL2 分解全 0.397。
- Phase 2（纯离线）：变体×身份数字表（见预注册 §4，判据原样引用 §4.2）。
- Phase 3（仅定罪后）：四路一致性修复 + 编译 + 全验收序列。

---
## 5. Phase 1 结论（2026-08-21，数据全在 `b29_slot_table.json` / `b29_m1_state_probe.json`）

### 5.1 门全部通过
- 自洽门 [g1_Tmix, T8 关]：relL2 = 0.39701093640754337（与 B26 route_V1 **逐位一致**，PASS）。
- 加性门：max|sum−full| = 3.4e-21，rel = 7.0e-17（机器精度，PASS）。
- 迁移率门 2.19e-15、V 门 2.72e-14（PASS）。
- T8 实测贡献：9.6e-18（理论恒零**实测确认**；b26a"跳过"注释保留有效，T8 非嫌疑）。

### 5.2 全 38.6%/39.7% 失配分解（leave-one-out ΔrelL2，全开基线 [g1_TIn, T8 开] relL2=0.38617）
- g1 口径（Tmix vs TIn）：ΔrelL2=0.01084 —— **caliber**（B26 用错公式，b29 已修正为生产 TIn）。
- T3（内部 U 折叠，g0 载体）：Δ=+0.5867，占 b 范数 **93.5%** —— 绝对主导。
- T6（边界 U 折叠）：Δ=+0.0266，占 27.6%。
- T2（δH^T）：Δ=+0.0115，占 3.5%。
- T1（hA 内部）：Δ=+0.0114，占 3.3%。
- 其余（T4/H7s1/T5/H7s1b/T7/H7s2/H7s2b/T8）：Δ<1e-3，relL2 层面 **innocent**。
- 逐槽标签（预注册 §2.1.6）：suspect = {T3, T6, T2, T1}；innocent = 其余。

### 5.3 M1 收缩预览（route 全开 vs 导出；lhs=b^T wH7，rhs=FD−C−Gx 为 B26 冻结件）
| 行 | route lhs | route ratio | 导出 ratio |
|---|---|---|---|
| D1 | +0.034728 | **+0.4665** | **−0.6561**（符号相反） |
| D2 | −0.006863 | −10.7446 | −10.0836 |
| D3 | +0.002241 | **+0.2049** | **+1.0573**（route 低 5 倍） |

### 5.4 M1 状态探针（决定性）
8 种状态（导出 / state-mix / all-1 / all-B2 / 单场切换 T、phiTh、phi_bvals、Tb→/1 代理）的 M1 矩阵：
**无任何状态复现导出的 M1 行**。D1 在所有 route 状态下为正（+0.17..+0.47）vs 导出 −0.656；
D3 为 −1.02..+0.20 vs 导出 +1.057。all-B2 cos 最高（0.9877）但 M1 最差（D1=+0.1733, D3=−0.6748）。
→ **D1 符号翻转与 D3 振幅缺口非状态可切换**，属于结构性（槽折叠）或输入级（Tb）。

### 5.5 g1 排除（口径核验，本会话完成）
- 生产 `computeObjective.H` L152-172（maximizeTotalHeatTransfer）：`dJ/dphi_f = -(T_out - T_in)/(M_frozen*Tref)`，`T_in = gAverage(T[coldInlet])`。
- 案例 inlet T = fixedValue uniform 600 → **T_in = 600.0**，与 instrument TIn 硬编码精确相等；
  M_frozen = sum(phiOutlet) snapshot，`wstate_outletMeta.mtx` Mflow=4.07337e-3 匹配。
- → g1 槽排除（非 D1 符号来源）。注意：g1@B2 存在 2.8% 校准差（111.855 vs log 108.816），留作记录，非符号源。

### 5.6 尺度消除证明（承接上一会话，已归档）
- 导出端 prodRhs=raw/scale，写回乘 prodMomentumScale/prodAreaScale → 尺度精确相消，**导出即 raw 生产 rhs**；
  1.82×U / 1.18×P 振幅缺口为真实缺口。

### 5.7 Tb 不可得（输入级）
- Stage B2 全盘无 Tb 导出（搜 `wstate*Tb*` 无果）；B26 M1 instrument L237 同样用 /1 Tb。
- Tb 是唯一不可状态切换的 g0 输入（adjDown = Tb[nei|own]，直乘 jumpT，经 T3 折叠占 b 范数 93.5%）。
- Tb@B2 与 /1 Tb 之差是 M1 D1/D3 缺口的**输入级候选**，本阶段无法直接测试。

### 5.8 Phase 1 结论
全失配已分解：T3 主导（93.5% 范数），T6/T2/T1 次级，其余 innocent；g1 为 caliber。
导出 vs route 的 M1 分歧（D1 符号、D3 振幅 5 倍）在任何可切换状态下均不可消除 →
**结构性（槽折叠符号/尺度）或输入级（Tb@B2 缺失）**。Phase 2 嫌疑槽 = {H7s2（B26 主犯）, T4（B26 次犯）} ∪ Phase 1 发散槽 {T3, T6, T2, T1}。

## 6. Phase 2 D3 判据解释（数据前决定，绑定）

- 预注册 §4.2 的 D1/D3 判据是**相对 B26 rhs 的绝对带**（rhs = FD_J(Q)−C−Gx，冻结，不重算）：
  - **D1_v ∈ [0.85, 1.15]**（B26 基线 −0.656 → 变体需把 D1 拉进带内）；
  - **D3_v ∈ [0.90, 1.20]**（"D3 当前 1.057" 指**导出的** D3 = 生产真值，**不是** route 的 0.2049）。
- route 基线 D3=0.2049 **不在带内** → 定罪变体必须同时：
  (a) 把 D1 从 +0.4665 提到带内（lhs 需 ~+0.0744，差 +0.0397，约为当前 2 倍）；
  (b) 把 D3 从 +0.2049 提到带内（lhs 需 ~+0.01094，差 +0.00870，约为当前 5 倍）。
- 这是预注册的**绑定读法**。不允许按 route 可达带重缩放判据（= 拟合因子）。
- 若变体扫描不能同时满足 (a)(b) → 按预注册 §4.2 排除规则：诚实报告"候选排除"，写结论，STOP，不进入 Phase 3。
- 提示（记录）：route D1=+0.4665（正）已比导出 −0.656 更接近目标 —— 按生产源码正确符号的 route 折叠本身已将 D1 符号修正；若变体仍不能入带，剩余缺口指向输入级（Tb@B2 缺失）。

---
## 7. Phase 2 结论（2026-08-21，数据在 `b29_variant_matrix.json` + `b29_crossscan.json`）

### 7.1 变体矩阵（15 变体，b29_phase2_variants.py，单/双切换）
V0_base = route 全开。对每个嫌疑槽（{H7s2, T4, T3, T6, T2, T1}）做 Z_i（归零=离开除）、F_i（翻号），
另加 X_Z_H7s2+T4、X_F_H7s2+T4。判据原样（D1∈[0.85,1.15] 且 D3∈[0.90,1.20]，rhs 冻结）：

| 变体 | D1 | D3 | D1ok | D3ok |
|---|---|---|---|---|
| V0_base | +0.4665 | +0.2049 | n | n |
| Z_T4 | **+0.8899** | −0.5607 | **Y** | n |
| Z_H7s2 | +1.8617 | −24.8031 | n | n |
| Z_T1 | +0.5675 | −0.3777 | n | n |
| Z_T3 | +4.4290 | −14.4970 | n | n |
| F_T3 | +8.3914 | −29.1990 | n | n |（relL2→1.61, cos→−0.87，结构崩坏） |
| Z_T6 / F_T6 | −4.4857 / −9.4380 | +17.86 / +35.51 | n | n |
| 其余（F_H7s2 +3.26/−49.8, Z_T2 +0.46/−1.21, F_T2 +0.46/−2.63, F_T1 +0.67/−0.96, X_Z_H7s2+T4 +2.29/−25.6, X_F_H7s2+T4 +4.10/−51.3） | 均出带 | 均出带 | n | n |

**单开关唯一进带者**：Z_T4（D1=+0.8899 入带）但 **D3=−0.5607 被彻底破坏**。

### 7.2 全双切换交叉扫描（42 变体，b29_phase2_crossscan.py，补全"双变体交叉"）
嫌疑槽 ∪ H7s2b（H7 族近相消对）上枚举**全部** C(7,2)=21 对 × {Z,F} = 42 变体。**无一定罪**：
- 仅 2 个达 D1ok：Z_T4+T2（D1=+0.8870, **D3=−1.9764**）、Z_T4+T1（D1=+0.9909, **D3=−1.1434**）；
- 其余 40 个 D1/D3 至少一个出带（幅度常达 O(10..100)，如 Z_H7s2+T3 D1=+5.82 D3=−39.5）。
- 结论：**D1 拉进带与 D3 进 [0.90,1.20] 在所有 零/翻号 组合下互斥**。D3 方向被近相消对主导
  （H7s2 D1=−0.104/D3=+0.274，H7s2b D1=+0.095/D3=−0.260；T3 +0.161 vs T6 −0.193），
  槽切换只能在 O(0.1..0.6) 步长内移动 D3，而 D3 需要 ~5 倍放大（目标窗口 rhs×[0.90,1.20]=[0.00985,0.01313]，
  route lhs=+0.00224）——槽符号/归零自由度**不覆盖**该方向。

### 7.3 逐槽 M1 贡献（delta = 全开 − 除 S，b29_variant_matrix.json per_slot_contrib）
T3 D1=−0.295/D3=+0.161；T6 D1=+0.369/D3=−0.193（主导）；H7s2 D1=−0.104/D3=+0.274；
H7s2b D1=+0.095/D3=−0.260（近相消对）；T4 D1=−0.0315/D3=+0.0084；T8≈0。

### 7.4 定罪判定（预注册 §4.2 绑定判据）
**无一变体通过 → 候选排除（H7s2/T4/Phase1 发散槽均无法使 M1 一致）**。
结构性（槽折叠符号/归零）假设在**全部** 单/双 零翻号 组合上被穷尽否定；
剩余缺口指向输入级（Tb@B2 缺失，§5.7）与 g1@B2 2.8% 校准差（§5.5），均非本阶段可修复。
**STOP：不进入 Phase 3（不触碰生产源码 solveDiscreteFlowAdjointProduction.H）**。

### 7.5 已知限制（如实记录，不拟合）
- g1@B2 校准差 2.8%（111.855 vs log 108.816）：B2 无 g1 源快照，不可离线判定其是否为 D3 缺口来源。
- Tb@B2 缺失：唯一不可状态切换的 g0 输入；需 B2 运行端补导出才能最终定位 D1/D3 缺口。
- 无 FD 重算（Phase 2 纪律）：rhs 为 B26 冻结件；若 B26 rhs 本身含与 B2 相同的 Tb/g1 输入，缺口被共同携带。
