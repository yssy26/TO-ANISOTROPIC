# BFINAL-028 — c 段参照系对照轮 REFCROSS_REPORT

日期：2026-08-21。分支 agent/dsH-stage-b-validation @ 88ee077。
纪律执行：零编译、零求解器运行、零生产代码改动、无 git 提交。
对照对象：我方 b_TC 折叠（T-elimination：dJ/dT → b_TC = M^T·g）vs 外部参照系。
D1 反号指纹（待解释形态）：ADJ_J=−0.0358 vs FD_J=+0.0425（O(1) 反号）；M1 lhs=−0.0488 vs rhs=+0.0744，ratio=−0.656。

---

## 1. 我方 b_TC 折叠完整项清单（预注册，来自 NOTEBOOK.md §2）

### g 层（面泛函 g_f，AdjNS_HT.H discrete 分支）

| # | 代数形式 | 符号 | 源码 file:line |
|---|---|---|---|
| g0 | g_f = −mask_f·Tb_down·(T_nei−T_own)，Tb_down=upwind(Tb) | 显式 MINUS | AdjNS_HT.H:11-26 |
| g1 | 出口面 g_f += +dJ/dφ_out = −(T_out−T_in)/(M_frozen·Tref) | +dJ/dφ | AdjNS_HT.H:28-41 |
| g2 | 其余边界 g_f = 0 | — | 无赋值 |

### 路由层——内部面（Production L351-386）

| # | 代数形式 | 符号 | 源码 file:line | 槽位 |
|---|---|---|---|---|
| T1 | HA[own] += αRel·rAU_o·w·(Sf·g_f)；HA[nei] += αRel·rAU_n·(1−w)·(Sf·g_f) | + | Production:359-362 | αRel·rAU·dH hA 累积 |
| T3 | U(own,c) += (1−αRel)·w·(Sf·g_f)_c；U(nei,c) += (1−αRel)·(1−w)·(Sf·g_f)_c | + | Production:363-369 | (1−αRel) 直接内插转置 |
| T4 | kf=interp(mob)·δ·|Sf|；P(own) += kf·g_f；P(nei) −= kf·g_f | own+/nei− | Production:370-378 | −kf_f(dp_n−dp_o) 转置 |
| H7s1 | g_c[own] += mob_o·w·g_f·Sf；g_c[nei] += mob_n·(1−w)·g_f·Sf | + | Production:379-385 | H7 stage1 g_c 累积 |

### 路由层——边界（Production L387-422；U-fixed 跳过 L389-392）

| # | 代数形式 | 符号 | 源码 file:line |
|---|---|---|---|
| T5 | HA[celli] += αRel·rAU_c·(Sf_b·g_b) | + | Production:403-404 |
| T6 | U(cell,c) += (1−αRel)·(Sf_b·g_b)_c | + | Production:405-409 |
| H7s1b | g_c[celli] += mob_c·g_b·Sf_b | + | Production:410-412 |
| T7 | P(cell) += mob_c·δ_b·|Sf_b|·g_b（仅 pressureFixed 出口） | + | Production:413-420 |

### deltaH^T 展开（Production L428-465）

| # | 代数形式 | 符号 | 源码 file:line |
|---|---|---|---|
| T2 | U(nei) −= prodUpper·(1/V_o)·HA[own]；U(own) −= prodLower·(1/V_n)·HA[nei] | − 两侧 | Production:428-441 |
| T8 | U(cell) += (−ic_c+cav)·(1/V_c)·HA[cell] | + | Production:442-465 |

### H7 stage 2（Production L466-497）

| # | 代数形式 | 符号 | 源码 file:line |
|---|---|---|---|
| H7s2 | q=Sf&(g_o/V_o−g_n/V_n)；P(own) −= w·q；P(nei) −= (1−w)·q | − 两侧 | Production:469-481 |
| H7s2b | P(cell) −= (g_c&Sf_b)/V_c | − | Production:482-497 |

### 装配（Production L499-508）+ 缩放（L834-845）

| # | 内容 | 源码 file:line |
|---|---|---|
| A1 | 汇总 U/P 槽进 discreteFlowVelocityRhs/DiscreteFlowPressureRhs | Production:499-508 |
| A2 | prodRhs = rhs/scale；导出 b18rhs_thermalCoupling.mtx（物理未缩放） | Production:834-845, L870-898 |

### 范围外（参与 M1 恒等式，不在 b_TC 内）

| # | 内容 | 源码 file:line |
|---|---|---|
| C | C=−Tb^T·R_T,x（thermalDiffusionDerivativeDTCell） | AdjHeatTransfer.H:189-328 |
| Gx | Gx=g^T·φ_x（通量直接 α 项，装配层授权外，生产缺） | B18 DERIVATION_FIX3 |

**共 17 项（g 层 3 + 路由内部 4 + 路由边界 4 + deltaH 2 + H7s2 2 + 装配 2），外挂 2 项（C/Gx）。**

---

## 2. 逐项参照对照表

参照系实测清单：
- **Fira thermalTopO**（https://github.com/Fira-Software/thermalTopO）— 连续伴随 CHT-TO，含显式热-动量折叠 gPhi。最接近 g0 的逐项参照。文件：`src/solvers/thermalAdjointSimple/thermalAdjointSimple.C`（gPhi L1280-1390，投影路由 L1578-1800）、`docs/derivation.md`（§3.1/3.2/6.2）、`docs/atc-t-open-channel.md`。克隆 /tmp/fira_thermalTopO。
- **Othmer continuous adjoint**（OpenFOAM-dev 随附 legacy）— 孔隙率伴随源符号。`adjointShapeOptimisationFoam.C`。克隆 /tmp/of_dev_ref。
- **SU2**（https://github.com/su2code/SU2）— 离散伴随 = AD/coDiPack，无手工 g_f 折叠（方法论差异）。`CHeatSolver.cpp:712/721-731`（AVG_TEMPERATURE=面积平均）、`CDiscAdjHeatIteration.cpp`。克隆 /tmp/su2_ref。
- **NTUA 双流体换热器 TO**（Crossref: 10.1108/hff-08-2024-0642, 10.1007/s00158-022-03330-w）— 连续伴随 TDDC，ATC 耦合 +T_a∇T。摘要级。

| 我方项 | 参照侧对应 | 符号对照结论 | 置信 |
|---|---|---|---|
| **g0** | Fira gPhiI = Tf·(AP−AN) = Ta_upwind·(T_own−T_nei)（C≡1, chi=0）。**与我 g0=Tb_down·(T_own−T_nei) 同形同号** | ✓ 同号，无翻转；迎风伴随×own−nei 温差结构一致 | 高 |
| **g1** | SU2 AVG_TEMPERATURE 为**面积平均**（非质量流量加权）；NTUA 出口目标经 T_a BC 进入（连续伴随） | ⚠ 功能形式差异：我方 Tmix 质量加权 vs SU2 面积平均；B0.2 已证 <Tb,dR_T/dφ_out> 恒零，RHS 恰为 dJ/dφ | 中 |
| **g2** | Fira 非目标边界 gPhi 走物面通量 BC | — | 中 |
| **T1/T2** | 无参照直接覆盖（SU2=AD；连续伴随无 αRel·rAU·dH^T 槽）。Fira atc-t 文档明言连续源单独不是 SIMPLE 映射的完备离散转置 → 此层正是完备性所在 | 无逐槽参照；依赖 B18 槽位镜像 + M1/M2 恒等式 | — |
| **T3** | 同上（(1−αRel) direct 为 relax 源项导数，连续伴随无对应槽） | 无逐槽参照 | — |
| **T4** | Othmer: Sp(α,U) 阻力伴随为 +Sp(α,Ua) 自伴随无翻号（对角）；Fira 面通量 kf 类槽无对应 | kf P-row 与算子自身 J_PP kf 槽同号（Production L339-342 注释，applyProdFlowJT 同文件镜像） | 中 |
| **H7s1/s1b/s2/s2b** | 无参照直接覆盖（H7=phiHbyA 压力通道，B24 推导 D5）；连续伴随用 lambda 投影替代（架构差异） | 无逐槽参照；镜像 applyProdFlowJT stage2 | — |
| **T5/T6/T7/T8** | 无参照直接覆盖；Othmer 伴随出口速度 phia·Sf/\|Sf\|²+切向（连续伴随出口范式） | 无逐槽参照 | — |
| **C / Gx** | NTUA dJ/dβ = J_β + ψ^T M_β 完整装配（含通量直接项）；我方生产缺 Gx（授权外） | ⚠ Gx 确认缺失（范围外项） | 高 |

**逐项对照规模：17 项我方 b_TC + 2 外挂，全部过了一遍；其中 g0 有高置信参照同号确认，g1 有功能形式差异确认，其余路由层项无直接参照（架构差异/AD 黑箱/连续伴随），依赖内部自洽证据。**

---

## 3. 候选错项排名（按 D1 反号解释力 + 双侧 file:line 证据）

D1 反号形态约束：O(1) 反号 → 符号/缺项类错误，非系数微差；D3 ratio=+1.057 干净 → 错误方向依赖（非全局翻号）；M1 |lhs|<|rhs|（−0.0488 vs +0.0744）→ 非纯整体符号翻（幅度也偏 ~35%，存在部分抵消）。

| 排名 | 候选（类/项） | 对 D1 的解释力 | 证据 | 置信 |
|---|---|---|---|---|
| **1** | **E2-路由层·压力行通道（H7s2 首选，T4 次之）** | O(1) 方向依赖反号 + 幅度偏移完全吻合：压力伴随 λ_p 符号翻转即 M1 lhs 翻号，且 H7s2/T4 仅被部分方向激活 → D3 干净可解释 | 双侧：我方 Production:370-378(T4)、469-481(H7s2) 无参照逐槽镜像（唯一可交叉的 Fira 用 lambda 投影替代，无此槽）；H7s2 为本轮最新加入项（B24），无 B18 级槽位复核 | **中高** |
| **2** | **E2-路由层·αRel·rAU·dH H^T 块（T1/T2）** | 方向依赖 off-diag U 行；B18 曾发现旧路由整体缺此块——残余符号/完备性问题在此层合理 | 双侧：Production:359-362, 428-441；B18 DERIVATION_FIX3（旧路由缺 αRel·rAU·dH H^T 与 (1−αRel) 因子）；Fira atc-t 明言连续源单独非 SIMPLE 完备转置 | 中 |
| **3** | **E3-g1 功能形式（质量加权 Tmix vs 面积平均）** | 出口泛函形式差异在回流出口方向依赖，可致局部反号；但功能形式差非 O(1) 翻号的典型因，且 D3 干净 | 双侧：AdjNS_HT.H:28-41（我方质量加权 Tmix）；SU2 CHeatSolver.cpp:712,721-731（AVG_TEMPERATURE 面积平均） | 低-中 |
| 4 | E1（g0/g1 整体符号） | 被否：Fira gPhi 与 g0 同形同号（高置信）；且若 g 层整体翻号则 D3 也应反（D3 干净） | 双侧：AdjNS_HT.H:11-26 vs Fira thermalAdjointSimple.C:1357（chi=0 时 Ta_upwind·(T_own−T_nei)） | 否 |
| 5 | E4（缺连续伴随源 +T_a∇T 整体项） | 被否：所有参照（Fira/NTUA/Othmer）在冻结物性、零浮力、非 buoyant 情形下均仅通过 φ/ATC 通道耦合，无直接 T→U/p 源；我方 discrete 分支形态正确 | Fira derivation.md §3.2 ATC-T；NTUA TDDC 摘要 | 否 |

**一行裁决（最高置信候选）：E2-路由层压力行通道——H7s2（H7 stage2 P-P 通道，Production:469-481）为首选，T4（kf P-row，Production:370-378）为次选——是唯一与"D1 O(1) 方向依赖反号 + D3 干净 + M1 幅度偏差"全部吻合且无参照覆盖的项；置信级：中高（类级判定可靠，但单项隔离需 M1 逐槽测试——超出本轮纯阅读授权）。**

---

## 4. 参照有而我方无（缺项清单）

| # | 参照项 | 参照位置 | 我方状态 | 分类 |
|---|---|---|---|---|
| M1 | **Gx = g^T·φ_x（通量直接 α 项）** | NTUA dJ/dβ = J_β+ψ^T M_β（HFF 2025 TDDC 摘要；Othmer λ·(Ua&U) 同族） | 生产缺（sensitivity.H 装配层，授权外） | **真缺项**（范围外） |
| M2 | **SU2 面积平均出口温度（AVG_TEMPERATURE）** | CHeatSolver.cpp:712,721-731 | 我方为质量流量加权 Tmix | 功能形式差异（g1 级） |
| M3 | Fira bounded 修正 chi·(AP·T[P]−AN·T[N]) | thermalAdjointSimple.C:1357 | 无此显式项 | 不适用：我方为对实际 bounded-Gauss-upwind 算子的精确离散转置（B0.2 恒等 1e-10），chi 修正已被转置吸收 |
| M4 | Fira lambda 投影（无散面通量投影 + w·Sf/V 路由） | thermalAdjointSimple.C:1658-1755 | 无投影 | 架构差异：连续伴随需投影强制无散；离散 M^T 精确转置无需 |
| M5 | Othmer 伴随出口速度范式（phia·Sf/\|Sf\|²+UtHat） | adjointOutletVelocityFvPatchVectorField.C | 我方 T5/T6/T7 槽位不同 | 架构差异（离散转置 vs 连续伴随出口条件） |

**关键缺项结论：参照系中唯一真正"有而我没"且会改变目标敏感度的是 M1（Gx），但其在 b_TC 链外（装配层，授权外），且 B26 已判 D1 反号对 C+Gx 偏差稳健（需 >180% 才翻号）。M2 为功能形式差异（影响 g1 的量值/回流方向，非 O(1) 翻号主因）。M3/M4/M5 为架构差异，不构成缺项。**

---

## 5. 开放项

1. **单项隔离需运行 M1/M2 逐槽测试**（H7s2 关/开、T4 号翻转）——本轮零运行授权，留待下游。
2. 真 ESI adjointOptimisation 模块不可达（auth/403，详见 refs/esi_adjointOptimisation_inaccessible.md）；离散路由层逐槽参照仍无外部镜像，依赖内部恒等式。
3. SU2 AVG_TEMPERATURE 与 Tmix 的差异是否影响 D1 方向——需构造回流出口数值检验（运行轮）。
4. NTUA HFF 2025 paywall，仅摘要级；若开放全文可补 ATC 逐式对照。

## 附：参照抓取物（refs/ 目录，均带来源 URL）

- `su2_avg_temperature.md` — SU2 面积平均出口温度（CHeatSolver.cpp:712/721-731），URL https://github.com/su2code/SU2
- `su2_disc_adj_methodology.md` — SU2 离散伴随=AD/coDiPack，无手工 g_f 折叠；CHT 接口无显式伴随方法
- `fira_thermalTopO_flux_sensitivity.md` — Fira gPhi 公式 + lambda 投影路由 + ATC-T open-channel 佐证，URL https://github.com/Fira-Software/thermalTopO
- `openfoam_othmer_adjoint.md` — Othmer Sp(α) 自伴随、UaEqn、敏感度 λ·(Ua&U)，URL https://develop.openfoam.com/Development/OpenFOAM-dev
- `esi_adjointOptimisation_inaccessible.md` — ESI 模块抓取失败记录与降级声明
- `ntua_bi-fluid_heat_exchanger.md` — NTUA 双流体换热器 TO（Crossref DOI 10.1108/hff-08-2024-0642 等）
