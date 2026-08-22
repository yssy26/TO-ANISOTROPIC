# BFINAL-031 — 装配链统一审判轮 — NOTEBOOK（执行磁盘笔记本）

> 纯离线测量轮：零编译、零求解器、零生产改动。仓库 `/home/ys/dsH/TO-ANISOTROPIC`，
> 分支 `agent/dsH-stage-b-validation`，HEAD `00a7f27`。任务书 `NEXT_TASK_BFINAL031.md`。
> 预注册：`PREREGISTRATION.md`（先于对账数据写就）。仪器：`b31_recon.py`（复用 b26a/b27 组件）。

## 0. 标签与符号约定（动手前写下，B13 惯例核对）

- **Uc = pressureDrop 标签伴随速度；Ub = thermalCoupling 标签伴随速度**（B13 惯例，
  已由 sensitivity.H L65 `fsenshMeanT = -dAlphaDxh*(U & Ub)`（J 链）与 L89
  `gsenshPressureDropMomentum = -dAlphaDxh*(U & Uc)`（gDP 链）双重确认：J 用 Ub、gDP 用 Uc）。
- **pc = pressureDrop 标签伴随压力**（L100 注释：T = J_P^T pc）；**pb = thermalCoupling 标签
  伴随压力**（L193-196：RX_ADJ_PRESSURE=pb → rxPressureRowTb）。
- 装配（xh 层，链式规则前）：
  - gDP：`gsensh_PD_total = momentum + pressureRow`，
    momentum = `−dAlphaDxh·(U&Uc)·V`（L89-90，含 V），pressureRow = `−rxPressureRowT·dAlphaDxh`（L115，无 V）。
  - J：`fsenshMeanT = momentum(−dAlphaDxh·(U&Ub)·V, L65) + pressureRow(−rxPressureRowTb·dAlphaDxh, L239)
    + fluxDirect(+rxFluxDirectT·dAlphaDxh, L255)`，随后 L400-402 `+= thermalDiffusionDerivativeDTCell`（C）。
- 链式规则（filter_chainrule.H，L10-232）：×designMask → ×drho → eta5 体保持修正
  （`+= mask·V·(1−drho)·C_eta`，C_eta=Σ(mask·dProjDeta·g)/Σ(mask·dProjDeta·V)）→
  **Helmholtz 滤波伴随求解** `(L−diag(b))·post = −b·pre`（全局椭圆）→ ×designMask。
- 投影（validateStageB2GradientAmplitude.H L477-599）：方向 Dk 为物理坐标正弦方向，
  max=1 归一化，支撑=designMask∩(0.02<x<0.98)；`projX_k = Σ_active dgdx[k]·Dk`；
  `dgdx[1] = pressureGradientScale·gsensPressureDrop`，`dfdx = objectiveGradientScale·fsensMeanT`。
  **b25 optProperties：两个 scale 均 =1.0**（干净，无隐藏缩放）。
- 源恒等式口径（B26a/B30）：`b_L^T w_true`（b18rhs mtx，PHYSICAL 未缩放 Q 型），
  w_true = b26_wstar_h7.npz（H7，闭合 1e-11）。
- z 方向（stageB6_rxc_z_analytic.mtx）：z = drho·y + dProjEta·[Σmask(1−drho)V·y/Σmask·dProjEta·V]，
  y 为滤波切线解 (L−diag(b))y = −b·d（stageB6RxDesignOracle.H L966-998）——
  **数学恒等式：<gsensh_pre, z> ≡ <gsens_post, D>**（同一算子的转置对，滤波算子对称）。
  ⇒ "z 投影 vs D 投影"不是独立口径，M4 数量级差不可能由该选择解释（预注册判断，见下）。

## 1. 输入核验（开工即查）

| 输入 | 状态 |
|---|---|
| `/home/ys/dsH/b25_qgate/1/` 52 场 | 在位；含 Uc/Ub/pb/pc/fsenshMeanT/fsensMeanT/gsenshPressureDrop/gsensPressureDrop/dfdx/dAlphaDxh |
| `/home/ys/dsH/b8_verify_diag/rxpr_prod_gsensh_{momentum,pressurerow,total}.mtx` + `rxpr_prod_gsens.mtx` | 各 33600 行，08-18 23:12:37 写入 |
| `b26_wstar_h7.npz` | 在位（闭合 1e-11 已验收） |
| `b25_qgate/b18rhs_{pressureDrop,thermalCoupling}.mtx` | 在位 |
| `stage_b2_fd_scan.tsv` + `stage_b2_summary.tsv` | 在位；ADJ_J=(−0.035835,−0.0249260,−0.0120030)，ADJ_gDP=(−5.23029,+1.0505250,+13.603423)，FD_gDP=(−2.536261@h=1e-5, +0.536366@h=1e-4, +6.480676@h=1e-4) |
| optProperties（b25_qgate） | objectiveGradientScale=1.0，pressureGradientScale=1.0，pressureDropMaxPa=25000，freezeColdFlowForValidation=false，stageB6RxDesignOracle=false |

## 2. 开工即发现（输入核验阶段的意外证据，先于预注册记录）

1. **b8_verify_diag/1 与 b25_qgate/1 的伴随场不一致**：
   - `Uc`：relL2(b8−b25)=**4.374e+02**（|Uc|L2: b8=4.698e+04 vs b25=1.074e+02，b8 大 ~437×）；
   - `pb`：b8 为 `uniform 0`（**b8 运行根本没解热标签压力伴随**），b25 为 nonuniform；
   - Ub/pc 差异待仪器量化（快测脚本在 Ub 处因 b8 侧尺寸/格式中断）。
   - B26 的逐位比对清单只含原始/设计场（alpha,U,p,phi,x,xh,xp,designMask,T,Tb,phiThermal,
     dAlphaDxh,dDTDxh），**未含任何伴随场**——"场状态与 b25 逐位相同"仅对原始态成立。
2. **时间戳**：rxpr mtx 与 b8/1/Uc、b8/1/gsenshPressureDrop 同秒写入（2026-08-18 23:12:37），
   b25/1/Uc 为 2026-08-21 12:56 ⇒ **rxpr 导出属于 b8 运行自己的伴随态，不是 state-B 伴随态**。
3. 推论（预注册为 M4 第一假设）：rxpr 导出 = b8 态 lambda 的分段；ADJ_gDP（tsv）= B2 门轮
   （b25 态）lambda 的投影。两者不是同一 lambda 的测量——M4 数量级差候选解释=**导出态错位（口径差）**，
   而非装配链真差。待仪器逐 cell 复算裁决。

## 3. 执行时间轴

- 2026-08-22：读任务书/定向包/仪器件（b27_m23_instrument.py、b26a_m1_instrument.py）；
  读 sensitivity.H（装配+导出+置零控制流）、filter_chainrule.H（链式规则全貌）、
  validateStageB2GradientAmplitude.H（方向+投影）、stageB6RxDesignOracle.H（z 构造）；
  输入核验 + §2 意外发现；写 NOTEBOOK+PREREGISTRATION。
- 下一步：b31_recon.py（M0 门控 → M1 三口径对齐 → M2 逐项对账矩阵）。

## 4. 未决问题

- rxpr_momentum 是否逐 cell 等于 b8 态复算（−dAlphaDxh·(U&Uc_b8)·V）？（判"导出=自身态自洽"）
- b25 态复算分段与源收缩 b_L^T w_true 是否闭合（gDP 三方向 0.82/0.37/1.04 缺口归因）？
- D3 28×：生产 dfdx 投影 vs 离线四段和的缺口落在哪一段/哪一层？

## 5. M0M1 结果（2026-08-22，运行日志 b31_m0m1_run.log）

门控：V 金字塔 relL2=2.719e-14 PASS。b8 伴随态盘点：|Uc8|L2=4.698e+04、|pc8|L2=1.983e+07、
Ub8=uniform(0,0,0)、pb8=uniform 0 —— **b8_verify_diag 运行只解了 pressureDrop 标签伴随，
且其 Uc/pc 与 state-B 完全不同量级**。

### 逐 cell 比较（b31_m0m1_compare.tsv）
| 比较 | relL2 | 结论 |
|---|---|---|
| rxpr_momentum vs −dAlphaDxh(U&Uc_b8)V | **2.844e-12** | 导出≡b8 态复算（机器级） |
| rxpr_momentum vs 同式(b25 Uc) | 1.009（2756 符号翻转） | 与 state-B 无关 |
| rxpr_total vs gsenshPressureDrop(b25 盘面) | 1.002 | 非 state-B |
| rxpr_gsens(后链,b8run) vs gsensPressureDrop(b25) | 1.041 | 非 state-B |

### 投影三口径表（b31_m4_projections.tsv）
| dir | rxprT·z(B26a 复现) | b25盘面前链·z | 后链盘面·Dk==tsv projDP | dfdx·Dk==tsv projJ |
|---|---|---|---|---|
| D1 | −5.624e-02 ✓复现 | **−16.116** | −5.230294 ✓ 全位 | −0.03583504 ✓ 全位 |
| D2 | +2.977e-01 ✓复现 | **−3.525（与后链反号！）** | +1.050525 ✓ | −0.02492599 ✓ |
| D3 | −3.634e-02 ✓复现 | **+80.410** | +13.60342 ✓ | −0.01200301 ✓ |

### 中期裁决
1. **M4 第一道闸通过（预注册 H1 成立）**：rxpr 系列 = b8 运行自身 lambda 的分段导出
   （机器级证实），非 state-B 装配段；B26a M4 数量级差定性为**导出态错位（口径差）**，
   rxpr 系列从证据链降级。M4 开放差异关闭。
2. **生产 tsv 腿闭合**：projDP/projJ 与盘面终场投影逐位一致 ⇒ 生产投影端无仪器问题。
3. **新开放项（升级为下一闸）**：<fsensh_PD_b25, z_file> ≠ <post, Dk> 且 D2 反号。
   z_file 与生产链不自洽（z 过期 or 生产链真缺陷）。判别 = 离线自建链式复算
   （drho/eta5 二分/Helmholtz b=1）作用于盘面前链场：若复算后链投影 == 盘面后链投影
   ⇒ z_file 失效（H3，仍属口径）；若不等 ⇒ 生产链真缺陷（H2，病灶证据）。

下一步：m2 阶段（explicitJT_H7 压力行段重算 + C/flux 段重算 + 链式复算 + 源收缩 + 对账矩阵）。

## 协调人中止注记（attempt 1，18:45）

executor-temp 首战失败：17:56 在 m2 形状崩溃（prow_gdp = -T_pd*dAlphaDxh，134400 vs 33600）
后未修复、空转至耗尽，无 SUMMARY/DONE/提交。协调人诊断：崩溃行不止形状错，
**逐元素乘的代数本身错误**——压力行项须按 stageB6 钉扎走 assembleWeightedRPa 加权装配
（w=dAlphaDxh·方向；anRPd=div(dphiHbyA_w − dflux_w)），构件在 b29_slot_instrument.py
（interp/drAU/dHbyA 路由件）与 B26 attempt-1 NOTEBOOK 的源码钉扎节均有。
已 salvage 的有效产出：M4 跨代伪影销案（b8 重算 2.8e-12 vs b25 场 1.01+2756 翻转——
b8 伴随场是 B18/H7 之前老求解器解的）、生产投影磁盘锚（dfdx/fsensMeanT/gsensPD 与
tsv 逐位）、m0m1 口径表、stage cache。attempt 2 以 b31b_ 前缀续作。

---

## 6. b31b_ 压力行项代数推导与符号约定（attempt 2 动手前强制先行）

### 6.1 任务要点（协调人钉扎）
- 禁止逐元素乘捷径：`prow = -T*dAlphaDxh` 中 T 若为 NUNK(134400) 维向量，形状错（134400 vs 33600）
  且代数本身错。
- 正确机制 = assembleWeightedRPa（stageB6 L426-519）：w（N 维逐 cell 权重）→ 面通量差 → div → N 维
  压力残差 anRPd；压力行贡献 = ⟨−λ_P, anRPd⟩。
- 转置 T = R_P^T λ_P 是 **N 维逐 cell 场**（面循环转置，非显式矩阵乘积）；逐 cell 乘 dAlphaDxh 才合法。

### 6.2 线性化点（与生产一致，关键）
- 生产 rxPressureRowTranspose.H 用 **当前 b25 态**（非 stageB2 基线）重建减化 SIMPLE（L150-188）：
  - rxUEqn = fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U) == −fvc::grad(p)
    （若 ransFlowModel 加 −div(nuEffFrozen·dev2(T(gradU))))；relax；fvOptions.constrain。
  - rxrAU = 1/A(U)，rxHbyA = constrainHbyA(rxrAU·H, U, p)。
- 实测（本会话）：b25 当前态 vs stageB2 基线 relL2 —— U=0.136、p=0.106、phi=0.134、alpha=1.6e-16。
  ⇒ 压力行必须用 b25/1/U,p,phi,alpha + nuEffFrozen（uniform 5.19009e-05），**不得用 wstate_baseline_***。
- 微分系数（stageB6 L407-425，ALPHAREL=0.4）：
  - drAU = −rAU_rel²/ALPHAREL（逐 cell）
  - dHbyA = (rAU_rel/ALPHAREL)·((1−ALPHAREL)·U − HbyA)（逐 cell，3 分量）

### 6.3 正向算子 R_P（assembleWeightedRPa，stageB6 L426-519）
对逐 cell 权重 w（N 维）：
- dHbyA_w = dHbyA·w（逐 cell 3 分量）；drAU_w = drAU·w（逐 cell）
- 面插值：dphiHbyA_w[f] = Sf&(wf·dHbyA_w[own] + (1−wf)·dHbyA_w[nei])
- dflux_w[f] = g0[f]·(wf·drAU_w[own] + (1−wf)·drAU_w[nei])，g0 = flux(fvm::laplacian(one,p))
  （scheme 精确面梯度，系数线性 ⇒ 上式对 drAU_w 线性成立）
- anRPd[cell] = div 收缩（owner +/ neighbour −，面循环）

### 6.4 转置 T = R_P^T λ_P（rxPressureRowTranspose.H 面循环，N 维）
- λ_P = pc（gDP 标签）或 pb（J 标签）；λdiff[f] = λ_P[own] − λ_P[nei]
- 通道 T1（HbyA）：T1[own] += wf·(Sf&dHbyA[own])·λdiff；T1[nei] += (1−wf)·(Sf&dHbyA[nei])·λdiff
- T1 边界（仅 U.boundaryField()[patchi].assignable() 的 patch）：T1[b] += (Sf_b&dHbyA[cell])·λ_P[cell]
- 通道 T2（g0·drAU）：T2[own] += −wf·g0[f]·λdiff·drAU[own]；T2[nei] += −(1−wf)·g0[f]·λdiff·drAU[nei]；
  T2[b] += 0
- T = T1 + T2（N 维逐 cell 场）
- 灵敏度：gsenshPressureDropPressureRow = −T·dAlphaDxh（N 维逐 cell 乘，合法）

### 6.5 代数论证（为何逐元素捷径错、面循环转置对）
- 方向 Dk 压力行灵敏度 = ⟨λ_P, anRPd⟩，anRPd = R_P·w，w = dAlphaDxh·Dk。
- ⟨λ_P, R_P·w⟩ = ⟨R_P^T λ_P, w⟩ = Σ_cells T[cell]·w[cell]。
- R_P 是"权重→div 通量"算子：其转置把 λ_P 沿同一面通量结构反向传播回 cell，形状天然 N 维；
  显式 explicitJT_H7（134400×134400）乘积把速度/压力分块混叠且维度错位，不可用。

### 6.6 动量项（B13 核对：J 用 Ub、gDP 用 Uc）
- mom_gDP = −dAlphaDxh·(U&Uc)·V（sensitivity.H L89-90，pressureDrop 标签）
- mom_J   = −dAlphaDxh·(U&Ub)·V（sensitivity.H L65，thermalCoupling 标签）
- 缓存交叉校验目标（前链逐 cell）：mom_J + prow_J + Gx + C ≈ fsenshMeanT（relL2<1e-6）；
  mom_gDP + prow_gDP ≈ gsenshPD（relL2<1e-6）。

### 6.7 J 链附加项（仅 J 标签）
- Gx = rxFluxDirectT·dAlphaDxh（L255，正号；仅 RX_FLUX_DIRECT_SOURCE=thermalCouplingFaceFunctional 定义时）：
  internal own += wf·((Sf&dHbyA[own]) − g0[f]·drAU[own])·g[f]；
  nei += (1−wf)·((Sf&dHbyA[nei]) − g0[f]·drAU[nei])·g[f]；
  边界 assignable += (Sf_b&dHbyA[cell])·g[b]
- C = thermalDiffusionDerivativeDTCell（AdjHeatTransfer.H 面公式，**不乘 V**）：
  C[own] += −w·dd_o·base；C[nei] += −(1−w)·dd_n·base；
  base = (Tb[own]−Tb[nei])·(T[own]−T[nei])·dc·magSfI

### 6.8 对账矩阵（2 标签 × 3 方向 Dk）
- 腿 1 生产 tsv：dot(dfdx,Dk)（ADJ_J 逐位锚 −0.0358350/−0.02492599/−0.01200301），
  dot(gsensPD,Dk)（ADJ_gDP 逐位锚 −5.230294/+1.050525/+13.60342）。
- 腿 2 离线分段：Σ_items chain_apply(item_i)·Dk，item ∈ {mom, prow, (J 另 +Gx, +C)}；
  chain_apply 必须过 G0/G1/G2 门（forward filter relL2、eta5 二分、xh 投影复现）。
- 腿 3 源收缩：b_PD^T w_true、b_TC^T w_true（B26a/B30 闭合 1e-11）。
- 目标：三腿闭合 ⇒ H2/H3 裁决；defcheck 缺口（b_PD^T w_true/ADJ_gDP = 0.8178/0.3660/1.0362）
  与 D3 28× 定量归因（缺口落在哪一段/哪一层）。

### 6.9 新增 pinning（attempt 2，2026-08-22）：reduced-SIMPLE 重建机制 + b8 pc8 门设计
**A. reduced-SIMPLE 重建（rxPressureRowTranspose.H L150-188，LOCKED 形式）**
- rxUEqn = fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U) == −fvc::grad(p) + fvOptions(U)
- ransFlowModel=TRUE（flowModel incompressibleRANSFrozen）→ rxUEqn −= fvc::div(nuEffFrozen·dev2(T(grad(U))))。
  注意 fvMatrix::operator-= 作用于 fvc 场 = source() −= div(dev2)，故源含 −div(dev2)（与 B18 模板一致）。
- nuEffFrozen = 5.19009e-05 全均匀（nutFrozen ≡ 0 于 b25 与 b8，已验证）。
- fvOptions 文件不存在 → fvOptions(U)=0。
- rxUEqn.relax()：D_rel = D0/alphaRel；A() = D_rel/V；rAU = 1/A() = V/D_rel（携带 V）。
- rxAlphaRel = mesh.equationRelaxationFactor("U") = 0.4（b25 fvSolution relaxationFactors U=0.4）。
- drAU = −rAU²/alphaRel（**NO extra V**）；dHbyA = (rAU/alphaRel)·((1−alphaRel)·U − HbyA)。
- 离线矩阵实现 = B18 build_dHbyA_drAU 模板（relax 语义：src_vec += (1−alphaRel)·D_rel·U；offU 经
  upper/lower 面系数；含 dev2 均匀-nu 项；边界 inflow −phi_b·U_b、fixed-U nuEff·magSf·dc·U_b）。

**B. T = R_P^T λ_P 面循环转置（L223-257；T1 边界项仅 assignable 面）**
- internal：pcDiff = λ[own]−λ[nei]；T1[own]+= wf·(Sf&dHbyA[own])·pcDiff；
  T1[nei]+= (1−wf)·(Sf&dHbyA[nei])·pcDiff；T2[own]+= −wf·g0[f]·pcDiff·drAU[own]；
  T2[nei]+= −(1−wf)·g0[f]·pcDiff·drAU[nei]。
- boundary：仅 U.assignable() 面：T1[celli]+= (Sf_b&dHbyA[celli])·λ[celli]；T2[b] += 0。
- g0 = flux(fvm::laplacian(one,p)) 全面；g0_int = magSf·dc·(p[nei]−p[own])（B18 符号，已验证）。
- T = T1+T2，N 维。生产项：gsenshPressureDropPressureRow = −T·dAlphaDxh（NO extra V，sensitivity.H L115）。

**C. b8 pc8 机器精度门（唯一直接锚；b25_qgate stageB6RxDesignOracle=false 不导出）**
- b8_verify_diag/rxpr_prod_gsensh_pressurerow.mtx（33600-dim，md5 f06cd716ad842f8fc2b19333730deec2，
  L2=0.0145680567705，max=±0.003418029543172264）= −T(pc8)·dAlphaDxh（非裸 T）。
  5 个拷贝 md5 相同（b7_scratch_pre/on、b8_scratch_prod/v1、/home/ys/dsH/）。
- stageB6_lambda.mtx（134400-dim，L2=19831720）布局：[:100800]=Uc8、[100800:]=pc8（rel 1.1e-12/1.3e-12
  vs b8 磁盘场）。pc8 L2=19831664（≠ b25 pc25，b8 门必须用 pc8）。
- b8 原始态与 b25 逐位一致（U/p/phi/alpha/nuEffFrozen/designMask rel=0），dAlphaDxh b8=b25。
- 运行点锚（Log.verify_diag.txt L2080）：sum(T*w)=−2.696485451749351e-08，|T|L2=0.000256656999292，
  |J_P*w|L2=6.22344977639e-10，|g0_bnd|max=0.0650620425279，w-support=5040，relErr=2.58e-15。
- w = sin(0.271·(celli+1)) 于 support {designMask>0.5 AND alpha > alphamin·(1+1e-10)}，
  alphamin=1e-4（transportProperties L17）。support 规模 5040 即第三个 gate 数。
- 门关闭：dot-T 数值、|T|L2、−T(pc8)·dAlphaDxh vs 生产 mtx（relL2<1e-8），三者全过。
- **pb 侧无直接生产锚**：b25 不导出 rxpr（false）、b24_diag rxj_*.mtx 为 NaN、b8 日志仅一条 [pc] 点测
  线（无 [pb]）。J 侧压力行只能经对账矩阵隐式一致性验证：fsensh − C − Fx − mom_J ≈ prow_J。

**D. assignable() 语义更正（2026-08-22 源码级核验，修正此前 "ALL patches EXCEPT outlet" 的旧注）**
- 旧注（§6.9B 及 6.4）说 assignable = 除 outlet 外全部 patch —— **错误，已被源码否决**。
- 源码事实（OpenFOAM 7）：
  - 基类 `fvPatchField::assignable()` **返回 TRUE**（fvPatchField.H L308-311，`return true;`）。
  - `fixedValueFvPatchField` 覆写为 **FALSE**（fixedValueFvPatchField.H L154-157）。
  - `slipFvPatchField` 覆写 FALSE（slipFvPatchField.H L136-139）；`sliced` 覆写 FALSE（slicedFvPatchField.H L142-145）。
  - `zeroGradient` **无覆写 → 继承基类 TRUE**。
  - `noSlip` 继承 `fixedValueFvPatchVectorField` → FALSE。
- b25 U 边界类型实测：inlet=fixedValue、outlet=zeroGradient、hotInlet=fixedValue、hotOutlet=fixedValue、
  solidEndWalls/bottomWall/topWall/sideWalls=noSlip。
- ⇒ **assignable = 仅 outlet（84 面）**。交叉验证：b16mesh_boundary.mtx col7 Ufix 非零计数 7796，
  7880−7796 = 84 = outlet 面数（inlet 84）。Ufix==1 ⇔ 非 assignable；Ufix==0 ⇔ assignable。
- 影响：T1 边界项（rxPressureRowTranspose.H L230-236）只在 **outlet** 面加
  `(Sf_b&dHbyA[celli])·λ[celli]`；T2[b]=0 不变。dot-test 正向复现的边界 dphiHbyA 同理仅 outlet。
- 与 b29 路由件注释一致（"boundary (assignable-U only = outlet)"），与 b31_recon.py 的
  `am = np.where(~Ufix)[0]` 选取一致（Ufix=1→非 assignable，取反→outlet）。

### 6.10 G0 门根因定案（2026-08-22 源码级，attempt-1 relL2=2.62e-3 已解释）
- createFields.H L404-410：`dimensionedScalar b("b", dims(0,-2,0,0,0,0,0), 1.0)` —— **b=1.0 均匀量**。
- filter_x.H：`Eqn4 += x*b;`（逐字段 add，非 fvm::Sp）。
- fvmSup.C L167-168（fvm::Sp(b, xp)）：`fvm.diag() += mesh.V()*sp.value()` = **V·1.0 加入对角**。
- fvMatrix.C L1081（operator+=(DimensionedField)）：`source() -= su.mesh().V()*su.field()` =
  **源减去 V·x**。
- ⇒ 正确离散滤波方程：**(L_int − V·I)·xp = −V·x**（L_int = 内部 Laplacian，
  negSumDiag 对角，upper=kfv；maskF 边界 0 → 无边界对角贡献）。
- attempt-1 的 b31_recon.py 用了 (L−I)xp=−x —— 缺 V 缩放 → relL2=2.62e-3 失败。
- **修正**：A = (L_int − V·I)；正向门 RHS = −V·x_full；chain_apply 的 RHS = −V·mid，
  解出 post = lu.solve(−V·mid)·designMask。
- fvMatrix::flux()（fvMatrix.C L863-930）= faceH(psi)：内部 = upper·(p[nei]−p[own])，
  边界 = internalCoeffs·patchInternal − coupled·patchNeighbour → 无对角贡献。
  证实 g0_int = magSf·dc·(p[nei]−p[own])，g0_bnd 仅诊断。
- 运行点确认：b8 pc8 门 anchors（Log L2080）sum(T*w)=−2.696485451749351e-08、
  |T|L2=0.000256656999292、w-support=5040；rxpr_pressurerow.mtx md5
  f06cd716ad842f8fc2b19333730deec2（= −T(pc8)·dAlphaDxh，L2=0.0145680567705）。
- stageB6_lambda.mtx（134400-dim）布局确认：[:100800]=Uc8、[100800:]=pc8（pc8 L2=19831664）。

### 6.11 dflux_w 精确离线复现 + b31b 运行计划（2026-08-22，源码级）
- fvMatrix::flux()（fvMatrix.C L863-930）：内部 faceH = upper·(psi[nei]−psi[own])；
  边界 = internalCoeffs·patchInternal − (coupled? boundaryCoeffs·patchNeighbour:0)；
  faceFluxCorrectionPtr 仅当 mesh.fluxRequired(vf.name()) 时设置 —— **p 非 fluxRequired，
  非正交修正不进 flux()**（b25 mesh 近正交，b18 已验证修正 1e-12 量级）。
- fvmLaplacianUncorrected（gaussLaplacianSchemes.C L33-87）：fvm.upper() = deltaCoeffs·gammaMagSf，
  gammaMagSf = gamma·magSf（L45-48）；negSumDiag 对角；边界 internalCoeffs =
  pGamma·gradientInternalCoeffs、boundaryCoeffs = −pGamma·gradientBoundaryCoeffs
  （gaussLaplacianScheme.C L63-84）。
- fvmLaplacian(volScalarField gamma) → laplacianScheme.C L97-103 调 interpolate(gamma)
  → 面插值 surfaceScalarField（默认 linear，mesh weights）。fixedValue p 的
  gradientInternalCoeffs=−deltaCoeffs、gradientBoundaryCoeffs=deltaCoeffs·value
  （fixedValueFvPatchField.C L128-145）。
- **dflux_w 离线公式**（drAUwField 内部 uniform 0、边界 calculated 默认 0）：
  - dflux_w[fi]（内部）= magSf·dc·(w·drAU[own]·wf[own] + (1−w)·drAU[nei]·wf[nei])·(p[nei]−p[own])
  - dflux_w[边界] = 0（gammaMagSf_bnd=0）
  - 与 b18 g0 捷径同形 ⇒ **dflux_w = g0_int·interp(drAU·w)** 为精确形式
    （Gx 捷径与正向 J_P*w dot-test 均可用此式）。
- b31b 运行计划（对账矩阵三大列）：
  1. 生产 ADJ tsv（b31_m4_projections.tsv：ADJ_gDP/D1-3、ADJ_J/D1-3）。
  2. 离线分段重算：mom_J/mom_gDP（b18 L246-299 当前态动量矩阵，ALPHAREL=0.4，
     phi_bvals 走 b29 read_field 路由）；prow_J/prow_gDP = −T·dAlphaDxh
     （T=面循环转置，rxPressureRowTranspose.H L223-257，T1/T2 + 边界仅 assignable=outlet）；
     C 场（b31_recon L434-444 配方）；Fx 场（b18 Gx_of，自门 GX_ANCHOR）；
     fsensh = mom + prow + C + Fx（J 再 + g_b 直连经 chain）。
  3. 源收缩：b_PD^T w_true、b_TC^T w_true（load_b 布局 [0:NV:3]=mx,[1:NV:3]=my,
     [2:NV:3]=mz,[NV:]=mp；w_true=b26_wstar_h7.npz 4N 全长）。
- chain 算子（G1/G2/G0 门）：A=(L_int−V·I)；chain_apply=lu.solve(−V·mid)·designMask；
  正向门 RHS=−V·x_full；每段 chain_apply 后投影 Dk。
- 三段源收缩自洽检查：b_TC^T w_true 应 ≈ ADJ_J_tsv（若 fsensh 组装自洽且 b 正确）；
  b_PD^T w_true vs ADJ_gDP_tsv 差距 0.82/0.37/1.04 为 defcheck 已知缺口，需定量解释或如实报断点。

### 6.12 强制收尾：当前状态、断点与续作路径（2026-08-22 21:0x）
协调人强制收尾令（上下文压缩频发征兆，NOTEBOOK 20:38 后停更）。先笔记后数字，立即执行。

**已完成（可核验成果）：**
1. 压力行项代数推导 §6.9（面循环转置 T=T1+T2，T1[own]+=wf·(Sf&dHbyA[own])·pcDiff、
   T1[nei]+=(1-wf)·(Sf&dHbyA[nei])·pcDiff、T2[own]+=-wf·g0·pcDiff·drAU[own]、nei 对称，
   边界仅 uAssignable=outlet：T1[cell]+=(Sf_b&dHbyA[cell])·λ[cell]、T2[b]=0；
   生产 prow=−T·dAlphaDxh；符号约定 T=J_P^T λ、pcDiff=λ[own]−λ[nei]）。
2. dflux_w 精确离线复现 §6.11（dflux_w[fi]=g0_int·interp(drAU·w)，边界=0，
   与 b18 g0 捷径同形）。
3. G0 门根因 §6.10：生产是 (L_int−V·I)·xp=−V·x；attempt-1 的 (L−I)xp=−x 错误
   （relL2=2.620e-03，chain fieldrel 5.85/10.1）。正确：diag−V、RHS=−V·x_full、
   chain_apply=lu.solve(−V·mid)·designMask。
4. 数据核验（全部 VERIFIED）：b8==b25 状态字段（U/p/phi/alpha relL2=0，dAlphaDxh 全同）；
   pc8 L2=19831663.927（stageB6_lambda.mtx[100800:]），Uc8=lam[:100800].reshape(33600,3)；
   rxprP L2=0.0145680567705、rxprM L2=0.0373267050442、rxprT=rxprP+rxprM（9.9e-18）；
   w-support=5040；b8 门 anchors（Log L2080）：sum(T*w)=−2.696485451749351e-08、
   relErr=2.57680498181e-15、|T|L2=0.000256656999292。
5. H1 裁决（b31_m0m1_compare.tsv）：rxpr_momentum vs mom_recomp_b8 relL2=2.844e-12（MATCH），
   vs mom_recomp_b25 relL2=1.009（MISMATCH）⇒ rxpr 导出与 b8 运行 lambda 自洽，非 state-B。
6. M4 三口径投影（b31_m4_projections.tsv）：D1 dot_gsensPDdisk_D=−5.2302938463e+00（=ADJ_gDP）、
   dot_dfdx_D=−3.5835040049e-02（=ADJ_J）；D2 dot_gsenshDisk_z=−3.5247766837e+00（与 ADJ_gDP 反号）；
   D3 dot_gsenshDisk_z=+8.0409993473e+01 vs dot_gsensPDdisk_D=+1.3603422582e+01。
7. attempt-1 跑日志（b31_m2_run.log）：G1 eta5=0.7480051615837531、G2 xh relL2=1.301e-10
   （投影复制门通过）；G0 门失败 relL2=2.620e-03（已定位根因，见 §6.10）；
   chain D1 gdp recomp=−1.935047e+01 vs 盘 −5.230294e+00（fieldrel 5.85e+00）、
   J recomp=−2.787855e-01 vs −3.583504e-02（fieldrel 1.01e+01）；随后崩溃于 b31_recon.py L430
   broadcast (134400,)/(33600,)（压力行 −T_pd·dAlphaDxh 用错显式 JT_H7 张量）。

**未完成（断点精确位置）：**
- b31b_m2.py **未写出**（强制收尾先于实现；attempt-2 全程做了源码钉扎与数据核验，
  未进入代码编写）。
- 两标签三方向对账矩阵 b31b_matrix.tsv **未构建**。
- b8 pc8 机器精度门（当前态重建 T 后 sum(T*w) vs −2.696485451749351e-08）**未执行**。
- defcheck 缺口 0.82/0.37/1.04 与 D3 28× 的定量解释 **未完成**（已知数值：
  b_PD^T w_true/ADJ_gDP = 0.8178/0.3660/1.0362；FD_J D3≈−0.000445 vs ADJ_J −0.0120）。

**续作路径（任何人可从笔记接手，b31b_m2.py 实现蓝图已定）：**
1. 内联 b29 read_field（含 boundary 解析，约 80 行）读当前态 b25/1 的 U/p/phi/alpha、
   phi_bvals/p_bvals/U_fixed/U_bvals（b18 L204-225 路由模式）；mesh 用 b31_recon.load_mesh
   （V 门 relL2<1e-12）。
2. 动量矩阵当前态（b18 L246-299，ALPHAREL=0.4）：nuEff_f=w·nuEff[o]+(1-w)·nuEff[n]、
   a_f=nuEff_f·magSf·dc、po/pn_、diag_u+=alpha·V−divphi、iC_b 边界、
   sumOff、D_rel=max(|D0|,sumOff)/ALPHAREL、mob=V/D_rel；upper=pn_+a_f、lower=−po+a_f。
3. build_dHbyA_drAU 当前态（b18 L462-513）：offU、gradp（outlet p_bvals=0）、dev2、
   src_vec=(−gradp−divdev)·V+src_bcomp+(1−ALPHAREL)·D_rel·U、HbyA=mob·(src_vec−offU)/V、
   drAU=−mob²/ALPHAREL、dHbyA=(mob/ALPHAREL)·((1−ALPHAREL)·U−HbyA)。
4. T=面循环转置（rxPressureRowTranspose.H L223-257）分别对 pc25/pb25/pc8：
   T1/T2 内部循环 + 边界仅 ~Ufix（am=np.where(~Ufix)[0]，outlet 84 面）；
   g0_int=magSf·dc·(p[nei]−p[owner])；prow=−T·dAlphaDxh。
5. b8 pc8 门：w=sin(0.271·(celli+1)) on {designMask>0.5 AND alpha>1e-4·(1+1e-10)}
   （support=5040）；门 sum(T_pc8·w)、|T_pc8|L2、prow_b8 vs rxpr_prod_gsensh_pressurerow.mtx
   relL2<1e-8。
6. mom_seg=−dAlphaDxh·einsum(U,Uc25)·V（gDP）、−dAlphaDxh·einsum(U,Ub25)·V（J）；
   C 场（b31_recon L434-444）、Gx 场（b31_recon L536-556 配方但 dHbyA/drAU 用当前态；
   自门 GX_ANCHOR D1 −1.5704339057e-03/D2 4.0431596817e-03/D3 8.4067076831e-03）；
   g 泛函（b31_recon L453-477，outlet flux 解析 + wstate_outletMeta.mtx）。
7. chain（§6.10 修正版）：A=(L_int−V·I)、chain_apply=lu.solve(−V·mid)·designMask、
   正向门 RHS=−V·x_full；源收缩 b_PD^T w_true/b_TC^T w_true（load_b 布局
   [0:NV:3]=mx,[1:NV:3]=my,[2:NV:3]=mz,[NV:]=mp；w_true=b26_wstar_h7.npz）。
8. 对账矩阵三列（生产 ADJ tsv | 离线分段重算 | 源收缩）写 b31b_matrix.tsv；
   defcheck 缺口定量解释或如实报断点；收尾 EXECUTOR_SUMMARY + 空 DONE + 选择性 git commit。

### 6.13 强制收尾补充：可救数字（2026-08-22 21:0x，直接复用缓存，无新开发）
跑了一个纯复用式的秒级 salvage（b31b_salvage.tsv/.json），给出对账矩阵的两条腿（源收缩 + 动量段）：
- 源收缩 b_PD^T w_true（xADJ_gDP）：D1 −4.27746124e+00（0.8178）、D2 +3.84462806e-01（0.3660）、
  D3 +1.40963533e+01（1.0362）——与 B30 defcheck 已知缺口 0.82/0.37/1.04 完全一致，确认缺口是
  真实存在且可复现（压力行段未算入时 b_PD 与 ADJ_gDP 的差）。
- 源收缩 b_TC^T w_true（xADJ_J）：D1 −4.88401870e-02（1.3629）、D2 −6.44099895e-03（0.2584）、
  D3 +1.15666175e-02（−0.9636）——与 ADJ_J 不同量级且符号不齐（D3 反号），说明 J 标签的
  对账在 mom+prow+flux+C 全段齐备前不能用单腿源收缩闭合（B13 已注明 b 为流中介源）。
- 动量段投影（pre，未过 chain）：mom_gDP preD = −3.128852e+00/+1.409900e+00/+2.193815e+00
  （D1/D2/D3）；mom_J preD = +2.044486e-02/−7.990468e-03/−1.206465e-02。
- D3 28× 复核：FD_J(D3)=−4.4495118952e-04 vs ADJ_J(D3)=−1.2003012269e-02，比值 26.97×
  （≈28× 与先前期一致，FD 与 ADJ 在此方向差一量级多）。
以上数字未过 chain（G0 门修正后未重跑），仅作 salvage；对账矩阵 b31b_matrix.tsv 未构建（断点见 §6.12）。
