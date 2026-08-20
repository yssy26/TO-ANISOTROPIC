# BFINAL-021 cycle-1 — FINAL REPORT（P 行表示级差异定位轮）

任务：BFINAL-020 停止后的跟进——用零/近零代码实验把「J 的 P 行切线与真实
系统响应之间的表示级差异」定位到具体机制，或证明它是量具伪影。四实验
（E1 构造敏感性审计 / E2 eps 阶梯 / E3 T1 级外部锚 / E4 表示一致性审计）。
分支 `agent/dsH-stage-b-validation`，起点 `45084c7`，执行 2026-08-21。

---

## 0. 结果速览

| 项 | 结果 |
|---|---|
| E1(a) 通量表示 | 扰动态**从未导出通量**（B20 的 dφ_state 必为算子重建）；基线导出 phi = 最终压力修正通量（\|contRes\|max=1.527e-14） |
| E1(b) 伪影份额 | ≈0：探针通量 div-free 至 1e-13；w(h) 阶梯自收敛 0.04–1.9%；恒等式 eps 平坦（<0.6%） |
| E1 意外主发现 | **A-B 线性化点/物理系统失配**：U 13.57%/p 10.63%/phi 13.41%（量化其对算子系数影响小） |
| E2 阶梯+对照 | 恒等式（精确 J 行）\|r_P\|/\|rxd_P\|=1.66/2.14/1.51 eps 平坦；干净 w* 解以 3e-11 满足恒等式 → 失配 ≡ J_P(A)(w_state−w*(A))，relP 0.76/1.02/1.17 |
| E3 外部锚 | 重跑（23 min，FD 逐位复现）+ 新通量/粘性导出；面级总通量切线：D1/D3 达仪器地板 3.4%（A 基）、8.7%（B 基）；D2 45%（两基） |
| E4 槽位清单 | 同点表示全部一致（5 项核查通过）；失配=跨点量具（SLOT-1）、stale kf 基（SLOT-2）、离线工具 da 链 bug（SLOT-3）、闭环缺失（SLOT-4）、面级切线残差（SLOT-5） |
| 代码 | 唯一新增 = stageB15StateExport 门控导出（+73 行，默认关）；零语义改动；1 次求解器重跑（标准流程） |

---

## 1. E1(a) — 恒等式构造敏感性：导出通量的表示

1. **扰动态无通量导出**：`validateStageB2GradientAmplitude.H` 的
   stageB15StateExport 块只写 U/p/T（B15/B16 轮新增）；`b13_probe3/stageB2/`
   的 wstate_* 文件清单证实。⇒ BFINAL-020 的「dφ_state」必然是算子切线
   （Psi 图）作用于状态差分 (dU,dp) 的重建，而非真实通量差分——审阅者
   C1(a) 的「最终 phi vs phiHbyA」歧义对扰动态**不存在**。
2. **基线 phi = 最终压力修正通量**：`b16_states/1/phi` 的 |contRes|max =
   1.527e-14（b20_derivation_check.log 既有数字，本轮复核）——phiHbyA
   的散度应为动量残差水平（O(1e-3)），不可能达到机器零。E3 重跑的新导出
   进一步证实（§4）。
3. b16_states/1 与 b15_export/1 逐位相同（U/p/alpha maxabs=0）——两者是
   同一 A 态（B17 已证，本轮复核）。

## 2. E1 数据考古 — A-B 线性化点/系统失配（本轮主发现）

### 2.1 事实链

- `wstate_baseline_*`（b13_probe3/stageB2，FD 差分参考）与 b16_states/1
  （B15–B20 全部离线量具的线性化点）相差 **U relL2 13.57%、p relL2
  10.63%（p 均值偏移 −4800 Pa）、phi relL2 13.41%**。
- b13_probe3/0（初始态）≈ b15_export/1（relL2 1.1e-7/8.3e-9 = ASCII 写
  盘精度）——探针算例从 A 态拷贝启动；b13_probe3/1 与 b16_states/1 逐位
  相同（主循环在 B2 模块前写盘）。
- 机理（源码级）：
  - **A 态**：主循环 `freezeTurbulenceForValidation=true`（b13_probe3/
    constant/optProperties:52）→ NS.H 跳过 turbulence->correct 且
    nutFrozen≡0、nuEffFrozen=分子常数 5.19009e-05 → NS.H 44 个
    corrector 后写盘（分子粘性系统、欠收敛）。
  - **B 态**：stageB2 模块 phase-1 `solveFullSST(1000,1e-6)`（freeze=
    false → nutFrozen=chi·nut 每步更新，冷区 nut 至 1.912e-3 = 37×
    分子粘性）246 次收敛后冻结 → 全部 FD 探针在 B（湍流粘性系统）。
- ⇒ **B15–B20 的全部离线仲裁（w′ vs w_true、λ-direct 1.77/0.86/2.27、
  P 行恒等式、b_TC^T w_true）都是「A 点算子 vs B 点 FD」的混合比较。**
  这解释了「为何七个算子耦合候选全灭」的量具学背景。

### 2.2 该失配的量化强度（E3 后获得）

用新导出的 nuEffFrozen(B)/phi(B) 重建 B 基并与 A 基对比：
- rAU_rel(B)/rAU_rel(A)：中位 1.0000，p10 0.990，min 0.560/max 1.572，
  仅 3.1% 面超 ±10%，无面超 [0.5,2]；
- kf(B)/kf(A) 同分布（3.3% 面超 ±10%）。
机理：Brinkman α·V（固体区 D 主导，α~1e8）与对流 |φ|（流体区 sumOff
主导）共同压制了 nut 的贡献 ⇒ **点/系统失配是真实的量具缺陷，但对算子
系数（进而 P 行切线）的直接影响小，单独不能解释 O(1) 恒等式失败。**

## 3. E1(b)/E2 — 伪影排除与 eps 阶梯

| 检验 | 数字 | 裁决 |
|---|---|---|
| 探针通量连续性（重跑直接导出） | \|div(φ+−φ−)\|L2 = 1.8e-14…1.2e-13（各方向各 h）；基线 \|div(φ_B)\|L2=8.3e-14 | 通量层伪影 ~5e-11 ≪ \|rxd_P\|≈2e-4（7 个数量级） |
| w(h) 阶梯自收敛 | \|w(3e-4)−w(1e-3)\|/\|w\|：U 2.3e-3/1.9e-2/4.0e-4，p 1.8e-4/1.0e-3/4.2e-4（D1/D2/D3） | FD 噪声 ≤1.9% |
| 恒等式 eps 平坦（精确导出 J 的 P 行） | \|r_P\|/\|rxd_P\| = 1.665/2.137/1.515（3e-4）→ 1.662/2.136/1.514（1e-3）→ 1.655（D1,3e-3）；cos = −0.17/−0.63/−0.06 | **系统性真失配签名**（噪声应随 1/h 移动） |
| 干净未钉扎 w*（SuperLU 1769 s） | r_P(w*) = 1.8/1.3/3.0 ×10⁻¹¹ | 机制验证：失配 ≡ J_P(A)·(w_state−w*(A))；relP 0.759/1.024/1.172, relU 0.180/0.526/0.190 |
| b15 wprime_cache 对照 | r_P(w′_cached) = 0.70/0.04/1.56 ≠ 0 | b15 缓存解带钉扎行伪影（B16 已知问题的 w′ 侧）；历史 w′ 数字小幅受染 |

**E1/E2 裁决：恒等式失败中伪影份额 ≈ 0（窗口收敛/连续性/eps 噪声全部
定量排除）；失败是系统性的 O(1) 真失配，且精确等于「真实系统响应与
A 点算子自身响应之差」在 P 行的投影。**

另发现（E2 副产物，SLOT-3）：b20_augmented_check.py 的 PsiA/Psi_alpha 用
`da = dAlphaDxh·rawD`，而生产 rxc 用 `da = dAlphaDxh·z`（z = 滤波+保体积
投影切线，|z|/|rawD| = 5.64/5.91/5.11，stageB6RxDesignOracle.H:1056-1077）
→ 其 alpha 通道量级小 ~7×（0.135/0.132/0.149，cos 0.93-0.96）。
BFINAL-020 的「含闭环修正 rxd_P_c（14% 幅度，cos 0.93）也不能闭合」等
alpha 通道类结论**作废**；closure_rel 4–7% 因 alpha 份额仅 4–18% 而近似
保持。

## 4. E3 — T1 级外部锚（重跑 + 面级通量响应量具）

### 4.1 重跑完整性

- 唯一代码新增：`validateStageB2GradientAmplitude.H` 两处
  stageB15StateExport 门控导出块（基线 + 探针态的 phi（内面+边界）、
  nutFrozen、nuEffFrozen、alpha；+73 行；默认关零影响）。
- 重跑 `b21_anchor`（b13_probe3 cp -a，新编译 build/bin/MTO_HF，标准
  流程），23 min，MTO_RC=0，FROZEN_GRADIENT_STATUS=FAIL（历史一致）。
- **FD 侧逐位复现**：stage_b2_fd_scan.tsv 的全部 FD_*/nIter 与存档相同
  （仅 ADJ_J 列不同——B18 fix3 源码演进影响伴随收缩，与 primal 无关）；
  wstate U/p 全部 bit-identical ⇒ 新通量导出精确对应历史量具的状态。

### 4.2 面级通量锚（本轮核心新量具）

真通量响应 dφ_state-FD = (φ+−φ−)/2h（重跑直接导出，无算子参与）vs
算子总切线（Psi_U dU + Psi_p dp + PsiA da，da 用正确的 z 链）：

| 方向/h=1e-3 | A 基（分子@A） | B 基（湍流@B） |
|---|---|---|
| D1 | cos 0.9994，relL2 **0.034** | cos 0.9972，relL2 0.087 |
| D2 | cos 0.9151，relL2 **0.451** | cos 0.9099，relL2 0.449 |
| D3 | cos 0.9994，relL2 **0.036** | cos 0.9979，relL2 0.090 |

- 离线重建仪器地板：绝对通量重建保真度内面 3.35%（A）/3.63%（B），
  出口边界面 40%/34% ⇒ **D1/D3 的 A 基匹配已达到仪器地板**；换到
  「正确」的 B 基不改善（反而 8.7-9.0%）——在仪器分辨率内 A/B 基不可
  区分，点失配对通量切线的影响被 (1−α)dUbar 主通道 + 系数不敏感进一步
  压制。
- **D2 在两个基下都有 45% 面级失配 = 真缺陷信号**（与 D2 的 relU 0.53、
  阶梯噪声 1.9% 最大、cos(w*,w_state)=0.982 最低一致）。
- div 层恒等式经我的重建在 B 基不闭（|r_P|/|rxd_re| = 3.14/2.29/3.51
  vs A 基重建 1.87/2.34/1.67）——但重建仪器（3.5% 面级）经 div 收缩
  放大后的读数不可分辨 A/B。
- **放大机制（定量自洽）**：面级 3.4% 的切线残差，其 div 分量即可达到
  ~1.7×|rxd_P|（|J_P w|≈|rxd_P| 的近对消结构；同 B20 在 U 行发现的
  1.5e-4 近对消放大）；压力响应经 Poisson 反解进一步放大到 relP O(1)。
  ⇒ **恒等式 O(1.5–2) 失败与 p 响应 O(0.76–1.17) 失配可以由 ≤3.5% 的
  面级通量切线残差完全承载**——该残差的来源（算子公式残差 / 重建地板
  / D2 状态质量）在零代码+单导出权限内不可再分解。

## 5. E4 — 表示一致性审计（槽位清单）

**一致（无缺陷）槽位**（算子假设表示 = 前向实际表示）：
1. P 行公式：`deltaPhiFacePU = Sf&(α·rAU_u·dH + (1−α)·dUbar) − kf·∇dp`
   （solveDiscreteFlowAdjoint.H:707-726）与 NS.H/validateCommon.H 的
   phiHbyA−pEqn.flux() 同点切线一致；rAUAdj=1/A(自重组未松弛矩阵)+α
   折算 = 前向 relax 后 rAU；kf 基=primalPressureMobility=前向 rAtU
   （B8 P1：J_PP kf-only vs 实际 pEqn.flux FD 1.44e-8）。
2. 边界：constrainHbyA 钉扎（固定 U 边界通量切线=0）与出口（唯一
   assignable）外推均已正确建模（:744-813）。
3. adjustPhi：出口 p=fixedValue → 无 fixedFluxPressure 面 → no-op。
4. 非正交修正：nNonOrthogonalCorrectors=1 + 正交网格 → faceFluxCorrection
   =0，单次公式=多次合成（B8 旁证）。
5. R_P 的 div：面散度（fvc::div(surfaceScalarField) 无 scheme 依赖）；
   stageB6 rxd（assembleWeightedRPa，真实 fvm::laplacian(drAU·z·w,p).flux()
   ）与算子 P 行同源同构。

**失配槽位**（表示级差异的具体形态）：
- **SLOT-1（量具·跨点）**：离线算子/rxd/rhs@A vs FD 真值@B
  （§2；真实但系数影响小，单独不充分）。
- **SLOT-2（量具·stale 成分）**：primalPressureMobility 仅在 NS.H:129
  捕获（A 点、分子基）；stageB2 模块 round-2 伴随在 B 点运行但 kf/J_PP
  通道沿用 A 基（mobB/mobA：中位 1.0，3.3% 面超 ±10%，极值 0.56/1.57）
  → 「生产 ADJ/FD=2.13×」量具受染（量级似不足主因，待修后重测）。
- **SLOT-3（离线工具 bug）**：b20_augmented 的 da 设计链缺失（§3）。
- **SLOT-4（既有量化）**：(I−Ψφ) 闭环缺失 4–7%。
- **SLOT-5（真残差·未定位）**：面级通量切线残差 ≤3.5%（D1/D3，=仪器
  地板）至 45%（D2），经 div 收缩/压力反解放大成 O(1)。

## 6. 联合裁决与缺陷地图更新

**裁决：混合。**
- 伪影 ≈ 0%（E1(b)/E2 定量排除）；
- 量具链缺陷（SLOT-1/2/3）= 本轮定位的新事实，量化后**单独不充分**；
- 真失配残余 = SLOT-5：百分比级面级通量切线残差 + 近对消放大机制，
  「表示级差异」从抽象概念收敛为可测的具体形态。

**缺陷地图（更新）**：
| 缺陷 | 形态 | 干净量具 | 状态 |
|---|---|---|---|
| 原「P 行 O(1) 失配」 | ≡ J_P(A)(w_state−w*(A))，由 ≤3.5–45% 面级切线残差经 div/压力反解放大 | 本轮面级锚 + 精确 J 恒等式 | 已定位到 SLOT-5 |
| 量具点失配 | A vs B（13.6%/10.6% + 分子vs湍流） | A-B 态比对 | 本轮坐实、量化 |
| 模块伴随 stale kf | round-2 伴随 kf 基@A | mobB/mobA 分布 | 本轮坐实、量化 |
| b20 离线工具 | da 链缺失（alpha 通道 7×） | PsiA vs rxc | 本轮发现、修正 |
| 闭环缺失 | (I−Ψφ)⁻¹ 4–7% | b20（近似成立） | 维持 |
| D2 特异 | 45% 面级失配 + 最大阶梯噪声 + 最低 cos | 面级锚 | 新增，需专项 |

**下一轮建议（按信息价值）**：
1. **B 点门控全 J 导出**（复用模块 round-2 伴随处的 explicitJT 导出路径
   加开关）→ w*(B) 精确 vs w_state：一举消除 SLOT-1 与仪器地板，直接
   裁决 SLOT-5 的算子侧份额；
2. 或把离线重建修到 <1%（出口边界通量公式 34–40% 误差 + D_rel clamp
   语义与 OF7 fvMatrix::relax 逐字核对）；
3. D2 专项：更深收敛重跑一对 D2 探针（其 1.9% 阶梯噪声 + 45% 面级失配
   提示状态质量）；
4. 一行修复 SLOT-2（模块伴随前重捕 primalPressureMobility=rAtU(B)）后
   重测 2.13× 量具。

## 7. 工件清单

| 文件 | 说明 |
|---|---|
| `b21_e2_identity_ladder.py/.log/.json` | E2 A 基恒等式阶梯 + PsiA vs rxd 对照 |
| `b21_e2b_exportedJ_identity.py/.log/.json` | E2 精确导出 J 行恒等式 + 钉扎对照 |
| `b21_wstar_clean.py/.log/.json` + `b21_wstar_cache.npz` | 干净未钉扎 SuperLU 解（1769 s） |
| `b21_e3_anchor_analysis.py/.log/.json` | E3 面级锚 + B 基重建 + mobility 量化 |
| `Log.verify_b021_anchor.txt` | 重跑日志（已归档） |
| `b21_anchor_fd_scan.tsv` | 重跑 FD scan（与存档 FD 侧逐位一致） |
| `/home/ys/dsH/b21_anchor/`（运行目录，只读保留） | 重跑算例与新导出 |
| `src/validateStageB2GradientAmplitude.H`（+73 行） | 唯一代码新增：门控导出 |

## 8. 纪律与不变量

- 零源码语义改动；唯一新增 = 开关门控导出（git diff 核实仅 +73 行导出）；
- 重跑走标准流程（绝对路径 build/bin/MTO_HF、算例 cp -a、日志归档）；
- 禁改清单（算子/装配/源/链/MMA/阈值/锁定开关）未触碰；无拟合因子；
- 工作树中既有的他人未跟踪文件未动；git 选择性提交，未 push。
