# BFINAL-029 — R1 预注册（Phase 1 + Phase 2 判据，数据之前）

> 任务书：`/home/ys/dsH/TO-ANISOTROPIC/NEXT_TASK_BFINAL029.md`。
> 本文档在**任何 Phase 1 数据计算之前**写就。判据原样引用任务书。
> Phase 3 修复授权：仅限定罪槽位，且只在离线定罪后。零拟合因子。

## 1. 目标与总判据（任务书 §2/§4 原样）

- Phase 1：c 段 b_TC 逐槽离线重建，**全 39.7% 失配必须分解到具体槽位**，逐槽判定（离线误差 / 生产误差 / 口径差异）。
- Phase 2：对嫌疑槽（主 H7s2、次 T4，加 Phase 1 发散槽）跑变体（清零/翻号）。
  定罪判据（任务书原文，预注册为绑定判据）：
  > **某个变体使 D1 行进入 |ratio−1| ≤ ~15% 且不破坏 D3 行（D3 当前 1.057）**。
  另测双变体交叉（如 H7s2+T4）。
  **若无一变体通过 → 诚实报告"候选排除"，停止，无 Phase 3。**
- Phase 3（仅定罪后）：修复定罪槽位；诊断同源槽位四路一致性（生产/诊断 J/诊断 J^T/CSR，B24 约定）；编译；验收序列。

## 2. Phase 1 槽位表（17+2）与逐槽可切换设计

每个槽位是 route 的**独立可切换加性分量**。开关 `{T1,T2,T3,T4,T5,T6,T7,T8,H7s1,H7s1b,H7s2,H7s2b}` ∈ {0,1}，
g0 恒开，g1 公式变体 ∈ {TIn（生产）, Tmix（B26 对照）}。A1/A2 为汇布局（恒开）。
C、Gx 不在 route 中（来自 b26a M1 件，Phase 1 不重算）。

| # | 槽 | 生产出处 | 内容（c 段折叠） | 符号 | 默认 |
|---|---|---|---|---|---|
| 1 | g0 | AdjNS_HT.H L23-25 | 内部面 g_f = −coldmask·adjDown·jumpT | − | 开 |
| 2 | g1_TIn | computeObjective.H L152-172 | 出口 dJ/dphi = −(T_out−T_in)/(M·Tref)，T_in=600 | − | **开（生产）** |
| 3 | g1_Tmix | b26a L289-291 | 出口 dJ/dphi = −(T_out−Tmix)/(M·Tref) | − | 变体 |
| 4 | T1 | Production:359-362 | hA += mob·w·(Sf·g)（αRel·rAU==mob） | + | 开 |
| 5 | T3 | Production:363-369 | U += (1−αRel)·w·(Sf·g) | + | 开 |
| 6 | T4 | Production:370-378 | P(own)+=kf·g；P(nei)−=kf·g | +/− | 开 |
| 7 | H7s1 | Production:379-385 | h7G += mob·w·g·Sf | + | 开 |
| 8 | T5 | Production:403-404 | 边界 hA += mob·(Sf_b·g_b)（assignable-U） | + | 开 |
| 9 | T6 | Production:405-409 | 边界 U += (1−αRel)·(Sf_b·g_b) | + | 开 |
| 10 | H7s1b | Production:410-412 | 边界 h7G += mob·g_b·Sf_b | + | 开 |
| 11 | T7 | Production:413-420 | 边界 P += mob·δ_b·\|Sf_b\|·g_b（pressureFixed=outlet） | + | 开 |
| 12 | T2 | Production:428-441 | δH^T：U(nei)−=upper/V_o·hA_o；U(own)−=lower/V_n·hA_n | − | 开 |
| 13 | T8 | Production:442-465 | 边界 H 对角 (−ic+cav)·Vinv·hA(c)；**b26a 跳过，本轮实测**。预期≡0：noSlip/fixedValue-U 的 (−ic+cav)=0，outlet zeroGradient-U 的 ic=0 | + | 开 |
| 14 | H7s2 | Production:469-481 | q=Sf&(h7G_o/V_o−h7G_n/V_n)；P(own)−=w·q；P(nei)−=(1−w)·q | − | 开 |
| 15 | H7s2b | Production:482-497 | P(c)−=(h7G&Sf_b)/V_c（仅 zeroGradient-p：inlet/hotInlet/hotOutlet/solidEndWalls/bottomWall/topWall/sideWalls） | − | 开 |
| 16 | A1 | Production:499-508 | 汇总布局 | | 开 |
| 17 | A2 | stageB18RhsExport L870-898 | bexp 排列（b[0:NV:3]=m[:,0]…b[NV:]=m[:,3]） | | 开 |
| 18 | C | sensitivity.H L396-399 | thermalDiffusionDerivativeDTCell（不乘 V） | | 外挂 |
| 19 | Gx | — | g^T·φ_x（真缺失，超范围） | | 外挂 |

### 2.1 Phase 1 判定协议（预注册）

**关键定义修正（写于数据前）**：b26a 的 0.397 是 [g1_Tmix, T8 关] 的 route；b29 默认是生产镜像
[g1_TIn, T8 开]。二者不是同一 route，0.397 是"旧 g 层口径 + 缺 T8"的失配，不是生产镜像的失配。

1. **仪器自洽门（首要）**：b29 重构 route 在 **B26 同设置（g1_Tmix，T8 关，其余全开）**下必须
   **精确复现 B26 的 0.397**：`relL2 ≈ 0.3970109364`（cos≈0.9582656）。这是对 b29 逐槽重构与 b26a
   route_V1 逐位一致的验证。若不成立 → 重构有误，BLOCKED 汇报（不强行分解）。
2. **加性门**：b = 各输入槽（T1,T3,T4,H7s1,T5,T6,H7s1b,T7）单独置开（仅该输入槽开、下游
   T2/T8/H7s2/H7s2b 恒处理其输入）之贡献之和，须在机器精度内等于全开 route：
   `Σ_j b(input_j on) == b(full)`（验证折叠管线的线性与实现正确性）。
3. **逐槽敏感性（leave-one-out）**：对全部 12 个可切换槽（含下游 T2/T8/H7s2/H7s2b）：
   `ΔrelL2(j) = relL2(full minus slot j) − relL2(full)`，生产镜像设置 [g1_TIn, T8 开] 为 full 基线。
   Δ 大（>0.01）→ 该槽是剩余失配的主要成分；Δ~0 → 该槽与失配无关。
4. **g1 口径贡献**：同掩码下 g1_Tmix vs g1_TIn 的 `ΔrelL2(g1) = relL2(Tmix) − relL2(TIn)` =
   g 层出口公式对 0.397 的贡献（B26 用了 Tmix，生产用 TIn；若 Δ 大则该贡献为 **caliber** 而非生产错误）。
5. **T8 实测**：T8 开关置 0 vs 置 1 的 ΔrelL2 必须实测；同时报告 ‖T8 贡献‖/‖b‖（理论：
   动量算子为标量各向同性 → (−ic+cav)≡0 → T8≡0）。若 >0.005 → b26a 的"恒零"注释作废并升级为嫌疑槽。
6. **逐槽判定标签**（每槽一个）：
   - `offline-error`：槽在 b26a 中实现有误（对照生产源码 L#）→ 修正后 route 失配下降；
   - `production-error`：槽实现正确但导出 b_TC 与生产源码不符 → Phase 3 嫌疑（需 w 求解确认）；
   - `caliber`：口径差异（如 g1 Tmix vs TIn、T8 恒零 no-op）→ 记录，不进入 Phase 2；
   - `innocent`：贡献 ~0 且实现正确。
7. **全表**（自洽门、加性门、逐槽 ΔrelL2、g1 口径、T8 实测、标签）写入 `b29_slot_table.json`。

### 2.2 Phase 1 通过标准

叠加门成立（全开 relL2 与 B26 0.397 一致到 ≥3 位），且逐槽贡献表给出各槽 ΔrelL2 与标签。
Phase 1 不判定生产有罪/无罪，只锁定嫌疑槽供 Phase 2 变体测试。

## 3. 预期数字（来自 B26/B27，非本次计算）

- M1（B26a，用 bexp 收缩，未受 Tmix 缺陷污染）：
  D1 lhs=−0.0488402 / rhs=+0.0744368 / ratio=−0.65613（sign_lhs=−, sign_rhs=+）；
  D2 lhs=−0.0064410 / rhs=+0.0006388 / ratio=−10.0836；
  D3 lhs=+0.0115666 / rhs=+0.0109403 / ratio=+1.05725（source-clean）。
  三段判定（B26）：a/b INNOCENT；c D1 CONVICTED（H7s2 主、T4 次）。
- FD 判据行（stage_b2_fd_scan.tsv，h=1e-3，正式判据）：
  D1 +0.04252469623637622（signJ=0）；D2 −0.003528972517971574（signJ=1）；D3 −0.0004449511895182612（signJ=1）。
  FD 主参考（h=1e-5/1e-4）：D1 +0.04244654087520727；D2 −0.003699503574594587；D3 −0.0004211702528400529。
- 目标变体收敛目标：D1 ratio→+1（|ratio−1|≤0.15），D3 保持 1.057（±小偏差）。
- PRODH7SHARE=0.701947825646（H7 相关槽占 c 段折叠的主导份额）。
- 案例参数：N=33600，NV=100800，NUNK=134400，ALPHA_REL=0.4，NU=5.19009e-05，
  Mflow=4.07337e-03，Tref=600.0，TIn=600.0，Tmix=687.840。

## 4. Phase 2 变体×身份数字表（预注册，判据绑定）

变体空间：对 Phase 1 锁定的嫌疑槽 S ∈ {H7s2（主）, T4（次）, Phase 1 发散槽}：

| 变体 | 定义 |
|---|---|
| V0 | 基线：现有折叠（Phase 1 全开 route） |
| Z_S | 槽 S 置 0（从折叠中移除） |
| F_S | 槽 S 翻号（±1 交换） |
| X | 双变体交叉：H7s2 与 T4 同时 Z 或同时 F（两种都测） |

每个变体计算 `b_TC_variant^T w_true(H7)` → 3 行（D1/D2/D3）数字。
身份矩阵（预注册字段）：`variant × identity = (dir, lhs, rhs, ratio, |ratio−1|, D3_ratio, verdict)`。

### 4.1 不变量（所有变体必须满足）

- w_true 与 rhs 来自 B26（**不重算 FD**）；
- C 与 Gx 固定（外挂件）；Gx=−0.0015704/−0.0040432…（D1/D2 已含）…见 b26a_m1_results.json；
- lhs 收缩用 `b_TC_variant^T w_true`，张量对齐（Ux Uy Uz P）与 B26 逐字节相同。

### 4.2 定罪判据（任务书原文，绑定）

对每个变体 v，检查 3 行：
1. **D1 行**：ratio_v,D1 ∈ [0.85, 1.15]（即 |ratio−1| ≤ 0.15）。B26 基线 −0.656 → 需变体把 D1 拉进带内。
2. **D3 行**：ratio_v,D3 ∈ [0.90, 1.20]（D3 当前 1.057；**不破坏**意味着仍在带内）。
3. D2 行记录但不作为定罪条件（B26 D2 自身被判定与 C/Gx 精度相关，见 M1 json reldev 11.08）。
定罪当且仅当存在变体 v 满足 1 与 2。

**明确否定规则**：若无一变体通过 1&2 → 诚实报告"候选排除（H7s2/T4/Phase1 发散槽均无法使 M1 一致）"，写结论，STOP，不进入 Phase 3，最终报告照常提交。

### 4.3 Phase 2 输出

`b29_variant_matrix.json`：逐变体 3 行数字 + 判定；NOTEBOOK 写 Phase 2 结论。
若定罪：Phase 3 计划（四路一致性修改 + 编译 + 验收）写入 NOTEBOOK 后再动源码。

## 5. Phase 3 计划（占位，仅定罪后激活）

- 修改 `src/solveDiscreteFlowAdjointProduction.H` 的定罪槽；若 `src/solveDiscreteFlowAdjoint.H` 同源槽一致（已核 L346-508 镜像逐槽同号），实现四路一致性（生产/诊断 J/诊断 J^T/CSR）。
- 编译命令（B28 惯例，禁止 set -e/set -u/tee/FD_SIGFPE 默认）：
  `export PATH=...; source /opt/openfoam7/etc/bashrc; export FOAM_USER_APPBIN=.../build/bin; export PATH=.../ccache-shim:$PATH; unset FOAM_SIGFPE; cd src; wclean >/dev/null 2>&1; wmake > ../wmake_current.log 2>&1`
- 验收：① b25_qgate FD 门全重跑（D1 修复：符号转正+振幅入带；D3 仍过冲=装配层 R2 事务；gDP ~2.05× 已知态，记录不判）；② b24_gauge 回归：FD 侧三列逐位一致，ADJ_J/gDP 列允许变化（记录逐列变化）；③ 双标签伴随收敛 ≤1e-12 + GRADSTABLE，等价门 1e-16；④ 静态测试全过。

## 6. 纪律

- 证据只追加（`evidence/agent-group/BFINAL-029/cycle-1/`，b29_ 前缀）；旧证据不动。
- 零拟合因子；不删/不禁用失败测试；不改验收阈值；不设 stageB4JacobianProbe=true；不建 FROZEN_GRADIENT_UNLOCK.txt。
- 结束：EXECUTOR_SUMMARY.md + 空 EXECUTOR_DONE + 一个选择性 git commit（841df95 风格），**永不 push**。
- 单步卡 >15 min 上报；歧义在安全点 BLOCKED。
