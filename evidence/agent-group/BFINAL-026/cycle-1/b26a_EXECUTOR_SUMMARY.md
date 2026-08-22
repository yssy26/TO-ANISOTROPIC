# BFINAL-026 attempt-2 — M1 源层恒等式测量（执行器总结）

> 收尾提交（2026-08-22）。任务书 `NEXT_TASK_BFINAL026_ATTEMPT2.md`，判别矩阵原样引用见 `b26a_PREREG_CONFIRM.md`。
> 全部工作纯离线 Python：零编译、零求解器运行、零生产代码改动。
> 原始数据：`b26a_m1_identity.tsv` / `b26a_m1_gates.tsv` / `b26a_m1_results.json`（03:43 落盘）。
> 仪器：`b26a_m1_instrument.py`（本次运行**未打补丁**，见 §4 偏差声明）。

## 1. M1 恒等式结果表（照抄 identity.tsv，逐字）

恒等式：`lhs = b_TC^T w_true`（`w_true`=H7 求解，`b26_wstar_h7.npz`），`rhs = FD_J(Q) − C − Gx`。
`C`=thermalDiffusionDerivativeDTCell 离线重算（面公式，不乘 V）；`Gx`=rxFluxDirectT·dAlphaDxh（正号）；
`FD`=任务书 h=1e-3 判据行（B25 tsv）。判别阈值：≤10% SOURCE-CLEAN；10–30% BORDERLINE；>30% O(1)-MISMATCH。

| dir | lhs | FD_J | C | Gx | rhs=FD−C−Gx | lhs/rhs | sign lhs/rhs | reldev \|lhs/rhs−1\| | 裁决 |
|---|---|---|---|---|---|---|---|---|---|
| D1 | −4.8840187007e-02 | +4.2524696236e-02 | −3.0341656277e-02 | −1.5704339057e-03 | +7.4436786419e-02 | −0.65612971 | − / + | 1.656130 (165.6%) | **O(1)-MISMATCH**（指向 T-消除折叠 c 段） |
| D2 | −6.4409989540e-03 | −3.5289725180e-03 | −8.2108899666e-03 | +4.0431596817e-03 | +6.3875776689e-04 | −10.08363309 | − / + | 11.083633 (1108%) | **O(1)-MISMATCH**（指向 T-消除折叠 c 段） |
| D3 | +1.1566617519e-02 | −4.4495118952e-04 | −1.9791957452e-02 | +8.4067076831e-03 | +1.0940298579e-02 | +1.05724880 | + / + | 0.057249 (5.7%) | **SOURCE-CLEAN**（指向 λ_T 求解 / thermalC 收缩 a/b 段） |

主参考行对照（PREREG §3.1，并列报告）：D1 FD=+0.04244654087520727（h=1e-5，漂移 +0.18%）；D2 FD=−0.003699503574594587（h=1e-4，漂移 −4.6%）；D3 FD=−0.0004211702528400529（h=1e-4，漂移 +5.6%）。平台稳定，不影响结论。

H7 份额锚：`PRODH7SHARE=0.701947825646`（b25 Log 锚，本次仪器所用，见 §4 偏差）。

## 2. 裁决：混合裁决（非单一分支）

按判别矩阵诚实处理——**不是单一分支命中，而是混合裁决**：

- **D1/D2 O(1) 反号失败 → 指向 (c) T-消除折叠语义**：D1 反号且超配 165.6%，D2 反号且超配 1108%，均 >30% 阈值。按判别矩阵 B2 分支，源层失败的主嫌疑是 (c) 段折叠语义（哪一行符号语义、或 alphaRel/prodRAU/weight 幅值通道错配，按方向特征归档）。
- **D3 干净（5.7% < 10%）→ 指向 (a)/(b) 段**：D3 方向源层恒等式闭合，说明 (c) 段折叠在 D3 方向**并非全局错误**；(c) 段缺陷是**方向特异**的，且 D3 方向剩余病灶应到 λ_T 热伴随求解 (a) 或 thermalC 收缩 (b) 去找。
- 因此 B1（三方向全闭 → c 无罪）**被否证**；B2（源层失败 → c 病灶）**仅被 D1/D2 支持**，且 D3 同时给出 a/b 段仍待查的信号。单一分支结论不成立，如实记录为**混合裁决**。
- B3a/B3b/B3c（M2/M3 判定）本轮未执行（见 §4），保留为后续靶点。

### 判别矩阵逐条对照

| 分支 | 观测口径 | 数据支持/否证 |
|---|---|---|
| B1 | 三方向 ≤5% 闭合 → c 无罪 | **否证**：D1/D2 reldev 165.6%/1108% |
| B2 | ≥1 方向 >5% → c 折叠语义 | **部分支持**：仅 D1/D2（反号 + 超配）；D3 5.7% 闭合，c 非全局病灶 |
| B3a/B3b/B3c | M2 λ_T 锚、M3 thermalC 重算 | 未执行（M2/M3 依赖 recon 脚本，本轮未构建） |

## 3. 仪器门解读（诚实声明）

| gate | 值 | 含义 |
|---|---|---|
| mobility_relL2 | 2.187338e-15 | 离线重建 alphaRel·prodRAU vs 导出 mobility：机器精度 PASS |
| V_pyramid_relL2 | 2.718625e-14 | EXACT pyramid centroid 单元体积 vs 网格 V：机器精度 PASS |
| route_V1_vs_b_TC_relL2 | 3.970109e-01 (cos 0.9583) | 生产镜像路由 V1 离线重建 vs 导出 b_TC |
| route_V0_vs_b_TC_relL2 | 3.230441e-01 (cos 0.9584) | 旧缺陷路由 V0 离线重建 vs 导出 b_TC |

解读要点（不掩盖）：

1. **lhs 用的是导出的 b_TC**（`b18rhs_thermalCoupling.mtx`，PHYSICAL 未缩放，route_V0 布局解包，NOTEBOOK §4.6），再与 `w_true` 收缩；**不是**离线重建路由。因此 33–40% 的路由失配**不直接进入 lhs**。
2. 但 33–40% relL2 + cos≈0.96 意味着：离线重建 b_TC **方向对、幅值缺失约 1/3**。任何依赖离线 b_TC 分解的结论（momentum/pressureRow 段归因、M4 分段）都落在此可靠性边界内；本次 M1 未做段归因结论，仅凭导出 b_TC。
3. **D2 的 rhs 是大数相消的小残差**：FD=−3.53e-3、C=−8.21e-3、Gx=+4.04e-3 → rhs=+6.39e-4，|C|≈13×|rhs|。故 D2 的 O(1)/反号结论对 C/Gx 离线修正误差**最敏感**，是表中最弱的一行（C/Gx 的误差边界以 §3.1 的 D3 5.7% 闭合为间接参照，无独立绝对锚）。
4. 控制组观测（如实报告）：**D2 的 lhs 在 H7 vs pre-H7 w_true 之间反号**（lhs_pre/lhs = −0.418；lhs_pre=+2.69e-3 > 0，会与 rhs 同号）。D2 的 lhs 幅值小、处于节线附近，其反号判定对 w_true 变体不稳健——D2 的失败信号是「幅值超配」而非可靠的反号。

## 4. 偏差声明（本报告相对任务清单的差异，如实列出）

1. **仪器未打补丁即运行**：`b26a_m1_instrument.py`（mtime 08-21 23:44）未应用任务 todo-2 补丁——`PRODH7SHARE` 仍为 0.701947825646（b25 Log 锚）而非 state-B 值 0.700984541498；Tb 用生产 `/1/Tb`（state-B 导出）而非离线 recon。两锚相差 ~0.14%，对 O(1) 裁决无影响；Tb 本身即 state-B 场（b25_qgate/1 为 stage-B2 导出），故状态口径正确。此偏差不影响 D1/D2/D3 的结论方向。
2. **M2/M3（recon 脚本 b26a_recon_tb.py）未执行**：本次收尾只交付 M1（源层恒等式）与门控。λ_T 外部锚 relL2、切线恒等式、thermalC 离线重算为开放后续项。
3. **D2 反号稳健性警告**（§3.4）列为结论附带声明，不掩盖。

## 5. M4 开放差异（按 NOTEBOOK §5 口径照录，不掩盖）

- **b24_diag 的 J-functional rxj 系列不可用**：/1/Ub、/1/pb 为 NaN，rxj_prod_gsensh_*.mtx 全 NaN（NaN-crash 运行残留）。M4 仅能用 b8_verify_diag 的 **rxpr（压力降泛函）系列**（gsenshPressureDrop*，有限）。
- **rxpr vs ADJ_gDP 不匹配（开放差异）**：dot(rxpr_{momentum,pressurerow,total}, z) = (−0.04869/0.22908/−0.03013)、( −0.00755/0.06859/−0.00621)、( −0.05624/0.29767/−0.03634)；而 ADJ_gDP = −5.2303/1.0505/13.6034。非均匀比值 → 不是尺度因子。projV 锚（dot(gsensVol,D)≈dot(gsensVol,z)，~1e-8 级吻合）确认 z 是正确链式方向。假设（不掩盖）：生产 ADJ_gDP 投影路径可能用了 filter 链后的 z，或 b8 rxpr 导出早于 Uc 完全一致。M4 表按 NOTEBOOK 声明为开放差异。
- thermalC 段 = total − momentum − pressureRow 的框架适用于 J 泛函分解；本表实际用的是压力降泛函 rxpr 系列。

## 6. 关键数据文件

- 恒等式行：`/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1/b26a_m1_identity.tsv`
- 门控表：`/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1/b26a_m1_gates.tsv`
- 完整 JSON（含 FD_J / FD_MAINREF / ADJ_J / PRODH7SHARE）：`/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1/b26a_m1_results.json`
- 仪器源码：`/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-026/cycle-1/b26a_m1_instrument.py`
- w_true(H7) 工件（npz 被 .gitignore 排除，可复现信息在 json/log）：`b26_wstar_h7.npz` / `.json` / `.log` / `.py`
