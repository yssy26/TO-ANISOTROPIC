# BFINAL-022 cycle-1 — FINAL REPORT（量具修复与复测轮）

任务：修复 BFINAL-021 挖出的三个量具缺陷（SLOT-1/2/3），建立 B 点一致量具，
复测全部头条数字（gDP 因子、J 符号表、恒等式、λ 一致性），并做 D2 专项。
分支 `agent/dsH-stage-b-validation`，起点 `8a8e10b`，执行 2026-08-21。
预期管理：本轮结果决定项目走向——**结果：三个缺陷全部修复后，头条数字
基本不动（gDP 因子 2.13×→2.15×，J 符号 2/3 不变）⇒ O(1) 幅度失配与
D1 符号翻转是同点、同基的真失配（SLOT-5 坐实为唯一残余解释）。**

---

## 0. 结果速览

| 项 | 结果 |
|---|---|
| W1 SLOT-2 修复 | `validateStageB2GradientAmplitude.H` 量具层 +59 行：full-SST 再基线（phase-1）+ 冻结后、round-2 伴随前，用与冻结求解器完全相同的装配/relax/constrain 序列重建动量矩阵并刷新 `primalPressureMobility = rAtU(B)`；真实 A→B kf 基变化 **relL2 0.1742** |
| W4.1 不变性 | 静态 **13/13**（新增 `test_bfinal022_b2_module_refreshes_kf_basis_at_state_b`）；WMAKE_EXIT=0；FD 侧 13 行、基线四值、ADJ_gV、round-1 全部（GRADPROXY/等价性门四值/残差）**逐位不变**；ADJ_J/ADJ_gDP 按预期变化 |
| W4.2 头条 | gDP 因子 ADJ/FD = **2.146/2.159/2.181**（修前 2.130/2.113/2.168；微幅变差）；J 符号 **D1 仍反（2/3）**；FROZEN_GRADIENT_STATUS=FAIL |
| W2 SLOT-3 修正 | 修正版 `b22_b20augmented_fixed.py`：closure_rel 4-7% → **9.7/28.0/8.1%**；rxd_P_c「14%/cos0.93」→ **0.92/0.12/0.99 幅度、cos≈−0.35（反号）** |
| W4.3 B 点量具 | 系数层仪器误差 **12.9% relL2**（新实锤，见 §3）；面锚@B D1 8.7%/D2 **44.9%**/D3 9.0%；恒等式@B O(1)（3.5/2.8/4.1 混合基）；rhs@B/ mobility@B 已导出 |
| W4.4 伴随 | 双标签 ≤9.8e-10（算例容差 1e-9；历史 b18/b19 的 1e-12 配置见 §5 注）；GRADSTABLE 4/4（各 1 次 continue 后稳定）；λ 身份：round-2 代理移动 0.3-0.5%，round-1↔round-2 差 7.6% = 纯 SLOT-1 点差 |
| W5 D2 | 失配能量 **98% 在 top-10% 面、97% 在 x-法向面**（真值 92%）；x 通道 relL2≈0.46 vs y 0.21/z 0.30；最差 x∈[0.1,0.2] |
| 代码 | +126 行（+76 量具模块、+50 静态测试）；生产路径零触碰；一次重编译 + 一次求解器运行 |

---

## 1. W1 — SLOT-2 修复（字段清单与核查）

**改动**：`src/validateStageB2GradientAmplitude.H`（仅量具层）：
1. phase-1 `solveFullSST` + `updateFrozenTurbulenceFields.H` 之后、
   `baselineState` 快照与 `runFullB2Chain()`（round-2 伴随）之前，插入
   kf 基刷新块：按 NS.H/solveFrozenPrimalConverged 的同一装配序列
   （`fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U)
   == −fvc::grad(p) + fvOptions(U)`，ransFlowModel 时减显式偏应力项，
   `relax()` + `fvOptions.constrain()`）重建动量矩阵，`rAU=1/A`，
   `simple.consistent()` 时 `rAtU=1/(1/rAU−H1())`，
   **`primalPressureMobility = rAtU().primitiveField()`**，并打印与陈旧
   A 基的 relL2（实测 **0.174181**）。
2. stageB15StateExport 门控块内新增 `wstate_baseline_primalPressureMobility.mtx`
   导出（W3：B 态精确 kf 基，供离线量具）。

**round-2 消费字段逐一核查**（全 B 态 ✓）：
- U/p/phi/alpha/k/omega/nut：phase-1 full-SST 终态 ✓
- nutFrozen/nuEffFrozen/alphaTurbulent：`updateFrozenTurbulenceFields.H`
  （模块内 line 114）✓
- primalPressureMobility：**修复对象**（此前唯一陈旧项，NS.H:129 之外无
  写点——静态测试第 13 项断言此唯一性）
- pressureConstraintDerivative：纯几何（createFrozenHotRegionFields.H 装配
  入口/出口面元），与流动状态无关 ✓
- thermalObjectiveDerivative / thermalObjectiveFluxDerivative：
  runFullB2Chain 内 computeObjective.H 在 B 态重算 ✓
- DT/DTEffective/x/xp/xh/eta5：设计链场，A/B 同值（flow solve 不改设计）✓

**等价性**：FD 探针（solveFrozenPrimalConverged）不读该场 ⇒ FD 侧逐位不变
（§4.1 实证）；模块出口 restoreStateInto(entryState) 不恢复该场，但 endTime=1
单次外循环 + NS.H 下轮重捕，且 validateSSTDirection/Gate6 不消费该场（grep
核实）⇒ 主循环无泄漏影响。

## 2. W4.1/W4.2 — 不变性门与头条复测

### 2.1 逐位不变性（全部通过）

| 项 | 证据 |
|---|---|
| FD_J/FD_gDP/FD_gV/Valid/nIter±（13 行） | b22_fd_scan vs b21_anchor 逐字段字符串相同 |
| wstate_* 导出 122 个文件 | 逐字节相同 |
| 基线 J/PD/gDP/gV | −1.16745040316938 / 48353.78151317623 / 0.934151260527049 / −7.887278274942844e-09 逐位 |
| round-1（A 点）GRADPROXY-FINAL | TC +1198.85085218、PD −205194.141124 逐位 |
| round-1 等价性门四值（两标签） | 1.13142348768e-16 / 8.25782114104e-17 / 1.42213008543e-16 / 1.34302013134e-16 逐位 |
| round-1 伴随迭代/残差 | TC 694/9.9016e-10、PD 902/9.9282e-10 逐位 |
| ADJ_gV（体积链不经流伴随） | 逐位 |
| round-2 等价性门四值 | 1.138e-16/1.074e-16/1.531e-16/1.334e-16（~1e-16，随新 kf 微移，符合预期） |

### 2.2 头条对照表（h=1e-3；完整表见 b22_fd_scan.tsv vs b21_anchor）

| 量 | D1 | D2 | D3 |
|---|---|---|---|
| FD_gDP（逐位不变） | −2.537459771851047 | 0.5234094521578969 | 6.477215547073478 |
| ADJ_gDP 前→后 | −5.4050587339 → −5.4440967986 | 1.1059269805 → 1.1300780777 | 14.0447142088 → 14.1276593948 |
| ADJ_gDP 变化 | −0.72% | +2.18% | +0.59% |
| **gDP 因子 ADJ/FD 前** | **2.130** | **2.113** | **2.168** |
| **gDP 因子 ADJ/FD 后** | **2.146** | **2.159** | **2.181** |
| FD_J（逐位不变） | +0.0428559144 | −0.0034689264 | −0.0002826582 |
| ADJ_J 前→后 | −0.0347544 → −0.0349244 | −0.0238141 → −0.0233185 | −0.0136414 → −0.0139741 |
| J 符号 | D1 **仍反** | 同号 | 同号 |

**裁决：修复 SLOT-1（比较双方同在 B 态）+ SLOT-2（kf 同基）后，gDP 幅度
失配与 D1 J 符号翻转原样复现且微幅变差。** 任务预设的「若恶化如实报告并
停止定位」条款触发：已如实报告，未调参。

## 3. W4.3 — B 点一致量具与离线仪器解剖

选项取舍：精确 J(B) 导出需 B4 诊断路径，被 validateMmaUnlockGate 禁止与
B2 同跑 ⇒ 走「门控导出算子数据」路线的可行子集：mobility@B（新导出）、
rhs@B（既有 stageB18RhsExport 开关，零代码）、FD 态@B（既有导出）；
rxd@B 用同基离线重建。

### 3.1 系数层仪器误差（本轮新实锤）

精确 mobility(B) vs 离线 D_rel 重建（B15–B21 全部离线量具的 rAU/kf 通道）：

| 指标 | 值 |
|---|---|
| relL2 | **12.93%** |
| maxRel / 胞>1% / 胞>0.1% | 39.6% / 21.1% / 36.6% |
| 结构 | 固体区（α>1e6）精确（2.6e-5）；偏差集中于 nut>1e-6 胞（均值 11.8%、98%>1%）与 clamp 绑定胞（sumOff>|D|，均值 13.0%） |

OF7 语义排查（b22_mobility_forensics.py）：
- S1 边界对角：relax 尾部 `-cmptMin(ic)` 与 `D()=diag+addCmptAvBoundaryDiag`
  对均匀分量 ic 恰好抵消 ⇒ 旧公式（clamp 前加 |ic|、无移除）**正确**；
- S2 bounded Gauss upwind 的 `-fvm::Sp(surfaceIntegrate(phi))` 带 V 因子
  （fvmSup.C）vs 脚本无 V——影响 2.8e-11（V=1.2e-10，可忽略）；
- S3 laplacian γ 面插值 vs 单侧——改善 0.9e-2 relL2（12.93→12.52%）。
- **残余 12.5% 未定位**（疑在 nut 梯度带的面系数更细语义）。
- 推论：B21「D1/D3 达仪器地板 3.35%」的评估需弱化——该地板是通量层读数，
  系数层带 13% 系统差；离线 <10% 的绝对读数不可分辨。

### 3.2 SLOT-2 真幅度（真值对真值）

mobility(B 真值)/mobility(A 重建)：中位 1.0000、**9.3% 胞超 ±10%**
（min 0.47/max 1.36）——B21 用「重建对重建」估计的 3.3% 低估了三倍。

### 3.3 面级锚与恒等式@B（b22_bpoint_gauge.py）

| 量 | D1 | D2 | D3 |
|---|---|---|---|
| 面锚@B 重建基 relL2 | 8.68% | **44.9%** | 8.99% |
| 面锚@B exact-kf 混合基 | 11.99% | 39.8% | 12.43% |
| 恒等式@B 混合基 \|rP\|/\|rxdP\| | 3.52 | 2.82 | 4.10 |
| （参考）A 点精确 J 恒等式 | 1.66 | 2.14 | 1.51 |

D2 的 45% 在两基下都复现；D1/D3 的 B 基 8.7-9.0% 在 13% 系数不确定度内
不可分辨。混合基（精确 kf + 重建 D_rel）内部不一致，仅作括号读数。
w′@B 未解：无精确 J(B)（§3 开头），离线 J(B) 13% 系数误差下无信息量——
由 FD 门的精确同点比较替代（§2.2）。

### 3.4 rhs@B 导出（W3，零代码开关）

PD：p 块与 A 逐位相同（几何量）、U 块两者皆零；TC：U 块差 28%、p 块差
164%（B 态 rhs 与 A 实质不同，可供未来 J(B) 精确导出后使用）。

## 4. W4.4 — 伴随双标签、GRADSTABLE、λ 身份

| 项 | round-1（A，主循环） | round-2（B，模块） |
|---|---|---|
| TC 伴随 | 694 it / 9.9016e-10（逐位同 b21） | 691 it / 9.7328e-10（b21：687/9.5518e-10） |
| PD 伴随 | 902 it / 9.9282e-10（逐位） | 675 it / 9.6057e-10（b21：679/9.8016e-10） |
| GRADSTABLE | 4/4 | 4/4（TC 一次 relChange 3.9e-6 越限后继续迭代至稳定，与 b21 同型） |
| λ 代理 | TC +1198.85 / PD −205194.14 | TC −322.81（−0.31%）/ PD −220817.60（−0.50%） |
| 模块指纹 | — | TC −479.346（+0.31%）/ PD −327895.71（−0.50%） |

注：本算例（b13_probe3 谱系）容差 1e-9，双标签收敛至 ≤9.9e-10；b18/b19
历史「≤1e-12」来自其 1e-12 容差算例配置。**λ 身份一致性的 W1 直接检验**：
修复使 round-2 代理仅移动 0.3-0.5% ⇒ round-1↔round-2 的 7.6% 持久差是
SLOT-1（两点差），非 kf 陈旧——与 B18 定位一致并加严。

## 5. W2 — B20 结论作废/复算清单（b22_b20augmented_fixed）

| B20 结论 | 状态 | 修正后 |
|---|---|---|
| (I−Ψφ)⁻¹ 修正 relL2 4.2–7.1% | **作废→复算** | **9.7% / 28.0% / 8.1%**（D1/D2/D3） |
| rxd_P_c「14% 幅度、cos 0.93」 | **作废→复算** | 幅度 **0.92/0.12/0.99 ×\|rxd_P\|**、cos **−0.35/−0.36/−0.29**（反号）——闭环修正是 O(1) 且反向的，非小而同向 |
| P 行恒等式（重建基，A 点） | 复算参考 | \|rP\|/\|rxdP\| 1.47/1.53/1.43，cos −0.71/−0.80/−0.73（与 B21 一致） |
| 七候选仲裁、rxd FD 验证（5.6%）、w_true 测试 | **不受影响** | 用求解器内导出 rxc（z 链正确） |
| B21「closure 4–7% 近似保持」 | 修正 | 8–28% |

## 6. W5 — D2 专项定位

D2（h=1e-3，重建基）面级失配 totR−dphi_true 的结构：

- **98.0% 失配能量在 top-10% 面**——高度集中，非均匀噪声；
- **97.1% 在 x-法向（流向）面**（真值能量 92.3% 为 x-法向）⇒ 分类 relL2
  ≈ x 0.46 / y 0.21 / z 0.30——**流向通量切线通道是坏通道**；
- x-十分位：0.37/0.76/0.35/0.53/0.37/0.34/0.46/0.26/0.54/0.18——最差
  x∈[0.1,0.2]，最好出口段 [0.9,1.0]；
- 固侧（owner α>1e3）面 0.458 vs 流侧 0.218；corr(失配, z_D2) = 0.52。

与 D2 方向构造（sin(4πx)+0.3cos(2πx)，纯 x 变化）一致：D2 的设计摄动
本身集中在 x 向通量上，其切线失配是**集中的、方向选择性的通道缺陷**。

## 7. 缺陷地图（更新）与下一轮建议

| 缺陷 | 本轮状态 | 量化 |
|---|---|---|
| SLOT-1 A-B 跨点混合 | 头条比较中消除（双侧@B） | round-1↔2 代理差 7.6%（纯点差） |
| SLOT-2 stale kf | **已修复**（W1） | 真幅度 kf relL2 17.4%；头条投影 ≤2.4% |
| SLOT-3 工具 da 链 | **已修复**（W2） | closure 8-28%、corrP 反号 |
| SLOT-4 闭环缺失 | 幅度修订 | 4-7% → 8-28%（D2 最差） |
| **SLOT-5 真失配** | **唯一残余解释（坐实）** | 同点同基下 gDP 2.15×、J 符号 2/3、恒等式 O(1) 全部复现 |
| 新增：离线系数层仪器差 | 实锤、未定位 | 12.9% relL2（nut 活跃+clamp 胞） |

**下一轮建议（按信息价值）**：
1. **SLOT-5 专项**：以 D2 为载体定位流向通量切线通道（W5 已给出空间/方向
   指纹：top-10% 面、x-法向、x∈[0.1,0.2]）——对照该区间算子 P 行与
   Rx-B 探针的逐面差；
2. 离线系数层修复（fvMatrix 装配逐语义核对至 <1%）——解锁恒等式/锚的
   离线绝对读数；
3. FD 门本身作为主量具常态化（本轮证明其精确、可复现、免疫离线仪器差）。

## 8. 工件清单

| 文件 | 说明 |
|---|---|
| `b22_bpoint_gauge.py/.log/.json` | B 点一致量具（mobility 对照/保真度/面锚/恒等式/D2 定位） |
| `b22_mobility_forensics.py/.log` | 系数层 12.9% 误差的结构解剖与三个 OF7 语义嫌疑排查 |
| `b22_b20augmented_fixed.py/.log/.json` | W2 修正版 augmented 工具（da=z 链） |
| `Log.verify_b022_gauge.txt` | 求解器运行日志（已归档） |
| `b22_fd_scan.tsv` / `b22_summary.tsv` / `b22_adjoint_identity.tsv` | FD 门全套输出 |
| `wmake_b22.log`（仓库根） | 全量重编译日志 WMAKE_EXIT=0 |
| `/home/ys/dsH/b22_gauge/`（运行目录，只读保留） | 含 stageB2 全部导出（含 mobility@B、b18rhs@B） |

## 9. 纪律与不变量

- 源码：仅 `validateStageB2GradientAmplitude.H`（+76）与
  `src/tests/test_stage_b_safety_gates.py`（+50）；git diff 核实；
  生产算子/装配/源/链/MMA/阈值/锁定开关未触碰；无拟合因子；
- `frozenGradientValidated=false`、`mmaUpdateEnabled=false` 维持；
- 重跑标准流程（cp -a b13_probe3 → b22_gauge，绝对路径 build/bin/MTO_HF，
  日志归档）；他人未跟踪文件未动；git 选择性提交，未 push。
