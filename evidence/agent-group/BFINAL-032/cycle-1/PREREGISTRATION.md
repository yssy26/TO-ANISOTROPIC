# BFINAL-032 预注册（先于数据）

> 写就时间：2026-08-22（编译/运行任何数据产出之前）。本轮动机与依据：
> B31 协调人会话（NOTEBOOK §7，commit 4999855）——`/1/` 场目录疑似混合不同伴随轮次，
> 使一切离线场级分解被跨轮伪影污染（`gsenshPD` 76.5% 范数残差之谜）。唯一可靠出路：
> **一次带开关门控导出的重跑**，让生产自己写出全部分解场（同轮次、同点）。
> 本文件锁定恒等式、预期与判定树（审阅判定树见 NEXT_TASK_BFINAL032.md）。

## 0. 基线事实（运行前已核）

- b25 运行恰好 **2 次** `sensitivity.H` 调用（Log.verify_b25_qgate.txt 实证）：
  - **轮 1**：主循环 MTO_HF.C L110（"MTO_HF: Starting/Finished sensitivity.H" 包裹，
    log 行 4784-4801），momentum L2 = 0.00351314475249。
  - **轮 2**：`validateStageB2GradientAmplitude.H` L460 `runFullB2Chain` lambda 内
    （裸 "sensitivity analysis"，log 行 10408），momentum L2 = 0.00226699400158。
  - 两轮 thermalC L2 = 0.00111307150907 相同。
- 混合机制（已锁定）：`writeOptimizationState.H`（MTO_HF.C L107，先于 sensitivity.H）
  `runTime.write()` 写 `/1/Uc /1/Ub /1/pc /1/pb /1/Tb`（AUTO_WRITE）= **轮 1 值**；
  轮 2 的 sensitivity.H L530-563 写块（finalSolvedIteration 真）覆写 on-disk 灵敏度族
  `gsensh*/gsens*/fsensh*` = **轮 2 值**。B31 §7 离线观测与此完全一致。
- 环境：N=33600，NV=100800，NUNK=134400，ALPHAREL=0.4，NU=5.19009e-05。
- 本轮开关 `stageB31WriteDecomposition`（optProperties，默认 false）。b32_decomp 克隆仅加
  `stageB31WriteDecomposition true;`，其余逐字同构 b25。

## 1. 硬闸（阶段 2，先于一切分解结论）

**FD 侧三列逐位一致**：b32_decomp 产出的 `stageB2/stage_b2_summary.tsv`
（projJ_D1/D2/D3、projDP_D1/D2/D3、projV_D1/D2/D3）与
`stageB2/stage_b2_fd_scan.tsv`（FD_J/FD_gDP/FD_gV 三列）必须与 b25 对应文件
**逐位一致**（`cmp` 字节级 或 sha256 一致）。任何非逐位 ⇒ 导出污染计算路径 ⇒
立即停止报告（实现污染，回炉）。

判定数据：b32_decomp 与 b25_qgate 两个 tsv 的 sha256 对照表 + `cmp` 结果。

## 2. 阶段 3 恒等式与预期

数据在手的同轮次分解（全部用 b32 写出的同轮场，不再碰 `/1/` 跨轮目录）：

- **门 1（同轮闭合门，应逐位）**：`gsenshPD_r? == gsenshPressureDropMomentum +
  gsenshPressureDropPressureRow`（生产两段直接求和）。relL2 ≤ 1e-12 预期。
- **门 2（T 环镜像）**：生产 `rxPressureRowT(pc)` vs 协调人离线 T 环（b31b_m2.py
  `T_row(lam)`）。**预期 relL2 < 2%**。若失败：差异即 T 输入的生产实现细节
  （协调人重建用其自有 dHbyA/drAU 通道，2% 为校准先例）。
- **门 3（跨轮伪影消解）**：用同轮 `Uc/pc` 重算 mom/prow 段投影，与 B31 §7
  跨轮数字对照。**预期 76.5% 残差大幅消失，同轮闭合残差 < 5%**。消失后剩余部分
  即为**真病灶段**（或证无）。判据（审阅树）：残差消解且 T 环镜像通过 →
  m2 对账用同轮数据完成（回 B31 收尾）；残差残留 → 真病灶现形，进入修复轮规格。
- **门 4（HbyA 校准）**：生产 `rxDHbyA/rxrAU` vs 协调人重建（HbyA 2% 误差校准）。

## 3. 预期数值（预注册）

| 量 | 预期 |
|---|---|
| FD 三列逐位 | sha256(b32 tsv) == sha256(b25 tsv) |
| 门 1 同轮闭合 relL2 | ≤ 1e-12（逐位） |
| 门 2 T 环镜像 relL2 | < 2% |
| 门 3 跨轮伪影消解 | 同轮闭合残差 < 5% |
| ADJ_J / ADJ_gDP 复现 | b25 tsv 值逐位（-0.0358350400485765 等，见 NOTEBOOK §2） |

## 4. 判定树（审阅视角，见任务书）

1. FD 不逐位 → 实现污染，回炉（硬闸，不可绕过）。
2. FD 逐位 且 门 1 闭合 且 门 2 通过 → 门 3 数字即答案：
   - 残差 < 5%：跨轮伪影消解确认，无真病灶（或病灶远小于 76.5% 表象）；
   - 残差仍大：真病灶现形，进入修复轮规格。
3. 门 2 失败：T 输入差异为生产实现细节，门 3 数字仍有效（T 输入不影响
   mom/prow 段的同轮闭合，因为两者都用生产自己的 T）。

## 5. 纪律声明

- 本文件先于一切分解数据/编译完成。
- 所有产物 `b32_` 前缀入本目录。
- 不修改任何计算路径：仅新增 `stageB31WriteDecomposition`（默认 false）门控的
  `.write()` 与日志打印；默认关闭时 b32 == b25 计算路径逐位。
