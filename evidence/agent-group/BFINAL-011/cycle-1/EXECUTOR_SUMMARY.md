# BFINAL-011 cycle-1 — EXECUTOR SUMMARY (2026-08-19)

**结论：PASS。pressureDrop GAMG 内层 NaN 的根因已定位并修复；两标签首次在生产迭代路径下同时收敛，完整程序 MTO_RC=0。**

```
thermalCoupling: FGMRES 750 iters  trueRelRes=9.4956818492e-10  (<=1e-9)
pressureDrop:    FGMRES 902 iters  trueRelRes=9.92818139401e-10 (<=1e-9)
GRADPROXY 平台:   TC +2217.48004 / PD -205194.141
PRODPRECGAMGCHECK 四组门值逐位不变 (1e-16 量级，两标签一致)
静态测试 7/7 OK；运行 120 s；MTO_RC=0（sensitivity.H 正常完成）
```

## 根因（一句话）

BFINAL-010 等价性门对非 const 矩阵调用 `upper()/lower()`，OF7 非 const
`lduMatrix::lower()` 会把对称存储的 lower 物化为独立数组（lduMatrix.C:168-182），
矩阵被判为非对称：PCG 被拒、GAMG 走非对称分支（`scaleCorrection_ =
matrix.symmetric()` → false，GAMGSolver.C:76），在该损坏路径上对
pressureDrop rhs 首次内层求解即发散（33600/33600 单元非有限）。
**修复 = 门内改用 `const lduMatrix&` 别名读取（const 版返回 upper 别名、
无副作用）**；门数值逐位不变，求解路径恢复对称 GAMG（scaleCorrection=true）。

## 实验矩阵与判别路径

| run | 状态 | 结果 |
|---|---|---|
| E0（修复前，默认） | 旧路径（物化非对称） | 完整复现 b010：门值逐位一致；TC 750/9.58e-10 收敛；PD 在 apply 1 被「非有限 psi」快检 FatalError（33600 单元非有限）|
| E1（修复前，`solveThermalCouplingFlowAdjoint false`） | 旧路径 | HT 完全不在日志中，PD 依旧 apply 1 NaN；PRODRHS PD sumRhsP=95.1198507991 不变 → 排除跨标签泄漏（H2）|
| E2（修复前，`discreteProdPressurePrecSolver PCG`） | 旧路径 | `Unknown asymmetric matrix solver PCG` → **发现矩阵被物化为非对称**（E2 日志已存档）|
| **E0'（修复后，纯默认）** | **对称存储恢复** | **TC 750/9.496e-10、PD 902/9.928e-10 双收敛，MTO_RC=0**；TC 轨迹与 BFINAL-009 b8 部分数据逐位一致（@80 3.04473267442e-4，@240 4.90659573665e-5）|

决策树在 E0' 解决：既不需要 PCG，也不需要 GAMG 参数排障
（scaleCorrection/nPreSweeps/directSolveCoarsest 覆盖项保留为稳健性保险）。
E1'/E2'（修复后变体）无需运行，留作可选检查。

## 对既有结论的修正

- BFINAL-010 FINAL_REPORT 中「pressureDrop NaN 属 GAMG 内层求解稳健性问题（BFINAL-011 范畴）」修正为：**触发者正是 BFINAL-010 门自身的存储副作用**，非 GAMG×rhs 固有失稳。BFINAL-010 的等价性结论（四组 1e-16）不受影响——物化不改变数值，只改变存储/求解路径。
- 计划中的 H1/H2 二分不完备：真因为 H3（矩阵存储被检查代码物化 → 求解器路径损坏）。E1 的「单标签仍 NaN」与 H3 同样自洽。

## 证据与诚实性说明

- `Log.verify_b011_E0p.txt`（胜局运行，sha256 已记录）、两份 checkpoint tsv、`optProperties.b11_diag.txt`、`Log.verify_b011_E2.txt`（副作用发现工件）已归档于本目录。
- **E0/E1（修复前）日志因运行脚本初版 `rm -f Log*` 被误删**（脚本已改为只删当前运行的日志）；其关键行已在上方表格与本目录 FINAL_REPORT.md 中原样引用。如需逐字节复核 E0/E1，可用 git 历史中的任一旧二进制路径重放（配置已记录）。

## 代码改动（默认全部=基线）

1. `createFields.H` + `MTO_HF.C`：`solveThermalCouplingFlowAdjoint` /
   `solvePressureDropFlowAdjoint`（默认 true）标签隔离开关。
2. `solveDiscreteFlowAdjointProduction.H`：`discreteProdPressurePrecSolver`
   （GAMG 默认 / PCG+DIC 备选，非法值 FatalError）；GAMG 稳健系数
   optProperties 覆盖（scaleCorrection/nPreSweeps/directSolveCoarsest，
   默认=OF7 现状）；内层 psi 非有限快检 FatalError（标签+apply 序号+提示）；
   **门改 const lduMatrix& 访问（本任务的核心修复）**；PRODPRECGAMGSETUP
   加打 solver 词。
3. `test_stage_b_safety_gates.py`：新增
   `test_bfinal011_label_isolation_switches_and_inner_solver_controls`，
   含「禁止非 const `.lower()` 物化访问」防回归断言（7 tests OK）。

## 下一步（BFINAL-011 后续 / BFINAL-012）

- 两标签生产伴随已收敛 → 冻结梯度验证链（BFINAL-012 FD 幅值门）具备输入。
- 可选：E2'（对称矩阵下的 PCG 行为对比）、并行伴随（P2，见 §8.5 顺序要求）。
