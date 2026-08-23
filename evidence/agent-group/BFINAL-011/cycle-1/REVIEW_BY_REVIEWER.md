# BFINAL-011 cycle-1 — Independent Review (REVIEW_BY_REVIEWER)

- 生成时间：2026-08-19 21:47 (+0800)，由「每30分钟检查BFINAL-011执行完成并自动审阅」自动化触发（应用户要求提前手动执行）
- 依据的 EXECUTOR_SUMMARY.md 修改时间：2026-08-19 21:00:21
- 审阅基线：`b01999c`；审阅时 HEAD：`2aaf649` "BFINAL-011 cycle-1: root-cause and fix pressureDrop GAMG NaN; both labels converge"
- 审阅方式：只读复核（OF7 源码、git diff、E0p/E2 归档日志、静态测试实跑），除本文件与 NEXT_TASK_BFINAL011.md 外未修改任何文件

## 1. 复核结论表（声明 vs 证据 vs 判定）

| # | 执行声明 | 证据 | 判定 |
|---|---|---|---|
| 1 | 根因 = BFINAL-010 门对非 const 矩阵调 `upper()/lower()`，OF7 非 const `lduMatrix::lower()` 物化 lower 为独立副本 → 矩阵被判非对称 → PCG 被拒、GAMG `scaleCorrection_=false` 走非对称分支 → PD rhs 首次内层解发散 | **OF7 源码独立核实**：`lduMatrix.C:167-182` 非 const `lower()` 确实 `new scalarField(*upperPtr_)` 物化；`:262-281` const 版返回 upper 别名不物化；`GAMGSolver.C:76` `scaleCorrection_(matrix.symmetric())`（本审阅者此前 grep 已见）。E2 日志 `Log.verify_b011_E2.txt:950-951` 实存 `FOAM FATAL IO ERROR: Unknown asymmetric matrix solver PCG` —— 物化发生的直接证据 | ✅ 成立 |
| 2 | 修复 = 门内经 `const lduMatrix&` 别名读取；门值逐位不变 | `git diff b01999c..2aaf649` 生产文件：`const lduMatrix& prodPrecLduMatrix = prodPressurePrecMatrix;` 后 const 访问 diag/upper/lower；E0p 日志 945-946/2375-2376 行四组门值 `1.13142348768e-16 / 8.25782114104e-17 / 1.42213008543e-16 / 1.34302013134e-16` 与 b010 **逐位一致**（两标签相同） | ✅ 成立 |
| 3 | TC 收敛 750 迭代 9.4956818492e-10 | E0p 日志 2365/2369 行；GRADPROXY-FINAL 2217.48004367（与 b010 的 2217.48004323 同平台，8 位一致） | ✅ 成立 |
| 4 | PD 收敛 902 迭代 9.92818139401e-10，全程 119.82 s，MTO_RC=0，零 FATAL | E0p 日志 3953/3955-3957 行；`ExecutionTime = 119.82 s`、`MTO_RC=0` 在日志尾；`grep -c "FOAM FATAL\|non-finite"` = 0；PRODRESID PD 有限（\|rU\|L2=1.02e-8，\|rP\|L2=1.47e-9） | ✅ 成立 |
| 5 | 对称路径恢复的独立佐证：TC 轨迹与 BFINAL-009 b8 部分数据逐位一致 | E0p @80=3.04473267442e-4、@240=4.90659573665e-5，与 `b8_verify_gamg/Log.verify_gamg.txt`（门存在之前的对称路径）逐位相同。b010（非对称路径）则为 750/751 迭代、9.489e-10 与 GRADPROXY 第 8 位后微差——与「路径改变但数值等价」自洽 | ✅ 成立（强佐证） |
| 6 | PD rhs 特征不变（U 块恒零） | E0p 2372 行 `sumRhsU=0 sumRhsP=95.1198507991 \|RhsP\|L2=10.3784265805`，与 b010/E1 引用值一致 | ✅ 成立 |
| 7 | 代码改动范围（隔离开关/内层求解器词/稳健系数覆盖/有限性快检/静态测试），默认全=基线，阈值零改动 | diff 逐块核对：`solveThermalCouplingFlowAdjoint`/`solvePressureDropFlowAdjoint`（默认 true）守两个 include；`discreteProdPressurePrecSolver`（默认 GAMG，非法值 FatalError）；scaleCorrection 默认 true / nPreSweeps 默认 0 / directSolveCoarsest 默认 false；psi 有限性快检 FatalError（标签+apply 序号+提示）；静态测试 7/7 实跑通过（含禁止非 const `.lower()` 物化访问的防回归断言）。外层容差 1e-9、等价性门 1e-8 未动 | ✅ 成立 |
| 8 | E1（PD 单独先跑）仍 NaN → 排除跨标签泄漏 H2 | **仅引用，原始日志已丢失**（见 §3 缺口）。判别逻辑本身正确：H3 下门在 PD 的 include 内照样物化，与「跳过 HT 无关」自洽 | ⚠️ 逻辑成立、原始证据缺失 |
| 9 | E0（修复前默认）复现 b010 基线 + 快检触发（33600 单元非有限） | 同上，仅引用。可由 E0p 门值逐位一致 + E2 物化证据间接支撑 | ⚠️ 同上 |

## 2. 结论分级

**已验证事实**
- F1 根因链完整成立（OF7 源码 + E2 归档日志 + 修复后行为，三方独立印证）：BFINAL-010 门自身的非 const 矩阵访问物化了对称存储，损坏求解器路径；TC 恰好能在损坏路径收敛、PD 不能，造成「标签特异性稳健性问题」的假象。
- F2 修复正确且最小（const 别名读取）；四组等价性门值逐位不变——**BFINAL-010 的等价性结论不受影响**（物化不改变数值，只改变存储/求解路径）。
- F3 两标签首次在生产迭代路径同时收敛 ≤1e-9（750/902 迭代，120 s，MTO_RC=0），**P0（生产伴随可扩展性与稳健性）实质闭合**：预条件器 setup O(nnz)、单次求解 2 分钟量级。
- F4 过程无掩盖：默认参数即胜局运行，PCG 备选与 GAMG 系数覆盖均未被使用；阈值未放宽。
- F5 静态测试 7/7（含物化防回归断言）。

**强假设**
- H1 「BFINAL-009 部分日志（240 迭代被杀）为健康数据」——现由 TC 轨迹逐位复现升级为近乎事实，但仍属单算例单设计点证据。
- H2 E1 的 H2 排除（跨标签无泄漏）——判别逻辑正确，原始日志缺失（见 F 缺口）；鉴于 H3 已完全解释现象，此项重要性下降。

**推测**
- C1 GAMG 非对称分支（无 scaleCorrection）对该 rhs 类必发散——单一样本观察（PD 发散、TC 幸存），机理层面（scaleCorrection=false 时 V-cycle 校正无缩放保护）合理但未做参数化验证。executor 保留的稳健系数覆盖项即为此兜底。

## 3. 与声明不符处 / 证据缺口（均不阻塞验收）

1. **E0/E1 原始日志丢失**（运行脚本初版 `rm -f Log*` 误删，脚本已修）：E0/E1 的载荷行仅存在于报告引用中。缓解：E2（物化直接证据）与 E0p（修复证据）均已归档，结论不依赖 E0/E1 的独有信息。**要求**：后续轮次日志一律先归档再清理（executor 已改脚本，遵守即可）。
2. E2'（对称矩阵下 PCG 行为）未跑——executor 声明为可选，同意：不影响验收，留作后续对照。
3. BFINAL-010 FINAL_REPORT 的「PD NaN 属 GAMG 稳健性问题」表述已被 executor 显式修正——处理得当，无需回改历史报告（保留修正记录即可）。

## 4. 审阅判定

**PASS。** 根因从 OF7 源码级别闭合、修复最小且经机器精度门值 + 独立轨迹复现双重验证、双标签生产收敛 + 干净退出、无掩盖操作、证据归档符合惯例（除已声明的 E0/E1 缺口）。计划决策树进入**分支 A**：下一任务 = **BFINAL-012 当前源 FD 幅值门**（见 `/home/ys/dsH/TO-ANISOTROPIC/NEXT_TASK_BFINAL011.md`）。
