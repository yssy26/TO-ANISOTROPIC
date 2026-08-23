# BFINAL-013 cycle-1 — Independent Review (REVIEW_BY_REVIEWER)

- 生成时间：2026-08-20 01:02 (+0800)，应用户告知提前执行（定时自动化同款流程）
- 依据的 EXECUTOR_SUMMARY.md 修改时间：2026-08-20 00:43
- 审阅基线：`f874b71`；审阅时 HEAD：`5aaae31` "BFINAL-013 cycle-1: localize gradient-assembly defects to sensitivity.H contraction terms"
- 审阅方式：只读复核（mtx 原始导出 L2 重算、三份日志 grep、复现 tsv 与 B12 逐位 diff、静态测试实跑、P0 行号抽查、git diff）

## 1. 复核结论表（声明 vs 证据 vs 判定）

| # | 执行声明 | 证据 | 判定 |
|---|---|---|---|
| 1 | gDP 缺陷定位于 `sensitivity.H:89-90` 动量行收缩项，原始梯度占比 97.6% | **本人从 `rxpr_prod_gsensh_*.mtx` 重算**：\|momentum\|L2=0.28155、\|pressurerow\|L2=0.00839、\|total\|L2=0.28864，momentum 占 97.5%；行号与我此前通读 sensitivity.H 的记录一致（`:89-90` `gsenshPressureDropMomentum = -dAlphaDxh*(U & Uc); *= mesh.V()`） | ✅ |
| 2 | J 缺陷定位于 `sensitivity.H:65-66` Brinkman 收缩项（热扩散项可忽略） | P2 分项：D1 brinkman +4.59e-3 vs thermal −1.4e-4；`:65-66` 行号核实；与 gDP 缺陷同构（同为 `−(U & U_adj)·V·dAlphaDxh` 收缩） | ✅（附注：见 §3-1 推断性说明） |
| 3 | P1 源点测试（新仪器）：源全部正确 | `Log.verify_b013_p2.txt:3991-4013` 实存：P1a rel 2.8e-8→8.8e-7（5 档全 PASS）、P1b rel 1.7e-10→1.2e-8（全 PASS）、P1c FD 序列 +1.42→+16.6→+266.9→−372.8→−307.95 渐近收敛至 ADJ=−301.39（1e-4 处 rel 2.2%） | ✅ |
| 4 | P1 首跑空转（方向支撑不覆盖源，0==0）被如实保留 | `Log.verify_b013.txt` 保留原始空转行，未选择性引用；修正方向改为源支撑上的确定性场 | ✅ 过程诚实 |
| 5 | P2 装配自洽（xTotal 逐位复现 B12 ADJ） | `Log.verify_b013.txt:8779-8782` D 行：D1 xTotal=−5.4051 / D2 +1.1059 / D3 +14.0447，与 B12 `stage_b2_summary.tsv` projDP 逐位一致；并明确记录 `filter_chainrule.H:148-186` 就地改写导致「pressureRow(diff)」列为混合量、原始分项以 P3 导出为准（防误读注释到位） | ✅ |
| 6 | P3 R_x 内部一致性复跑 relErr 1.7e-15 | `Log.verify_b013_b6oracle.txt:3971` 实存：relErr=1.737e-15，alphaRel=0.4 | ✅ |
| 7 | 探针零扰动（B12 全套 FD 扫描逐位复现） | **diff B12/B13 复现 tsv：逐位一致**；主日志中伴随 750/902 与 b11/b12 逐位一致；MTO_RC=0 | ✅ |
| 8 | 范围纪律（定位即停，未修复；仅新增开关门控诊断代码） | `git diff f874b71..5aaae31`：新增 `stageB13GradientProbe.H`（260 行）/`stageB13ContractionProbe.H`（115 行）+ 接线（MTO_HF.C +3、validateStageB2 +14，均为 `stageB13GradientProbe` 默认 false 门控）+ 静态测试 +42 + 证据；**sensitivity.H / AdjHeatTransfer.H / 算子层零改动**；静态测试 8/8 实跑通过 | ✅ |
| 9 | P0 双侧因子对照表 | 抽查两处行号属实：`computeObjective.H:49-56`（`−φ/(M·Tref)`）、`createFrozenHotRegionFields.H:640+`（`+ρ·A_f/(PDmax·A_tot)` 入口、出口按 `!fixesValue()` 门控——与 fixedValue 出口 p=0 的状态导数语义一致）；表内未见单侧 ≈2.13 因子 | ✅ |

## 2. 结论分级

**已验证事实**
- F1 两个缺陷的**载体**定位成立：目标定义层（P0）、目标导数源层（P1 机器精度）、滤波/投影链（gV 1e-7 穿同一链）、R_x 压力行（2.9% 占比 + 1.7e-15 内部一致）、求解器/矩阵（BFINAL-010/011 门 + 逐位复现）——全部被定量证据排除；剩余载体 = `−(U & U_adj)·V·dAlphaDxh` 类动量行收缩（gDP 侧占原始梯度 97.5%）。
- F2 探针零扰动：B12 FD 扫描与伴随在探针二进制上逐位复现。
- F3 执行纪律良好：定位即停、无修复、空转探针如实保留、证据先归档。

**强假设**
- H1 gDP 的 ≈2.13 恒定因子源于该收缩项的**数学定义不完整/不一致**（H-F1 松弛语义 / H-F2 算子中 Sp(α,U) 之外的 α 依赖 / H-F3 FD 不动点与线性化不动点的系统差异）——载体已锁定，机理待修复轮推导。注意原始占比 97.5% 是**链前**份额；链式层是线性但空间再分布的，严格来说修复轮若仅修动量项而比值未归一，则需回头查链后分摊（预期不会，因为错误是干净的方向恒定标量）。
- H2 J 与 gDP 共享同一收缩缺陷（同构性 + Brinkman 项主导 + 源正确）——**最干净的确认实验就是修复轮：单一收缩修正应同时使 J 与 gDP 归一**；若只归一其一，说明 J 另有叠加因素。

**推测**
- C1 修复不需要触碰已验证算子层（J_PU/J_PP/R_x/J^T），仅收缩层即可闭合——若推导表明需要动算子，属重大范围变更，须单独授权。

## 3. 注记（不阻塞）

1. **J 的定位是「强但推断性」的**：直接证据（恒定标量、97.5% 占比）只在 gDP 侧；J 侧依赖同构论证（源正确 + Brinkman 主导 + 同一收缩公式）。修复轮把「J 与 gDP 同时归一」列为验收项即可闭环。
2. P1c 的中间 eps 偏差解读（有理函数曲率）物理上合理，且渐近收敛是硬证据；无需进一步验证。
3. P1 源点测试已具备「常驻验证链」价值，建议修复轮后纳入静态测试清单引用。

## 4. 审阅判定

**定位成立（PASS），进入分支 A → BFINAL-014 修复轮。** 修复轮的核心纪律：**先推导后实现**——把 FD 链实际收敛的离散不动点残差写出来、对其求导、与现行收缩逐项对照；禁止用拟合 2.13 的方式猜因子（executor 自己已声明 2.5/1.92 都不是 2.13）。任务书见 `/home/ys/dsH/TO-ANISOTROPIC/NEXT_TASK_BFINAL014.md`。
