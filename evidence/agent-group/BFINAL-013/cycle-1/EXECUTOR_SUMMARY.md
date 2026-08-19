# BFINAL-013 cycle-1 — EXECUTOR SUMMARY (2026-08-20)

**结论：定位完成，未做任何修复（定位轮纪律）。两个 BFINAL-012 缺陷均定位到
`sensitivity.H` 的具体收缩项；其余所有层（目标定义、导数源、滤波/投影链、
R_x 压力行、算子、求解器）被定量证据逐一排除。**

## 定位结论（到文件、到项、附证据）

| 缺陷 | 定位 | 定量证据 |
|---|---|---|
| **gDP 恒定 ≈2.13× 偏大** | `src/sensitivity.H:89-90` `gsenshPressureDropMomentum = −dAlphaDxh·(U & Uc)`（×mesh.V()）——动量行收缩项 | 原始梯度 L2 分解：momentum 0.28155 / total 0.28864（**97.6%**）；ADJ/FD = 2.113/2.130/2.169（三方向恒定）；源点测试 PASS |
| **J 符号翻转 + 6.7–11.2× 放大** | `src/sensitivity.H:65-66` `fsenshMeanT = −dAlphaDxh·(U & Ub)`（Brinkman 项；热扩散项可忽略：D1 上 1.4e-4 vs 4.6e-3） | 源点测试 P1b/P1c PASS（源正确）；错误原始场经（已验证的）投影链再分布后表现为翻转/放大 |

两项缺陷同构：都是流场伴随速度（Uc/Ub）对 `Sp(alpha,U)` 的 α-导数收缩。
链式层无罪的证明：gV 穿过同一 drho/eta/滤波链与 FD 吻合到 1e-7。

## 探针结果速览

- **P0（读码双侧审计）**：FD 侧 `evaluateObjective/evaluatePressureDrop/costfunction.H` 与 ADJ 侧源/缩放公式逐项一致；无单侧 ≈2.13 因子。因子对照表已归档。
- **P1（源点测试，新永久仪器）**：gDP p 源 rel 2.8e-8→8.8e-7、J T 源 rel 1.7e-10→1.2e-8（全 eps 阶梯机器精度 PASS）；J φ 源渐近吻合（eps→0 收敛到 −301.39，1e-4 处 2.2%；中间偏差为有理函数曲率，非缺陷）。首跑 P1a/b 因方向支撑不覆盖源而空转（0==0），已如实保留并把方向修正为源支撑上的确定性场。
- **P2（B2 基线态分项分解）**：xTotal 列逐位复现 B12 的 ADJ（装配自洽）；`filter_chainrule.H` 就地改写 xh 场的事实已记录（防止把「链式值−原始动量项」误读为原始压力行项——原始分项以 P3 导出为准）。
- **P3（R_x oracle 复跑）**：当前 HEAD 内部一致性 relErr 1.7e-15（alphaRel=0.4）；伴随与 b11/b12 逐位一致。
- **确定性**：主探针运行完整重跑 B2/B3 FD 扫描，所有数字与 B12 逐位相同，探针零扰动。

## 修复轮假设（已标注，未验证）

- H-F1：FD 所微分的不动点携带 relaxed-SIMPLE 结构（alphaRel=0.4），而收缩按未松弛动量行处理——恒定幅度因子（符号不变）正是此类失配的特征；注意 1/0.4=2.5、1/(1−0.4·0.6)=1.92 都不等于 2.13，须推导而非拟合。
- H-F2：动量算子中 `Sp(alpha,U)` 之外的 α 依赖（冻结 laplacian 系数、fvOptions）未进收缩。
- H-F3：FD 侧 B2 探针原始解与伴随所线性化的不动点存在细微差异（p 松弛、最终非正交迭代、bounded-div Sp 路由）。
- 建议首选探针：用松弛语义重算动量行收缩；小网格 FD-of-FD 交叉验证。

## 仪器与代码

新增（开关 `stageB13GradientProbe`，默认关，只读）：`stageB13GradientProbe.H`（P1，主循环态）+ `stageB13ContractionProbe.H`（P2，B2 基线态）+ 接线 + 静态测试（8/8 OK）。两次编译均 WMAKE_EXIT=0。禁改清单未触碰；阈值/开关均未动。

## 证据

`evidence/agent-group/BFINAL-013/cycle-1/`：三份完整日志（主探针 / 修正方向 P1 / B6 oracle）、原始逐单元导出（momentum/pressureRow/total，链式前后）、B12 复现 tsv、算例快照、P0 因子表与决策表（FINAL_REPORT.md）。
