# BFINAL-016 cycle-1 — 第〇阶段（量具对齐）FINAL REPORT

**结论：第〇阶段完成。A1 裁决「产/导差异」为钉扎伪影——生产与导出是同一算子；
算子对真值的方向依赖失配（1.77 / 0.86 / 2.27）成立且量具已清洁；B1 确认 b_TC
缺陷并给出清洁幅度表；另发现一处新的装配层内部矛盾（sensitivity 装配 vs
λ-direct，D2 相差 147%）。按任务「先报告再继续」条款，A2/B2 修复以本报告的
修正后靶点为准进入下一 cycle。**

## A1 — 偏应力对齐（裁决 C1）

`b16_a1_dev.py`（J_dev 按 production.H:105-231 公式 numpy 复现，网格来自
`stageB16MeshExport` 一次性导出，96860 内面/33600 单元/7880 边界面）：

```
加入 J_dev 前后 b_PD^T w′: -4.3260030677 → -4.3260185401  (Δrel ≈ 3.5e-6)
与生产 ADJ 的偏差保持: D1 20.0% | D2 61.5% | D3 1.5%
```

**偏应力块数值上可忽略**（nuEffFrozen 量级 ~1e-2），不解释产/导差异。

## A1b/A1c — 「第三处差异」的定位与裁决（重要修正）

`b16_a1b_lambda.py` + `b16_a1c_clean.py`：
- λ_prod（状态运行写盘的 Uc/pc，生产 FGMRES 解）与 λ_exp（导出系统转置解）
  一致：cos 0.99992、relL2 3.6%（3.6% 来自离线解的 P0 行钉扎伪影）；
- **干净块差分：`M_unpinned @ λ_prod − bPD` 的 |r| = 8.02e-11（P 行 5.8e-16）**
  ——生产伴随以机器精度满足导出系统方程。

**裁决：不存在第三处算子差异；生产 matrix-free 与导出显式矩阵是同一算子。**
BFINAL-015 报告中「两个彼此不同、各自错误的算子」([A] 检查) 是离线钉扎行在
转置求解中的副作用所致，予以更正。产/导对齐完成（以 λ_prod 满足 M·λ=b 为证）。

## 修正后的算子缺陷量具（干净口径）

用未钉扎恒等式（λ 满足 M λ = b、rxc 外部验证 ~1e-7）：

```
λ-direct 收缩  −λ_prod^T rxd:  D1 −4.4859 | D2 +0.4479 | D3 +14.7303
真值           b^T w_true:      D1 −2.5375 | D2 +0.5234 | D3 +6.4772
比值                            1.769      | 0.856      | 2.274
```

方向依赖失配成立（非标量）；失配块签名（BFINAL-015 的 P 行集中度系钉扎口径
下测得，**完整算子+未钉扎口径的复测留待 A2 修复轮**，本报告不沿用旧签名）。

## 新发现 — 装配层内部矛盾（超出本轮预注册框架，如实报告）

```
sensitivity.H 装配（B12 ADJ）:  D1 −5.4051 | D2 +1.1059 | D3 +14.0447
λ-direct（同一 λ、同一 rxc）:   D1 −4.4859 | D2 +0.4479 | D3 +14.7303
相对差:                          20.5%      | 147%       | 4.9%
```

数学上 −λ^T·(R_x d) 与「先解后收缩、链式后置」应严格相等（标量的转置），
除非 sensitivity.H 的分项装配与 rxc 的 R_x 构造在**动量行或压力行的转置语义**
上不一致。D2 的 147% 差异远超算子缺陷本身——**gsenshPressureDrop 的装配链
（gsenshMomentum + rxPressureRowT + filter_chainrule 的组合）存在独立于算子
的缺陷**。该矛盾使 BFINAL-013「定位到 −(U&Uc)V·dAlphaDxh 单项」的结论需要
修订为：该项公式在给定 λ 与给定 R_x 时正确，但装配链整体与 λ-direct 不等。

## B1 — 热介导量化（缺陷 2 清洁化）

`b16_b1_thermal.py`（T 状态设计 FD + P1b 验证过的 dJ/dT 源；恒等式两种符号
约定均不闭合，任务判据触发——缺陷 2 有未建模成分）：

```
                  D1         D2          D3
FD_J（真值）      0.04286    −0.003469   −0.000283
热介导（实测）    −0.16687   +0.08230    −0.05665
真流介导=FD−热    +0.20973   −0.08577    +0.05637
b_TC^T w_true     +0.11251   −0.02559    −0.04481   ← 与真流介导差 86%/70%/符号翻转
```

**b_TC ≠ T-消除后的 dJ/d(U,p)**，缺陷 2 确认且独立量化（所用仪器全部外部验证：
w_true 真值、dJ/dT 机器精度、FD_J 真值）。

## 修正后的缺陷地图（A2/B2 的输入）

| 缺陷 | 修正后定位 | 干净量具 |
|---|---|---|
| 算子（原 gDP 2.13×） | 单一共享算子的 P 行/通量耦合语义（产=导，无第三差异） | −λ^T rxd / b^T w_true = 1.769/0.856/2.274 |
| 装配（新） | sensitivity.H 收缩装配链与 λ-direct 不一致（D2 147%） | sensitivity 值 vs −λ^T rxd |
| 源（原 J 目标） | b_TC 热消除折叠（AdjNS_HT/AdjHeatTransfer 链） | b_TC^T w_true vs (FD_J − 热介导) |

## 下一 cycle 的修复顺序建议

1. 先审计装配矛盾（最小代价、可能揭示算子缺陷的真实形状）：逐项对比
   gsenshMomentum/gsenshPressureDropPressureRow 与 rxc 的 anRUa/assembleWeightedRPa
   构造；
2. A2 算子修正（以 1 的发现为准）；
3. B2 源修正。

## 过程记录（诚实性）

- 钉扎伪影两次干扰读数（BFINAL-015 的 [A]「双算子」结论、A1b 的 [3] 块残差），
  均已在证据中保留原始输出并更正；
- 全部代码新增为开关门控量具（stageB16MeshExport、T/dJdT/出口元数据导出），
  默认关，零数值影响；静态测试 9/9；未触碰任何算子/源/装配语义。

## 工件

| file | note |
|---|---|
| `b16_a1_dev.py` + `a1_console.log` + `a1_dev_alignment.json` | A1 偏应力对齐 |
| `b16_a1b_lambda.py` + `a1b_console.log` + `a1b_lambda_localization.json` | λ 定位（含钉扎伪影原始读数） |
| `b16_a1c_clean.py` + `a1c_clean.json` | 干净块差分与 λ-direct 收缩 |
| `b16_b1_thermal.py` + `b1_thermal.json` | B1 热介导量化 |
| `../b16_mesh/Log.verify_b016_mesh.txt`, `../b16_states/Log.verify_b016_states.txt`（运行目录） | 量具运行日志 |
