# BFINAL-019 cycle-1 — 执行摘要（结论先行）

**任务**：J 装配完成轮——补 `−pb^T·R_P,x`（热伴随压力 pb 的压力行收缩）与
`Gx = g^T·φ_x`（面通量直接 α 项）两项进生产 J 梯度（功能口径：J(D1) 符号恢复）。

## 一句话结论

**两项按 DERIVATION_FIX3 精确落地（遗留改动审计结论 A：方向正确、经修正零处直接
沿用；全部不变性门逐位通过：等价性门四值、gDP/gV 投影与 FD 列、双标签 1e-12、
GRADSTABLE、静态 12/12、WMAKE_EXIT=0）；但 FD 门 J 符号仍 2/3——D1 仍反
（FD 平台 +0.0429 vs ADJ −0.03475）。两项的 D1 净贡献 +0.00274，仅为翻号所需
（0.0348）的 7.9%。按任务条款 5 停止并完成三分量定位：装配层已闭合且无罪，
残余全部落在缺陷①（算子通量图语义）在近对消 TC 泛函上的作用 + 已实证的状态
敏感性；gDP 侧（符号 3/3、单因子 2.130/2.113/2.169）与 gV 与 B18 逐位相同。**

## 验证序列结果对照

| 步骤 | 要求 | 实测 | 判定 |
|---|---|---|---|
| 1 静态+编译 | 12/12 + WMAKE_EXIT=0 | **12/12**（原 11 + 新 1 `test_bfinal019_j_assembly_pressure_row_and_flux_direct`）；WMAKE_EXIT=0（wmake_b19.log，~37 min 全量）；二进制含 `BFINAL-019 J assembly decomposition` 标记 | ✅ |
| 2 等价性门 | 四值逐位不变；gDP 投影与 B18 逐位相同；gV 逐位相同 | 四值两标签 **1.13142348768e-16 / 8.25782114104e-17 / 1.42213008543e-16 / 1.34302013134e-16 逐位相同**（b19_main 全日志 diff：除时间戳/新增 J 分解行/fsensMeanT 过滤残差/dJ 活动范围外逐行相同；dV、dDP 范围逐位相同）；FD 门 projDP/projV/ADJ_gDP/ADJ_gV 及全部 FD 列 **13 行字符串级相同** | ✅ |
| 3 双标签 1e-12 + GRADSTABLE | ≤1e-12 + 稳定 | TC 898 iter 9.80383374526e-13 / PD 1130 iter 9.68035507374e-13（主循环轮）；TC 894 iter 9.82e-13 / PD 855 iter 9.36e-13（B2 基线轮）——**与 B18 逐位相同**；GRADSTABLE 4/4（5.14e-11/3.53e-13/3.65e-10/2.82e-12）；0 次 NOT stable；MTO_RC=0 ×2 | ✅ |
| 4 FD 门（功能口径） | **J 符号 3/3**；gDP/gV 符号不回归 | **J 符号 2/3：D1 仍反（FD +0.0416…+0.0437 恒正 vs ADJ −0.034754）**；D2/D3 保持（ADJ −0.023814/−0.013641 vs FD −0.003469/−0.000283）；gDP 符号 3/3 且比值 2.130/2.113/2.169 逐位复现；gV 符号 3/3、bestRelV 1.58e-6/3.60e-6/6.83e-6 逐位复现 | ⚠️ J(D1) 未过 |
| 5 J(D1) 仍反 | 停止 + 三分量定位，不重试不调参 | 已停止（无任何重试/调参/拟合）；定位见 FINAL_REPORT §4 | ✅ 按条款执行 |

## 交付内容

1. `src/rxPressureRowTranspose.H`：宏参数化（RX_ADJ_PRESSURE/RX_OUTPUT_FIELD/
   RX_OUTPUT_NAME/RX_LABEL/RX_FLUX_DIRECT_SOURCE；无宏定义时展开与历史文本
   token 级一致——`audit_token_compat.txt`，仅诊断 Info 标签字符串与一处换行
   不同，数值零影响）+ Gx 直接项基场 rxFluxDirectT（逐槽位镜像算子
   dphiHbyA_w − dflux_w 通道，U-assignable 边界通道与 b18_folding_lab.py
   `Gx_of` 完全一致）。
2. `src/sensitivity.H`：热耦合面泛函 g 的逐值复刻（与 AdjNS_HT.H 5–41 行同构）、
   pb 参数化二次收缩、`gsenshMeanTPressureRow = −rxPressureRowTb·dAlphaDxh`、
   `gsenshMeanTFluxDirect = +rxFluxDirectT·dAlphaDxh`、接入 fsenshMeanT 链
   （C 项之后、freeze 零化语义保持）、分解 L2 打印 + 门控 mtx 导出。
3. `src/tests/test_stage_b_safety_gates.py`：第 12 项静态测试。
4. 证据：本目录（两算例日志、FD 门 tsv、对比表、审计工件、编译日志）。

## 关键数字

- 新两项运行内分解（B2 基线态）：|pressureRow(pb)|L2=4.1494e-4、
  |fluxDirect(Gx)|L2=3.2261e-4、|momentum|L2=2.2712e-3、|thermalC|L2=1.1131e-3。
- 新两项 D1 净投影（同 λ、in-pass）：**+0.002735**（B18 冻结态 λ 估计 +0.004714，
  两状态符号一致——装配正确性的旁证）。
- D1 翻号所需：≥ +0.0348（ADJ 归零）；两项提供 7.9%；占 ADJ↔FD 全距 0.0776 的 3.5%。

## 纪律

未触碰：J 算子/预条件、AdjNS_HT/AdjHeatTransfer rhs 折叠、
filter_chainrule/过滤/投影/MMA/目标定义、等价性门阈值、
`frozenGradientValidated=false`、`mmaUpdateEnabled=false`；无拟合因子；
FD 门 J(D1) 未达功能口径即按条款 5 停止（详见 FINAL_REPORT.md）。
