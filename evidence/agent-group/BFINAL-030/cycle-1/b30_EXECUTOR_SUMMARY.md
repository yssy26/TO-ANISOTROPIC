# BFINAL-030 — T3/(1-alphaRel) 语义审判轮（执行器总结）

> 收尾提交（2026-08-22）。任务书 `NEXT_TASK_BFINAL030.md`，预注册判据
> `b30_PREREGISTRATION.md`，执行日志与推导 `b30_NOTEBOOK.md`。
> 全部工作纯离线 Python（路由仪器 / explicitJT_H7 系数级补丁 / RHS 文本文件）：
> **零编译、零求解器运行、零生产代码改动**；`src/` 未触碰。
> 仓库 `/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`，
> HEAD `cd9bd7a`（BFINAL-029 提交）——本轮无任何提交之外的变化。

## 1. 一句话裁决

**部分（partial）**：T3 槽的 (1-alphaRel)=0.6 直接松弛系数**被证实为一个真实贡献者**
（校正 0.6→1.0 使 fp2 U 行幅值 0.696→0.959、fp3 D1 −0.656→+1.126、gDP D1 因子
2.05→1.013），但**单一 T3 槽修正不能闭合指纹 2/3**（fp2 relL2 地板 25.8% vs
目标 <5%；fp1 主判据三方向全 FAIL；fp3 D3 从 +1.057 恶化到 −1.763）。
统一单槽假说**未定罪**；T3 是至少两个病灶之一，不是唯一病灶。

## 2. 推导核心结论（Stage 0）

- OF7 relax() 语义（fvMatrix.C 521-668，源码逐行验证）：松弛对流通量的切线为
  `dphi_f = Sf&( alphaRel*(w*rAU_o*dH_o+(1-w)*rAU_n*dH_n)
  + (1-alphaRel)*(w*dU_o+(1-w)*dU_n) ) − kf_f*(dp_n−dp_o)`。
- 生产算子实现的正是该码径：T3 槽 (1-alphaRel)=0.6（Production:363-369），
  B18 DERIVATION_FIX3 亦同——两者一致说明**码径本身就是 0.6**。
- 判定按**不动点语义**采纳主修正常数 **c*=1.0**（SIMPLE 松弛因子在外循环收敛处
  是单步预条件子，(1-alphaRel)*dU 项在定点消失；正确稳态通量切线为 Sf&dU）。
  备用灵敏度：c*=1.09（0.6·1.82）、c*=1.23（0.6·2.05）。
- B24 GateOracle 佐证：码径算子切线 relL2=1.36 vs 真原始手 FD（1.65e-11 闭合），
  导出 explicitJT_H7 == 隐式算子（5.48e-16）⇒ 失配在算子语义本身。

## 3. 三指纹实测 vs 预注册（主常数 c*=1.00）

| 指纹 | 预注册目标 | 实测（c*=1.00） | 判定 |
|---|---|---|---|
| ② 路由 vs b_TC relL2 | <5% | 0.397→0.296（c*=1.0）；**扫掠最低 0.2582 @ c*=0.85**；U 幅值 0.696→0.959 | **否证阈值**（地板 25.8%）；方向成立 |
| ③ b_TC_corr^T w_true | 三方向 \|ratio−1\|≤15% | D1 **+1.1265 过**；D2 −32.8（大数相消，预注册非担责）；D3 **−1.763**（生产 +1.057 反被破坏） | **仅 D1 成立**；D3 给出反证信号 |
| ① 补丁重解（主判据） | 三方向 ≤15% | D1 −0.502、D2 40.5、D3 13.5（T6=corr）；T6=keep 更差 | **全 FAIL** |
| ① gDP 因子预测 | ADJ/FD 2.05 → 1.0±0.15 方向均匀 | D1 **1.0134 过**；D2 0.161、D3 0.752 | **仅 D1**；方向均匀性未复现 |

定罪条件（四项同用 c*=1.00 同时满足）**不满足** ⇒ 按预注册映射为**部分**，
非无罪（多通道方向性移动确凿）。

## 4. 关键新发现（本轮）

1. **对角 U 行/P 列 T3 位置被 hA/hB Stage-2 路径污染**（solveDiscreteFlowAdjoint.H
   2722-2807 写入相同 (row,col)：`-upperF*Vinvo*hBx`、`-lowerF*Vinn*hBx`、hA
   cellFacesOf 循环）。'nn' 位原值 ≈ −0.6·vv（反号，偏差达 vmax）；'oo' 位
   orig/vv≈1.66。⇒ 离线补丁只能打**非对角槽**。
2. 非对角槽也非全纯：193,720 个显著位中 183,812（94.9%）恰为 0.6·vv
   （中位 rel 1.05e-14）；9,908（5.11%）在边界邻胞处复合（maxrel 0.398）。
   ⇒ 补丁门控必须用**位置级掩码 |rel|<1e-6**，排除位保持生产复合值并记录。
3. 边界 U 行/P 列块为**复合块**（hB `Bdiag*Vinvc*(alphaRel*rAUb)*patchSf` 比例
   ~0.3 + 直接 (1-alphaRel)*patchSf + 内面累积），非纯 T3 槽，补丁排除。
4. **MMD_AT_PLUS_A 在本近奇异矩阵上永不收敛**；scipy 1.10.1 `splu` 传
   `diag_pivot_thresh` 必须用直接关键字（`options=dict(...)` 报非法关键字）；
   唯一验证路径 = 默认 COLAMD `splu`（b26 已验证，~22 分钟）。
5. **defcheck（必须报告）**：`b_PD^T w_true / ADJ_gDP = 0.8178/0.3660/1.0362`
   —— b_PD 连 w_true 都不能复现 ADJ_gDP（B13 分解本身有缺口）；gDP 因子预测
   承载此模型误差。ADJ_gDP 锚 D1 −5.2303 / D2 1.0505 / D3 13.6034；b24 实测
   因子 D1 2.062 / D2 2.007 / D3 2.100。B13 已将 2.05× 的主体定位在动量行通道，
   非 T3 通量槽。

## 5. 偏差声明（相对任务清单，如实列出）

1. **fp1 补丁范围 ≠ 预注册配方全文**：预注册配方含对角 4 项 + 边界项，实测
   证明对角被 hA/hB 污染、边界为复合块，故补丁实际仅打内部**非对角 2 项**
   （[U(own),P(nei)]/−w·sf 与 [U(nei),P(own)]/+(1−w)·sf）。这是对预注册配方
   的必要收窄，已在 NOTEBOOK §1.3 与 resolve JSON meta 记录。
2. **灵敏度 c*=1.09 / 1.23 未做矩阵重解**（每个需 ~22 分钟 COLAMD splu）。
   协调人指示以现有数据收尾、不扩展扫描；路线级 fp2/fp3 已覆盖 c*=1.09/1.23
   （方向性一致：D1 超配、D2/D3 发散），不影响「部分」裁决。
3. **B31 生产修复规格未写**：预注册约定定罪才写 B31；裁决为部分，未达阈值。
   `src/` 零改动。
4. **fp2 细扫 T8=on/off 结果相同**（T8 槽贡献为零）——作为仪器自洽性观察记录，
   不构成对 b29 排除结论的推翻。
5. 预算实际超过 180 分钟（长 splu 1339.8s + 两次背景轮询），为任务书允许的
   长计算路径，未产生额外风险。

## 6. 产物清单（`evidence/agent-group/BFINAL-030/cycle-1/`）

- 推导与日志：`b30_NOTEBOOK.md`；预注册：`b30_PREREGISTRATION.md`
- fp2 路由重建：`b30_fp2_route_rebuild.py` / `b30_fp2_results.json` / `b30_fp2_route.log`
- fp2 细扫：`b30_fp2b_fine_sweep.py` / `b30_fp2b_results.json`
- fp3 M1 比值：`b30_fp3_m1_ratios.py` / `b30_fp3_results.json`
- fp1 补丁重解：`b30_fp1_patch_resolve.py` / `b30_fp1_resolve_c1.00.json|log`
  / `b30_fp1_wcorr_c1.00.npz`（npz 被 .gitignore 排除；可复现信息在 json）
- fp1 分析：`b30_fp1_analyze.py` / `b30_fp1_analyze_c1.00.json`
- 本总结：`b30_EXECUTOR_SUMMARY.md` + `b30_EXECUTOR_DONE`
