# BFINAL-021 cycle-1 — EXECUTOR_SUMMARY

任务：P 行表示级差异定位轮（零代码优先）。分支 `agent/dsH-stage-b-validation`，
起点 `45084c7`。执行 2026-08-21。**结果：表示级差异被定位到五个具体槽位；
「窗口收敛伪影」被定量排除（份额≈0）；三个量具链缺陷为本轮新发现并量化；
残余真失配收敛为「≤3.5%（D1/D3，=仪器地板）至 45%（D2）的面级通量切线
残差，经 div 收缩/压力反解放大成 O(1)」。** 唯一代码新增 = 开关门控导出
（E3 需要），无任何语义改动。

## 一行裁决

E1/E2/E3 联合：**混合**——伪影≈0%；量具链缺陷（SLOT-1/2/3，真实、已量化、
单独不充分）；真失配残余（SLOT-5：面级通量切线百分比级残差 + 放大机制）。

## 关键数字

- 恒等式（精确导出 J 的 P 行，h=1e-3）：|r_P|/|rxd_P| = 1.662/2.136/1.514，
  跨 eps（3e-4→3e-3）变化 <0.6% → **eps 平坦 = 系统性真失配签名**。
- 探针终态通量 div-free：|div(φ± 差分)| = 1.8e-14…1.2e-13 → 伪影上限
  ~5e-11 ≪ |rxd_P|≈2e-4（7 个数量级）→ **E1(b) 排除窗口收敛伪影**。
- 干净未钉扎 SuperLU（1769 s）：r_P(w*) = 1.8/1.3/3.0 e-11（机器精度）→
  失配 ≡ J_P(A)·(w_state − w*(A))；|w_state−w*|：relP 0.759/1.024/1.172
  （relU 0.180/0.526/0.190）——p 块 O(1) 失配是硬事实。
- **A-B 线性化点/系统失配（本轮主发现）**：FD 真值@B（模块 phase-1 246 次
  收敛、冻结湍流）vs 算子/rxd/rhs@A（主循环 44 corrector、分子粘性）：
  U 13.57%/p 10.63%（−4800 Pa 均值偏移）/phi 13.41%。但其对算子系数影响
  小（rAU、kf 的 B/A 比中位 1.0000，仅 3.1-3.3% 面超 ±10%）→ 单独不充分。
- 面级总通量切线 vs 真 flux-FD（重跑新量具）：A 基 relL2 0.034/0.451/0.036，
  B 基 0.087/0.449/0.090；仪器地板 3.3-3.6%（内面）——D1/D3 达地板，
  **D2 两基均 45% = 真缺陷信号**。
- b20-augmented 工具 bug（本轮发现）：da 缺设计链（rawD vs z=filter+投影
  切线，|z|/|d|≈5.1-5.9）→ alpha 通道量级小 7×——B20 的「rxd_P_c 14%」类
  闭环结论作废。
- 模块 round-2 伴随 stale kf 基：primalPressureMobility 只在 NS.H 捕获
  （A 点分子基），B 点运行时未刷新 → 污染「ADJ/FD=2.13×」量具。

## E3 重跑

`b21_anchor`（b13_probe3 配置 cp -a + 新导出），23 min，MTO_RC=0；
FD 侧与存档逐位一致（fd_scan FD 值/nIter 全同、wstate U/p bit-identical；
仅 ADJ_J 列因 B18 fix3 源码演进不同）。新导出：基线+14 探针的
phi/phiB/nutFrozen/nuEffFrozen/alpha。

## 纪律

- 源码唯一改动 = `validateStageB2GradientAmplitude.H` 两处 stageB15StateExport
  门控导出块（+73 行，默认关零影响；git diff 已核）；
- 未触碰算子/装配/源/链/MMA/阈值/锁定开关；无拟合因子；未 push。
