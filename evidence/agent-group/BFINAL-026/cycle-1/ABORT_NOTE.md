# BFINAL-026 cycle-1 尝试 1 — 中止记录（ABORT_NOTE）

- 记录时间：2026-08-21 22:45 (+080)，协调人（审阅者）
- 执行者：executor-working-time（DeepSeek v4 Flash），运行 8.5 小时、2209 次工具调用、85.9M token
- 结局：**空完成**（最终消息为空；无 EXECUTOR_SUMMARY、无 EXECUTOR_DONE、无测量数据、无 git 提交）

## 事实

1. 有效产出仅两件：`PREREGISTRATION.md`（14:09，判别矩阵，质量合格）与
   `NOTEBOOK.md`（17:59 后再无更新）——后者含高价值源码钉扎（金字塔质心公式、
   D1/D2/D3 精确构造、filter/projection 切线、drAU/dHbyA 语义、模板缺陷清单、
   b8_verify_diag ≡ b25_qgate 场状态逐位相同的发现）。
2. 17:59–22:30 约 4.5 小时零磁盘产出。死因：w_true 所需的 SuperLU 分解
   （523MB / 15.09M nnz 显式矩阵）被以 280s 超时 × 线程数扫描的方式反复
   重试（该求解实测需 ~25 分钟，任何 280s 超时必然失败），叠加上下文压缩
   循环（token 消耗为 B25 的 7 倍），直至耗尽。协调人 21:57 的解堵消息
   （后台长超时 + MMD 配方 + H7 矩阵纠正）送达时已过深，未能挽救。
3. 无生产代码改动（git 零修改），纪律未破——是能力/结构问题，不是违规。

## 教训（写入后续派发实践）

1. 单轮范围必须匹配执行者上下文容量：本轮打包了「源码钉扎 + 模板考古 +
   自建离线实验室 + 长时 LU 求解 + 四项测量 + 报告」，过重。
2. 长杆数值任务（LU/GMRES 类）要么由协调人预先算好作为输入交付，要么
   在任务书里给出精确到命令级的配方 + 后台运行模式，不给执行者留下
   「自行摸索求解策略」的空间。
3. 已潜伏 ≥2 小时零产出（NOTEBOOK 停更）应触发协调人主动干预，而不是
   等到自然结束。

## Salvage 路径

- 协调人直接执行 `b26_wstar_h7.py`（b21 模板 + H7 矩阵 + MMD 配方 +
  b25_qgate 状态验证），产物 `b26_wstar_h7.npz/json` 作为仪器工件供
  尝试 2 使用；闭合验收 |r_P(w*)|/|rxdP| ~1e-11（>1e-9 则 BLOCK）。
- 尝试 2（瘦身派发）：M1 源层恒等式为硬交付，M2–M4 为弹性目标；
  输入全部就位（w_true、b25_qgate/b18rhs_thermalCoupling.mtx（Q 型 b_TC）、
  b8_verify_diag 的 rxc/rpa/gsensh 系列导出、B25 FD tsv、NOTEBOOK 源码钉扎）。

本目录中 attempt-1 的 PREREGISTRATION.md 与 NOTEBOOK.md 保留原样，
尝试 2 的新文件一律加 `b26a_` 前缀避免覆盖。
