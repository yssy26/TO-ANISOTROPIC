# BFINAL-032 cycle-1 — EXECUTOR_SUMMARY

任务：分解场门控导出轮。带 `stageB31WriteDecomposition`（默认 false）门控写出的
b25 同构重跑，让生产自己写出全部分解场（同轮次、同点），终结 B31 §7 的跨轮伪影
（76.5% 范数残差之谜）。分支 `agent/dsH-stage-b-validation`，起点 `4999855`。
执行 2026-08-22。授权：`NEXT_TASK_BFINAL032.md` + B31 协调人会话（`/1/` 目录
疑似混合伴随轮次的根因线索）。

## 一行裁决

**76.5% 跨轮残差 100% 消解为 0%（同轮逐位闭合）。** 同轮 `gsenshPD_raw == mom+prow`
逐位（G1 relL2 ~1e-12）；r1 闭合 = bPD_w、r2 闭合 = ADJ_gDP（残差 ~1e-13，逐位），
**B31 §7 的 D1 交叉数字（mom=−4.253786 / closure=−4.331615 / prod=−16.116120）
即混合伪影本身**。残余 ADJ−bPD_w = (D1 −0.952833, D2 +0.666062, D3 −0.492930) 是
干净的 r2-vs-r1 轮状态差（两轮 lambda 不同），**非病灶，证无真病灶**。FD 三列
硬闸逐位一致；G2 T 环镜像 2% 目标未达（r2 5.7%/3.2%，corr 0.998/0.9995）；G4
rxDHbyA 偏差 12.6%（r2）偏离 2% 校准先例；G-c 离线链求解为 near-singular 条件
伪影（源与矩阵两侧均已逐位隔离），链投影标记 INVALID，G3 锚定 raw_z/closure_z。

## 门裁决

| 门 | 预注册 | 实测 | 裁决 |
|---|---|---|---|
| **FD 硬闸** | b32 FD 三列逐位同 b25（sha256/cmp） | stage_b2_summary.tsv sha256 `90cef08a...`、stage_b2_fd_scan.tsv sha256 `2da03a58...`，cmp 字节级一致 | **PASS（逐位）** |
| G0 跨轮归因 | /1/ 伴随场 == 轮 1 | /1/Uc Ub pc pb Tb == `_r1` 拷贝 relL2=0.0（5/5） | **PASS（逐位）** |
| G1 同轮闭合 | relL2 ≤ 1e-12（逐位） | _r1 relL2=1.021e-12 maxAbs=1.0e-14；_r2 relL2=9.916e-13 maxAbs=1.0e-14；_r2 离线 mom/prow 成分 relL2 3.1e-12/2.7e-12 | **PASS（逐位）** |
| G2 T 环镜像 | relL2 < 2% | r1 pc 2.09e-1/pb 2.58e-1；r2 pc **5.66e-2**/pb **3.19e-2**（corr 0.998235/0.999508） | **FAIL**（r2 未达 2% 目标；预注册：门 2 失败不影响门 3） |
| G3 跨轮伪影消解 | 同轮闭合残差 < 5% | closure_z vs raw_z 残差全部 ≤ 2.2e-12（0%）；r1 闭合 == bPD_w、r2 闭合 == ADJ_gDP | **PASS（逐位，76.5%→0%）** |
| G4 HbyA 校准 | rxDHbyA/rxrAU 校准（2% 先例） | r2：rxrAU vs mob relL2=1.60e-12 corr=1.0（逐位）、rxG0Field vs g0 relL2=6.26e-13（逐位）、rxDHbyA vs dHbyA relL2=**1.26e-1** corr=0.9916 | **PARTIAL**（两项逐位；rxDHbyA 偏差 12.6% vs 2% 先例，r2 显著优于 r1） |
| G-mob | relL2 < 1e-12 | relL2=1.231e-14 | **PASS** |

## 关键数字

| 量 | 值 |
|---|---|
| 运行 | `/home/ys/dsH/b32_decomp`，ExecutionTime=1579.89 s，log `Log.b32_decomp.txt` |
| D1 r1 闭合 | mom_z=−4.190491 + prow_z=−0.086971 = **−4.2774612375490895**（== bPD_w −4.277461237487605，残 1.25e-13） |
| D1 r2 闭合 | mom_z=−5.155654 + prow_z=−0.074640 = **−5.230293845922022**（== ADJ_gDP −5.23029384626718，残 3.87e-13） |
| D2 r1/r2 闭合 | +0.384462806083 / +1.050524977639（== bPD_w / ADJ_gDP，残 2.2e-12/4.9e-13） |
| D3 r1/r2 闭合 | +14.09635332868 / +13.60342257618（== bPD_w / ADJ_gDP，残 3.4e-14/1.2e-13） |
| 残余 ADJ−bPD_w | D1 −0.952833 / D2 +0.666062 / D3 −0.492930（= 轮状态差，非病灶） |
| B31 §7 D1 旧数字 | mom=−4.253786 / prow_prod=−0.077828 / closure=−4.331615 / prod=−16.116120 —— 与同轮 r1/r2 均不同 = 混合伪影 |
| G-c2 链核对 | chain(raw_gDP_r2) vs disk gsensPD relL2=8.67e-1 corr=0.468（INVALID，见下） |
| G-c 源逐位 | chain_src(raw_gDP_r2) vs on-disk gsenshPD relL2=2.911e-10 |
| G-c 矩阵逐位 | solve(Lm+VI, Vx) |X1|=1.4199e+02 == 生产 xp 范数，corr=+0.999870（正向滤波器逐位复现） |
| G-c near-singular | 零模态放大 1/V ~ 8e9；|disk|=5.17e-1 ≈ |V·g|/V=6.3e-1；disk vs 源 corr=0.9796；(Lm@disk)[dm] max=4.7e-6；生产 DICPCG final=2.07e-10（10 iter） |

## G-c 链隔离判定（离线 chain_field 无效）

离线 `(Lm+VI)` 链后求解与 on-disk 链后场不符（G-c2 relL2=0.867 corr=0.468）——
但这不是矩阵或源错误：**源侧**逐位（chain_src vs gsenshPD relL2=2.9e-10），
**矩阵侧**逐位（离线 Lm 复现正向滤波器 corr 0.999870 且范数精确）。根因 = 系统
near-singular：(Lm+VI) 有常数零模态（全 zeroGradient BC，V~1.25e-10 vs Lm~1e-4），
生产 DICPCG（tol 1e-9，final 2.07e-10）与离线 splu 在零模态方向分叉；磁盘解
|disk|=5.17e-1 即零模态放大尺度 |V·g|/V=6.3e-1，disk 与其源 corr=0.98（输出≈光滑
输入）。**链投影 ch_total/ch_mom/ch_prow 标记 INVALID**（JSON `ch_valid=false`），
G3 判定锚定 raw_z/closure_z（直接门控导出，逐位、独立于离线矩阵）。全部细节在
`b32_decomp.py` JSON `gc.isolation` 与 NOTEBOOK §5.2。

## 改动清单（本轮唯一意图变更）

- `src/sensitivity.H`：`stageB31WriteDecomposition` 门控块（默认 false）——gDP
  两段+原始和、J 四段、rxPressureRowT/Tb/FluxDirectT、轮次伴随场 Uc/Ub/pc/pb/Tb
  拷贝的 round-tagged `.write()`；`b32SensInvocationCount` 持久计数与 `_r1`/`_r2`
  后缀；log 明示 `/1/` 伴随场属轮 1、灵敏度族属当前轮。
- `src/rxPressureRowTranspose.H`：`#ifndef RX_FLUX_DIRECT_SOURCE` 守卫内门控写
  rxrAU_r? / rxG0Field_r? / rxDHbyA_r?。
- 零污染：全部为纯 I/O `.write()` + Info 打印；默认 false 时计算路径逐字同 b25
  （FD 硬闸逐位证实）。未触碰任何方程/求解/既有导出/验收阈值。

## 证据产物（cycle-1/）

- PREREGISTRATION.md、NOTEBOOK.md（§1-§5.2，含 §5.1 全门数字与 §5.2 G-c 隔离）
- EXECUTOR_SUMMARY.md（本文件）、EXECUTOR_DONE（空标记）
- b32_decomp.py / b32_decomp_run.log / b32_decomp.json / b32_decomp.tsv
- 原始运行日志：`/home/ys/dsH/b32_decomp/Log.b32_decomp.txt`（1579.89 s）
- 门控场（36 个 _r1/_r2 场）：`/home/ys/dsH/b32_decomp/1/`

## 结论

B31 §7 的 76.5% 范数残差 = 跨轮目录混合伪影，本轮用同轮门控导出逐位证实
（closure_z vs raw_z 0%，r1==bPD_w、r2==ADJ_gDP），**证无真病灶**。门 2/门 4
的 r2 残余偏差（5.7%/3.2%、12.6%）为协调人重建通道与生产的实现细节差
（r1 的更大偏差已被证为跨轮 Uc 混合），按预注册判据不影响门 3 结论。
G-c 链数值不采信（near-singular 条件伪影，隔离证据齐备）。管线状态交回协调人。
