# BFINAL-032 NOTEBOOK — 分解场门控导出轮

> 轮次目标（任务书 NEXT_TASK_BFINAL032.md）：一次带开关门控导出的重跑，让生产自己写出
> 全部分解场（同轮次、同点），终结 B31 §7 的跨轮伪影（76.5% 残差之谜）。三阶段：
> 门控写出实现（默认 off，一次编译）→ ccache 编译 + b32_decomp 同构运行（FD 三列
> 逐位一致硬闸）→ 同轮场直接分解。
> 执行者：executor-nonworkingtime（DeepSeek v4 Flash）。分支 `agent/dsH-stage-b-validation`
> @ 4999855（本轮开头）。

## §1 轮结构判定（运行前，log 实证）

b25 运行恰好 **2 次** `sensitivity.H` 调用：

| 调用 | log 行 | 出处 | 包裹标记 | momentum L2 |
|---|---|---|---|---|
| 轮 1 | 4784-4801 | MTO_HF.C L110 主循环 | "MTO_HF: Starting/Finished sensitivity.H" | 0.00351314475249 |
| 轮 2 | 10408+ | validateStageB2GradientAmplitude.H L460 `runFullB2Chain` lambda | 裸 "sensitivity analysis" | 0.00226699400158 |

（第三潜在点 validateStageBRepeatability.H L455 被 `stageBEnabled=false` 门掉，b25 未触发。）

**混合机制（B31 §7 根因确认）**：
- `writeOptimizationState.H`（MTO_HF.C L107）`runTime.write()` 写 `/1/Uc /1/Ub /1/pc /1/pb /1/Tb`
  （AUTO_WRITE）= **轮 1 值**；
- 轮 2 sensitivity.H L530-563 写块（finalSolvedIteration 真）覆写 on-disk 灵敏度族 =
  **轮 2 值**。

因此 `/1/` 目录的伴随场（轮 1）与 on-disk 灵敏度族（轮 2）**非同轮** —— 这正是
B31 §7 离线 mom_J（/1/Ub）匹配 r1 打印而非 r2、gsenshPD 76.5% 残差的机制。

## §2 钉死基线数字（b25 门轮）

- ADJ_J（projJ）：D1=-0.0358350400485765，D2=-0.02492599038981423，D3=-0.01200301226991504
- ADJ_gDP（projDP）：D1=-5.23029384626718，D2=+1.0505249811003，D3=+13.60342258150915
- ADJ_gV（projV）：D1=0.1555186511144686，D2=-0.08449667741658891，D3=0.06396237684432686
- FD 侧见 `stage_b2_fd_scan.tsv`（FD_J/FD_gDP/FD_gV），逐位对照基准 = b25 文件字节级。
- 关键 log 锚：轮 2 前有 `[FGMRES-PROD pressureDrop] iter=1102`、`[GRADSTABLE pressureDrop]
  cycles=14`；thermalC L2 两轮相同 = 0.00111307150907。

## §3 实现设计（阶段 1）

`sensitivity.H` / `rxPressureRowTranspose.H` 新增 optProperties 开关
`stageB31WriteDecomposition`（默认 false）。true 时：
- 生产两段 + J 四段分解场 `.write()`：gsenshPressureDropMomentum / gsenshPressureDropPressureRow /
  gsenshMeanTMomentum / gsenshMeanTPressureRow / gsenshMeanTFluxDirect / gsenshMeanTThermalC
  （thermalC = thermalDiffusionDerivativeDTCell 包装场）；
- T 环与中间量：rxPressureRowT（pc 基）/ rxPressureRowTb（pb 基）/ rxFluxDirectT（Gx 基）/
  rxrAU / rxDHbyA / rxG0Field（内部面段，仅在 pc 首次包含时写）；
- 轮次标注：optProperties 运行时计数（两包含点共享同一 IOdictionary），后缀 `_r1`/`_r2` 写
  Uc/Ub/pc/pb/Tb；log 打印明示 `/1/` 伴随场属轮 1、on-disk 灵敏度族属当前轮。

**零污染保证**：全部为纯 I/O `.write()` + Info 打印，不触碰任何方程/求解/既有导出。
默认 false 时计算路径逐字同 b25。

## §4 阶段 2 进度

**阶段 1 实现完成（2026-08-22，两次 Edit 于 sensitivity.H / rxPressureRowTranspose.H）**：

`sensitivity.H` 新增（L1-26 顶部区）：
- optProperties 持久计数 `b32SensInvocationCount`（避免 per-include-site static 陷阱），
  调用号 `_r1`/`_r2` 由 `b32RoundTag = word("_r") + Foam::name(b32CallNumber)` 生成；
- 开关 `stageB31WriteDecomposition`（lookupOrDefault<Switch>，默认 false）；
- log 明示声明：`/1/` Uc Ub pc pb Tb 为轮 1 值（writeOptimizationState.H 所写），
  本调用写出的灵敏度族属轮 `b32CallNumber`。

`sensitivity.H` BFINAL-019 打印块后（freeze 块前）新增门控块 `if (b32WriteDecomposition)`：
- gDP 两段 + 原始和：gsenshPressureDropMomentum_r? / _PressureRow_r? / _Raw_r?（拷贝带名写入，
  捕获 filter_chainrule.H 之前的真 pre-filter 值）；
- T 环基：rxPressureRowT_r?（pc）/ rxPressureRowTb_r?（pb）/ rxFluxDirectT_r?（Gx）；
- J 四段：gsenshMeanTMomentum_r? / _TPressureRow_r? / _TFluxDirect_r? /
  gsenshMeanTThermalC_r?（thermalDiffusionDerivativeDTCell 包装场）；
- 轮次伴随场拷贝：Uc_r? / Ub_r? / pc_r? / pb_r? / Tb_r?（GeometricField copy-with-IOobject ctor，
  NO_READ 不读盘，捕获本调用内存值）。

`rxPressureRowTranspose.H` 生产赋值后新增门控（`#ifndef RX_FLUX_DIRECT_SOURCE` 守卫 = 仅 pc 首次包含）：
- rxrAU_r? / rxG0Field_r?（tmp 构造场，带名拷贝后 write）；
- rxDHbyA_r?（raw vectorField 包装为 volVectorField 后 write）。

**污染保证**：全部为带名拷贝 `.write()` + Info 打印；不触碰任何方程/求解/既有导出路径。
默认 false 时计算路径逐字同 b25（静态测试 15/15 OK 于编译前验证）。

**编译**：`wmake` 后台（ccache-shim），预计 ~40 min（MTO_HF.C 单 TU）。

**Phase 2 准备（2026-08-22，编译等待期间完成）**：
- b32_decomp 克隆就绪：`rsync -a --exclude=1/ --exclude=Log* --exclude=optimization_* --exclude=adjointCheckpoint_* --exclude=b18rhs_*.mtx --exclude=stageB2/ b25_qgate/ b32_decomp/`；
  diff 校验：`system/` 与 `0/` 逐字节一致（diff exit 0），`constant/` 仅 optProperties 新增
  `stageB31WriteDecomposition true;`（5 行，含注释）—— b25 同构成立。
- FD 门控基线 sha256（b25 原盘）：
  - stage_b2_summary.tsv = `90cef08a2c09e8bfbef0a09f3eeb9e8fb35d939b22eda69ed42ce861fa27c7c0`
  - stage_b2_fd_scan.tsv = `2da03a58c6f32745f7a9ef7cf4af593ff7b36c2ecaf2a2de793462bf5eb19bd6`
- b25 固定参考数：projJ=(−0.0358350400485765, −0.02492599038981423, −0.01200301226991504);
  projDP=(−5.23029384626718, +1.0505249811003, +13.60342258150915);
  projV=(0.1555186511144686, −0.08449667741658891, 0.06396237684432686)。
- b25 配置复核：controlDict endTime=1, deltaT=1, writeControl timeStep, writeFormat ascii,
  writePrecision 12; optProperties stageB2Enabled true, mmaUpdateEnabled false,
  freezeColdFlowForValidation false, discreteFlowAdjointMaxIter 4000。

**运行计划**：编译完成后 `bash run_bfinal025.sh /home/ys/dsH/b32_decomp Log.b32_decomp.txt`
→ 完成后立即 FD 三列（J/gDP/gV）逐位对比（sha256 对照上述两 tsv）→ 硬门限。

## §5 阶段 3 进度

（数据到手后立即分解，勿拖延 —— B26 教训）

### §5.0 运行完成 + FD 硬闸判定（2026-08-22 23:2x）

- b32_decomp 运行正常结束：`ps` 无 MTO_HF 进程；Log.b32_decomp.txt 尾部
  `--- Stage B2/B3 complete (state restored) ---`、`ExecutionTime = 1579.89 s`。
- 18 个 _r1 + 18 个 _r2 门控场全部落盘 `/home/ys/dsH/b32_decomp/1/`（详见运行前清单）。
- **FD 硬闸 PASSED（逐位一致）**：`cmp` 字节级验证，且 sha256 复核：

| 文件 | b32 sha256 | b25 sha256 | 结论 |
|---|---|---|---|
| stage_b2_summary.tsv | `90cef08a2c09e8bfbef0a09f3eeb9e8fb35d939b22eda69ed42ce861fa27c7c0` | 同左 | 逐位一致 |
| stage_b2_fd_scan.tsv | `2da03a58c6f32745f7a9ef7cf4af593ff7b36c2ecaf2a2de793462bf5eb19bd6` | 同左 | 逐位一致 |

→ 导出实现未污染计算路径，硬闸清空，进入同轮直接分解（本行）。

### §5.1 同轮分解结果（2026-08-22 23:4x，b32_decomp.py 重跑入环）

全部门数字来自 `b32_decomp.py`（2026-08-22 23:45 版，含 G-c 隔离注记）输出
`b32_decomp_run.log` / `b32_decomp.json` / `b32_decomp.tsv`。

**G0 跨轮归因（逐位）**：`/1/Uc /1/Ub /1/pc /1/pb /1/Tb` == `_r1` 拷贝 relL2=0.0
（5/5）—— `/1/` 伴随场 = 轮 1 值，与 §1 混合机制判定一致。

**G1 同轮闭合门（逐位级）**：
- `_r1`：raw == mom+prow relL2=1.021e-12，maxAbs=1.0e-14；
- `_r2`：relL2=9.916e-13，maxAbs=1.0e-14；
- 离线成分复核：`_r2` mom=-dAlphaDxh(U·Uc_r2)V relL2=3.1e-12、prow=-rxT·dAlphaDxh
  relL2=2.7e-12（`_r1` mom relL2=7.8e-2 仍 = 旧 /1/Uc 跨轮混合，r2 同轮逐位）。

**G2 T 环镜像（生产 rxPressureRowT vs 离线 T_row，2% 目标）**：

| 轮 | 标签 | relL2 | corr | 判定 |
|---|---|---|---|---|
| r1 | pc | 2.09e-1 | 0.990996 | FAIL |
| r1 | pb | 2.58e-1 | 0.994459 | FAIL |
| r2 | pc | 5.66e-2 | 0.998235 | FAIL（<2% 未达） |
| r2 | pb | 3.19e-2 | 0.999508 | FAIL（<2% 未达） |

r2 显著优于 r1（5.7% vs 21%，prow 重建 relL2 2.7e-12 逐位）——r1 的差是跨轮 Uc 混合，
r2 剩余差为协调人重建通道（dHbyA/drAU）与生产的实现细节差。预注册判据：门 2 失败不影响
门 3 有效（mom/prow 两段同轮闭合用生产自己的 T）。

**G4 HbyA 校准**：
- `_r2`：rxrAU vs mob relL2=1.60e-12 corr=1.0（逐位）；rxG0Field vs g0 relL2=6.26e-13
  （逐位）；rxDHbyA vs dHbyA relL2=1.26e-1 corr=0.9916（偏差 12.6%，前轮先例 2%）。
- `_r1` 对应 1.97e-1 / 8.56e-2 / 2.16e-1 —— 与 G2 同因（跨轮 Uc）。

**G3 跨轮伪影消解（本轮核心判定）**：

| D | 轮 | mom_z | prow_z | closure_z | raw_z | vs bPD_w / ADJ_gDP |
|---|---|---|---|---|---|---|
| D1 | r1 | -4.190491 | -0.086971 | **-4.2774612375491** | -4.2774612375496 | == bPD_w（残 1.25e-13） |
| D1 | r2 | -5.155654 | -0.074640 | **-5.2302938459220** | -5.2302938459240 | == ADJ_gDP（残 3.87e-13） |
| D2 | r1 | +0.366399 | +0.018064 | **+0.384462806083** | 同 | == bPD_w（残 2.2e-12） |
| D2 | r2 | +1.026292 | +0.024233 | **+1.050524977639** | 同 | == ADJ_gDP（残 4.9e-13） |
| D3 | r1 | +14.245498 | -0.149145 | **+14.09635332868** | 同 | == bPD_w（残 3.4e-14） |
| D3 | r2 | +13.832100 | -0.228678 | **+13.60342257618** | 同 | == ADJ_gDP（残 1.2e-13） |

**判定（FINAL）**：同轮 closure_z 与 raw_z 逐位闭合（全部 relL2 ≤ 2.2e-12，target ≤ 5%）。
r1 闭合 = bPD_w（主线期望约束）、r2 闭合 = ADJ_gDP（B31 验证标签）—— B31 §7 的 76.5%
范数残差 **100% 消解为 0%（逐位）**。残余 ADJ−bPD_w = (D1 −0.952833, D2 +0.666062,
D3 −0.492930) 为干净的 r2-vs-r1 轮状态差（两轮 lambda 不同），**非病灶**。D1 跨轮对照：
B31 §7 的 mom=−4.253786/prod=−16.116120 与同轮 r1 mom=−4.190491 / r2 mom=−5.155654
均不同——旧数字即混合伪影本身。

### §5.2 G-c 链隔离判定（离线 chain_field 无效 → 锚 raw_z）

`b32_decomp.py` 的 `chain_field` 离线 (Lm+VI) 求解与 on-disk 链后场不符
（G-c2 relL2=8.67e-1 corr=0.468；G-c1 relL2=9.98e-1），链上投影 ch_total/ch_mom/ch_prow
标记 **INVALID**（JSON `ch_valid=false`）。隔离证据（/tmp 调查，2026-08-22）：

- **源逐位**：`chain_src(raw_gDP_r2)`（链反向作用在未链源上）vs on-disk gsenshPD
  relL2=2.911e-10（which_round.py）—— 链输入侧逐位；
- **矩阵逐位**：离线 Lm 复现**正向**滤波器逐位（scale_test.py：solve(Lm+VI, Vx) 的
  |X1|=1.4199e+02 == 生产 xp 范数，corr=+0.999870）；
- **near-singular 条件伪影**：(Lm+VI) 有常数零模态（全 zeroGradient BC → Lm 常数模态
  特征值 0；V~1.25e-10 对角 vs Lm 尺度 ~1e-4 → 零模态放大 1/V ~ 8e9）。生产 DICPCG
  （Log.b32_decomp.txt L10420-10423：init=0.0879 final=2.07e-10，10 iter）与离线 splu
  在零模态方向分叉；磁盘解 |disk|=5.17e-1 ≈ 零模态尺度 |V·g|/V=6.3e-1，disk 与其源
  corr=0.9796（输出≈光滑输入），(Lm@disk)[dm] max=4.7e-6（disk 是近零向量）。
- **结论**：G3 判定锚定 raw_z/closure_z（直接门控导出，逐位），不锚链投影。

## §4b 编译进展日志（2026-08-22）

| 时刻 | 状态 |
|---|---|
| 22:10 | wmake 后台启动（ccache-shim），wclean 先行 |
| 22:20 | 9/10 .o 完成（harmonicSymmTensor, 8 adjoint* + MMA）；MTO_HF.C cc1plus 进行中 |
| 22:27 | MTO_HF.C cc1plus CPU 17:13（108%），无 error:，无 WMEXIT 哨兵 |
| - | 阶段 3 脚本 b32_decomp.py 已就绪（复用 b31b_m2.py 骨架，读 round-tagged 场，G0-G4 门） |
