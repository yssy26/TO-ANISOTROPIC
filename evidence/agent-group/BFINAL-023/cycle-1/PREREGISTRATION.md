# BFINAL-023 cycle-1 — PREREGISTRATION（先于本轮任何新数据落盘）

日期：2026-08-20。分支 `agent/dsH-stage-b-validation`，起点 `c920212`。
任务：SLOT-5 专项轮（算子通量切线 vs 真实系统响应）。零生产代码改动；
全部检验离线（既有 B 态导出 + b16 网格），无重编译、无求解器运行（除非
判别表明确要求且不影响生产）。

## 0. 本轮开始前已在源码层确认的事实（非本轮数据）

以下为源码阅读结论，构成本轮假设的机理基础（引用行号为 c920212）：

- F1（前向钳制语义）：`src/NS.H:26-45` 组装 `fvm::div(phi,U) −
  fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U)`，随后 `UEqn.relax()`
  （OF7 `fvMatrix.C:521-670`：边界对角并入 → `D = max(|D|, Σ|off|)` 钳制
  → `D /= alphaRel` → 边界对角回减 → 源补偿 `S += (D−D0)·psi`）。
  `rAU = 1/UEqn.A()`（NS.H:113）与 `primalPressureMobility = rAtU`
  （NS.H:129；SIMPLE 无 `consistent` → rAtU = rAU）都取自**钳制+松弛后**
  的对角。
- F2（算子无松弛语义）：`src/solveDiscreteFlowAdjointProduction.H:144-173`
  与 `src/solveDiscreteFlowAdjoint.H:22-81` 组装同一矩阵但**从不调用
  .relax()**；`prodRAU = 1/A_unrelaxed`。P 行通量切线的 U 通道用
  `alphaRel*prodRAU`（Production:357/394/557/594），kf 通道用精确钳制
  `primalPressureMobility`（Production:367-375）。**结构性前提成立：
  U 通道 mobility 因子与松弛源系数在钳制胞上与真前向语义不一致。**
- F3（OF7 装配符号，来自 /opt/openfoam7 源码）：`fvm::div` 上风
  （gaussConvectionScheme.C）：`lower = −w·phi, upper = (1−w)·phi,
  diag = negSumDiag`；`fvm::laplacian`（gaussLaplacianScheme.C
  fvmLaplacianUncorrected）：`upper = lower = +γ_f·|Sf|·δ,
  diag = negSumDiag`（即矩阵代表 −∇·(γ∇)·V，进入动量方程的
  `− fvm::laplacian` 后：**diag += kfv_f, upper/lower −= kfv_f，
  kfv_f = 线性插值 γ_f·|Sf|·δ，上下对称**）。
- F4（离线仪器符号缺陷，本轮新发现，W-B 候选根因）：b16 谱系离线重建
  （b22_bpoint_gauge.py build_basis / b22_mobility_forensics.py rebuild）
  装配 `upper = qn·phi + kfv_n, lower = −qp·phi + kfv_o, diag −= kfv`，
  即**扩散符号相对生产方程（F3）整体翻转，且 γ 用单侧 owner/nei 值而非
  对称插值**。固体区（α·V 支配、kfv 均匀）符号不可见 → 与 B22 观测
  「固体精确 2.6e-5、偏差集中 nut 活跃 + 钳制绑定胞」一致。B22 已测
  S3（仅 γ 插值、符号不变）仅移动 0.4%——与「符号才是主导」不矛盾。

## 1. 预注册假设

### H-CLAMP（审阅者主假设，本轮主检验对象）

真前向通量映射（frozen nut，SIMPLE 非一致）：
`HbyA = rAtU·H`，`H = [lduH(dU) + relaxSrc(dU,dα)]/V`，
`phi = interp(HbyA)·Sf − interp(rAtU)·δ|Sf|·(p_n−p_o)`，
其中 `rAtU = V/max(|D0+bc|, Σ|off|)/alphaRel`（钳制胞取 Σ|off|）。

在**非钳制**胞：`rAtU·A_u = alphaRel` 恒等 → 生产切线
`alphaRel·rAU_u·dH + (1−alphaRel)·dU` 恰为真切线。
在**钳制**胞（sumOff > |D0|）三处失配：
  (a) mobility 因子：真 = `rAtU = alphaRel·V/sumOff`，
      生产 = `alphaRel·V/D0`，高估 `sumOff/|D0| ≥ 1` 倍；
  (b) 松弛源直通系数：真 = `1 − rAtU·A_u = 1 − alphaRel·D0/sumOff > 1−alphaRel`，
      生产 = `1−alphaRel`；
  (c) 设计通道（da）：真 `drAtU = 0`（钳制值 Σ|off| 与设计无关，nut 冻结、
      对流冻结于 B），`d(relaxSrc)/da = −V·da·U`；生产/离线仪器
      `drAU = −mob²·da/alphaRel ≠ 0`、`dHsrc = (1/alphaRel−1)·da·U`。

判别预测（预注册）：若 H-CLAMP 携带 D2 的 45% 面级失配，则把
(a)+(b)+(c) 全部换成真语义的变体 VC 应使 D2 面锚 relL2 大幅回落
（阈值见 §3）；且修正量应集中在钳制绑定面邻域。

### H-WB（W-B 根因，执行者在 F4 基础上提出）

b16 谱系离线基的 12.9% mobility 仪器误差主因是扩散装配符号翻转（F4）。
判别预测：仅修正符号+对称插值 γ_f 后，精确 mobility vs 重建 relL2
从 12.9% 跌破 1%，且非钳制胞恒等式 `mob_exact = alphaRel·V/|D0_corr|`
在非钳制胞集上残差 <1%（双通道自证）。

## 2. 判别实验（顺序固定）

- **T0 复现控制**：逐位复现 B22 仪器读数（mobility 12.9%；D1/D2/D3
  面锚 8.68%/44.9%/8.99%，重建基，h=1e-3）。不符则本轮作废排查仪器。
- **T1（H-WB）**：修正基（生产符号 + 对称 γ_f + 其余同 b22）重建
  mobility；报 relL2、非钳制恒等式残差、钳制分数变化
  （b22 基 14.1% → 修正基 ?）。真钳制胞集合与 A_relaxed/A_unrelaxed
  比值分布由修正基 + mob_exact 双重确定（钳制胞 mob_exact·sumOff/
  (alphaRel·V) ≈ 1 自检）。
- **T2（重合度）**：D2 失配面集合（B22 W5 定义：rebuild 基 h=1e-3，
  top-10% 失配能量内部面，x-法向 |nx|>0.9）与真钳制胞的重合度：
  top 面中 owner/nei 至少一个钳制胞的份额；钳制面 vs 非钳制面的
  relL2 分解；修正基下的 D2 失配面-钳制胞相关。
- **T3（H-CLAMP 变体面锚，主实验）**：同一 dphi_true（既有 FD 导出，
  h=1e-3 主读数、3e-4 桥），五个变体（全部用精确 kf 通道 =
  interp(mob_exact)，与生产 W1 后一致）：
   - V0old：B22 仪器原样（旧基，符号翻转）——复现控制；
   - V0c：B22 形式（mob·dH + 0.6·dU + 旧 da 通道），修正基 + mob_exact；
   - VP：生产镜像（alphaRel·rAU_u(未钳)·dH + 0.6·dU + 旧 da 通道），
     修正基——直接估计生产算子的面级切线误差；
   - VC：H-CLAMP 真语义（mob·dH + (1−mob·A_u)·dU + 钳制 da 通道
     [钳制胞 drAU=0, dHsrc=−da·U]）；
   - VC1/VC2：VC 只修 (b) / 只修 (c)——通道归因（W-C）。
   读数：全体内部面 + D2 top-10% x-法向面子集的 relL2/cos；
   W-C 三通道（dU/dp/da）幅度份额与 VC−V0c 逐通道差在 D2 载体面上
   的分布。
- **T4（裁决，见 §3）**。不满足证实条件即按退出坡道报告，不扩展新
  假设轮（除非 W-C 归因明确指向单一可修槽位且判别代价为零代码）。

边界/杂项语义（b22 ic 模型、outlet 处理）在所有变体间保持一致，
其残余属 S1 级（bnd 邻胞 err 7e-3），不参与本轮归因。

## 3. 预注册判读规则

- **H-WB 证实**：T1 relL2 < 1% 且非钳制恒等式残差 < 1%。
  部分证实：1–3%；否证：> 3%（此时回到 γ 插值/边界项排查，不阻塞
  H-CLAMP 主线——VC 用 mob_exact + 修正基 A_u，对 A_u 残差做 ±15%
  敏感性包络）。
- **H-CLAMP 证实**：VC 使 D2 全场面锚 relL2 相对 V0c 下降 ≥ 1.5×，
  且 D2 载体（top-10% x-法向）面 relL2 下降 ≥ 1.5×，且改善集中
  （钳制邻域面改善 > 非钳制面改善 2×）。
  VC1/VC2 给出通道归因（预注册预测：(b) 松弛源系数为主、(c) 次之、
  (a) 已被 B22 W1 部分吸收）。
- **H-CLAMP 否证**：VC 相对 V0c 改善 < 20%（相对），或改善不集中于
  钳制邻域 → 如实报告 SLOT-5 仍无明确可修机制，执行退出坡道
  （项目转用户战略决策），附 W-C 通道分解与 E4 式逐槽审计表。
- 无拟合因子：所有系数（alphaRel=0.4、mob_exact、修正基）来自导出
  数据与 OF7 源码语义，不调参。

## 4. 退出坡道

若 H-CLAMP 否证且 W-C 无单一通道携带 >30% 失配能量：本轮停止，报告
「SLOT-5 无明确可修机制」+ 完整判别表；不进入实现轮，不提出新假设
（除非用户指示）。

## 5. 数据与产物

输入：`/home/ys/dsH/b22_gauge/stageB2/wstate_*`（B22 导出，FD 真值 +
baseline + mob_exact）、`/home/ys/dsH/b16_mesh/*`、
`/home/ys/dsH/b16_states/1/{U,p,phi,alpha,dAlphaDxh}`、
`/home/ys/dsH/b15_export/stageB6_rxc_z_analytic.mtx`。
产物：`b23_hclamp_gauge.py`、`b23_hclamp_gauge.log/json`、本文件、
判别表 `VERDICT_TABLE.md`、EXECUTOR_SUMMARY.md、FINAL_REPORT.md、
空 EXECUTOR_DONE。先归档再清理；git 选择性提交，不 push。
