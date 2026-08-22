# TO-ANISOTROPIC 代码状态报告（供外部专家审阅）

> 日期：2026-08-23 ｜ 分支：`agent/dsH-stage-b-validation` ｜ 本文档随分支一并推送。
> 目的：请专家独立审查当前梯度链实现与我们的诊断结论。文档力求自包含；所有断言的证据
> 位置均在 `evidence/agent-group/` 下可复核。

## 1. 一句话现状

OpenFOAM-7 两流湍流共轭传热拓扑优化求解器：**求解层与装配层已机器精度级验证干净，
但端到端 FD 梯度门仍失败**（热目标 J 的 D1 方向反号、D2/D3 幅值超配 ~6.7×/28×；
压降约束 gDP 符号全对但均匀 2.05× 超配）。诊断刚刚经历一次重大重置（§6），当前头号
假设是"双伴随轮次的线性化态与 FD 基态一致性"问题。

## 2. 物理与优化问题

- 冷/热双流 + 中间固体，冷侧设计域做密度法拓扑优化（Brinkman 阻化 + 各向异性导热插值）；
- 热侧流场冻结（U_h/p_h/k-ω 冻结，T_h 活）；冷侧每外层态全 SST-RANS，敏感度阶段冻结湍流；
- 目标：换热总量 Q（新增，`maximizeTotalHeatTransfer`，J_Q=−Q/(Tref·M_frozen)）或冷侧出口
  混合温度（旧型）；约束：冷侧压降、固相体积分数；优化器 MMA（门控锁定中）；
- 链条：x → Helmholtz 滤波 → 体积保持 Heaviside 投影 → xh → α(xh)/DT(xh) →
  冻结湍流 SIMPLE 原始解 → 伴随（热标签 λ_T + 生产离散流伴随 λ_TC/λ_PD）→ 装配 → FD 门。

## 3. 构建/运行/复现

```bash
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
export PATH=/home/ys/dsH/tools/ccache-shim:$PATH   # ccache（可选）
unset FOAM_SIGFPE    # 勿用 set -e/-u、勿用 tee
cd src && wclean && wmake   # 单一大 TU ~40 分钟；ccache 命中时 1s
```
验证算例（外部路径，不在仓库）：`/home/ys/dsH/b25_qgate`（Q 目标 FD 门）、
`/home/ys/dsH/b32_decomp`（+同轮分解场写出）。运行 ~26 分钟，产出
`stageB2/stage_b2_{fd_scan,summary}.tsv`。详见 `OPENFOAM_USAGE_AND_FILES.md` 与
`ORIENTATION_FOR_EXECUTORS.md`。

## 4. 代码架构地图（src/）

| 文件 | 职责 | 状态 |
|---|---|---|
| `MTO_HF.C` + 49 个 `.H` 碎片 | 主程序（单 TU） | — |
| `NS.H` | 原始 SIMPLE/RANS + 冻结湍流 | 验证干净 |
| `HeatTransfer.H` / `AdjHeatTransfer.H` | 热前向 / 热伴随 λ_T(=Tb) | λ_T 已机器精度外部锚定（B27） |
| `solveDiscreteFlowAdjointProduction.H` | **生产离散流伴随**（J^T、b_TC/b_PD 折叠、H7、pressureGAMG） | 求解层验证干净；折叠=源收缩已机器级证（B32） |
| `solveDiscreteFlowAdjoint.H` | 诊断 oracle（显式 CSR/COO、点测试） | 干净 |
| `sensitivity.H` | 梯度装配（动量/压力行/通量直接/thermalC 分解 + B19 打印） | **同轮闭合逐位**（B32） |
| `rxPressureRowTranspose.H` | 压力行转置 T=J_P^T λ（BFINAL-005） | 机器级精确（B32） |
| `costfunction.H` / `computeObjective.H` | 目标值与 dJ/dT、dJ/dphi（三类型） | Q 型新加，旧型逐位不变（B25） |
| `filter_x.H` / `filter_chainrule.H` | 滤波/投影与前向链/转置链 | 转置链的 eta 修正非对称性已推导（B31 §7） |
| `validateStageB2GradientAmplitude.H` | **FD 门**（17 列 tsv、三判据） | 口径的最终裁判 |
| `validateMmaUnlockGate.H` | MMA 硬门 | 锁定（勿绕） |

## 5. 验证状态总表

**已机器精度级验证干净（证据轮）**：
1. 生产伴随求解（pressureGAMG、双标签 1e-12 + GRADSTABLE）（B10/011/018）；
2. H7 压力通道修复：phiHbyA 含 −interp(rAtU·∇p)·Sf，四处一致（B23 推导/B24 实现）；
3. λ_T 热伴随 vs 直接解：relL2 1.15e-12、零符号翻转（B27）；
4. thermalC 收缩重算：比值全 1.000000（B27）；
5. **装配链 = 源收缩**：同轮 gsenshPD=mom+prow 逐位；r1 闭包==b_PD^T w_true 1.25e-13（B32）；
6. rxPressureRowT（T 环）与 rxrAU/rxG0 逐位（B32）；T(pb) 镜像 1.6%（B31/B19 打印）；
7. Q 目标实现（B25，旧型逐位回归）；gV 体积链（1e-7 级 FD）；静态测试 15/15。

**仍失败（同轮自洽数据，B25 tsv 为最终口径）**：
| 量 | D1 | D2 | D3 |
|---|---|---|---|
| J：FD / ADJ | **+0.0425 / −0.0358（反号）** | −0.0036 / −0.0249（6.9×） | −0.00044 / −0.0120（27×） |
| gDP：ADJ/FD 因子 | 2.06 | 1.96 | 2.10（均匀 ~2.05×，符号 3/3 对） |
| gV | 机器级 ✓ | ✓ | ✓ |

**真缺项**：Gx=g^T·φ_x 生产装配层缺失（B28 参照对照确认；对 D1 裁决影响 <180% 鲁棒界内）。

## 6. 关键上下文：双轮结构与 2026-08-23 的战略重置

B2 模块每轮验证跑**两个伴随轮**：r1（主循环态）与 r2（全 SST 再基线态，B18 已知两轮
线性化态不同 J=−1.146/−1.167）。最终梯度取 r2。**B32 发现 `/1/` 场目录混合两轮的场**
（未标注的 Uc/Ub/pc/pb/Tb 与 _r1 副本逐位相同）——因此：
- B26 的 M1 源层恒等式判词（"D1 稳健反号定罪折叠层"）、B29 的槽位排除、B30 的
  T3 三指纹，**凡跨用 `/1/` 场与轮特定值者均降级为疑似污染**；
- B30 的 defcheck 缺口（0.82/0.37/1.04）已证为纯 r1/r2 轮态差（B32）；
- 新头号假设：**线性化态与 FD 基态的一致性**（状态记账类问题）。注意纯态差解释不了
  gDP 的均匀 2.05×（r1 轮因子非均匀 1.68/0.72/2.18 vs r2 轮均匀）——需与之共存的
  均匀尺度机制尚无候选。
- b32_decomp 算例已产出 36 个**同轮**标注场（`*_r1/*_r2`），B33（轮清洁 J 恒等式 +
  态一致性裁决，纯离线）为下一步。

## 7. 证据索引（`evidence/agent-group/`，每轮含 FINAL/EXECUTOR_SUMMARY + REVIEW_BY_REVIEWER）

里程碑：B3(J_PU)/B5(R_x)/B8(J_PP) → B10/011(GAMG+生产收敛) → B12(FD 门首次正式失败)
→ B18(容差/判据+折叠算子一致化) → B23(三仪器修正+H7 发现) → B24(H7 落地，gDP 2.13→2.06)
→ B25(Q 化否证泛函结构假设) → B26(M1 源恒等式首测) → B27(λ_T/thermalC 双清白)
→ B28(参照系对照：g 层洗清) → B29/B30(槽位/T3，证据后被降级) → B31(跨轮伪影发现)
→ B32(同轮分解：装配链无罪 + 战略重置)。全局：`DEFECT_MAP_2026-08-22.md`（注意其
B26-B30 部分已被 B32 重置部分推翻）、`STRATEGIC_REASSESSMENT_2026-08-22.md`。

## 8. 安全不变量（审阅时请勿建议绕过）

MMA 锁定（`frozenGradientValidated=false`）；`validateMmaUnlockGate.H` 无方向豁免；
锁定里程碑 B3/B5/B8(+H7/B24 升级)未有新证据不得改写；诊断与生产分离；零拟合因子。

## 9. 想请教专家的三个问题

1. 在装配层/源层机器级干净、λ_T 干净的前提下，您对 **J 的 D1 反号 + 方向依赖超配**
   与 **gDP 均匀 2.05×** 的并存模式有何独立解释？（我们当前最优假设：线性化态一致性
   + 未知均匀尺度机制）
2. 双伴随轮（主循环态 vs 全 SST 再基线态）取哪个作为 FD 门的线性化态才是良定义的？
   两个基态的差（B32: Δ=(−0.95,+0.67,−0.49)）是否应在门内显式对账？
3. 您是否见过 relaxed-SIMPLE 离散伴随里类似的**均匀因子**现象（文献或经验）？
   例如 αRel/矩阵 relax 语义在 H() 通道 vs 通量通道的不一致放大。

---
*审阅环境提示：全部历史证据、日志、逐位对照表均在仓库内；FD 真值在
`evidence/agent-group/BFINAL-025/cycle-1/stage_b2_fd_scan.tsv`；同轮分解场清单见
`BFINAL-032/cycle-1/EXECUTOR_SUMMARY.md`。*
