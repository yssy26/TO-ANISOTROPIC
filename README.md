# TO-ANISOTROPIC

OpenFOAM-7 实现的分层多目标拓扑优化求解器：各向异性热-流-固体耦合，frozen-RANS 冻结湍流框架，reduced discrete flow adjoint 离散伴随灵敏度，MMA 优化器。面向共轭散热结构（导热-对流联合）的形状/材料分布优化。

## 目录结构

```
src/
├── MTO_HF.C              # 主程序（优化循环 + 阶段验证编排）
├── NS.H                  # 不可压 RANS 流动求解（SIMPLE，frozen 湍流场）
├── HeatTransfer.H        # 各向异性热传导求解
├── AdjHeatTransfer.H     # 热伴随求解
├── AdjNS_PD.H / AdjNS_HT.H / AdjNS_FF.H   # 各目标伴随（压降/热耦合）
├── solveDiscreteFlowAdjoint.H  # reduced 离散流动伴随（J/J^T 算子、GMRES、CSR 导出）
├── computeObjective.H / costfunction.H    # 目标与约束
├── sensitivity.H         # 灵敏度组装（伴随投影）
├── MMA/                  # MMA 优化器
├── adjointInletPressure* / adjointOutlet* # 伴随边界条件
├── validate*.H           # 阶段验证（B1 可重复性 / B2-B3 梯度幅值 / Gate6 / 链式法则）
├── tests/                # 单元/分区验证
└── Make/                 # wmake 构建
```

## 核心流程

1. **流动**：`NS.H` 求解稳态不可压 RANS（SIMPLE），冻结 k-omega 湍流场（`nuEffFrozen`）。
2. **热**：`HeatTransfer.H` 求解各向异性热传导（含扩散系数对角化与有效张量合成）。
3. **灵敏度**：`sensitivity.H` 组装 pressureDrop / thermalCoupling 等目标与体积约束的梯度（frozen-RANS 离散伴随）。
4. **优化**：`MMA` 更新设计变量（过滤-链式法则 `filter_x.H` / `filter_chainrule.H`）。

## 离散伴随（reduced deltaPhi 路径）

`solveDiscreteFlowAdjoint.H` 实现正向算子 `applyDiscreteFlowJ` 与转置算子 `applyDiscreteFlowJT`（reduced 连续性路径：`δU → δH → δphi → div(δphi)`），配套：

- 同基态重装配 Gate 验证（H2 系数 / H3-H4 offdiag / FD-H oracle / P 行多 h reduced FD）；
- J/JT 块级与全向量 dot test（12 组随机向量，max 相对误差 1.1e-13）；
- 显式 CSR 组装与 `explicitJT.mtx` 导出（vs matrix-free 4.0e-16）；
- GMRES + ILU0 预条件求解伴随（当前限制：求解层数值未收敛，见下）。

## 构建

依赖 OpenFOAM-7（`/opt/openfoam7`），在 `src/` 下：

```bash
source /opt/openfoam7/etc/bashrc
wmake
```

## 运行

case 目录需包含 `constant/optProperties`（目标、缩放、求解器开关）、`system/`、`0/`（初值场）与 frozen 湍流场输入。主程序：

```bash
MTO_HF -case <caseDir>
```

关键开关（`constant/optProperties`）：

| 开关 | 作用 |
|---|---|
| `adjointMode discrete` | 使用 reduced 离散伴随 |
| `stageB4JacobianProbe true` | B4 Jacobian 诊断（Gate 输出，仅诊断） |
| `discreteExportOnly true` | 仅导出显式矩阵/右端，跳过 GMRES |
| `stageB2Enabled true` | B2/B3 梯度幅值验证（FD vs adjoint） |

## 验证链

- **Stage A**：各向异性实现的分区/单元级验证。
- **Stage B1**：frozen 湍流可重复性与噪声底。
- **Stage B2/B3**：D1/D2/D3 梯度幅值（FD vs TLM vs ADJ，sign + best-step 门槛）。
- **Stage B4**：reduced 离散 Jacobian 一致性（Gate H1-H4、P 行 FD、dot test、explicitJT oracle）。

## 已知限制

- `solveDiscreteFlowAdjoint` 仅支持串行（processor-patch 转置项未实现）。
- 伴随 GMRES 求解层尚未收敛（Krylov 向量数值爆炸，疑似 ILU0 预条件问题）；矩阵本身已通过 dot test 与显式 oracle 验证至机器精度，待修复预条件/缩放后解锁 D1 验收。
- 当前阶段未启用 MMA 更新与 objective J 梯度修复。

## 说明文档

- `src/SRC_CODE_GUIDE.md` — 源码级阅读指南
- `src/ANISOTROPIC_CONDUCTIVITY_PLAN.md` / `src/ANISOTROPIC_STAGE_A_VALIDATION.md` — 各向异性实现与验证说明
