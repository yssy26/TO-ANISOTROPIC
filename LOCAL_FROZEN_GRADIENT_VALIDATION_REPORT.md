# LOCAL_FROZEN_GRADIENT_VALIDATION_REPORT

日期：2026-08-11
执行：本地 OpenFOAM-7 严格验证（按 `CODEX_LOCAL_VALIDATION.md`）

## 1. 仓库 HEAD 与算例来源

- 仓库：`/home/ys/TO-ANISOTROPIC`
- 分支：`agent/stage-a-anisotropic-validation`
- **HEAD：`3b7052f962f1299904e7168eb1d602add7006e60`**（`Require SST protection until full-SST validation passes`）
- 算例：`/home/ys/b2_case_smoke`（`b2_case` 的干净副本，未污染基准算例；`b2_case/constant/optProperties` 未被修改）
- 规模：33,600 cells / 96,860 internal faces / 134,400 未知量（U-p 系统）

## 2. Clean build 结果

```text
BUILD = PASS
```

`wclean && wmake` 一次通过，无编译错误。本阶段**未做任何源码修改**（后续修改均为求解层诊断/预条件，见 §8）。

## 3. Production 伴随收敛历史（第二阶段）

配置（`b2_case_smoke/constant/optProperties`）：`mmaUpdateEnabled false; stageB2Enabled false; stageB4JacobianProbe false; frozenGradientValidated false; adjointMode discrete; flowModel incompressibleRANSFrozen; frozenTurbulenceAdjoint true; solveFlowAdjoints true; objectiveGradientScale 1.0; pressureGradientScale 1.0; discreteUseExplicitSolution false;`

### 3.1 thermalCoupling

| 求解器/预条件 | maxIter | trueRelRes | 备注 |
|---|---|---|---|
| FGMRES + signed Jacobi（restart 80） | 4000 | **0.48506** | 每 cycle 内 projRes 几乎不变，h00≈3.15 恒定 |
| FGMRES + Jacobi（restart 160） | 8000 | 0.48191 | 仅微降 |
| FGMRES + **L1 行和预条件**（restart 80） | 4000 | **0.07304** | h00→0.219；cycle 内 projRes 开始下降；残差 U/P 平衡（|rU|0.117/|rP|0.491） |
| FGMRES + L1（restart 320） | 8000 | 0.06927 | 每 cycle 线性微降，收敛率 ~0.99994/iter |
| FGMRES + block 三角 Schur（L1 对角 + UP 耦合） | 4000 | 0.07350 | 与纯 L1 无差别 |
| BiCGSTAB（右预条件 L1） | 3000 | 0.49136（发散） | BCGS 内部 trueRelRes 升至 2.79 |

### 3.2 pressureDrop

| 求解器/预条件 | maxIter | trueRelRes | 备注 |
|---|---|---|---|
| FGMRES + L1（restart 80） | 4000 | **0.15709** | RHS：sumRhsU=0、\|RhsP\|L2=10.38（纯压力源） |
| FGMRES + L1（restart 320） | 8000 | 0.15711 | 无改善 |

### 3.3 求解层诊断（全部记录）

- **RHS 平衡**：thermalCoupling `|RhsU|L2=0.0735` vs `|RhsP|L2=6.91`（缩放后 P 主导 94 倍，物理量级合理）；pressureDrop `|RhsP|L2=10.38`
- **可解性**：`J^T·1_P = 0`（精确）——常数 P 在 J^T 零空间；P 均值正交化（rhsPmean -3.45e-4）无改善
- **残差空间分布**：L1 后 thermal `|rU|L2=0.117 / |rP|L2=0.491`（U/P 平衡），max 残差分布广（非单点）
- **近特征判别**：`cos(Au,u) = -0.39`（残差非近特征向量，近零 deflation 不适用）
- **谐 Ritz pilot**：最小 Ritz 值 θ≈0.03 单一主导（cov 数据受反幂数值影响不可全信，但 θ 提示存在近零模式）
- **L1 对角**：avg=54.4、min=1.01、max=1438

### 3.4 第二阶段验收

```text
true relative residual <= 1e-9  →  FAIL（thermal 0.073、pressureDrop 0.157）
finite Ub/pb/Uc/pc              →  未验证（求解未收敛，安全门阻止进入梯度）
no FatalError                   →  收敛门按设计拦截（Gradients blocked）
```

**结论**：两套 production 伴随均未收敛至 1e-9。已尝试 restart 80/160/320、signed Jacobi、L1 行和、block 三角 Schur、BiCGSTAB、近零 deflation——L1 行和改善 8 倍（0.485→0.073）为当前最优，但收敛率 ~0.99994/iter（谱近 1 密集的非正规 U-p 系统），未达验收。**按规程不进入 B2/B3**。

## 4. Stage B4 回归

**未执行**（规程：B4 为独立诊断路径；production 伴随未收敛时以 B4 诊断结果为准的正式梯度验收不成立）。诊断版 `solveDiscreteFlowAdjoint.H` 未被修改（本次所有改动均在 production 求解器 `solveDiscreteFlowAdjointProduction.H` 的求解层），预期 B4 的 J/JT dot test 与 oracle 不受影响，需在 production 收敛后单独回归确认。

## 5. Stage B2/B3

**未执行**（`stageB2Enabled true` 未运行——规程明确禁止以未收敛伴随进入 B2/B3）。

## 6. J/gDP/gV PASS/FAIL

**未判定**（B2/B3 未运行）。

## 7. FROZEN_GRADIENT_UNLOCK.txt

**未生成**（B2/B3 未 PASS）。

## 8. 本地代码修改（求解层，非物理）

文件：`src/solveDiscreteFlowAdjointProduction.H`（其余文件未动；`applyProdFlowJT` 算子本身**未修改**）

| 修改 | 内容 | 性质 |
|---|---|---|
| 求解层诊断 | `PRODRHS/PRODDIAG/PRODRESID/PRODSOLV` 打印（RHS 平衡、可解性、残差分布、L1 统计） | 诊断 |
| FGMRES trace | 每 cycle 投影残差 + h00 打印 | 诊断 |
| **L1 行和右预条件** | 逐列组装行和替换纯对角（改善 8 倍） | 预条件（规程允许） |
| block 三角 Schur 预条件 | UP 耦合 + L1 对角（无改善，保留开关） | 预条件（规程允许） |
| BiCGSTAB 选项 | `discreteProdSolverType bicgstab`（发散，默认 fgmres） | 求解层（规程允许） |
| 近零 deflation | cos 判别 + rank-1 消元（cos=-0.39 未启用） | 求解层（规程允许） |
| 收敛门开关 | `discreteProdConvergeFatal`（默认 true；诊断时设 false 以判别跨目标收敛） | 求解层（默认恢复拦截） |

未改变：目标/约束、`objectiveGradientScale`/`pressureGradientScale`（保持 1.0）、B2/B3 阈值、J/JT 数学定义、`applyProdFlowJT`、MMA 开关。未添加任何推测性 `alpha -> phi` 设计导数。

## 9. 最终状态

```text
BLOCK_MMA
```

依据：production 伴随未收敛（thermalCoupling 0.073、pressureDrop 0.157，均 > 1e-9），按规程禁止解锁梯度与 MMA。`frozenGradientValidated` 保持 `false`（未手工设置）。

## 10. 下一步建议（信息增益排序）

1. **谐 Ritz GCRO-DR（GMRES-DR）**：pilot 显示最小 Ritz 值 θ≈0.03 单一主导，deflation 保留 k≈8-32 个谐 Ritz 对是当前最可能突破点（GPT 专家建议首选）；需可靠的 H 特征分解（QR）替代失效的反幂法。
2. **two-level coarse correction**：若最终残差为低 graph-frequency 模式（GPT 建议先做 graph-Laplacian 商 qP 判别），构造粗网格 U-p 算子做长波校正（保持 applyProdFlowJT 不变）。
3. **face-Vanka / 两 cell 8×8 block smoother**：若 qP 判据显示高频 face-to-face 振荡主导，实现 face 块局部平滑。
4. production 收敛后依次执行：Stage B4 回归（dot test / 块测试 / Gate H1·H4 / P-row FD / explicit-core oracle）→ Stage B2/B3 正式验证（D1/D2/D3 逐方向）→ 解锁 → MMA smoke。
