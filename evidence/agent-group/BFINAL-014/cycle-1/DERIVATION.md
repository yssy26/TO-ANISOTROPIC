# BFINAL-014 cycle-1 — 第〇阶段推导文档（含停止裁决）

日期：2026-08-20。基线 `5aaae31`。本文件按任务书要求：从 NS.H 实际前向路径出发推导，
逐项对照现行收缩，裁决 H-F1/H-F2/H-F3。**结论：推导不支持任何 sensitivity.H 收缩项修改；
缺陷定位收敛到「J 作为不动点系统 Jacobi 矩阵」这一已验证算子层 —— 按任务停止规则停止。**

## 1. NS.H 实际前向路径（逐行）

每个 corrector（`NS.H:23-1192`，冻结湍流、`fvOptions` 为空——日志确认 "No finite volume options present"）：

| 步 | 行号 | 操作 |
|---|---|---|
| 1 | 26-43 | 组装 `UEqn = fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U)`，RHS `−fvc::grad(p) + fvOptions(U)`；ransFlowModel 时 `UEqn −= fvc::div(nuEffFrozen·dev2(T(grad U)))` |
| 2 | 45 | `UEqn.relax()`（equationRelaxationFactor("U") = 0.4） |
| 3 | 52 | `solve(UEqn == −grad(p))`（动量预测） |
| 4 | 113-116 | `rAU = 1/UEqn.A()`（**松弛后**对角）；`HbyA = constrainHbyA(rAU·UEqn.H(), U, p)`；`phiHbyA = flux(HbyA)` + adjustPhi |
| 5 | 1018-1034 | `pEqn: fvm::laplacian(rAtU(), p) == div(phiHbyA)`（非正交循环；**rAtU=rAU=松弛后**）；`phi = phiHbyA − pEqn.flux()` |
| 6 | 1038-1041 | `p.relax()`；`U = HbyA − rAtU·grad(p)`；correctBCs |

α（设计）依赖项枚举：
- `fvm::Sp(alpha, U)`（NS.H:30）——动量对角，**唯一**的动量行 α 依赖；
- `nuEffFrozen = nu + nutFrozen`（updateFrozenTurbulenceFields.H:60）——**与 α 无关**（H-F2 动量侧排除）；
- `fvOptions` 空（H-F2 fvOptions 侧排除）；
- 压力行：rAtU = 1/A_relaxed，A_relaxed 含 α（经 Sp）→ ∂rAtU/∂α = −rAU²/α_rel（未钳制单元）；rxPressureRowTranspose.H:15 的 `drAU = −rAU²/alphaRel` 与此一致 ✓（该项仅占原始梯度 2.9%，BFINAL-013 实测）；
- 对流系数经 phi 依赖状态（phi 是被消除的状态；J^T 已有 velocityJump/lambdaPdiff 路由——见 §3 讨论）。

## 2. OF7 `fvMatrix::relax()` 精确形式（fvMatrix.C:521-670）

```
D0 = diag();                         // 含边界 internalCoeffs
sumOff = Σ|offdiag|（+coupled 边界 |boundaryCoeffs|）
D = max(|D|, sumOff);                // 支配性钳制（非线性、依赖系数）
D /= alpha;                          // ÷ α_rel
S += (D − D0)·psi();                 // 源项用当前 ψ
```

**不动点代数**（ψ_prevIter = U）：松弛系统 `D_rel·U = S + (D_rel−D0)·U` ⇒ `D0·U = S` ——
**不动点残差 = 未松弛动量方程，与松弛因子无关**（钳制同理被源项精确抵消）。
因此：
- `∂R_U/∂α = +U`（每单位体积；fvMatrix 体积分 semantics，diag 含 α·V）→
  现行收缩 `−dAlphaDxh·(U & U_adj)·V`（sensitivity.H:65-66/89-90）**公式正确**；
- H-F1（纯标量松弛失配）**被拒**：不动点与 α_rel 无关，任何 1/α_rel 类纯因子都无法同时满足
  「恒定 2.13×」与「不动点代数」；且实测比值跨方向有 1.3% 离散（2.113–2.169），非纯标量特征。

**钳制机制**（H-F3 的具体载体之一）：被钳制单元 A_rel = sumOff/(α_rel·V) ≠ A_unrel/α_rel，
破坏 BFINAL-003 P←U 切线的 α_rel-比例假设。**定量排除**：本算例 |alpha|₂≈1.46e10
（StageB6 fingerprint，alpha RMS ≈ 659）支配对角 2–3 个量级，钳制单元占比可忽略。

## 3. 对照现行收缩与裁决

| 项 | 现行 | 推导 | 判定 |
|---|---|---|---|
| 动量行 `−dAlphaDxh·(U&U_adj)·V` | sensitivity.H:65/89 | `−λ_U^T(∂R_U/∂α)=−λ_U·U·V·dα/dxh` | **公式一致** |
| 压力行 `−rxPressureRowT·dAlphaDxh` | rxPressureRowTranspose.H | `−λ_P^T(∂R_P/∂α)`（drAU=−rAU²/α_rel） | 一致（2.9% 份额） |
| J 动量行 | production.H:57-62 未松弛矩阵 | 不动点 R_U 的未松弛导数 | 一致 |
| J 压力行 mobility | NS.H:129 实际 rAtU（松弛后） | 不动点 R_P 的实际 mobility | 一致 |
| **J 作为系统 Jacobi** | 冻结 phi 系数 + relaxed-SIMPLE P←U 切线（BFINAL-003） | 消去 phi 后的真实不动点 Jacobi | **从未被外部验证**（见 §4） |

三个假设的裁决：
- **H-F1（松弛语义纯标量）：拒**（§2 不动点代数）；
- **H-F2（Sp 之外 α 依赖）：拒**（nuEffFrozen/fvOptions 均 α-无关，读码确证）；
- **H-F3（FD 不动点 vs 线性化系统的系统差异）：存活，且收窄为「J 的系统-Jacobi 语义」**
  ——具体即消去 phi 后的完整耦合（对流系数对状态的敏感性）与 BFINAL-003 P←U 切线的组合是否等于
  真实不动点映射的 Jacobi。该问题存在于**已验证算子层**（J_PU/J_PP/J_U,phi 的系统组合），
  不在本轮授权范围内。

## 4. 本轮新增裁决实验（零编译，B4 诊断求解器 + B6 oracle，`Log.verify_b014_b4b6.txt`）

1. **RX-A（外部 FD，冻结态）**：`∂R_U/∂α` 解析 vs 设计-FD relL2 = 9.2e-8/3.2e-7/9.1e-7/3.0e-6
   （eps 1e-3→3e-5），cos=1；`∂R_P/∂α` 同样 ~1e-6。**R_x 算子外部正确。**
2. **StageB6 加权 D1/D2/D3（历史「D_mom-prod ≈ 2e-9」验证的复跑）**：本次诊断路径 λ 退化
   （`lambda_U^T R_U,x d ≈ 7.4e-47`，随后 MMA 非有限门 abort，MTO_RC=134）——**该历史验证的构造被
   揭穿为 λ-相对内部一致性**（同一 λ 自比），从未构成 λ 的外部验证。附带发现：**B4 诊断求解器路径
   在当前 HEAD 产生退化/非有限伴随**（生产路径不受影响——BFINAL-010/011 已修）。
3. 外部验证覆盖表（终版）：
   - 源（∂Obj/∂w）：BFINAL-013 P1 机器精度 ✓
   - R_x（∂R/∂α）：本轮 RX-A 外部 FD ~1e-7 ✓
   - 链（xh→x）：gV 端到端 1e-7 ✓
   - **λ（= J⁻ᵀ 源）：从未外部验证** ✗ —— 历史 BFINAL-006 任务书的 T2/T3（切线解 vs 状态 FD）
     被明确搁置且从未运行；P8 只证明了 J·w′=−R_x·d 可精确求解（6e-14），未对照状态真值。

## 5. 结论与停止

端到端失恒 2.13×（及 J 的结构失配）在给定正确 λ 下无法由任何收缩项公式解释；
所有其他环节已被外部证据封闭。**缺陷收敛于 λ 本身，即 J 作为消去-phi 不动点系统 Jacobi 的语义
——属已验证算子层，按任务书「推导表明必须修改已验证算子层 → 立即停止并报告」执行停止。**

下一轮建议（按性价比排序）：
1. **T2/T3 补课（首选，零新数学）**：用 BFINAL-008 P8 既有机制解 `J·w′ = −R_x·d`
   （D1/D2/D3），与「x±h·d 全链重收敛后的 (U,p) 状态差分」直接对照——一次运行即可测得
   λ-coded 与真值之比的方向依赖谱，直接判定 J 系统语义缺陷的具体块（J_U,phi / J_PU / 消去项）。
2. J_PU 切线的独立外部锚（phiHbyA 的 U/p-FD vs 解析切线）。
3. B4 诊断求解器路径修复（独立小任务，恢复 λ-相对 oracle 的可用性）。
