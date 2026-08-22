# TO-ANISOTROPIC — 下一阶段任务书（BFINAL-014：梯度收缩修复轮）

> 用法：把下面「任务指令」整段复制给执行 agent。分支 `agent/dsH-stage-b-validation` @ `5aaae31`（BFINAL-013 已由独立审阅确认：两缺陷定位到 `sensitivity.H` 的 `−(U & U_adj)·V·dAlphaDxh` 类动量行收缩项；审阅见 `evidence/agent-group/BFINAL-013/cycle-1/REVIEW_BY_REVIEWER.md`）。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC（OpenFOAM-7 拓扑优化求解器，`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）的**执行 agent**。当前任务：**BFINAL-014 — 梯度收缩修复轮**。BFINAL-013 已把 gDP（恒定 ≈2.13× 偏大）与 J（符号翻转 + 6.7–11.2×）定位到 `sensitivity.H:65-66/89-90` 的动量行收缩项；本轮推导正确形式、实现修复、并完整重验。

### 第〇阶段：先推导，后动手（本轮最重要的纪律）

**禁止用拟合因子的方式猜 2.13**（已知 1/α_rel=2.5、1/(1−α_rel(1−α_rel))≈1.92 都不等于 2.13）。必须：

1. 从 `NS.H` 的**实际前向路径**写出 FD 链所收敛的离散不动点残差 R_U(x; U, p)：包括方程松弛（`equationRelaxationFactor("U")`=0.4、p 松弛）、`bounded Gauss upwind` 的 bounded 路由、非正交修正、fvOptions、压力修正闭环——逐项列出哪些项依赖设计 α（直接或通过系数场）；
2. 对该不动点做符号微分，得到正确的动量行设计导数收缩 `−λ_U^T (∂R_U/∂xh)`，逐项与现行 `−dAlphaDxh·(U & U_adj)·V` 对照，指出缺失/错标项；
3. 把推导写成文档（进证据目录），标明每一项的来源行号（NS.H / sensitivity.H / fvSolution）。

三个待判假设（BFINAL-013 遗留）：H-F1 松弛语义失配、H-F2 `Sp(alpha,U)` 之外的 α 依赖（冻结 laplacian 系数、fvOptions）、H-F3 FD 不动点与伴随线性化不动点的系统差异。推导应能唯一裁决或排除三者。

### 实现阶段

- 授权修改范围：`sensitivity.H` 的收缩项（:65-66 / :89-90 及直接相关装配），以及推导证明必需的最小新增场/系数读取；
- **如果推导表明必须修改已验证算子层（J_PU / J_PP / R_x / matrix-free J^T / 目标定义 / 过滤投影）→ 立即停止并报告**，那是单独授权的范畴，且需全部既有门重跑；
- 修复保持 BFINAL-005 R_x 压力行项不动（它只占 2.9% 且内部一致 1.7e-15）；
- 每个独立的修复项单独提交（commit 粒度 = 可独立回滚）。

### 验证阶段（顺序执行，全部通过才算完成）

1. **小网格快检**（若已有小网格复现台则用之，否则在现有网格上跳过此步直接进 2）：FD-of-FD 交叉验证修复后的收缩；
2. **P1 源点测试复跑**：必须维持机器精度 PASS（源未受修复影响）；
3. **完整 BFINAL-012 重跑**：三方向 × 三目标 × 全 eps 阶梯。验收：
   - gDP：ADJ/FD ∈ 1±0.05（有清晰平台支撑的稳定 1±0.10 可附条件接受，须明示）；
   - J：符号全对、幅值误差 <5%（同上 5–10% 平台条款）；
   - **J 与 gDP 必须由同一收缩修正同时归一**——若只归一其一，如实报告并停下（说明 J 另有叠加因素，回到定位逻辑）；
   - gV：维持既有 1e-7 量级（不回归）；
4. **无回归三件套**：PRODPRECGAMGCHECK 四值逐位不变；伴随收敛不变（≤1e-9，迭代数允许因梯度装配无关而逐位一致）；静态测试全过（8+，修复若改变相关字符串断言须同步更新测试）；
5. 修复后的梯度语义变化写入文档：`AI_AGENT_HANDOFF.md` 加日期注记、`OPENFOAM_USAGE_AND_FILES.md` 相应条目更新。

### 保持不变

- `mmaUpdateEnabled=false`、`frozenGradientValidated=false`——**即使全部通过，解锁（置 frozenGradientValidated=true）也是用户在审阅你的报告后的决定**；你只需准备好解锁所需的证据摘要（哪个 commit、哪些门、哪些数字）；
- 外层容差 1e-9、等价性门 1e-8、gV formal 阈值 1e-6 均不动；
- 停止条件：任一验收检查点失败 → 立即停止并定位报告，不得调参/放宽/重试掩盖。

### 证据与提交

- 证据目录：`evidence/agent-group/BFINAL-014/cycle-1/`：推导文档（第〇阶段产物）、修复前后对比（分项投影 + FD 对照）、完整验证日志、EXECUTOR_SUMMARY.md、FINAL_REPORT.md；
- **先归档再清理**任何日志；完成后写 `EXECUTOR_SUMMARY.md` 并创建空标记 `EXECUTOR_DONE`；
- git 选择性提交（修复 + 测试 + 证据 + 文档），**一律不 push**。

---

## 审阅者备注（不需复制给执行 agent）

- 审阅重点：第〇阶段推导是否真的从 NS.H 实际路径出发（而非从伴随代码反推）；修复是否最小；J 与 gDP 是否同时归一；PRODPRECGAMGCHECK/伴随是否逐位不变；若修的是「松弛语义」，注意与 BFINAL-003 已验证的 P 行松弛语义的一致性论证；
- 若推导走向 H-F3（FD 不动点差异），修复可能在 FD 侧收敛规程（更深的重收敛）而非梯度侧——同样合法，但属「验证方法修正」，须在报告中明确区分「改了物理」还是「改了验证」；
- 修复通过后，下一站即 C0 串行 MMA smoke 准备（解锁决定留给用户）。
