# TO-ANISOTROPIC — 下一阶段任务书（BFINAL-025：目标函数 Q 化 + FD 门复测轮）

> 用法：把下面「任务指令」整段复制给执行 agent。**派发本任务书 = 用户授权**（2026-08-21，用户在 BFINAL-024 审阅后的战略菜单中选定「选项 1：目标函数 Q 化 + 带保护烟雾测试」；本轮是其第一半）。授权范围：**新增**一个 `thermalObjectiveType` 目标类型 `maximizeTotalHeatTransfer`（换热总量 Q），及其在目标求值/导数/伴随源/候选评估各挂钩的实现与验证。**两种既有目标类型（`legacyMeanTemperature`、`maximizeColdOutletTemperature`）的代码路径必须逐位不变**——它们连同算子层（`solveDiscreteFlowAdjointProduction.H`/`solveDiscreteFlowAdjoint.H` 的 J/J^T/H7 折叠）、`rxPressureRowTranspose.H`、`filter_chainrule.H`/过滤/投影、MMA、`validateMmaUnlockGate.H`、一切验收阈值**全部仍在禁改清单**。烟雾测试与解锁**不在本轮**（等审阅后单独派发）。分支 `agent/dsH-stage-b-validation` @ `f0034c0`（审阅见 `evidence/agent-group/BFINAL-024/cycle-1/REVIEW_BY_REVIEWER.md`）。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC（OpenFOAM-7 拓扑优化求解器，`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）的**执行 agent**。当前任务：**BFINAL-025 — 目标函数 Q 化与 FD 门复测轮**。背景：15 轮调查后，J（冷侧出口混合温度 `maximizeColdOutletTemperature`）的 FD 门仍失败（D1 反号、幅值病理敏感），已归因于该泛函的**近对消结构**（B16 实测：J = −N/(M·Tref)，M=sum(phi_out)、N=sum(phi_out·T_out)，dJ/dphi=−(T_f−Tmix)/(M·Tref) 中 (T_f−Tmix)≈0 且 1/M 小，流/热两路几乎抵消）；而 handoff §1 本来就首选「换热总量 Q」。本轮把目标换成 Q，检验其 FD 门能否用当前算子精度直接通过。

### 第〇阶段：推导（先于一切实现，写入 PREREGISTRATION.md）

1. **Q 的定义**（冷侧能量平衡形式，只用边界量）：Q = N − T_in·M，其中 N、M 同上，T_in = 冷侧入口 T 的 fixedValue 值（从算例 BC 读取，勿硬编码）。稳态下即冷流吸热量 ∝ m_dot·(Tmix−T_in)。在报告中给出与热侧/固体界面传热的一致性论证（能量守恒叙述即可，不要求新增界面积分器）。
2. **归一化（关键设计决定）**：J_Q = −Q/Qref，Qref = thermalReferenceTemperature × M_frozen，其中 M_frozen = sum(phi_out) 在**每个外层优化状态的开头快照一次**，该状态内的**所有** J 求值（含 FD ±h 扰动态、候选评估）都用同一快照。硬性要求：同一次梯度/FD 评估内 Qref 必须是常数——**绝不允许 J 求值时用当前 phi 重除**（那会重新引入商结构/近对消病理）。等价形式：在快照口径下 J_Q = −(Tmix−T_in)/Tref。快照的落点与刷新时机在实现里显式化并在报告说明。
3. **导数**（相对当前 J 逐项对照）：
   - dJ_Q/dT[cell] = −phi_f/(Tref·M_frozen)（冷侧出口贴壁胞，支承同现有 dJ/dT，仅分母换成冻结常数）；
   - dJ_Q/dphi[f] = −(T_f−T_in)/(Tref·M_frozen)（出口面）——**这是本轮的核心动机**：(T_f−T_in) 是 O(温升) 的非小量，近对消结构消失；
   - 审计 dJ/dT、dJ/dphi、以及 sensitivity.H/伴随链中**每一个消费目标泛函的槽位**，逐槽标注「泛型（自动正确）/需新增分支/旧类型专用假设（列出）」。特别核查 BFINAL-018/024 建立的 b_TC 面泛函折叠与 `discreteExternalFaceFluxAdjoint` 路径是否完全经 `thermalObjectiveDerivative`/`thermalObjectiveFluxDerivative` 泛型驱动。
4. **预注册预测**（先于数据）：J_Q 的 FD 门符号 3/3、幅值因子 ADJ/FD ∈ [0.9,1.1]（预期机理：源不再近对消 → 不再触发病态放大/符号翻转）；gV 逐位不变；gDP 与 B24 相同（~2.0×，已知状态，非本轮失败判据）。

### 第一阶段：实现（最小侵入，新增分支，不碰旧类型路径）

挂钩清单（逐个加 `maximizeTotalHeatTransfer` 分支）：
- `src/createFrozenHotRegionFields.H`：白名单（L136-154 一带）+ Qref 快照变量的创建与读取；
- `src/costfunction.H`：Q 值求值 + 归一化 + M_frozen 快照刷新点；
- `src/computeObjective.H`：dJ/dT 与 dJ/dphi；
- `src/evaluateCandidate.H`：候选评估分支（L66/L141 一带）；
- `src/AdjHeatTransfer.H` 及 b_TC 折叠路径：仅当第〇阶段审计发现旧类型专用假设时才动，逐条引用审计条目；
- `src/tests/test_stage_b_safety_gates.py`：新增第 15 项静态测试（Q 类型存在性/白名单/Qref 冻结语义的静态断言）；
- 新算例 `/home/ys/dsH/b25_qgate`：b24_gauge 同构克隆，仅改 `thermalObjectiveType maximizeTotalHeatTransfer;`（保持 tolerance 1e-12 等全部其他控制）。

### 第二阶段：验证（顺序执行，任何一步意外失败即停下定位报告，禁止调参掩盖）

1. `wclean; wmake` → WMAKE_EXIT=0，存日志；
2. 静态测试 15/15；
3. **旧目标回归**：重跑 `b24_gauge`（未改任何输入），FD 侧三列（J/gDP/gV）、ADJ_gV、gDP 因子、等价性门四值必须与 B24 证据**逐位相同**（证明旧路径零扰动）；
4. **新目标 FD 门**：跑 `b25_qgate` 完整 BFINAL-012 口径 FD 门（ε 阶梯 + 三方向），产出与 B24 同构的 tsv（FD_J/ADJ_J/relJ/signJ、gDP、gV、双标签伴随 trueRelRes、GRADSTABLE、P 行恒等式指纹）；
5. F1 热通道硬门（validateFrozenGradient.H：median≤1%、max≤5%、零符号失配）在新目标下必须通过；
6. 对照预注册预测逐项裁决（过/不过/意外，意外要归因）。

### 纪律与不变量

- `frozenGradientValidated=false`、`mmaUpdateEnabled=false` 保持不动；**不写 FROZEN_GRADIENT_UNLOCK.txt**；
- 零拟合因子；不改任何验收阈值；
- 先 PREREGISTRATION.md 再跑数据；证据进 `evidence/agent-group/BFINAL-025/cycle-1/`（追加式，不动 B24 及更早）；
- 结束时：EXECUTOR_SUMMARY.md（含改动清单、关键数字表、与预注册对照）+ 空标记文件 EXECUTOR_DONE；
- 选择性 git 提交一次（信息风格同 B24），**不 push**；
- 诚实报告：若 J_Q 的 FD 门不过，如实报数字并归因，不得为过门而改口径。

---

## 审阅者备注（不需复制给执行 agent）

- 审阅重点：①Qref 冻结语义是否真的在全链路成立（重点盯 FD ±h 求值点有没有偷偷用当前 M 重除——这是本轮最大的实现风险）；②第〇阶段的槽位审计是否逐槽完成、旧类型专用假设是否全部列出并处理；③旧目标回归是否逐位（不是「差不多」）；④预注册是否先于数据落盘；⑤J_Q 符号表与幅值因子逐方向复核（本人从 tsv 重算）。
- 判定树：J_Q 门过（3/3 符号 + 幅值 ≤10% 稳定平台 + 旧路径逐位回归）→ 派发 BFINAL-026 带保护烟雾测试轮（解锁证据包 + 5–10 MMA 迭代 + SST 候选验收/回滚层验证）；J_Q 门不过 → 回战略菜单（共同模态假设将获得新的强证据：泛函结构不是主因）。
- 风险提示：若执行者发现 Q 的 dJ/dphi 消费链在某处硬编码了 Tmix 形式，修正属本轮授权范围，但必须逐条引用审计条目并保持旧类型逐位不变。
