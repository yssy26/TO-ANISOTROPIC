# TO-ANISOTROPIC — 下一阶段任务书（BFINAL-012）

> 用法：把下面「任务指令」整段复制给执行 agent。分支 `agent/dsH-stage-b-validation` @ `2aaf649`（BFINAL-011 已由独立审阅判定 PASS，见 `evidence/agent-group/BFINAL-011/cycle-1/REVIEW_BY_REVIEWER.md`）。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC（OpenFOAM-7 拓扑优化求解器，`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）的**执行 agent**。当前任务：**BFINAL-012 — 当前源最终 FD 幅值门**（冻结湍流梯度验证链的最后一门）。

### 背景与地位

BFINAL-011 已闭合 P0：两标签（thermalCoupling 750 迭代 / pressureDrop 902 迭代）均以真残差 ≤1e-9 在生产迭代路径收敛，全程 ~120 s，MTO_RC=0；四组矩阵等价性门值 1e-16 量级逐位不变。**生产伴随已可信**。BFINAL-012 要回答的是：用这套伴随算出的梯度，与「对 raw 设计变量做有限差分、完整重跑过滤/投影/物性/冻结原始链」得到的导数，方向与幅值是否一致到足以驱动 MMA。

### 目标

在**当前源**（`2aaf649`）与当前验证算例上，对三个目标完成定向有限差分验证：

1. 换热目标（thermalCoupling 链）；
2. 冷侧压降约束（pressureDrop 链，含 BFINAL-005 R_x 压力行贡献）;
3. 固相体积分数（解析导数，应远紧于流场项）。

### 硬性要求（规格来源：`AI_AGENT_HANDOFF.md` §13）

- **扰动施加在 raw 设计变量 x 上**（x ± h·d），每个扰动点必须完整重跑同一条 filter → projection → material → frozen-primal 链，不做任何「只重解部分」的捷径；
- **多个确定性方向**（可复用 Stage-B2 的 D1 类物理方向构造，如 `validateStageB2GradientAmplitude.H` 既有机制与 `src/validation/stage_b/stage_b_b2b3_*` 工具链），不要随机方向；
- **eps 阶梯**（建议 1e-2 → 1e-5 跨四个量级）建立稳定平台，报告每个方向的平台区间；
- 流场相关梯度（换热、压降）验收口径：**符号一致、方向稳定、幅值误差目标 <5%；有清晰 eps 平台支撑的稳定 5–10% 带可接受**；幅值差 ≥2 倍或符号翻转 = FAIL；
- 体积分数导数应达到远高于流场项的精度（解析项，预期 ~1e-10 量级或更好）；
- **重复性**：关键方向至少完整重跑一次（噪声底口径，沿用 Stage-B1 方法论），报告 run-to-run 差异；
- 每个 FD 点的伴随求解必须保持 BFINAL-011 的收敛状态（真残差 ≤1e-9）——若任何 FD 点伴随不收敛，停下报告，不得用未收敛伴随继续。

### 范围禁改清单（违反即中止）

- J_PU（BFINAL-003）、J_PP（BFINAL-008）、R_x / `rxPressureRowTranspose.H`（BFINAL-005）、matrix-free 物理 `J^T`；
- 目标与约束定义、过滤/投影链、MMA；
- 外层容差 1e-9 与等价性门 1e-8 阈值；
- `mmaUpdateEnabled=false`、`frozenGradientValidated=false` 保持不变——**即使全部通过，`frozenGradientValidated=true` 也只能由用户在审阅你的报告后决定设置，不是你的权限**。

### 验收与停止条件

- 三目标 × 多方向 × eps 阶梯的完整幅值对照表（adjoint FD vs analytic/FD，含 cos、relL2、幅值比）；
- 检查点失败（某方向幅值差 ≥2 倍 / 符号翻转 / 伴随不收敛）→ **立即停止并定位**，不得调参掩盖；
- 验收只看上述 FD 口径；GRADPROXY 平台（TC +2217.48 / PD −205194.14）只是伴随稳定性的佐证，**不能替代 FD 门**。

### 证据与提交

- 证据目录：`evidence/agent-group/BFINAL-012/cycle-1/`（FD 扫描数据、每个扰动点的运行日志摘要、eps 平台表、EXECUTOR_SUMMARY.md 与 FINAL_REPORT.md）；
- **先归档再清理**：任何日志删除前先拷贝入证据目录（BFINAL-011 的 E0/E1 教训）；
- 完成后写 `EXECUTOR_SUMMARY.md` 并创建空标记文件 `EXECUTOR_DONE`（独立审阅自动化依赖此协议）；
- git 选择性提交（代码如有 + 测试 + 证据 + 文档），**一律不 push**。

### 参考但不得盲信

- 历史 `STAGE_B_B2B3_GRADIENT_AMPLITUDE_REPORT.md` 及其数据**早于 J_PU/R_x/J_PP 修正与 BFINAL-010/011 修复，只能作方法学参考，数值结论一律以本轮重跑为准**；
- 运行/编译环境细节见 `OPENFOAM_USAGE_AND_FILES.md`（干净 PATH、`unset FOAM_SIGFPE`、绝对路径 `build/bin/MTO_HF`、单次全量编译 30–35 分钟）。

---

## 审阅者备注（不需复制给执行 agent）

- 若 BFINAL-012 全 PASS：下一步为 C0 串行 MMA smoke（5–10 迭代完整优化环），仍由独立审阅后解锁；
- 若压降链幅值超 10% 但方向稳定：先查 R_x 贡献占比与 eps 平台形状再下结论，勿直接判定 FAIL；
- 并行化（P2）继续后置于首个串行 MMA smoke。
