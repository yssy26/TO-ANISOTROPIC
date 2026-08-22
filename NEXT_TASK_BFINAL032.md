# TO-ANISOTROPIC — 任务书（BFINAL-032：分解场门控导出轮 — 一次性终结跨轮伪影）

> 背景：B31 协调人会话（`4999855` 提交，NOTEBOOK §7）发现 `/1/` 场目录疑似混合不同伴随轮次，
> 使一切离线场级分解被跨轮伪影污染（76.5% 残差之谜）。唯一可靠出路：**一次带开关门控导出的
> 重跑**，让生产自己写出全部分解场（同轮次、同点）。
> **授权范围**：仅限开关门控（默认关闭）的**诊断写出**——生产数学、求解逻辑、既有导出零改动。
> 分支 `agent/dsH-stage-b-validation` @ `4999855`。

---

## 任务指令

### 阶段 1：实现门控写出（一次编译）

在 `sensitivity.H` / `rxPressureRowTranspose.H` 加开关（optProperties 键
`stageB31WriteDecomposition`，默认 false；true 时 `.write()` 或 setWrite）：
- `gsenshPressureDropMomentum`、`gsenshPressureDropPressureRow`（生产两段，已有 volScalarField）；
- `rxPressureRowT`（pc 基）与 `rxPressureRowTb`（pb 基）、`rxDHbyA`、`rxrAU`、`rxG0Field`（内部面段）；
- `gsenshMeanT` 全分解四段（B19 已有：momentum/pressureRow/fluxDirect/thermalC——若当前仅打印，
  补成字段写出）；
- 按轮次标注：两轮伴随解出的 `Uc/Ub/pc/pb/Tb` 分别以后缀 `_r1`/`_r2` 写出（如已写出的
  `/1/` 值属于哪轮，在 log 打印明示）。
- 静态测试若需同步（不删断言）。

### 阶段 2：编译 + 运行

1. 编译用定向包 §1 ccache 命令块（改了 TU，全编 ~40 分钟属预期）；
2. b25 同构克隆 `b32_decomp`（仅加 `stageB31WriteDecomposition true;`），全流程运行（~26 分钟）；
3. **回归底线**：FD 侧三列（J/gDP/gV）与 b25 tsv **逐位一致**（导出不改变计算路径——若有
   任何非逐位，立即停止报告，那是实现污染了计算）。

### 阶段 3：直接分解（数据到手后立即做，勿拖延）

用写出的同轮次场做场级分解（复用 `BFINAL-031/cycle-1/b31b_m2.py` 的加载骨架）：
1. `gsenshPD_r? == momentum + pressureRow`（同轮闭合门，应逐位）；
2. 生产 `rxPressureRowT(pc)` vs 协调人离线 T 环（`b31b_m2.py` 的 T_row）——relL2 与 corr；
3. 用同轮 `Uc/pc` 重算 mom/prow 段投影与 B31 §7 的跨轮数字对照——76.5% 残差应大幅消失；
   消失后剩余部分即为**真病灶段**（或证无）；
4. 生产 `rxDHbyA/rxrAU` vs 协调人重建（HbyA 2% 误差校准）。

### 纪律

- NOTEBOOK + 预注册（预期：①跨轮伪影消解后残差 <5%；②T 环镜像 relL2<2%；若 ② 失败，
  差异即 T 输入的生产实现细节）先于数据；
- 产物 `b32_` 前缀入 `evidence/agent-group/BFINAL-032/cycle-1/`；
- 长编译/运行后台 + 轮询；单步卡 >15 分钟写 NOTEBOOK 报告；预算 240 分钟；
- 结束 EXECUTOR_SUMMARY + DONE + 一次选择性提交（不 push）。

---

## 审阅者备注（不复制）

- 审阅重点：FD 逐位回归（计算零污染的铁证）；写出场与打印值的抽读一致性；同轮分解闭合门；
  与 B31 §7 表的对照表（跨轮伪影消解量化的账要算清）。
- 判定树：残差消解且 T 环镜像通过 → m2 对账用同轮数据完成（回 B31 收尾）；残差残留 →
  真病灶现形，直接进入修复轮规格；FD 不逐位 → 实现污染，回炉。
