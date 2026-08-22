# TO-ANISOTROPIC — 下一阶段任务书（BFINAL-026：热伴随链定向诊断轮）

> 用法：把下面「任务指令」整段复制给执行 agent。**派发本任务书 = 协调人按判定树授权**（BFINAL-025 审阅后的最优判断路径；用户对 B/C 选项保留裁决权）。**本轮是纯诊断轮：零生产代码改动**。授权范围：开关门控的诊断导出（默认关闭）、离线量具脚本、诊断算例运行、证据与文档。**禁改清单不变**：两个伴随算子文件的 J/J^T 数学、`AdjHeatTransfer.H` 的生产数学、`rxPressureRowTranspose.H`、`filter_chainrule.H`、MMA、`validateMmaUnlockGate.H`、一切验收阈值、两种既有目标类型路径。分支 `agent/dsH-stage-b-validation` @ `055db8c`（审阅见 `evidence/agent-group/BFINAL-025/cycle-1/REVIEW_BY_REVIEWER.md`）。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC（OpenFOAM-7 拓扑优化求解器，`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）的**执行 agent**。当前任务：**BFINAL-026 — 热伴随链（dJ/dT 通道）定向诊断轮**。

背景（B25 已确立的事实）：J 的设计梯度**全部**由 dJ/dT 通道承载（dJ/dphi 通道贡献 ~1e-10 相对）；J 门失败（D1 反号、D2/D3 幅值超配 ~6.7×/28×）与目标泛函结构无关；自 B18 以来 J 的 ADJ 列对所有干预响应 ≤1e-10。病灶在三段之一：
- **(a) 热伴随求解**：λ_T（`AdjHeatTransfer.H` + adjointOutlet*Heat 边界条件）——λ_T 本身从未被外部验证过；
- **(b) thermalC 收缩**（`sensitivity.H` 的 BFINAL-019 装配分解第四项，L2=0.001113）；
- **(c) T-消除折叠**：dJ/dT → b_TC 的 (U,p) 行 → λ_TC → momentum/pressureRow/fluxDirect（L2 合计 ~0.0023，momentum 主导）。

核心机会：**源层恒等式 b_TC^T w_true vs FD_J 自 B15 之后从未用清洁仪器复测**（当时三项仪器缺陷全部在场：A-B 线性化点错配、stale kf、离线基符号翻转——B21–23 已全部修复）。B15 的旧读数 2.6×/7.4×/158× 不可信。

### 第〇阶段：预注册（先于一切数据）

1. 重建清洁量具协议：state-B 单点、B23 修正符号的离线基、`stageB18RhsExport` 机制、Q 泛函（b25_qgate，量规更干净：源无非对消结构）。
2. 写下判别矩阵（先于数据）：
   - 若 **源层恒等式闭合**（b_TC^T w_true ≈ FD_J(Q) − C − Gx，三方向 ≤5% 量级）→ 折叠 (c) 无罪，病灶在 (a) 或 (b)；
   - 若**源层失败** → 病灶在 (c) 的折叠语义，按方向特征归档；
   - (a) vs (b) 的判别：λ_T 外部锚（见下）与 thermalC 离线重算的偏差归属。
3. 预注册各分支的预期数字（幅值量级即可），含 w_true 的获取方式（清洁未钉扎 SuperLU 解，B21 惯例，闭合验证 ~1e-11 级）。

### 第一阶段：三段测量（全部离线/诊断，零生产数学改动）

**M1 — 源层恒等式（最高优先）**：在 b25_qgate（或同构克隆 b26_diag）state-B 导出：修正基下的 J 算子、b_TC（分 dJ/dT 部与 dJ/dphi 部两块导出）、直接项 C（热直接项）与 Gx（通量直接-alpha 项）。离线计算 b_TC^T w_true，对照 FD_J(Q) − C − Gx（FD 值直接取 B25 tsv，勿重跑 FD）。三方向 D1/D2/D3。
**M2 — λ_T 外部锚**：生产 λ_T（AdjHeatTransfer.H 解出的场）与「导出的热算子转置直接解」逐场对照（relL2）；再加一个切线恒等式锚：⟨λ_T, R_T,x·d⟩ vs 热通道 FD 响应（F1 型：冻结 U,p,phi 只重解 T 的 FD，可用开关门控在 B2 模块内实现或复用既有 Stage-A/B0 机器）。
**M3 — thermalC 收缩重算**：用导出的 λ_T 离线重算 thermalC 场，与生产 thermalC（B19 分解打印）对照；任何差异按符号/尺度/支承归类。
**M4 — 分段投影表**：对 D1/D2/D3 把 (a)热通道贡献、(c)流介导贡献（momentum/pressureRow/fluxDirect 经 λ_TC）、Gx 直接项分别投影，与 FD_J(Q) 真值并列——哪段贡献了超配的 6.7×/28×，哪段贡献了 D1 的反号，一目了然。

### 第二阶段：裁决与产出

1. 判别矩阵逐项裁决，写入 FINAL_REPORT（已验证/强假设/推测分级，B24 惯例）；
2. 若定位到**可实现的具体缺陷**（像 B23 之于 H7 那样）：写修复轮规格（不动手，等授权）；
3. 若三段全清洁 → 结论「热链装配层无罪，缺陷在更高链层/共同模态」，把 B24 的共同模态假设升级为最强残余解释，建议其检验路径。

### 纪律与不变量

- 纯诊断：生产数学零改动；新增导出一律开关门控默认关闭；静态测试 15/15 不得回归；
- `frozenGradientValidated=false`、`mmaUpdateEnabled=false` 不动；不写解锁文件；
- 零拟合因子；先 PREREGISTRATION.md（含判别矩阵与预期）再跑数据；
- 证据进 `evidence/agent-group/BFINAL-026/cycle-1/`（追加式）；EXECUTOR_SUMMARY.md + EXECUTOR_DONE；
- 诊断算例若复刻 b14/b24_diag 模板，注意其 legacy-ILU NaN 已知崩溃（B24 的 Run B 教训）：全部交付物要在该崩溃点之前完成或避开该求解器路径；
- 选择性 git 提交一次（信息风格同 B25），**不 push**。

---

## 审阅者备注（不需复制给执行 agent）

- 审阅重点：①M1 的 w_true 是否清洁未钉扎且闭合验证过；②FD 值是否直接取自 B25 tsv（不得重跑重定义口径）；③M2 的「导出热算子」是否用 B23 修正基（这是全部历史读数被污染的根源，重蹈即全废）；④判别矩阵是否真的先于数据落盘；⑤M4 分段投影的代数是否独立推导（各段之和 == 总 ADJ_J 的自洽校验）。
- 判定树：定位到具体缺陷 → 用户授权修复轮（B27）；三段全清 → 共同模态假设升级，回战略菜单（那时选项将是共同模态检验 vs 方法语义近似 vs 停止）。
- 风险提示：λ_T 的求解容差/收敛性本身先审计（AdjHeatTransfer 的解若未收敛，M2/M3 的对照会被污染）——B18 给流伴随加了 1e-12+GRADSTABLE，热伴随没有等价审计记录。
