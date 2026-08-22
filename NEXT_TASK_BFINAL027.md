# TO-ANISOTROPIC — 任务书（BFINAL-027：M2/M3 测量轮 — a/b 段判别）

> 背景与授权：BFINAL-026 审阅（`88ee077`）确立混合裁决——c 段（T-消除折叠）在 D1 定罪（稳健），
> D3 源干净但其全梯度 28× 超配必产于下游 **a（λ_T 求解）/ b（thermalC 收缩）** 段。本轮完成三段
> 判别的最后一步。**纯离线 Python，零编译、零求解器运行、零生产代码改动。**
> **本轮不提交 git**（与并行轮防冲突；协调人审阅后统一提交）。

---

## 任务指令

你是 TO-ANISOTROPIC 的执行 agent，任务 **BFINAL-027 — M2/M3 测量**。仓库
`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`（HEAD 88ee077）。

### 输入（全部已就位）

| 输入 | 用途 |
|---|---|
| `/home/ys/dsH/b25_qgate/1/Tb`（+同目录 phi/phiThermal/T/alpha/DTEffective 类场） | **生产 λ_T**（热伴随场，state-B） |
| `src/AdjHeatTransfer.H`（L25-100：转置方程构造 + 内置 ATa 点测试） | 生产 λ_T 求解语义（只读） |
| `src/sensitivity.H`（thermalC/装配分解，B19 打印值 0.001113） | b 段收缩公式（只读） |
| `evidence/agent-group/BFINAL-018/cycle-1/` 的 b18_folding_lab.py（热算子装配 L661-764，注意其 0.995-1.003 的精度边界与 DEFECTS 注记） | 离线热算子模板 |
| `evidence/agent-group/BFINAL-026/cycle-1/NOTEBOOK.md` + `b26a_NOTEBOOK.md` | 源码钉扎与缺陷清单 |
| `evidence/agent-group/BFINAL-026/cycle-1/b26a_m1_instrument.py`（read_field/网格/C/Gx 等可复用件） | 复用优先于重写 |

### 硬交付 1：M2 — λ_T 外部锚（判 a 段）

1. 从 b25_qgate state-B 场离线装配热算子 A_T（`fvc::div(phiThermal,T) − laplacian(DTEffective,T)`
   的离散语义，含边界系数；模板 b18 lab，修正其已知 DEFECTS）；
2. 装配 Q 型源 b_Q = dJ_Q/dT（公式：出口贴壁胞 −phi_f/(Tref·M_frozen)；M_frozen 用
   state-B 的 sum(phi_out)，与生产快照口径一致）；
3. 直接解 A_T^T λ_ref = b_Q（33600 未知元，规模小，splu 即可）；
4. 对照生产 `1/Tb`：relL2、符号翻转计数、边界带（出口/入口/壁面）分段误差。
   **判别**：relL2 O(1) 或结构性符号差异 → **a 段（λ_T 生产求解）有缺陷**；relL2 在
   离线装配精度界内（~1e-2 量级，受 b18 模板 0.5% 边界约束）→ a 段无罪。
5. 顺带审计：生产热伴随求解的收敛判据/容差（AdjHeatTransfer.H 的 solve 参数与日志证据）。

### 硬交付 2：M3 — thermalC 重算（判 b 段）

用生产 `1/Tb`（非 λ_ref）按 sensitivity.H 的 thermalC 公式离线重算场与投影
（对 D1/D2/D3），对照：①B19 分解打印的 L2=0.001113；②M1 已测的 rhs 链。
**判别**：重算 ≠ 生产值（超出离线精度界）→ **b 段（thermalC 收缩装配）有缺陷**；
一致 → b 段无罪，则 D3 的 28× 超配来源需上移重审（如实报告）。

### 纪律

- 开工建 `evidence/agent-group/BFINAL-027/cycle-1/NOTEBOOK.md`（每小步追加）；
- 预注册先行：`PREREGISTRATION.md` 写判别矩阵与预期（含 a/b 有罪/无罪四象限的预期形态）再跑数据；
- 所有中间数落 tsv；预算 **150 分钟**；单步卡 >15 分钟写 NOTEBOOK 报告，勿自行换路径；
- 结束：`EXECUTOR_SUMMARY.md` + 空 `EXECUTOR_DONE`；**不 git 提交**（协调人统一处理）；
- 数据齐即强制收尾（B26 教训：总结是交付物，勿为扩展分析耗尽上下文）。

---

## 审阅者备注（不复制）

- 重点核：离线 A_T 的边界系数语义（zeroGradient T 出口、fixedValue 入口对转置的影响）；
  b_Q 的 M_frozen 口径；四象限判别的诚实性（无罪结论必须给出精度界论证而非"差不多"）。
- 若 a 段有罪：修复候选集中在 AdjHeatTransfer 的 BC/转置细节，对照资料由并行轮 BFINAL-028 提供。
