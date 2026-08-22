# TO-ANISOTROPIC — 任务书（BFINAL-026 尝试 2：M1 源层恒等式测量，瘦身版）

> 本任务书取代尝试 1 的执行部分（中止记录见 `evidence/agent-group/BFINAL-026/cycle-1/ABORT_NOTE.md`）。
> **范围严格收窄**：纯离线 Python 测量，零编译、零求解器运行、零生产代码改动。
> 所有输入文件已就位（协调人已备好 w_true）。你的工作是"用现成仪器做测量并裁决"，不是"自建实验室"。

---

## 任务指令（直接复制这段）

你是 TO-ANISOTROPIC 的执行 agent，任务 **BFINAL-026 尝试 2 — M1 源层恒等式测量**。
仓库 `/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`（HEAD 055db8c，勿提交源码）。

### 输入清单（全部已存在，逐个核验后使用）

| 文件 | 含义 |
|---|---|
| `evidence/agent-group/BFINAL-026/cycle-1/b26_wstar_h7.npz` | **w_true（D1/D2/D3）**，H7 算子上的清洁未钉扎解（协调人执行，闭合 1e-11 级，验收记录在同目录 `.json/.log/.py`） |
| `evidence/agent-group/BFINAL-026/cycle-1/PREREGISTRATION.md` | 尝试 1 的预注册（判别矩阵与恒等式代数，**作为本轮代数依据**） |
| `evidence/agent-group/BFINAL-026/cycle-1/NOTEBOOK.md` | 尝试 1 的源码钉扎（模板缺陷清单、公式出处——**凡引用 b18_folding_lab.py / b23 模板必读其 DEFECTS 注记**） |
| `/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx` | Q 型、B 态的 b_TC（热耦合流伴随源，已折叠） |
| `/home/ys/dsH/b25_qgate/stageB2/stage_b2_fd_scan.tsv` | FD_J(Q) 真值（勿重跑 FD；取 h=1e-3 行并报平台稳定性） |
| `evidence/agent-group/BFINAL-021/cycle-1/b21_wstar_cache.npz` | pre-H7 w_true（对照控制组） |
| `/home/ys/dsH/b8_verify_diag/stageB6_*` 系列 | rxc/rpa/gsensh 分段导出（M4 弹性目标用） |

### 硬交付：M1 恒等式（判别矩阵的第一判据）

对 D1/D2/D3 各方向计算并并列报告：

```
lhs = b_TC^T · w_true(H7)
rhs = FD_J(Q) − C − Gx        （C/Gx 的代数按 PREREGISTRATION 钉扎的口径离线重算；
                               b18_folding_lab.py 的 C_field/Gx_of 是模板但含 DEFECTS，
                               修正后使用并在报告列出你修正了什么）
比值 lhs/rhs、符号、相对偏差；同时报对照控制组 b_TC^T·w_true(pre-H7)（H7 敏感度）。
```

裁决（对照判别矩阵，先于数据已在 PREREGISTRATION 写死）：
- 三方向 |lhs/rhs − 1| ≤ ~10% → **源层干净** → 病灶在 λ_T 求解或 thermalC 收缩（a/b 段）；
- O(1) 失败 → **病灶在 T-消除折叠（c 段）**，按方向特征归档（对照 D1 反号/D2/D3 超配 6.7×/28× 的指纹）。

### 弹性目标（硬交付完成后有余力才做，勿为此延长大轮）

- **M4 分段投影**：用 b8_verify_diag 的 rxpr_prod_gsensh_{momentum,pressurerow,total} 导出对 D1/D2/D3 投影，thermalC 段 = total − momentum − pressureRow（口径声明），并列 FD 真值——指认超配与反号各来自哪段。
- **M2 预备**：仅盘点 λ_T 生产场导出是否存在（Ta 字段等），不实施求解。

### 纪律

- 开工先建 `b26a_NOTEBOOK.md`（磁盘笔记本，每小步追加）；新产物一律 `b26a_` 前缀，写入 `evidence/agent-group/BFINAL-026/cycle-1/`；
- 先把 PREREGISTRATION 的判别矩阵与预期数字**原样引用**到 `b26a_PREREG_CONFIRM.md`（若你要改任何口径，先停下报告，不得静默改）；
- 零拟合因子；所有中间数（lhs/rhs 各分量）落 tsv 供审阅者重算；
- 结束：`b26a_EXECUTOR_SUMMARY.md`（三行裁决 + 关键数字表 + 与判别矩阵对照）+ 空文件 `b26a_EXECUTOR_DONE`；一次选择性提交（信息风格仿 08ff950，注明 attempt 2 + 协调人 w_true 工件），**不 push**；
- 预算：全部离线计算应在 **90 分钟内**完成；任何单步卡住 >15 分钟即写 NOTEBOOK 后报告，勿自行更换算法路径。

---

## 审阅者备注（不复制）

- 审阅重点：lhs 的符号约定（b_TC 导出的符号语义 vs w_true 求解的 −rxd 约定——b21 模板 trans='T' 与 M=J^T 的复合，最容易错一版符号）；C/Gx 修正是否逐条引用 NOTEBOOK 的 DEFECTS；FD 行选取（h=1e-3）是否与 B22/B24 口径一致；控制组（pre-H7）差值是否与 H7 份额 ~0.7 的既有认知自洽。
- 判定树不变：定位 (a)/(b)/(c) → 修复轮规格请授权；源层干净但 M4 显示 thermalC 段超配 → λ_T 外部验证轮。
