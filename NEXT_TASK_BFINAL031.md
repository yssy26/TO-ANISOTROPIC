# TO-ANISOTROPIC — 任务书（BFINAL-031：装配链统一审判轮）

> **战略背景**：`STRATEGIC_REASSESSMENT_2026-08-22.md` 暂停期重评后恢复执行；本轮是三条路线
> 共需的前置信息轮。三条独立证据线汇聚到 `sensitivity.H` 装配/收缩链：
> ①J-D3 源干净（M1 ratio 1.057）但全梯度 28× 超配（B27 排除法）；②gDP 源恒等式缺口
> `b_PD^T w_true / ADJ_gDP = 0.818/0.366/1.036`（B30）；③M4 开放差异：rxpr 导出投影
> (−0.056/−0.006/−0.036) vs 生产 ADJ_gDP (−5.230/1.051/13.603)（B26a，注意差两个数量级，
> 先校口径再谈真相）。**纯离线测量轮，零生产改动。** 分支 `agent/dsH-stage-b-validation`
> @ `00a7f27`。

---

## 任务指令

你是 TO-ANISOTROPIC 的执行 agent，任务 **BFINAL-031 — 装配链统一审判**：对 J 与 gDP
两条链的装配/收缩层做**逐项对账**，找出打破恒等式闭合的项。

### 输入（全部已就位）

| 输入 | 含义 |
|---|---|
| `/home/ys/dsH/b25_qgate/1/`（Ub、Uc、pb、fsensMeanT、fsenshMeanT、dfdx、gsensPressureDrop、dAlphaDxh、U、p、alpha…） | state-B 生产伴随场与装配结果场（Uc=pressureDrop 标签伴随速度、Ub=thermalCoupling 标签，B13 惯例） |
| `/home/ys/dsH/b8_verify_diag/rxpr_prod_gsensh_{momentum,pressurerow,total}.mtx` + `rxpr_prod_gsens.mtx` | 装配分段导出（场状态与 b25 逐位相同，B26 NOTEBOOK 已证） |
| `evidence/agent-group/BFINAL-026/cycle-1/b26_wstar_h7.npz` | w_true（闭合 1e-11 已验收） |
| `b25_qgate/b18rhs_{pressureDrop,thermalCoupling}.mtx` | 两标签源（b_PD / b_TC） |
| B27/b26a/b29 仪器件（read_field/网格/V 金字塔/D1-3 方向场/过滤器与投影切线） | 复用优先 |
| `b25_qgate/stage_b2_fd_scan.tsv` + `stage_b2_summary.tsv` | FD 真值与生产 ADJ 投影（勿重跑） |

### 第〇阶段：预注册（先于数据）

写下对账恒等式与预期闭合精度：
```
对每标签 L ∈ {J(thermalCoupling), gDP(pressureDrop)}，每方向 Dk：
  生产 ADJ 投影（tsv） == Σ_项 [momentum_L + pressureRow_L + fluxDirect/Gx_L + thermalC_L(仅J)]
  其中每项由导出场（Uc/Ub/pb/dAlphaDxh/U/V）离线重算
  且 Σ_项（流介导部分） ≈ b_L^T w_true（源恒等式，M1/B30 口径）
```
预注册：哪条腿预期机器级闭合（装配 vs 自身分解的重算）、哪条腿预期就是病灶（D3 的 28×、
gDP 缺口 0.82/0.37/1.04、M4 的数量级差先判口径还是真差）。

### 第一阶段：三口径对齐（先于真相）

三个"仪器"的归一化/投影约定必须先对齐：生产 tsv 投影、rxpr mtx 投影（z 链式方向 vs D 方向，
B26a 已验 z 为正确链向）、源恒等式。M4 的两个数量级差必须先解释（scale 因子？投影方向？
导出时序？）——**若是口径差，修正口径后重算；若是真差，升级为病灶证据。**

### 第二阶段：逐项对账与定位

对两标签三方向完成恒等式三腿对账（生产 / 离线重算分段 / 源收缩），产出对账矩阵。
**病灶判定**：打破闭合的具体项（哪个 term、哪个方向、差多少、什么形态——符号/系数/缺项/双计）。
特别注意：D3 的 28× 超配应能被某项的异常贡献解释；gDP 的 0.82/0.37/1.04 同理。

### 第三阶段：裁决

- 定位到具体项 → 写修复轮规格（不动生产；含 B13 分解式口径的修订声明——B13 曾定位
  gDP 2.13× 于动量行项，其对账口径若被本轮修正，需显式声明）；
- 未定位 → 如实报告断点在哪一腿，给下一判别设计。

### 纪律

- NOTEBOOK + 预注册先于数据；产物 `b31_` 前缀入 `evidence/agent-group/BFINAL-031/cycle-1/`；
- 读 `ORIENTATION_FOR_EXECUTORS.md`（含 ccache 与上下文纪律）；
- 单步卡 >15 分钟写 NOTEBOOK 报告；预算 180 分钟；数据齐即收尾（总结是交付物）；
- 结束 EXECUTOR_SUMMARY + DONE + 一次选择性提交（不 push）。

---

## 审阅者备注（不复制）

- 重点核：三口径对齐是否先于病灶判定（M4 数量级差的解释是本轮第一道闸）；D3 的 28× 是否
  被定量解释（不是"方向一致"级解释）；Uc/Ub 标签用反与否（B13 惯例 vs 实际导出核对）；
  分段重算是否用到 b27 的精确 V 金字塔。
- 判定树：定位到项 → 修复轮（B32）请授权；断点在源收缩腿 → 回看 B30 的 defcheck 口径；
  断点在生产 tsv 腿 → 仪器问题，回 B2 模块。
