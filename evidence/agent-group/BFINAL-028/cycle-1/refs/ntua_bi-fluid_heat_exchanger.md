# 文档级参照 — NTUA 连续伴随双流体换热器 TO（best-effort tier 3）

来源：Crossref API（https://api.crossref.org），抓取于 2026-08-21。
arXiv 检索（export.arxiv.org API）无命中（Papoutsis-Kiachagias 的双流体换热器工作未见 arXiv 预印本）；Semantic Scholar API 持续 429，弃用。

## 命中文献（NTUA 组，Galanos/Papoutsis-Kiachagias/Giannakoglou）

1. **"Synergistic use of adjoint-based topology and shape optimization for the design of Bi-fluid heat exchangers"**
   Struct Multidisc Optim (2022), DOI: 10.1007/s00158-022-03330-w
   （Crossref 记录无公开摘要文本；标题确认其框架为 adjoint-based TO + ShpO 双流体换热器设计）

2. **"A continuous adjoint cut-cell formulation for topology optimization of bi-fluid heat exchangers"**
   IJNMF Heat & Fluid Flow (2025-02-20), DOI: 10.1108/hff-08-2024-0642
   摘要关键句（Crossref jats 摘要）：
   - 采用 "Think Discrete – Do Continuous" (TDDC) adjoint methodology；
   - 设计变量为人工不可渗透性(impermeability)场；流固界面(FSI)逐周期计算，在 FSI 上施加精确边界条件求解 CHT；
   - 相比标准密度法(denTopO)，避免在流动方程里加 Brinkman 惩罚项（改用 cut-cell FSI）。

3. **"The Cut-Cell Method for the Conjugate Heat Transfer Topology Optimization of Turbulent Flows..."**
   Energies (2024), DOI: 10.3390/en17081817（开放获取）

## 对 BFINAL-028 的参照价值（结论性）

- NTUA 组 CHT-TO 全部采用**连续伴随**，其热→动量耦合即为 Fira 文档所引的 NTUA residual convention 下的 `+T_a∇T` (ATC-T) 源——与 Fira thermalTopO 同源同构。
- TDDC 方法论命名（"Think Discrete – Do Continuous"）与 Fira 的 ATC-T open-channel 结论一致：连续源本身需配合 fixed-point map transpose 才是分离式 SIMPLE 的完备离散转置。
- 对 D1 反号诊断的独立增量：**双流体换热器（冷热两股流）的 TO 参照确认，伴随温度场 T_a 携带的正是"上游/冷侧来源"的敏感性信息，其与原始温差在动量源中的耦合为正号 ATC 形式**；未提供任何支持"在热折叠上加整体负号"的证据。
- 局限：paywall，未取得全文逐式；其逐槽离散折叠（αRel·rAU·dH^T 等）仍需依赖 BFINAL-018 槽位镜像论证与 M1/M2 恒等式自洽检查。
