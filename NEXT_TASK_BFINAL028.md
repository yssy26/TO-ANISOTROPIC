# TO-ANISOTROPIC — 任务书（BFINAL-028：c 段参照系对照轮）

> 背景与授权：BFINAL-026 审阅（`88ee077`）把 c 段（T-消除折叠，dJ/dT → b_TC 的 (U,p) 行）
> 在 D1 上定罪（O(1) 反号，对 C/Gx 误差稳健）。本轮用**外部参照实现**逐项对照我们的折叠代数，
> 直接找 D1 反号的具体错项。**纯阅读/对照轮：零编译、零求解器、零生产代码改动。**
> **本轮不提交 git**（与并行轮防冲突；协调人审阅后统一提交）。

---

## 任务指令

你是 TO-ANISOTROPIC 的执行 agent，任务 **BFINAL-028 — c 段参照系对照**。仓库
`/home/ys/dsH/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`（HEAD 88ee077）。

### 我方被审对象（先读熟）

1. `src/solveDiscreteFlowAdjointProduction.H`：b_TC 折叠（B18 的"算子一致转置"实现 + B24 的 H7 面
   泛函通道）——**写出场/面级的完整项清单**（每一项：系数、符号、作用面/胞、语义出处）；
2. `src/AdjNS_HT.H`（若在）：discreteExternalFaceFluxAdjoint 的边界通量项；
3. `evidence/agent-group/BFINAL-018/cycle-1/` 的推导文档（当时"exact operator-consistent
   transpose"的论证）与 `BFINAL-026/cycle-1/b26a_NOTEBOOK.md` 的钉扎。

### 参照系（按可达性递补，至少完成前两个）

1. **SU2**（github.com/su2code/SU2，shallow clone 或 curl raw 文件）：多物理/CHT 离散伴随中
   热方程如何折叠进流伴随系统（找其 adjoint 的 multiphysics 耦合源装配——同结构问题的
   独立实现）；重点：热残差对 U/p 的偏导如何进入伴随源、边界泛函（出口温度类）的处理；
2. **OpenFOAM ESI adjointOptimisation**（github.com/OpenFOAM/Plus 或镜像，v2312 含官方拓扑优化）：
   官方伴随对 porosity(Brinkman)/多孔介质的敏感度代数与伴随边界类（adjointOutlet*，与我方同名
   谱系）——重点对照其动量方程伴随源中来自其他方程（能量/压力）的耦合项符号约定；
3. 文档级（PDF 可能抓不动，能抓则用）：FOAMacademy 伴随训练、Papoutsis-Kiachagias 双流体
   换热器 TO 论文（arXiv 版优先）、UPC 离散伴随论文。

### 硬交付：对照报告

`evidence/agent-group/BFINAL-028/cycle-1/REFCROSS_REPORT.md`，内容：
1. **我方折叠项清单表**（编号、代数形式、符号、源码 file:line）；
2. **逐项对照表**：每个我方项在参照系中的对应物（或"参照系无此项"），符号/系数一致性；
3. **候选错项排名**（按解释 D1 反号的可能性），每条给双方 file:line 证据与机理论证；
4. 参照实现里**有而我方没有的项**清单（缺项类候选——H7 就是这么找到的）。

### 纪律

- 开工建 NOTEBOOK；我方项清单先于读参照系完成（预注册式：先写我方理解再对照，防锚定）；
- 抓取的参照源码片段存 `evidence/agent-group/BFINAL-028/cycle-1/refs/`（带来源 URL）；
- 网络抓取失败不要死磕 >10 分钟，换源（clone/raw/镜像）或降级到下一参照；
- 预算 **150 分钟**；结束 EXECUTOR_SUMMARY + DONE；**不 git 提交**；
- 数据/清单齐即强制收尾。

---

## 审阅者备注（不复制）

- 重点核：项清单是否完整（对照 B18 推导文档查漏）；"参照系无此项"的论断是否真的翻过该代码路径；
  候选排名是否与 D1 反号形态（而非泛泛"可能有问题"）挂钩。
- 产出若给出高置信错项 → 与 BFINAL-027 的 a/b 裁决合成完整缺陷地图 → 修复轮规格（需用户授权）。
