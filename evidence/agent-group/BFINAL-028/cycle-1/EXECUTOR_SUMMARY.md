# BFINAL-028 — EXECUTOR_SUMMARY

日期：2026-08-21。分支 agent/dsH-stage-b-validation @ 88ee077。纯阅读/对照轮。

## 一行裁决

**最高置信候选错项 = E2-路由层压力行通道（首选 H7s2：Production:469-481 的 H7 stage2 P-P 通道；次选 T4：Production:370-378 的 kf P-row），置信级中高**——唯一与 D1 反号形态（O(1) 方向依赖反号、D3 干净、M1 |lhs|<|rhs| 幅度偏差）全部吻合且无任何参照镜像覆盖的项；g0（g 层显式 MINUS）被 Fira 参照**同形同号确认无罪**。

## 交付物（evidence/agent-group/BFINAL-028/cycle-1/）

- `NOTEBOOK.md` — 输入核验 + 预注册 b_TC 项清单（先于参照系）
- `REFCROSS_REPORT.md` — 4 节：我方 17+2 项清单表 / 逐项对照表 / 候选错项排名（5 级）/ 缺项清单（5 项）
- `refs/`（6 个文件，均带 URL）：su2_avg_temperature、su2_disc_adj_methodology、fira_thermalTopO_flux_sensitivity、openfoam_othmer_adjoint、esi_adjointOptimisation_inaccessible、ntua_bi-fluid_heat_exchanger
- 本 EXECUTOR_SUMMARY.md + 空 EXECUTOR_DONE

## 过程要点

- 顺序纪律：预注册项清单（T0→T1）→ 参照抓取 → 对照 → 收尾。
- 参照系 1（SU2）：离散伴随 = AD/coDiPack（无手工 g_f 折叠）→ 提供方法论差异 + 功能形式差异（AVG_TEMPERATURE 面积平均），不提供逐槽符号对照。
- 参照系 2（OpenFOAM ESI adjointOptimisation）：模块不可达（auth/403），按 tier 规则降级为 Othmer 连续伴随（Sp(α) 自伴随无翻号）+ Fira thermalTopO（连续伴随 CHT-TO，gPhi 显式热-动量折叠）。
- **关键正结果**：Fira gPhiI（chi=0）与 g0 同形同号（迎风伴随×own−nei 温差），官方 FD 门控通过——g0 显式 MINUS 得到参照系独立背书。
- **关键佐证**：Fira atc-t-open-channel.md 明言"连续伴随热源单独不是分离式 SIMPLE 映射的完备离散转置"——与 B18 结论同构，指向路由层。
- 文档级（tier 3）：NTUA 双流体换热器 TO 经 Crossref 定位（HFF 2025 TDDC 方法论，DOI 10.1108/hff-08-2024-0642）；arXiv/SemanticScholar 不可用或 429。
- 无编译、无求解器运行、无生产代码改动、无 git 提交。

## 开放项（转下游）

1. H7s2/T4 单项隔离需 M1/M2 逐槽测试（零运行授权外的运行轮）。
2. 离散路由层仍无外部逐槽参照（真 ESI 模块不可达；连续伴随/AD 均不镜像该层）。
3. SU2 面积平均 vs 我方质量加权 Tmix 对回流出口方向的影响待数值检验。
4. NTUA 全文 paywall，仅摘要级。

## 时间

- 启动→收尾在 150 分钟预算内；数据齐后即强制收尾。
