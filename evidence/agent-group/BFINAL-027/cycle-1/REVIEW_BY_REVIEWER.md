# BFINAL-027 cycle-1 — Independent Review (REVIEW_BY_REVIEWER)

- 生成时间：2026-08-22 09:40 (+080)，协调人
- 依据：EXECUTOR_SUMMARY + b27_* 全套产物抽查复核（gates/m2_compare/m3 tsv、run.log）+ 与 M1 锚交叉对照
- 状态：无 git 提交（按并行轮约定，协调人统一提交）

## 1. 复核结论

| # | 声明 | 复核 | 判定 |
|---|---|---|---|
| 1 | 仪器门：A_T·T−src−Q = 1.18e-11；A_T^T·Tb−dJdT = 1.91e-11 | b27_gates.tsv 直读 | ✅ |
| 2 | M2：λ_ref vs 生产 Tb relL2=1.15e-12、cos=1.0、零符号翻转、全边界带 ~1e-12 | b27_m2_compare.tsv 逐行：INTERNAL/inlet/outlet/hotInlet/hotOutlet/solidEndWalls/bottomWall/sideWalls 全部 0.7-2.9e-12；topWall 1.38e-8 但 maxabs 6.2e-30（物理零场的归一化伪影，声明属实） | ✅ |
| 3 | M3：重算 C_L2 与三方向投影 vs 生产 | b27_m3.tsv：比值全部 1.000000；与 B26 M1 的 C 锚（−0.0303/−0.0082/−0.0198）逐位吻合 | ✅ |
| 4 | 无罪论证带精度界（非"差不多"） | 两个判词都低于离线精度界（~1e-2）十个数量级；门 (b) 同时锁死装配语义与生产转置求解 | ✅ |

## 2. 审阅判定与战略后果

**PASS（判词双清白，机器精度级，证据链完整）。** 战略后果重大：

1. **a 段（λ_T 求解）与 b 段（thermalC 收缩）排除**——热伴随链的"老嫌疑"（B14 起挂账"λ_T
   从未被外部验证"）正式销案：λ_T 现在有直接解外部锚，逐场逐带 1e-12。
2. **D3 的 28× 全梯度超配被逼入 c 段/装配层**：源干净（M1）+ λ_T 无罪（M2）+ thermalC 无罪
   （M3）⇒ 超配只能产于流介导收缩/装配（momentum/pressureRow/fluxDirect）或折叠本身。
   这与 B26 的 M4 开放差异（rxpr 导出投影 vs 生产 ADJ_gDP 严重不匹配）指向同一层——
   该开放项从"仪器疑点"升级为"主嫌疑现场"。
3. 开放项 3（b26a read_field 的 uniform 值潜在 bug 已在 b27 修复）记录在案，不影本轮判词
   （其影响若存在只作用于 B26 的离线支路，M1 的 lhs 用导出 b_TC 不受染）。

## 3. 下一步输入

判别域外遗留两个精确目标：①c 段折叠的 H7s2/T4 逐槽隔离（B28 首选嫌疑）；②D3 装配层间隙
（M4 差异的正式归位）。见 defect map（agent-group 级合成文档）。
