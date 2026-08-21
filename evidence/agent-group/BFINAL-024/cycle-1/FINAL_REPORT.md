# BFINAL-024 cycle-1 — FINAL_REPORT（H7 实现轮）

## 0. 执行概况

预注册（PREREGISTRATION.md，先于全部运行数据）→ 实现（两算子文件四处一致
+ 预条件器口径）→ 一次编译（WMAKE_EXIT=0）+ 静态 14/14 → 两次求解器运行
（Run A 生产路径 b24_gauge，MTO_RC=0；Run B 诊断路径 b24_diag，交付件在
崩溃前全部完成，崩溃为 b14 模板固有）→ 离线修正仪器归因。全部门通过。

## 1. H7 推导条目与实现对照（每处变更引用）

| 推导条目（PREREGISTRATION §1） | 实现 | 验证 |
|---|---|---|
| D1 源级 −∇p 一次拷贝（operator== 复制矩阵、source+=V·su，fvMatrix.C:1428 本轮源码复核；第二份只进 solve 临时） | 前向通道存在性依据 | 与 B23 T7 条目互证；NS.H:26-34/52/113-115/1032 |
| D2 HbyA = rAU·H() ⊃ −rAU·∇p；SIMPLE 字典无 consistent → rAtU≡rAU≡mob=primalPressureMobility（=alphaRel·prodRAU，2.19e-15） | 通道系数 mob = primalPressureMobility（kf 同源） | b24_reattribution [0]：mobility 复现 2.185e-15 |
| D3 通道 = −interp(mob·∇dp)·Sf（P 行 p 切线的缺失槽，B23 H7 裁决） | deltaPhiFacePU 追加项 + 出口边界 −(S_b&mob_c·∇dp_c) | 转置点测试 6.4e-14；oracle 5.5e-16 |
| D4 ∇dp = Gauss linear（案例 gradSchemes）；zeroGradient-p → dp_b=dp_c、fixedValue-p(出口) → 0 | 梯度基逐面装配（内部 owner+/neighbour−，边界 per-BC） | 与 T7b grad_p 同构；CSR 条目含边界自项 |
| D5 转置槽位（stage1 g_c 累积 / stage2 G^T 展开；对角元逐项复核 −w²mob_o\|S\|²/V_o + w(1−w)mob_n\|S\|²/V_n = ∂phi/∂dp_o） | applyProdFlowJT / applyDiscreteFlowJT stage1+2 + 边界两段；COO 显式展开（distance-2 P-P + 自项 + λP(c) 列） | BlockDot PP 1.53e-14；ExplicitJToracle 5.48e-16（CSR +1.53M nnz） |

四处一致清单：生产 J^T（applyProdFlowJT）、诊断前向 J（applyDiscreteFlowJ）、
诊断 J^T（applyDiscreteFlowJT）、显式 CSR；另 b_TC 外部面泛函折叠（两文件）
同步镜像（face functional 代 λ_f）。

## 2. 预条件器/等价性门口径（审阅者 §3-2 注记的裁决）

**决定：pressureGAMG Poisson 矩阵保持 kf-Laplacian 近似不变。** 推导：H7
通道把 J_PP 变为 (kf Laplacian) + (G^T·diag(mob/V)·M·div) 的 distance-2
耦合；fvMatrix laplacian 是 distance-1 模板，结构上不能表示 H7 部分；
预条件器无需等于块本身。等价性门参照系声明更新为「J_PP 的 Laplacian(kf)
部分」，四值与 B22 **逐位相同**（1.1314/0.8258/1.4221/1.3430 与
1.1379/1.0744/1.5306/1.3337 e-16，r1/r2）——参照与矩阵均未动，满足
「若 J_PP 本体不变则四值必须逐位」的强条件在 Laplacian 部分上成立。
新增 PRODH7SHARE（方向 sin(0.173(c+1))）：|H7|/|kf| = 0.702/0.701
（r1/r2）——被排除份额的量化读数。Jacobi/L1 对角未动。

## 3. 验证序列结果

1. **编译**：WMAKE_EXIT=0（wmake_b24.log）；静态 14/14（新增
   test_bfinal024_h7_phihbya_pressure_channel）。
2. **转置自洽**（Run B，b24_diag）：12 随机点测试 max rel 6.409e-14
   （b14：6.864e-14）；BlockDot PU/PP/UU/UP/boundaryU =
   2.08e-13/**1.53e-14**/7.31e-14/1.80e-12/1.80e-13（全部 ≤1e-12）；
   ExplicitJToracle maxRelL2 5.48e-16、maxRelL_P 5.48e-16（显式 vs
   matrix-free，H7 显式条目被该 oracle 强制覆盖）。
3. **等价性门**：见 §2（四值逐位）。
4. **伴随双标签**：tolerance 1e-12 下 TC 9.81e-13（1123 迭代）、PD
   9.48e-13（1323）[r1]；9.25e-13（1112）、9.69e-13（1102）[r2]；
   GRADSTABLE 4/4（relChange ≤3.4e-12）。迭代数变化解释：容差收紧
   （B22 为 1e-9）+ J_PP 谱变化（H7 0.70×kf 幅度）。
5. **修正仪器重估归因**（b24_reattribution）：见 §4。
6. **FD 门**：见 §5。
7. **残差-FD 锚**：GateH4（0.00334657148198, relL2=0, cos=1）与
   GateOracle U 通道（1.6537e-11/1.36377032534/0.790384568764）与 b14
   **逐位相同**——H7 对 U-only 方向恒为零（推导预期），核心块语义未回归。
   RxPressureRowTranspose[pc] relErr=0（设计链未触碰）。

## 4. 归因表（修正仪器 + 新算子）

面级（b24_gauge wstate，FD 侧与 B22 逐位相同；du+kf 无 da = 生产语义）：

| 方向 | pre-H7 | post-H7 | +da 参考 | \|H7\|/\|FD\| | \|div(tangent)\| 变化 |
|---|---|---|---|---|---|
| D1 | 4.92% | 4.75% | 2.68% | 0.7% | −17.8% |
| D2 | **31.17%** | **22.12%** | 38.16%（da 有害） | 9.7% | −23.2% |
| D3 | 4.93% | 4.56% | 2.45% | 0.6% | −17.3% |

post-H7 数字与 B23 T8（0.2216）/T9 A-form（22.12%）精确互证 → 实现与
被证实的仪器同一数学。da 参考列复现 B23 T9 的方向特异结论（da 在 D2 净
有害、D1/D3 必需——生产 J 无 da 槽为 B18/B19 闭层现状）。

**gDP 2.15× 与 J(D1) 的分层归因：**
- **层 1（H7 修正，已落地）**：面级 D2 −9.0pp、D1/D3 −0.2~0.4pp；系统级
  ADJ_gDP 三方向一致向 FD 移动（−3.9%/−7.0%/−3.7%），因子
  2.145/2.159/2.181 → 2.061/2.007/2.100。
- **层 2（D2 集中缺陷放大）**：残余 D2 面级 22.1% 失配仍向 P 行注入
  4.8e-4 连续性残差（真值 6e-11），经压力反演放大（B21 机制）；H7 已
  削减该负担 23%。
- **层 3（J 链更高层）**：**决定性观察——H7 后三方向 gDP 因子方向均匀
  （2.01-2.10×，极差 0.094），而三方面级切线误差相差 4×（4.6% vs
  22.1%）**。若 O(1) 失配来自 P 行面级切线缺陷，D1/D3 的因子应显著小于
  D2 —— 数据否定该归因。剩余 O(1) 失配（以及 J(D1) 符号翻转：ADJ −0.0358
  vs FD +0.0429）属共同模态/更高链层（目标收缩、design-row、或标定），
  非算子 P 行切线。

## 5. FD 门（G6，如实报告）

FD 侧与 B22 逐位相同（FD 不读算子）。ADJ 列（h=1e-3 接受步）：

| 方向 | 指标 | FD | ADJ(B22) | ADJ(B24) | 符号 |
|---|---|---|---|---|---|
| D1 | J | +0.04286 | −0.03492 | −0.03584 | 反（不变） |
| D1 | gDP | −2.5374 | −5.4441 | −5.2303 | 同 |
| D2 | J | −0.00347 | −0.02332 | −0.02493 | 同 |
| D2 | gDP | +0.52341 | +1.13008 | +1.05052 | 同 |
| D3 | J | −0.00028 | −0.01397 | −0.01200 | 同 |
| D3 | gDP | +6.4772 | +14.1277 | +13.6034 | 同 |
| 全部 | gV | — | 逐位同 B22 | 逐位同 B22 | 同 |

**J 符号表不变（D1 反，2/3）**；gDP 因子头条 2.145/2.159/2.181 →
**2.061/2.007/2.100**。J(D1) 仍反号：其归因分解 = 层 3（方向均匀的
O(1) 共同失配叠加在 thermal 目标链上；D1 面级切线不劣于 D3，符号翻转
不能由 P 行切线解释）。

## 6. H8 决定（不实现，本轮）

推导支持存在精确的边界松弛源切线（exact_coef_c = mob_c·(clamp_c/alphaRel
− bnd_post_c − D0_int_c)/V_c 替代 (1−alphaRel)），但 B23 T8 实测其切线
与 D2 残差正交（corr 0.001-0.003，D2 0.2216→0.2224），值级 2.60→2.00%
为离线仪器映射性质（生产 mobility 本为精确导出）；且完整精确化牵涉
B23 未预注册的 preBnd-in-clamp 与 mob vs alphaRel·prodRAU 边界胞差。
按「不确定则只做 H7」执行。同列下轮候选：design-row（rxPressureRowTranspose）
的 −interp(drAU·∇p)·Sf 设计导数一致性（B18/B19 闭层，本轮禁改）。

## 7. Run B 崩溃归因（非 H7）

b24_diag 复刻 b14 模板：legacy ILU 诊断求解 thermalCoupling relRes=nan →
fsensMeanT nan → sensitivity.H:483 MMA 数据门 abort（MTO_RC=134）。
b14 模板（H7 前代码）同因同位同 nan（b14 日志 1282 行）。Run B 交付件
（点测试/BlockDot/oracle/GateH1-H4/PR/Rx-dot/CSR 导出 explicitJT_H7.mtx
523MB）全部在崩溃前完成并通过。

## 8. 里程碑升级声明

BFINAL-003/008 P 行通量切线语义正式升级（日期注记 2026-08-20）：
AI_AGENT_HANDOFF.md（BFINAL-008 锁定声明之后新增 BFINAL-024 升级节）+
OPENFOAM_USAGE_AND_FILES.md §10（状态速查新增 BFINAL-024 条目）。原
BFINAL-003/008 槽位在升级后的映射中不变（等价性门四值逐位为证）。

## 9. 产物清单

PREREGISTRATION.md、EXECUTOR_SUMMARY.md、本报告、Log.verify_b24_gauge.txt、
Log.verify_b24_diag.txt、b24_reattribution.{py,log,json}、
stage_b2_{summary,fd_scan,adjoint_identity}.tsv、wmake_b24.log、
sha256_artifacts.txt、EXECUTOR_DONE（空）。运行目录：/home/ys/dsH/b24_gauge
（含 stageB2 全部 wstate/b18rhs 导出）、/home/ys/dsH/b24_diag（含
explicitJT_H7.mtx）。b8/b14/b22/b23 归档未触碰。
