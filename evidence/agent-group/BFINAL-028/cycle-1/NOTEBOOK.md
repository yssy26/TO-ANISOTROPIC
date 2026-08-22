# BFINAL-028 cycle-1 — c 段参照系对照轮（纯阅读/对照）

日期：2026-08-21。分支 agent/dsH-stage-b-validation @ 88ee077。
任务书：NEXT_TASK_BFINAL028.md。纪律：零编译、零求解器运行、零生产代码改动、不 git 提交。

## 0. 输入核验表

| 输入 | 状态 |
|---|---|
| NEXT_TASK_BFINAL028.md | 已读 |
| ORIENTATION_FOR_EXECUTORS.md | 已读 |
| ABORT_NOTE.md（BFINAL-026 前车之鉴） | 已读 |
| solveDiscreteFlowAdjointProduction.H（b_TC 折叠 L346-508 + J^T L512-803 + solve L834-2210） | 已读全 |
| AdjNS_HT.H（g_f 公式 discrete 分支 L1-115） | 已读全 |
| AdjHeatTransfer.H（A_T^T + C 项 L189-328） | 已读全 |
| computeObjective.H（三目标导数） | 已读全 |
| B18 DERIVATION_FIX3 / FINAL_REPORT | 已读 |
| B23 PREREGISTRATION_T7_ADDENDUM | 已读 |
| B24 FINAL_REPORT | 已读 |
| B26 b26a_NOTEBOOK / REVIEW_BY_REVIEWER / ABORT_NOTE | 已读 |

## 1. 符号约定（贯穿本轮）

- φ_f：SIMPLE 收敛面通量；(U,p)：动量/压力解。
- g_f = dJ/dφ_f − (dR_T/dφ_f)^T·Tb：面泛函（B18 fix-3，AdjNS_HT.H discrete 分支装配）。
- M = ∂φ/∂(U,p)：relaxed-SIMPLE 通量切线（BFINAL-003 路由基准）。
- b_TC = M^T·g（= prodExternalTranspose，先 U 后 P 槽位）。
- 恒等式：b_TC^T·w_true == FD_J − C − Gx；C = −Tb^T·R_T,x（thermalDiffusionDerivativeDTCell）；Gx = g^T·φ_x。
- w_true：流伴随解（U 行 = λ_U，P 行 = λ_p），由 prodRhs 缩放前 = b_TC 的线性系统解出。
- D1 指纹：ADJ_J = −0.0358 vs FD_J = +0.0425（反号，O(1)）；M1 lhs=−0.0488/rhs=+0.0744，ratio=−0.656。
- 候选错项排名必须与 D1 反号形态挂钩：O(1) 反号最可能是符号/缺项类错误，非系数微差。

## 2. 预注册式：我方 b_TC 折叠完整项清单（先于任何参照系抓取）

### g 层（面泛函 g_f 的装配，AdjNS_HT.H discrete 分支）

| # | 代数形式 | 符号 | 源码 file:line | 语义出处 |
|---|---|---|---|---|
| g0 | g_f = −mask_f·Tb_down·(T_nei − T_own)；Tb_down = φ_T≥0 ? Tb[nei] : Tb[own]（upwind） | **显式 MINUS** | AdjNS_HT.H:11-26 | B18 DERIVATION_FIX3：g 公式无罪（内部面热残差转置） |
| g1 | 出口面 g_f += +dJ/dφ_out = −(T_out−T_in)/(M_frozen·Tref)（B25）或 −(T_out−T_mix)/(M·Tref)（旧型） | **+dJ/dφ** | AdjNS_HT.H:28-41; computeObjective.H:141-143(L164-171 B25) | B0.2：<Tb,dR_T/dφ_out> 在 bounded-Gauss-upwind 下恒零 → 合并 RHS 恰为 dJ/dφ |
| g2 | 入口/其余边界 g_f = 0 | — | AdjNS_HT.H 无赋值 | U-fixed 边界在路由层被跳过（见 T5 组） |

### 路由层——内部面（Production L351-386）：M^T 逐槽位转置

| # | 代数形式 | 符号 | 源码 file:line | 对应 relaxed-SIMPLE 通量切线槽位 |
|---|---|---|---|---|
| T1 | HA[own] += αRel·rAU_o·w·(Sf_f·g_f)；HA[nei] += αRel·rAU_n·(1−w)·(Sf_f·g_f) | + | Production:359-362 | αRel(w·rAU_o·dH_o+(1−w)·rAU_n·dH_n) 的 hA 累积 |
| T3 | U(own,c) += (1−αRel)·w·(Sf_f·g_f)_c；U(nei,c) += (1−αRel)·(1−w)·(Sf_f·g_f)_c | + | Production:363-369 | (1−αRel)(w·dU_o+(1−w)·dU_n) 直接内插转置（relax 源项导数） |
| T4 | kf = interp(mob)·δ_f·\|Sf_f\|；P(own) += kf·g_f；P(nei) −= kf·g_f | own+/nei− | Production:370-378 | −kf_f(dp_n−dp_o) 的转置；与算子自身 J^T kf 槽同号（applyProdFlowJT P(own)+=kf·λPdiff） |
| H7s1 | g_c[own] += mob_o·w·g_f·Sf_f；g_c[nei] += mob_n·(1−w)·g_f·Sf_f | + | Production:379-385 | H7 stage1：phiHbyA 压力通道 −interp(mob·∇dp)·Sf 的 g_c 累积（B24） |

### 路由层——边界（Production L387-422；U-fixed patch 整段跳过 L389-392）

| # | 代数形式 | 符号 | 源码 file:line | 说明 |
|---|---|---|---|---|
| T5 | HA[celli] += αRel·rAU_c·(Sf_b·g_b) | + | Production:403-404 | assignable-U 边界 hA 累积 |
| T6 | U(cell,c) += (1−αRel)·(Sf_b·g_b)_c | + | Production:405-409 | assignable-U 边界 direct |
| H7s1b | g_c[celli] += mob_c·g_b·Sf_b | + | Production:410-412 | assignable-U 边界 H7 stage1 |
| T7 | P(cell) += mob_c·δ_b·\|Sf_b\|·g_b | + | Production:413-420 | 仅 pressureFixed（出口 p=0）；入口 p zeroGradient 无此项 |
| T5skip | U-fixed patch 跳过 T5/T6/H7s1b | — | Production:389-392 | inlet/墙不计边界通量 g |

### deltaH^T 展开（Production L428-465）

| # | 代数形式 | 符号 | 源码 file:line | 说明 |
|---|---|---|---|---|
| T2 | U(nei,c) −= prodUpper·(1/V_o)·HA[own].c；U(own,c) −= prodLower·(1/V_n)·HA[nei].c | **−** 两侧 | Production:428-441 | H 切线转置 off-diag；与 applyProdFlowJT hA 路径同号 |
| T8 | U(cell,c) += (−ic_c+cav)·(1/V_c)·HA[cell].c；cav=mean(ic) | + | Production:442-465 | 边界 H 对角（各向同性时恰零，保留算子奇偶性） |

### H7 stage 2（Production L466-497）

| # | 代数形式 | 符号 | 源码 file:line | 说明 |
|---|---|---|---|---|
| H7s2 | q = Sf_f&(g_o/V_o − g_n/V_n)；P(own) −= w·q；P(nei) −= (1−w)·q | **−** 两侧 | Production:469-481 | G^T 展开 distance-2 P-P（B24 推导 D5）；镜像 applyProdFlowJT stage2 |
| H7s2b | P(cell) −= (g_c&Sf_b)/V_c | **−** | Production:482-497 | 非 fixedValue-p 边界自项（zeroGradient） |

### 装配（Production L499-508）+ 求解层缩放（L834-845）

| # | 内容 | 源码 file:line |
|---|---|---|
| A1 | discreteFlowVelocityRhs[celli] += prodExternalTranspose[U]；discreteFlowPressureRhs[celli] += prodExternalTranspose[P] | Production:499-508 |
| A2 | prodRhs[U] = rhs_U/prodMomentumScale；prodRhs[P] = rhs_P/prodAreaScale；导出 b18rhs_thermalCoupling.mtx 为物理未缩放 rhs（4 列 Ux Uy Uz P） | Production:834-845, L870-898 |

### 范围外但参与 M1 恒等式的项（不在 b_TC 内，供对照）

| # | 内容 | 源码 file:line |
|---|---|---|
| C | C = −Tb^T·R_T,x：thermalDiffusionDerivativeDTCell（内部面 baseDerivative=(Tb_o−Tb_n)(T_o−T_n)·δ·\|Sf\|，own −= w·dDnDxh_o·base，nei −= (1−w)·dDnDxh_n·base；物理边界 += dDnDxh·Tb_c·snGradT·\|Sf\|；coupled −= wL·dDnDxh·base） | AdjHeatTransfer.H:189-328 |
| Gx | Gx = g^T·φ_x（通量直接 α 项），sensitivity.H 装配层（授权外，生产缺） | B18 DERIVATION_FIX3 |

### 对 D1 反号形态的预注册候选错项机制（后验对照用）

按 B26 裁决，D1 反号对 C/Gx 误差稳健（需 C+Gx 偏差 >180% 才翻号）→ 候选错项必在
b_TC 链（g 层符号 / M^T 路由符号 / 缺项）：
- E1（g 层符号）：g0 的显式 MINUS 或 g1 的 +dJ/dφ 若是反号，全部方向同翻（D1/D2/D3 同号翻转——但 D3 ratio=+1.057 干净，若 E1 真则 D3 也应反）→ **E1 被判为 D3 不支持**（预注册，后验核对）。
- E2（路由层符号）：T4 kf own+/nei− 与 J_PP 同号 vs 反号；T2/T8 deltaH^T 号；H7s2 号。
- E3（缺项类）：参照实现有而我方无的通道（H7 即此类）。例：热方程对 U 的直接耦合源（如 buoyant/物性随 T 变）、R_T 对 φ 的入口贡献、出口温度泛函的其他边界项。
- E4（整链缺失）：参照系中热目标折叠进流伴随源的方式与我方整体不同（如连续伴随能量耦合源 +T·∇Tb 项——我方 discrete 分支无此项，continuous 分支有 AdjNS_HT.H:133）。

预注册锚点：若参照系确认热残差对 (U,p) 的耦合在目标为"出口温度/总换热量"时**仅通过 φ 通道**存在（冻结流场、零浮力、物性冻结），则我方 g 层结构（仅 φ 依赖）在形态上正确，嫌疑集中在 M^T 路由符号/缺项；若参照系有**直接 T→U/p 耦合**，则为 E4 类缺项。

## 3. 参照系抓取计划（按可达性递补）

1. SU2（github.com/su2code/SU2）：CHT/多物理离散伴随——热方程折叠进流伴随、热残差对 U/p 偏导进伴随源、出口温度类边界泛函。
2. OpenFOAM ESI adjointOptimisation（github.com/OpenFOAM/Plus）：porosity(Brinkman) 敏感度、adjointOutlet* 伴随边界、动量伴随源耦合项符号。
3. 文档级（尽力）：Papoutsis-Kiachagias 双流体换热器 TO（arXiv）、UPC 离散伴随论文、FOAMacademy。

抓取规则：git clone --depth 1 或 curl raw 优先；单源死磕 >10 分钟换源；片段存 refs/ 带 URL。

## 4. 时间轴

| 时刻 | 事件 |
|---|---|
| T0 | 目录建立；任务书/定向包/源码/历史文档读完 |
| T1 | 本预注册项清单写入（先于参照系） |
