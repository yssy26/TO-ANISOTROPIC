# BFINAL-024 cycle-1 — 预注册（实现轮：H7 phiHbyA 压力通道进生产算子）

先于本轮全部运行数据落盘。授权凭证 = 本派发（用户已授权算子层修复）。
授权范围：`solveDiscreteFlowAdjointProduction.H` 与 `solveDiscreteFlowAdjoint.H`
的算子装配、CSR 导出、预条件器一致性；禁改其他。

## 1. 推导条目（每处变更的依据）

**D1（源级 −∇p 一次拷贝）**：NS.H:26-34 构造
`fvm::div(phi,U) − fvm::laplacian(nuEffFrozen,U) + fvm::Sp(alpha,U) == −fvc::grad(p) + fvOptions(U)`；
OF7 `operator==(fvMatrix, Field)` 复制矩阵后 `source() += V·su`
（fvMatrix.C:1428-1436，本轮源码复核逐行确认）→ UEqn 自身源含 **−V·∇p 一次**；
`solve(UEqn == −fvc::grad(p))` 的第二份 −∇p 只进临时拷贝（B23 T7 预注册条目）。
依据：B23 PREREGISTRATION_T7_ADDENDUM §已落盘事实。

**D2（HbyA/phi 分解）**：`fvMatrix::H() = (diag·U − offdiag·U + source + bnd)/V`
（fvMatrix.C:760-812 复核）→ NS.H:113-115
`rAU=1/UEqn.A()`（relax 后）、`HbyA = rAU·UEqn.H()` 含 **−rAU·∇p**；
`phiHbyA = fvc::flux(HbyA)` 含 **−interp(rAU·∇p)·Sf**；
NS.H:1032 `phi = phiHbyA − pEqn.flux()`。案例 fvSolution SIMPLE 字典无
`consistent` 项 → `simple.consistent()==false` → NS.H:118-127 的 SIMPLEC
修正分支不执行，`rAtU ≡ rAU`；因此 **mob = rAU = rAtU = primalPressureMobility**
（B23 修正基重建 2.19e-15 机器精度互证；生产 `alphaRel·prodRAU ≡ rAtU` 同证）。

**D3（H7 通道 = P 行通量切线缺失槽）**：真前向 phi 的 p 切线 =
`−kf·(dp_n−dp_o) − interp(mob·∇dp)·Sf`；生产 J_P 仅有 kf 槽 → H7 缺失
（B23 VERDICT_TABLE H7：corr(修正,残差) −0.996/−0.879/−0.948，正确符号
值门 3.01→2.60%，反号 3.80% 恶化；D2 45.34→35.68）。

**D4（∇dp = Gauss linear）**：案例 gradSchemes `default Gauss linear`；
边界：zeroGradient-p（壁/入口）dp_b = dp_c；fixedValue-p（出口 p=0）dp_b = 0
（B23 T7 同构；T7b grad_p 复刻）。

**D5（转置槽位）**：通道 Ch_f(dp) = −[w·mob_o·(∇dp)_o + (1−w)·mob_n·(∇dp)_n]·S_f；
P 行 = div，故 J^T 用散度权 λ_f = λP(o)−λP(n)（纯权，无 convTerm/alphaRel 混合，
mob 已含 alphaRel）：
- stage 1：g_c += mob_c·κ(c,f)·λ_f·S_f（κ = w 当 c=owner，否则 1−w）；
  出口（assignable-U）边界面 g_c += mob_c·λP(c)·S_b；
- stage 2（G^T 展开）：q_f' = S_f'·(g_o'/V_o' − g_n'/V_n')，
  out[P(o')] −= w'·q_f'，out[P(n')] −= (1−w')·q_f'；
  zeroGradient-p 边界面 out[P(c)] −= (g_c·S_b)/V_c。
对角元逐项复核：out[P(o')] 对 λ_f' 的系数 =
−w'²mob_o|S|²/V_o' + w'(1−w')mob_n|S|²/V_n' = ∂phi_f'/∂dp_{o'} ✓。

## 2. 实现清单（三处一致 + 折叠 + 预条件器口径）

1. `applyDiscreteFlowJ`（诊断前向）：P 行通量切线 deltaPhiFacePU 增加
   `−(Sf & interp(mob·∇dp))`；出口边界通量增加 `−(S_b & mob_c·∇dp_c)`。
2. `applyDiscreteFlowJT` / `applyProdFlowJT`（诊断/生产 J^T）：D5 转置
   stage1+stage2。
3. COO→CSR 显式导出（explicitJT.mtx 谱系）：H7 全部槽位的显式条目
   （distance-2 P-P 模板 + zeroGradient-p 边界自项 + assignable-U λP(c) 列）。
4. 外部面泛函折叠（b_TC 路由，B018 fix-3 同槽位镜像）：两个文件各加
   H7 转置（face functional 代 λ_f）。
5. 预条件器口径（审阅者 §3-2 注记的裁决）：**pressureGAMG Poisson 矩阵
   保持 kf-Laplacian 近似不变**——fvMatrix laplacian 是 distance-1 模板，
   不能表示 H7 的 distance-2 耦合；预条件器无需等于块本身。
   PRODPRECGAMGCHECK 四值参照系 = J_PP 的 Laplacian(kf) 部分（声明更新）；
   预期四值逐位不变（参照与矩阵都没动）。新增 PRODH7SHARE Info：确定性方向
   dp=sin(0.173(c+1)) 上 |H7 dp 通道|/|kf dp 通道| L2 比（预条件器近似度的
   量化读数）。

**不实现（本轮）**：H8 松弛源边界后段。理由：H8 的生产算子相关部分只有切线
（bnd_post·dU·mob/V），B23 T8 实测其切线与 D2 残差正交（corr 0.001-0.003，
D2 0.2216→0.2224 不动）；值级 2.60→2.00% 是离线仪器映射性质（生产 mobility
本为精确导出，无值门）。完整精确边界松弛源切线还牵涉 B23 未预注册的
preBnd-in-clamp 与 mob vs alphaRel·prodRAU 两项边界胞差。按「不确定则只做
H7」执行；H8 与 design-row（rxPressureRowTranspose 的 −interp(drAU·∇p)·Sf
设计导数，B18/B19 闭层禁改）同列为下轮候选。

## 3. 预注册预测与门（先于数据）

| # | 项 | 预测/门 |
|---|---|---|
| G1 | wmake WMAKE_EXIT=0；静态 14/14 | 必过 |
| G2 | 转置自洽（诊断路径）：12 随机点测试 max rel ≤1e-12（历史 ~2e-13）；BlockDot 五块 ≤1e-12 | 必过（新增项转置必须精确） |
| G3 | ExplicitJToracle（显式 vs matrix-free）maxRelL2 ≤1e-14（历史 3.8e-16） | 必过 |
| G4 | 等价性门四值：与 B22 逐位相同（1.131/0.826/1.422/1.343 e-16 型） | 必过（参照系未动） |
| G5 | 伴随双标签 trueRelRes ≤1e-12 + GRADSTABLE（tolerance 已设 1e-12；迭代数变化须可解释：J_PP 谱变化） | 必过 |
| G6 | FD 门三方向三目标：FD 侧与 B22 逐位相同（FD 不读算子）；ADJ 列移动——预测方向：H7 改变 λ 求解 → ADJ_J/ADJ_gDP 移动量级 ~%；J 符号表与 gDP 因子如实报告（不设硬性通过线；若 J(D1) 仍反号给归因分解） | 如实 |
| G7 | 残差-FD 锚：Gate PR/H3/H4（U-only 方向）H7 通道恒为零 → 数字不变 | 必过（核心块未回归） |
| G8 | 离线重估（修正仪器 + H7）：面锚 D1/D2/D3 复现 T7b dp-扩展型读数（D1≈2.7%/D2≈35.7%/D3≈2.5%，du+kf+H7 无 da） | 预测 |

## 4. 运行清单

- Run B（b24_diag，stageB4JacobianProbe=true，stageB2 关）：点测试 +
  ExplicitJToracle + CSR 导出（explicitJT_H7.mtx）。
- Run A（b24_gauge，b22_gauge 同构 + tolerance 1e-12 + b18rhs 导出）：
  等价性门 + 双标签 + FD 门 + wstate 导出（FD 侧应与 B22 逐位相同）。
- 离线：修正仪器（b23 修正基）+ H7 通道的面锚复算与 gDP/J(D1) 归因分解。

## 5. 纪律声明

- 零拟合因子；无阈值/开关/链改动；B23 证据只附加不修改；
- git 选择性提交（H7 单独提交），不 push；
- 任一门失败即停并定位报告。
