# BFINAL-027 — M2/M3 测量轮（a/b 段判别）— 预注册

Repo: /home/ys/dsH/TO-ANISOTROPIC, branch agent/dsH-stage-b-validation (HEAD 88ee077)
Date: 2026-08-21 (run cycle-1)
Mode: 纯离线 Python; 零编译; 零求解器运行; 零生产代码改动; 不 git 提交

## 0. 任务定位

M1 (BFINAL-026 cycle-1) 已把病灶区间锁定在 λ_T solve / thermalC contraction (a/b 段):
- D3 (thermalCoupling source) SOURCE-CLEAN: lhs=+1.1566617519e-02 vs rhs=+1.0940298579e-02, reldev 0.057, verdict SOURCE-CLEAN
- D1/D2 O(1)-MISMATCH (T-elimination folding c 段), 与本轮无关
- B26 对 D3 的临时判词 "lesion in lambda_T solve / thermalC contraction a/b" — 本轮 M2/M3 正是要在此二段之间做判别

## 1. 判别对象（四象限）

- 段 a = A_T^T λ_ref = b_Q 生产热伴随求解（λ_ref 生产值 = Tb 场）
- 段 b = thermalC 收缩（sensitivity.H 公式由 Tb 场 → thermalCoupling C 场 → L2/投影）

判别矩阵：四象限 a-guilty/innocent × b-guilty/innocent，用 M2 与 M3 两个独立测量互相锁定。

## 2. M2 测量（判 a 段）

方法: 离线装配 A_T (b18 模板 L661-764 + 已确认 DEFECT 修正) → 直接解 A_T^T λ_ref = b_Q
  (33600 未知量, scipy splu, 后台运行, 预计秒级) → 与生产 Tb 场比较:
  - relL2 = ||λ_ref − Tb||_2 / ||Tb||_2
  - 符号翻转计数 (sign-flip count)
  - 边界带分段误差: outlet(84) / inlet(84) / hotInlet(196) / hotOutlet(196) /
    solidEndWalls(280) / bottomWall(1120) / topWall(1120) / sideWalls(4800) 逐 patch 的
    banded relL2 与 banded max-abs-err

b_Q 公式 (生产 snapshot 语义, computeObjective.H L86-103):
  dJ_Q/dT = −phi_f/(Tref·M_frozen), 仅 cold outlet 84 单元非零
  M_frozen = state-B sum(phi_out) = 4.0733699999641404e-03
  Tref = 600.0
  已验证: 我的 b_Q 与生产导出 wstate_dJdT.mtx 的 relL2 = 4.6e-10 (phiThermal == phi on outlet)

A_T 装配: div(phiThermal,T) bounded upwind − laplacian(DTEffective,T) 离散语义
  (b18 模板 + 边界系数: outflow conv diag+=phi_b / inflow src−=phi_b·T_b /
  fixedValue-T diag+=kfTb src+=kfTb·T_bvals / zeroGradient 无贡献; laplacian 方向投影
  DTEffective 6 分量)

## 3. 关键开放性项目: A_T 符号约定

b18 L685-686 注释 ("A += -lap: diag -= kfT, off += kfT") 与代码 (diag += kfT;
A_off_own -= kfT; A_off_nei -= kfT) 相反; b18_folding_lab.log 以 SyntaxError 收尾,
part-8 门从未跑通 → A_T 符号约定必须由本仪器经验仲裁:
  - 门 (a): A_T·T − src − Q ≈ 0 (Q=state-B uniform 0, T=state-B T 场, src 含边界项)
  - 门 (b): A_T^T·Tb − dJdT ≈ 0 (Tb=生产伴随场, dJdT=生产 wstate_dJdT.mtx)
两个门必须同时通过 (残差 L2 相对量级 ≪ 0.1); 若 (a)(b) 相互矛盾 → 报告 BLOCKED 型发现,
不臆断。门通过后其符号即被锁定, 用于 M2 主测量。

## 4. M3 测量（判 b 段）

方法: 由生产 Tb 场 (非 λ_ref) 按 AdjHeatTransfer.H L189-327 公式重算 thermalC 场
  - 内部面: C[own] -= w·dDnDxh_own·(Tb_own−Tb_nei)(T_own−T_nei)·deltaCoeffs·magSf;
            C[nei] -= (1−w)·dDnDxh_nei·same_base
  - 边界: C[cell] += dDnDxh·Tb_cell·snGrad(T)·magSf
    (snGrad 经 mixed BC valueFraction; zeroGradient 无贡献)
  再算 L2 与投影, 对照:
  (1) B19 分解锚点 |thermalC|L2 = 1.1131e-3 (state-B; 主循环态 3.5092e-3/4.3447e-4/
      2.8471e-4/1.1131e-3)
  (2) M1 实测 rhs 链 (D3 rhs=+1.0940298579e-02 等, 若 C 场投影可复现)

注意 M3 用生产 Tb 而非 λ_ref → 它测量的是 "thermalC 收缩" (段 b) 而非 λ_T solve (段 a)。
M3 与 M2 正交。

## 5. 四象限判别矩阵与预期形态

| 象限 | M2 (a 段, 离线 A_T^T λ vs Tb) | M3 (b 段, Tb→thermalC vs B19 0.001113) | 结论 |
|------|-------------------------------|----------------------------------------|------|
| a guilty / b innocent | relL2 O(1) 或结构符号差; 边界带结构性错 | L2 ≈ 1.1131e-3 (1e-2 内), 场形态吻合 | a 段缺陷: 生产 Tb solve 或 A_T 装配病灶 |
| a innocent / b guilty | relL2 ≲ 1e-2 (离线精度界, 由 b18 模板 0.5% 边界项约束) | L2 与 1.1131e-3 偏差 > O(1e-2), 场形态异 | b 段缺陷: thermalC 收缩病灶 |
| a guilty / b guilty | relL2 O(1) | L2 偏差大 | 双段缺陷 (需再细分, 报告诚实的上限结论) |
| a innocent / b innocent | relL2 ≲ 1e-2 | L2 ≈ 1.1131e-3 | 双清白 → D3 的 28× over-projection 源必须向上复审 (报告 B19/B26 分解锚点或目标导数结构问题) |

离线精度界 ~1e-2 的依据:
- b18 模板块回归系数 0.995/1.003/1.001 (hA=0.9947/direct=1.0032/kf=1.0010, 重建残差 6.2%)
- 生产 dJdT 本身 relL2 4.6e-10 级; b_Q 不是误差源
- 若 M2 relL2 落 (1e-2, 0.3) 灰色带 → 报告为 "弱判别", 附带逐 patch 分段误差定位
  最可疑 patch, 不下强判词

## 6. 通过判据 (预先承诺, 不事后修改)

- a 清白: M2 relL2 ≤ 1e-2 且无结构性符号翻转 (非零翻转仅孤立点) 且边界带无 O(1) 带状误差
- a 有罪: relL2 > 1e-2 (至少比离线精度界大一个量级) 或结构性符号差 (整 patch/整带翻转)
- b 清白: M3 重算 L2 与 1.1131e-3 的 reldev ≤ 1e-2 且投影方向吻合 (符号一致)
- b 有罪: reldev > 1e-2 (量级差异 ≥ O(1e-2)) 或符号相反
- 灰带 (relL2 1e-2..0.3) → 弱判别: 以边界带分布为次判据, 报告诚实

## 7. 已知 DEFECT 清单 (b18 模板修正, 用于 A_T 装配)

- 网格质心: 精确金字塔质心 (primitiveMeshCellCentresAndVols.C) 而非顶点均值 (b16 DEFECT)
- 边界系数: physical patches internalCoeffs/boundaryCoeffs = 0 (transpose 语义);
  coupled patches 交换; 本任务中 A_T 直接装配含边界行, 需按 OpenFOAM fvm 语义补边界贡献
- bounded Sp: diag += −divphiT (bounded upwind 的 fvm::div 隐式源)
- 生产 Tb solve 收敛证据 (已从 Log.verify_b25_qgate.txt L6810 确认):
  PBiCGStab/DILU tol 1e-12 relTol 0 maxIter 2000; state-B 末次 Tb final residual 3.68e-13
  (首迭代 6.58e-13), 66 次迭代; transpose dot-test 1.10e-13; 前向 T final 9.30e-13 (37 次)
  → 生产 Tb solve 本身收敛到机器精度, 不引入 1e-2 级误差 (a 段病灶不在迭代收敛,
  若有病灶则在 A_T 装配或求解器离散语义, 由 M2 的离线重建对照判定)

## 8. 运行纪律

- 所有中间数落 tsv; 每个步骤后追加 NOTEBOOK.md
- M2 splu 33600×33600 后台运行, 预计秒级; 若 >15 min 单步卡死 → NOTEBOOK 记录,
  不擅自换路径
- 数据齐即强制收尾: EXECUTOR_SUMMARY.md + 空 EXECUTOR_DONE, 不 git 提交
