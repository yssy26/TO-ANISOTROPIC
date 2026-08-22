# BFINAL-027 — M2/M3 测量轮（a/b 段判别）— NOTEBOOK

Repo: /home/ys/dsH/TO-ANISOTROPIC, branch agent/dsH-stage-b-validation (HEAD 88ee077)
Start: 2026-08-21 (session resume after context rollover)
Budget: 150 min; 15-min single-step rule; data-ready forced wrap-up.

## Step 0 — 输入固定与生产语义钉住 (completed)

| 项 | 值 | 来源 |
|---|---|---|
| N cells | 33600 | b16 mesh |
| Tref | 600.0 | b25_qgate/constant/optProperties |
| M_frozen | 4.0733699999641404e-03 | wstate_outletMeta.mtx line1 |
| b_Q 公式 | dJ_Q/dT = −phi_f/(Tref·M_frozen), 仅 outlet 84 单元 | computeObjective.H L86-103 |
| b_Q 验证 | 与生产 wstate_dJdT.mtx relL2 = 4.6e-10 (nz=84, maxabs=3.712e-05, L2=1.952e-04) | 离线复算 vs 导出 |
| Tb solver | PBiCGStab/DILU tol 1e-12 relTol 0 maxIter 2000 | fvSolution L48-55 |
| Tb 收敛 (state B) | initial 0.088024 → final 3.68e-13, 66 iter | Log L6810 |
| Tb 收敛 (首迭代) | initial 1.0 → final 6.58e-13, 66 iter | Log L937 |
| 前向 T (state B) | initial 0.0010358 → final 9.30e-13, 37 iter | Log L6807-6808 |
| transpose dot-test | max rel err = 1.09689091013e-13 | Log L6809 (state B) / L936 (首迭代) |
| B19 thermalC 锚 | \|thermalC\|L2 = 1.1131e-3 (state-B; 主循环态 3.5092e-3/4.3447e-4/2.8471e-4/1.1131e-3) | BFINAL-019 FINAL_REPORT.md L105 |
| M1 D3 锚 | lhs=+1.1566617519e-02 rhs=+1.0940298579e-02 reldev 0.057 SOURCE-CLEAN; C=-1.9791957452e-02 Gx=+8.4067076831e-03 | b26a_m1_identity.tsv |
| div schemes | div(phi,T) bounded Gauss upwind; div(-phi,Tb) Gauss upwind | fvSchemes L33,36 |
| laplacian/snGrad | Gauss linear corrected / corrected | fvSchemes L47,54 |
| DTEffective 分量 | 6-comp 对称张量 (方向投影经 weights) | b18 模板 L661-764 |
| patch 尺寸 | inlet84 outlet84 hotInlet196 hotOutlet196 solidEndWalls280 bottomWall1120 topWall1120 sideWalls4800 (nB=7880) | b16 mesh owner/startFace |
| outlet 单元 | [5679,5759,...,12319] step 80 (84) | owner array |
| A_T 符号约定 | UNRESOLVED → 由门 (a) A_T·T−src−Q≈0 与 (b) A_T^T·Tb−dJdT≈0 经验仲裁 | b18 L685-686 注释 vs 代码矛盾 |

生产伴随求解收敛审计结论: Tb solve 已收敛到 1e-13 级 (tol 1e-12, relTol 0), transpose
dot-test 1e-13 → 生产 λ_T 是可信参考; a 段病灶 (若有) 在 A_T 离散语义而非迭代收敛。

## Step 1 — 预注册 + 仪器构建 (completed)

PREREGISTRATION.md 已写 (四象限判别矩阵/预期形态/通过判据/灰带规则)。
构建 b27_m23_instrument.py (复用 b26a 组件 + A_T 装配 + 门仲裁 + M2/M3 测量)。
关键修正 (vs b18 模板 L661-764, 已知 DEFECT):
- 精确金字塔质心 (V gate relL2=2.719e-14 通过)
- DTEffective 先线性插值再方向投影 (gaussLaplacianScheme 语义)
- bounded upwind Sp: diag += −divphiT (含边界 phi)
- 混合 BC bottomWall (vf=0.00175607, refValue=600, refGrad=0): diag += vf·kfTb,
  src += vf·kfTb·600 + (kfTb/bdel)·(1−vf)·rg(0)
- 门 (a) 归一化改为 vs 算子尺度 (Q 场 uniform 0 → |Q|L2=0 除零是 b18 DEFECT)
- 读 field 支持 `internalField uniform 0;` 单行形式
运行: b27_run.log, 全程 ~14s (零编译/零求解器运行/零生产代码改动/不 git)。

## Step 2 — M2: A_T^T λ_ref = b_Q 离线直解 vs 生产 Tb (completed)

门仲裁 (符号约定经验锁定):
- 门 (a) A_T·T − src − Q: maxabs=6.8415e-13, relL2=1.184e-11 (vs 算子尺度 1.80)
- 门 (b) A_T^T·Tb − dJdT: maxabs=1.0894e-16, relL2=1.909e-11 (vs |dJdT|=1.952e-04)
→ A_T 装配符号/边界语义与生产 fvm 语义一致到机器精度。

M2 主测量:
- splu 33600×33600 分解 7.8s, 求解 <0.1s; 残差 |A_T^T·λ_ref−b_Q|/|b_Q|=7.006e-15
- λ_ref vs 生产 Tb: relL2=1.149494e-12, cos=1.00000000, sign_flips=0 (of 33600)
- λ_ref 与 Tb 同为 max=−2.467e-38, min=−4.3518e-01, L2=4.744213e+01 (逐位一致)
- 边界带分段 (b27_m2_compare.tsv):

| region | n | banded relL2 | banded maxabs |
|---|---|---|---|
| INTERNAL | 26208 | 1.150e-12 | 1.603e-12 |
| inlet | 84 | 1.130e-12 | 1.184e-12 |
| outlet | 84 | 9.750e-13 | 8.719e-13 |
| hotInlet | 196 | 2.877e-12 | 1.012e-12 |
| hotOutlet | 196 | 1.761e-12 | 4.544e-14 |
| solidEndWalls | 280 | 1.985e-12 | 1.622e-12 |
| bottomWall | 1120 | 7.086e-13 | 5.185e-13 |
| topWall | 1120 | 1.379e-08 | 6.224e-30 |
| sideWalls | 4800 | 1.229e-12 | 1.513e-12 |

topWall 的 relL2=1.38e-8 为归一化假象: |Tb|max=1.589e-22 (物理 ~0), 而 maxabs=6.2e-30,
绝对误差仍为机器零 → 非结构误差。全部 9 带 sign_flips=0。

## Step 3 — M3: 生产 Tb → thermalC 重算 vs B19 0.001113 / M1 rhs (completed)

按 AdjHeatTransfer.H L189-327 (sensitivity.H 公式, mixed-BC snGrad) 由生产 Tb 场重算
thermalC (C 场), 对照 b27_m3.tsv:

| qty | 重算值 | 锚点 | 相对偏差 |
|---|---|---|---|
| C_L2 | 1.1130715107e-03 | B19 1.11307150907e-03 | 2.397e-12 |
| proj_D1 | −3.0341656277e-02 | M1 C −3.0341656277e-02 | 1.4e-12 |
| proj_D2 | −8.2108899666e-03 | M1 C −8.2108899666e-03 | 2.9e-12 |
| proj_D3 | −1.9791957452e-02 | M1 C −1.9791957452e-02 | 8.2e-12 |

→ 段 b (thermalC 收缩) 复现 B19 分解锚点与 M1 C 投影到 ~1e-12。

## Step 4 — 四象限判词 (completed)

M2 门 (b) relL2=1.909e-11 同时证明 A_T 装配与生产 transpose solve 均机器精度正确;
M2 λ_ref vs Tb relL2=1.15e-12, sign_flips=0, 全带 ≤3e-12 (topWall 为归一化假象);
M3 C_L2 与 3 个投影 reldev ≤1e-11。

四象限: **a innocent (λ_T solve 清白), b innocent (thermalC 收缩清白)**
判据核对: a 清白需 relL2 ≤ 1e-2 且无结构符号翻转且边界带无 O(1) 带状误差 —
  实测 1.15e-12, flips=0, 全带 ≤3e-12 → 远优于离线精度界 (~1e-2, b18 模板 0.5% 边界)。
  b 清白需 C_L2 reldev ≤ 1e-2 且投影符号一致 — 实测 2.4e-12, 3 投影全符号一致。

诚实后果 (向上复审): b 清白 + a 清白 → D3 的 28× over-projection
(lhs/rhs=1.057, SOURCE-CLEAN) 的病灶不在 λ_T solve (a) 也不在 thermalC 收缩 (b),
不在本轮的判别域内 → 源必在 T-elimination folding (c 段) 或更上游 (目标导数结构/B19
分解锚点), 交由上级任务复审 (BFINAL-028 建议路径)。

收尾: EXECUTOR_SUMMARY.md + 空 EXECUTOR_DONE 已写; 不 git 提交 (协调员复核后提交)。
