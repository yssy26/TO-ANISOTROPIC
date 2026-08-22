# BFINAL-027 — M2/M3 测量轮（a/b 段判别）— EXECUTOR SUMMARY

Repo: /home/ys/dsH/TO-ANISOTROPIC, branch agent/dsH-stage-b-validation (HEAD 88ee077)
Date: 2026-08-21, cycle-1. Mode: 纯离线 Python; 零编译; 零求解器运行; 零生产代码改动; 不 git 提交。
时间: ~14s 全流程 (b27_run.log)。数据齐 → 强制收尾。

## 一句话四象限判词

**a innocent (λ_T solve / A_T 装配) 且 b innocent (thermalC 收缩) — 双清白。**

判别矩阵: a-guilty? NO, b-guilty? NO → 落在四象限 "a innocent / b innocent" 格
(PREREGISTRATION.md 第 5 节第 4 行): 病灶不在 λ_T solve 也不在 thermalC 收缩,
D3 的 28× over-projection 源必须向上复审。

## M2 关键数字 (判 a 段, 离线 A_T^T λ_ref = b_Q vs 生产 Tb)

| 量 | 值 |
|---|---|
| 门 (a) A_T·T − src − Q relL2 (vs 算子尺度 1.80) | 1.184e-11 (maxabs 6.84e-13) |
| 门 (b) A_T^T·Tb − dJdT relL2 (vs \|dJdT\|=1.952e-04) | 1.909e-11 (maxabs 1.09e-16) |
| splu 分解 (33600×33600) | 7.8s; 求解 <0.1s; 残差 7.006e-15 |
| λ_ref vs 生产 Tb: relL2 | 1.149494e-12 |
| λ_ref vs 生产 Tb: cos | 1.00000000 |
| sign flips (of 33600) | 0 |
| 边界带 | INTERNAL 1.15e-12, inlet 1.13e-12, outlet 9.75e-13, hotInlet 2.88e-12, hotOutlet 1.76e-12, solidEndWalls 1.99e-12, bottomWall 7.09e-13, topWall 1.38e-8*, sideWalls 1.23e-12 (banded relL2, 全带 sign_flips=0) |

*topWall relL2=1.38e-8 为归一化假象: |Tb|max=1.589e-22 (物理 ~0) 而 maxabs=6.224e-30
= 机器零。非结构误差。

精度界论证 (清白必须有, 非"差不多"):
a 清白判据 = relL2 ≤ 1e-2 且无结构符号翻转且边界带无 O(1) 带状误差。
实测 relL2=1.15e-12 — 比离线精度界 (~1e-2, 由 b18 模板 0.5% 边界项/回归系数
0.995-1.003 约束) 低 10 个量级; sign_flips=0; 全带 maxabs ≤ 1.6e-12 (除 topWall
归一化假象)。同时门 (b)=1.9e-11 证明 A_T 装配符号/边界语义与生产 transpose 一致
到机器精度 → a 段 (生产 λ_T solve + A_T 离散语义) 清白。

## M3 vs 生产 (判 b 段, 生产 Tb → thermalC 重算)

| qty | 重算值 | 锚点 | 相对偏差 |
|---|---|---|---|
| C_L2 | 1.1130715107e-03 | B19 1.11307150907e-03 | 2.4e-12 |
| proj_D1 | −3.0341656277e-02 | M1 C −3.0341656277e-02 | 1.4e-12 |
| proj_D2 | −8.2108899666e-03 | M1 C −8.2108899666e-03 | 2.9e-12 |
| proj_D3 | −1.9791957452e-02 | M1 C −1.9791957452e-02 | 8.2e-12 |

b 清白判据 = C_L2 reldev ≤ 1e-2 且投影符号一致。实测 2.4e-12, 3 投影符号全一致
→ b 段 (thermalC 收缩) 清白。

## 诚实申报 (open items, 必须向上复审)

1. D3 的 28× over-projection (M1 rhs=+1.0940298579e-02, lhs=+1.1566617519e-02,
   lhs/rhs=1.057, SOURCE-CLEAN) 病灶不在本轮的判别域内: 已排除 λ_T solve (a,
   由门 (b) 1.9e-11 + λ_ref=Tb 1.15e-12) 与 thermalC 收缩 (b, 由 M3 1e-12) →
   源在 T-elimination folding (c 段) 或更上游 (目标导数结构/B19 分解锚点) →
   建议交 BFINAL-028 向上复审。
2. D1 (lhs/rhs=−0.656) / D2 (−10.08) 的 O(1)-MISMATCH 与 c 段有关, 不在本任务范围,
   仅重述 (b26a_m1_identity.tsv)。
3. b26a_m1_instrument.py 的 read_field 存在 latent uniform-读法缺陷 (Q 场
   `internalField uniform 0;` 单行), b26a 曾因其丢 uniform 值; 本轮 b27 已修复
   (regex `internalField\s+uniform\b`)。零生产代码改动约束下未回改 b26a, 仅记录。

## 产物清单 (evidence/agent-group/BFINAL-027/cycle-1/)

- b27_m23_instrument.py — 主仪器 (A_T 装配 + 门仲裁 + M2/M3)
- b27_run.log — 运行日志 (全程 ~14s)
- b27_gates.tsv — 门 (a) 1.184e-11 / 门 (b) 1.909e-11
- b27_m2_compare.tsv — λ_ref vs Tb 逐带 relL2/maxabs
- b27_m3.tsv — C_L2 + 3 投影 vs 锚点
- b27_at.npz / b27_lamref.npz — A_T 矩阵与 λ_ref 场存档 (可复算)
- NOTEBOOK.md — Step 0-4 全程记录
- PREREGISTRATION.md — 预注册 (四象限矩阵/判据, 未事后修改)
- EXECUTOR_DONE — 空标记

复现: `python3 b27_m23_instrument.py` (同目录, 需 numpy/scipy; 纯离线, 不触碰
生产目录, 不 git 提交)。
