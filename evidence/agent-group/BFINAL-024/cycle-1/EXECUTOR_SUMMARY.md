# BFINAL-024 cycle-1 — EXECUTOR_SUMMARY

任务：H7 实现轮（全案第一个被证实的生产算子缺陷的修复）。分支
`agent/dsH-stage-b-validation`，起点 `2aca224`。执行 2026-08-20。
授权：本派发（用户已授权算子层修复；范围 = 两个伴随算子文件的装配、
CSR 导出、预条件器一致性）。

## 一行裁决

**H7（phiHbyA 压力通道）已按 B23 推导条目逐槽进入生产 J^T、诊断 J/J^T、
b_TC 折叠与显式 CSR（四处一致）；转置机器精度（点测试 6.4e-14 / oracle
5.5e-16）；等价性门四值与 B22 逐位相同（预条件器保持 kf-Laplacian 近似，
口径声明更新 + PRODH7SHARE=0.702 量化）；双标签 ≤9.8e-13 + GRADSTABLE；
FD 门：gDP 因子 2.145/2.159/2.181 → 2.061/2.007/2.100（三方向一致右移），
J 符号表不变（D1 反）；修正仪器归因：面锚 D2 31.2→22.1%、D1/D3 4.9→
4.6-4.8%，残余 O(1) 方向均匀 → 不能归因于 P 行面级切线缺陷。**
H8 未实现（推导支持但 B23 实测切线正交/仪器值级性质，按「不确定则只做
H7」纪律，连同 design-row H7 一致性列为下轮候选）。

## 关键数字

| 门 | 值 | 对照 |
|---|---|---|
| wmake | WMAKE_EXIT=0 | — |
| 静态测试 | 14/14（新增第 14 项） | 13/13 (B22) |
| 随机点测试（12×）max rel | 6.409e-14 | 6.864e-14 (b14) |
| BlockDot PP | 1.533e-14 | 2.293e-14 (b14) |
| ExplicitJToracle maxRelL2 | 5.477e-16 | 3.757e-16 (b14)；CSR nnz 13.56M→15.09M（+1.53M H7 项） |
| GateH4 / GateOracle（U 通道 FD 锚） | 逐位同 b14 | H7 对 U-only 方向恒零（推导预期） |
| 等价性门四值 r1/r2 | 1.1314/0.8258/1.4221/1.3430 与 1.1379/1.0744/1.5306/1.3337 e-16 | **与 B22 逐位相同** |
| PRODH7SHARE | 0.702 / 0.701（r1/r2） | 预条件器排除份额 |
| 双标签 trueRelRes | TC 9.81e-13 / PD 9.48e-13（r1）；9.25e-13 / 9.69e-13（r2）+ GRADSTABLE | B22 9.9e-10（tol 1e-9→1e-12） |
| FGMRES 迭代 | 1123/1323（r1）1112/1102（r2） | B22 694/902/691/675（容差收紧+谱变，可解释） |
| FD 侧（J/gDP/gV 全列） | 与 B22 逐位相同 | FD 不读算子 |
| ADJ_gDP（h=1e-3） | −5.2303/1.0505/13.6034 | B22 −5.4441/1.1301/14.1277 |
| **gDP 因子 ADJ/FD** | **2.061/2.007/2.100** | **B22 2.145/2.159/2.181** |
| ADJ_J | −0.0358/−0.0249/−0.0120 | B22 −0.0349/−0.0233/−0.0140；符号表不变（D1 反，2/3） |
| ADJ_gV | 逐位同 B22 | 体积链无流算子 |
| 面锚（修正仪器，du+kf±H7） | D1 4.92→4.75%、**D2 31.17→22.12%**、D3 4.93→4.56% | B23 T8/T9 du+dp+press 22.12% 复现 |
| +da 参考 | D1 2.68%、D2 38.16%（da 有害）、D3 2.45% | B23 T9 A/B-form 复现 |
| P 行连续性残差 \|div(tangent)\| | −17.8%（D1）/−23.2%（D2）/−17.3%（D3） | 真值 ~1e-11 |

## Run B 的 core-dump 说明

b24_diag 复刻 b14 模板（stageB4JacobianProbe=true）在点测试/oracle/
GateH1-H4/PR/Rx-dot 全部通过后，于 legacy ILU 诊断求解（GMRES-diag
relRes=nan）→ fsensMeanT nan → MMA 数据门 abort（MTO_RC=134）——
与 b14 模板逐位同因（b14 同处同 nan 同 abort）。非 H7 回归；Run B 的全部
交付件（点测试、oracle、CSR 导出 explicitJT_H7.mtx 523MB）在崩溃前完成。

## 改动清单

- `src/solveDiscreteFlowAdjointProduction.H`：applyProdFlowJT H7 转置
  （stage1/stage2/边界）+ 外部面泛函折叠 H7 + PRODH7SHARE 诊断 +
  等价性门口径声明（参照系 = J_PP 的 Laplacian 部分）；
- `src/solveDiscreteFlowAdjoint.H`：applyDiscreteFlowJ 前向 H7（Gauss-linear
  梯度基 + deltaPhiFacePU + 出口边界）+ applyDiscreteFlowJT 转置 + 折叠
  H7 + COO/CSR 显式 H7 条目（distance-2 P-P + zeroGradient-p 自项）；
- `src/tests/test_stage_b_safety_gates.py`：第 14 项静态测试；
- `AI_AGENT_HANDOFF.md` + `OPENFOAM_USAGE_AND_FILES.md`：BFINAL-003/008
  里程碑正式升级声明（日期注记 + usage 文档）；
- 证据（本目录）：PREREGISTRATION.md（先于数据）、两个运行日志、
  b24_reattribution.{py,log,json}、三个 TSV、wmake 日志、sha256 清单、
  EXECUTOR_DONE（空）。

## 纪律

- 零拟合因子；禁改清单（装配/源/链/MMA/目标定义/阈值/锁定开关）零触碰
  （rxPressureRowTranspose.H 未动；B23 证据未动）；
- git 选择性提交不 push；
- 案例：b24_gauge（=b22_gauge 同构 + tolerance 1e-12 + b18rhs 导出）、
  b24_diag（=b14 同构，explicit 矩阵导出到新路径，不覆盖 b8 归档）。
