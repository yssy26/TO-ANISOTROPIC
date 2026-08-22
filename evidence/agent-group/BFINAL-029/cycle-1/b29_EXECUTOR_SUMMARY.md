# BFINAL-029 — R1 修复轮：候选排除（执行器总结）

> 收尾提交（2026-08-22）。任务书 `NEXT_TASK_BFINAL029.md`，绑定判据原样引用 `b29_PREREGISTRATION.md` §4.2。
> 全部工作纯离线 Python：零编译、零求解器运行、零生产代码改动（`src/` 未触碰，Phase 3 未触发）。
> 阶段结论已写入 `b29_NOTEBOOK.md`（§5 Phase 1 / §6 判据解释 / §7 Phase 2）。

## 1. 结论一句话

**候选排除**：H7s2、T4 与 Phase 1 发散槽（T3/T6/T2/T1）的**全部** 零/翻号 单变体（15 个）与 双变体 交叉（42 个）均无法使 M1 恒等式同时满足预注册定罪判据（D1∈[0.85,1.15] 且 D3∈[0.90,1.20]）。
按预注册 §4.2 排除规则：**不进入 Phase 3**，生产源码保持不动；剩余 M1 D1/D3 缺口指向输入级（Tb@B2 缺失）与 g1@B2 2.8% 校准差，如实记录为已知限制。

## 2. 输入核验（全过，详见 NOTEBOOK §1）

| 项 | 值 |
|---|---|
| 网格 | N=33600, nF=104740, nIF=96860, nB=7880 |
| 导出 b_TC | `/home/ys/dsH/b25_qgate/b18rhs_thermalCoupling.mtx` (33600×4) |
| w_true / rhs | B26 冻结件 `b26_wstar_h7.npz` + `b26a_m1_results.json`（无 FD 重算） |
| 场 | `b25_qgate/1` 全 12 项（含 Tb、phiThermal、coldFaceMask） |
| 目标 | maximizeTotalHeatTransfer，TIn 硬编码=600.0（=生产 gAverage 精确值） |

## 3. Phase 1：门全过 + 全失配分解（数据：`b29_slot_table.json`）

| gate | 值 | 判定 |
|---|---|---|
| self-gate [g1_Tmix, T8 off] vs B26 route_V1 | relL2=0.39701093640754337 | 逐位匹配 PASS |
| additivity（逐槽贡献求和=全开） | max abs 3.4e-21 / rel 7.0e-17 | PASS |
| mobility gate | 2.19e-15 | PASS |
| V gate（精确金字塔体积） | 2.72e-14 | PASS |
| T8 贡献（理论零：各向同性动量算子） | 9.6e-18 | PASS |

全开 [g1_TIn, T8 on] relL2=0.38617 vs 导出。leave-one-out 分解：
**T3 Δ=+0.5867（范数份额 93.5%，主导）；T6 +0.0266（27.6%）；T2 +0.0115（3.5%）；T1 +0.0114（3.3%）**；
T5 +1.6e-4、T4 +5.4e-6、H7s1 +8.6e-7、H7s2 +6.3e-7、T7 +3.3e-7、H7s1b −3.6e-7、H7s2b −1.4e-8、T8≈0 全 innocent。

## 4. 三个假设逐一证伪（Phase 1，详见 NOTEBOOK §5.4–5.6）

1. **尺度/口径假设 → 证伪**：导出端 prodRhs=raw/scale、写回乘 scale → 精确相消，导出即 raw 生产 rhs；U 1.82× / P 1.18× 振幅缺口为真实缺口。
2. **g1 假设 → 证伪**：生产 `computeObjective.H` TIn=gAverage(T[coldInlet])=600.0 与 instrument 硬编码精确相等，M_frozen=4.07337e-3 匹配；g1 非 D1 符号源。
3. **状态假设 → 证伪（决定性）**：M1 状态探针（8 状态）无任何状态复现导出 M1 行；D1 在所有 route 状态为正（+0.17..+0.47）vs 导出 −0.656，D3 为 −1.02..+0.20 vs 导出 +1.057。

## 5. Phase 2：变体矩阵（15）+ 全双切换交叉（42）→ 无一定罪（数据：`b29_variant_matrix.json` + `b29_crossscan.json`）

判据（预注册 §4.2 绑定，rhs 冻结）：**D1_v∈[0.85,1.15] 且 D3_v∈[0.90,1.20]**。

| 关键变体 | D1 | D3 | D1ok | D3ok |
|---|---|---|---|---|
| V0_base（route 全开） | +0.4665 | +0.2049 | n | n |
| **Z_T4**（唯一单变体 D1 入带） | **+0.8899** | **−0.5607** | **Y** | **n** |
| Z_H7s2 | +1.8617 | −24.8031 | n | n |
| Z_T3 / F_T3 | +4.429 / +8.391 | −14.50 / −29.20 | n | n |
| Z_T6 / F_T6 | −4.486 / −9.438 | +17.86 / +35.51 | n | n |

全双切换交叉（嫌疑槽∪H7s2b 上 C(7,2)=21 对 × {Z,F}=42 变体，穷举）：
仅 2 个达 D1ok —— Z_T4+T2（D1=+0.8870，**D3=−1.9764**）、Z_T4+T1（D1=+0.9909，**D3=−1.1434**）；
其余 40 个至少一判据出带（幅度常达 O(10..100)）。

**不可达成结构**：D3 方向被近相消对主导（H7s2 D3=+0.274 vs H7s2b −0.260；T3 +0.161 vs T6 −0.193），
槽零/翻号自由度的 D3 移动步长仅 O(0.1..0.6)，而 D3 需 ~5 倍放大（route lhs=+0.00224 → 目标 ~+0.01094）；
**D1 入带与 D3 入带在所有零/翻号组合下互斥**。

## 6. 已知限制（如实记录，不拟合）

- **g1@B2 校准差 2.8%**（离线 111.855 vs b25 log 108.816）：B2 无 g1 源快照，不可离线判定是否 D3 缺口来源。
- **Tb@B2 缺失**：Stage B2 全盘无 Tb 导出；Tb 是唯一不可状态切换的 g0 输入（adjDown=Tb[nei|own]，经 T3 占 b 范数 93.5%）。需 B2 运行端补导出才能最终定位 D1/D3 缺口。
- rhs 为 B26 冻结件；若 B26 rhs 与 B2 共同携带 Tb/g1 输入差，缺口被共同携带。

## 7. 产出文件（`/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-029/cycle-1/`）

| 文件 | 内容 |
|---|---|
| `b29_NOTEBOOK.md` | 全程磁盘笔记本（§1-7） |
| `b29_PREREGISTRATION.md` | 预注册 + 绑定判据 §4.2 |
| `b29_slot_instrument.py` / `b29_slot_table.json` | Phase 1 逐槽分解 |
| `b29_state_test.py` / `b29_state_test.json` | 状态重建 |
| `b29_m1_state_probe.py` / `b29_m1_state_probe.json` | M1 状态探针（决定性） |
| `b29_phase2_variants.py` / `b29_variant_matrix.json` | Phase 2 变体矩阵（15） |
| `b29_phase2_crossscan.py` / `b29_crossscan.json` | Phase 2 全双切换交叉（42） |
| `b29_EXECUTOR_SUMMARY.md` / `b29_EXECUTOR_DONE` | 本总结 + 收尾标记 |

## 8. 偏差声明

无。按任务书三阶段执行：Phase 1（离线分解）→ Phase 2（变体×判据）→ 排除 → STOP；Phase 3 未触发，生产源码未改动。
