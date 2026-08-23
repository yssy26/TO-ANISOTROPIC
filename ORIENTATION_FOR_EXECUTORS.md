# 执行者定向包（ORIENTATION_FOR_EXECUTORS.md）

> 每轮执行 agent 开工时**先读本文件**（约 3 分钟），替代对仓库/历史/惯例的自行探索。
> 维护人：协调人（审阅者）。执行者发现过期内容：在 EXECUTOR_SUMMARY 里注明，不要自行改动本文件。
> v1（2026-08-21，BFINAL-026 轮起生效）

## 0. 上下文纪律（最重要——本文件存在的原因）

你的上下文有限且会被压缩。遵守以下规则可把压缩损失降到近零：

1. **磁盘笔记本**：开工即在证据目录建 `NOTEBOOK.md`，每完成一小步就追加（做了什么/关键数字/下一步/未决问题）。压缩后恢复 = 重读本文件 + NOTEBOOK，而不是重新探索。
2. **抽取式读日志**：**永远不要 Read 完整运行日志**（8 万行级）。用 `grep -n` 定位 → `sed -n 'X,Yp'` 抽取所需片段 → 直接存为证据文件。归档全日志用 `cp` + `sha256sum`（不经过你的上下文）。
3. **复用已归档数据**：FD 真值、门值、历史数字优先从 `evidence/agent-group/BFINAL-*/cycle-1/*.tsv` 取，**不要重跑求解器来重新获得已有数字**。
4. **一次编译原则**：单 TU 编译 30–35 分钟。把所有代码改动（含开关门控导出）一次性完成再 wmake；编译期间做离线脚本/文档工作。
5. **长跑管理**：求解器运行用后台 + 定时轮询（间隔 ≥60s），不要忙等。

## 1. 环境 / 编译 / 运行（三条铁律 + 现成命令块）

```bash
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"  # 干净 PATH 先行
source /opt/openfoam7/etc/bashrc
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
export PATH=/home/ys/dsH/tools/ccache-shim:$PATH   # ccache（g++/c++ shim；2026-08-22 验收：三次重建 sha 逐位一致，热重建 1s）
unset FOAM_SIGFPE        # 不 unset 则伴随 NaN 直接 abort
# 不要 set -e / set -u；不要用 tee（SIGPIPE 会带走 wmake）
cd /home/ys/dsH/TO-ANISOTROPIC/src
wclean > /dev/null 2>&1
wmake > ../wmake_current.log 2>&1; echo "WMAKE_EXIT=$?" >> ../wmake_current.log
```

ccache 说明：源码未变的 TU 走缓存（秒级）；改动过的 TU 照常全编（MTO_HF.C 单 TU ~40 分钟）。
版本切换（checkout 回旧版重建）会命中缓存。缓存统计：`/home/ys/dsH/tools/ccache/bin/ccache -s`。

## 2. 每轮标准流程（纪律清单）

PREREGISTRATION.md（含预期数字，先于数据）→ 实现/测量 → 验证 → 证据入
`evidence/agent-group/BFINAL-0XX/cycle-1/`（追加式）→ EXECUTOR_SUMMARY.md +
空文件 EXECUTOR_DONE → 一次选择性 git 提交（长信息风格见 `git log` 上一轮）
→ **绝不 push**。零拟合因子；失败如实报；禁改清单见当轮任务书。

## 3. 热点文件地图（省去全局搜索）

| 位置 | 内容 |
|---|---|
| `src/computeObjective.H` | dJ/dT 与 dJ/dphi（三种目标类型分支；M_frozen 快照在此顶部） |
| `src/costfunction.H` | 目标值 J（分支同上） |
| `src/sensitivity.H` L110-280 | BFINAL-019 J 装配分解：momentum/pressureRow/fluxDirect/thermalC |
| `src/validateStageB2GradientAmplitude.H` | B2 FD 门（tsv 产出、判据、17 列） |
| `src/solveDiscreteFlowAdjointProduction.H` | 生产伴随 + pressureGAMG + b_TC 折叠（H7 在此） |
| `src/solveDiscreteFlowAdjoint.H` | 诊断 oracle（显式 CSR/COO、点测试） |
| `src/AdjHeatTransfer.H` | 热伴随 λ_T（源 = thermalObjectiveDerivative） |
| `src/validateFrozenGradient.H` | Stage F——**结构性不触发**（L17 守门 + 量具算例 gradientValidated=true），别浪费时间找它的输出 |
| `src/Make/files` | 编译目标清单（MTO_HF.C 为单一大 TU） |

## 4. 关键数字速查（B24=f0034c0 / B25=08ff950，同态）

| 量 | D1 / D2 / D3 |
|---|---|
| ADJ_J（两轮同，1e-10 内） | −0.0358 / −0.0249 / −0.0120 |
| FD_J（B25=Q 型） | **+0.0425** / −0.0036 / −0.00043（D1 反号；ADJ 超配 ~6.7×/28×@D2/D3） |
| gDP 因子 ADJ/FD | 2.061 / 2.007 / 2.100（均匀，已知状态） |
| gV relV | 1e-5..1e-7（机器级，过） |
| state-B 装配分解 L2 | momentum 0.002267 / pressureRow 0.000409 / fluxDirect 0.000323 / thermalC 0.001113 |
| 双标签伴随 trueRelRes | ~1e-12 级 + GRADSTABLE 4/4（tol 1e-12） |

## 5. 已知陷阱（踩过一次的坑）

1. **旧符号离线基**：B23 之前的一切离线装配读数被污染（扩散符号翻转）。离线基必须用 B23 修正版。
2. **b14/b24_diag 模板**的 legacy-ILU 诊断求解在后期 NaN→abort（同位同因）——交付物要在该点前完成或避开该路径。
3. **tee/SIGPIPE、PATH 污染、set -e、FOAM_SIGFPE**：见 §1 铁律。
4. Stage F / F1 不触发（见 §3），别找它的日志。
5. FD 重收敛迭代数与目标泛函相关（外层收敛判据含 J 历史），跨轮比较 NIter 差异是良性的。

## 6. 可复用工件清单（调用优先于重写）

| 工件 | 用途 |
|---|---|
| `evidence/agent-group/BFINAL-024/cycle-1/b24_reattribution.py` | 归因/对照表生成模板 |
| `evidence/agent-group/BFINAL-025/cycle-1/` 各 tsv | FD 真值与门值（勿重跑） |
| `stageB18RhsExport`（代码开关） | b_TC/rhs 导出（mtx） |
| B23 修正基构造 | 见 `BFINAL-023/cycle-1/FINAL_REPORT.md` 的 W-B 节 |
| 各轮 PREREGISTRATION.md | 预注册格式模板 |
