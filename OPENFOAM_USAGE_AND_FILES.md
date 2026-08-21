# TO-ANISOTROPIC — OpenFOAM 使用说明与文件内容说明

> **用途：** 面向接手本仓库（`yssy26/TO-ANISOTROPIC`，分支 `agent/dsH-stage-b-validation`）代码修改、编译、运行与验证的协作者 / AI agent 的**操作手册**。本文档只讲「怎么编译、怎么跑、每个文件是什么、有哪些坑」；项目物理背景、阶段历史与当前任务请读仓库根目录的 `AI_AGENT_HANDOFF.md`（优先）以及 `src/validation/stage_b/` 下的阶段报告。
>
> 本文档信息以 **HEAD `6a0004b`**（2026-08-19，`Point README to current AI agent handoff`）为准。仓库在 `/home/ys/dsH/TO-ANISOTROPIC`（git 根 == 工作区）。若后续有更新的 commit，先 `git pull --ff-only` 再核对本文档是否有过期段落。

---

## 1. 环境速览

- 宿主：Windows + WSL2（Ubuntu 20.04.6 LTS），内存 62 GiB，磁盘充足。
- OpenFOAM：**官方 OpenFOAM Foundation v7**，安装在 `/opt/openfoam7`。
- 编译产物：`/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`（**本项目自建构建树**，`FOAM_USER_APPBIN` 指向这里）。
- 注意：`/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin/` 下**另有一个旧版 `MTO_HF`**（历史遗留，2026-08-15），不要用它。运行一律用绝对路径 `/home/ys/dsH/TO-ANISOTROPIC/build/bin/MTO_HF`。

### 1.1 环境变量（source bashrc 之后）

| 变量 | 值 |
|---|---|
| `WM_PROJECT_DIR` | `/opt/openfoam7` |
| `WM_PROJECT_VERSION` | `7` |
| `WM_OPTIONS` | `linux64GccDPInt32Opt`（64 位 / Gcc / 双精度 DP / 32 位标签 / Opt） |
| `FOAM_APPBIN` | `/opt/openfoam7/platforms/linux64GccDPInt32Opt/bin` |
| `FOAM_USER_APPBIN` | `/home/ys/OpenFOAM/ys-7/platforms/linux64GccDPInt32Opt/bin`（bashrc 默认值；**本项目编译/运行时要覆盖为 `/home/ys/dsH/TO-ANISOTROPIC/build/bin`**） |

### 1.2 环境加载的正确姿势（三条铁律）

1. **干净 PATH 优先**：先 `export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"` 再 `source /opt/openfoam7/etc/bashrc`。若 PATH 里残留 Windows 侧含空格的路径，`bashrc` 或 `wmake` 会静默失败（看起来没报错，二进制却没更新）。
2. **`unset FOAM_SIGFPE`**：关闭 OpenFOAM 浮点异常捕获，否则伴随求解中的 NaN/Inf 直接 abort。
3. **不要 `set -e`**：OpenFOAM 的 `bashrc` 有中间命令返回非零，`set -e` 会提前退出。

```bash
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
```

---

## 2. 编译（wmake）

```bash
cd /home/ys/dsH/TO-ANISOTROPIC
export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
source /opt/openfoam7/etc/bashrc
export FOAM_USER_APPBIN=/home/ys/dsH/TO-ANISOTROPIC/build/bin
unset FOAM_SIGFPE
cd src
wclean > /dev/null 2>&1
wmake > ../wmake_current.log 2>&1   # 直接重定向到文件，不要用 tee 管道
echo "WMAKE_EXIT=$?" >> ../wmake_current.log
```

### 2.1 编译要点

- **单入口大 TU**：`src/MTO_HF.C` 通过 `#include` 引入全部 `.H` 模块。改任何一个 `.H` 都会触发 `MTO_HF.C` 整文件重编译，**单次全量编译约 30–35 分钟 CPU 时间**（`cc1plus` 单进程约 109% CPU）。耐心等，不要中途杀。
- **不要用 `tee`**：本会话曾出现 `wmake 2>&1 | tee ../log` 时 tee 收到 SIGPIPE 导致 `WMAKE_EXIT=141`（后台作业 stdout 管道被捕获器关闭所致）。直接 `wmake > ../log 2>&1` 即可。
- 编译入口文件 `src/Make/files` 列出目标：`harmonicSymmTensor.C`、8 个伴随边界条件 `.C`、`MMA/MMA.C`、`MTO_HF.C`；`EXE = $(FOAM_USER_APPBIN)/MTO_HF`。
- `src/Make/options` 链接 `turbulenceModels`、`finiteVolume`、`meshTools`、`fvOptions` 等，带 `-DOMPI_SKIP_MPICXX`。
- 验证二进制含新代码：`strings build/bin/MTO_HF | grep <新标记>`（如 `PRODPRECGAMGSETUP`、`pressureGAMG`）。

### 2.2 编译前静态测试

```bash
cd /home/ys/dsH/TO-ANISOTROPIC
python3 src/tests/test_stage_b_safety_gates.py   # 期望 Ran 5 tests OK
```

CI（`.github/workflows/ci.yml`）在 GitHub-hosted runner 上只跑纯 Python 测试 + 仓库卫生检查（无 >50MB 跟踪文件等）；OpenFOAM 构建无法在 hosted runner 上跑，必须在本地。

---

## 3. 运行

### 3.1 运行器脚本

仓库根目录有封装好的运行器，**推荐直接使用**：

- `/home/ys/dsH/run_verify_prod.sh <caseDir> <logPath>` —— 生产路径验证运行（stageB4/B5/B6/B8 全关、走 `solveDiscreteFlowAdjointProduction.H`）。内部已用绝对路径调用 `build/bin/MTO_HF`，不会误用旧二进制。
- `/home/ys/dsH/run_verify_diag.sh <caseDir> <logPath>` —— 诊断路径验证运行（stageB4/B5/B6/B8 探针开启、export-only）。
- `/home/ys/dsH/run_mtohf.sh <caseDir>` —— 通用运行器（注意：它不覆盖 `FOAM_USER_APPBIN`，若直接用它可能解析到旧二进制；优先用上面两个 verify 脚本）。

生产路径调用（示例）：

```bash
bash /home/ys/dsH/run_verify_prod.sh /home/ys/dsH/b8_verify_prod /home/ys/dsH/b8_verify_prod/Log.verify_prod.txt
```

### 3.2 算例目录（都是独立副本，勿动真实算例）

| 目录 | 用途 |
|---|---|
| `/home/ys/b2_case_smoke` | **真实算例（只读参考）**：含 `explicitJT.mtx`、`explicitSol_*`、历史日志。不要改它。 |
| `/home/ys/dsH/b2_case_smoke` | 可写副本（诊断路径基础，stageB4/B5/B6 开、export-only） |
| `/home/ys/dsH/b8_verify_diag` | 诊断路径验证（stageB8 探针开） |
| `/home/ys/dsH/b8_verify_prod` | 生产路径验证（探针全关、`discreteProdPreconditionerSetup diagonal`） |
| `/home/ys/dsH/b8_verify_gamg` | 生产路径 GAMG 预条件器验证（`pressureGAMG`、1e-3/50）—— **BFINAL-009 停止点所在** |

新建验证算例的推荐姿势：`cp -a <基础算例> <新目录>`，然后只改 `constant/optProperties` 里需要改的开关，并在运行前 `rm -f` 旧的 `*.mtx` / `Log*` / `optimization_history.csv` / `optimization_log.dat` / `adjointCheckpoint_*.tsv`。

### 3.3 运行耗时

- 诊断路径（stageB4 export + 探针）：约 8–9 分钟。
- 生产路径（FGMRES 4000 迭代）：thermalCoupling + pressureDrop 合计约 30–50 分钟（取决于预条件器；GAMG 内层每步都解，可能更久）。
- 离线 Python 直接求解（scipy splu 134400² 稀疏矩阵）：约 11–15 分钟。

---

## 4. 算例结构与关键开关

算例目录结构：

```
<case>/
  0/         初场与所有场（U,p,T,alpha,x,伴随场 Uc/pc/Ua/pa,梯度场 gsensPressureDrop/gsensVol,掩码 designMask 等）
  constant/  optProperties（优化与验收开关总表）、polyMesh、thermalProperties/transportProperties/turbulenceProperties
  system/    blockMeshDict、controlDict、fvSchemes、fvSolution、decomposeParDict、topoSetDict
  根目录     运行产物：explicitJT.mtx、explicitRhs_*.mtx、explicitSol_*.mtx、log.*、adjointCheckpoint_*.tsv
```

网格：`blockMeshDict` 矩形多层板，`convertToMeters 0.001`，尺寸 40×7×14.8 mm，四层结构（底部固定固体 / 冷侧设计区 / 隔板 / 热流体区），**33600 单元**。离散伴随未知量 4×33600 = **134400**（每单元 3 个速度分量 + 1 压力）。

### 4.1 optProperties 关键开关（按功能分组）

**优化/验收模式（当前必须保持）**
```text
mmaUpdateEnabled                  false;   // MMA 锁定，不更新设计
frozenGradientValidated           false;   // 冻结梯度未解锁（B-final 门）
gradientValidated                 true;    // 跳过旧版梯度验证层
ransDirectionValidated            true;
solveFlowAdjoints                 true;
adjointMode                       discrete;
flowModel                         incompressibleRANSFrozen;
frozenTurbulenceAdjoint           true;
```

**Stage-B 探针（诊断路径用，生产路径全关）**
```text
stageBEnabled false;  stageB2Enabled false;
stageB4JacobianProbe        true/false   // 诊断/oracle 路径（与 stageB2Enabled 互斥，同时开会 FatalError）
stageB5BoundaryRelaxOracle  true/false
stageB6RxDesignOracle       true/false
stageB8JPPActualResidualFD  true/false   // BFINAL-008 实际残差 Jv-FD 探针
```

**离散伴随求解（生产路径）**
```text
discreteFlowAdjointTolerance 1e-9;
discreteFlowAdjointRestart   80;
discreteFlowAdjointMaxIter   4000;
discreteProdSolverType       fgmres;             // 或 bicgstab
discreteProdPreconditionerSetup pressureGAMG;    // pressureGAMG | diagonal | bruteForceL1
discreteProdPressurePrecTolerance 1e-3;          // GAMG 内层容差
discreteProdPressurePrecMaxIter   50;            // GAMG 内层最大迭代
discreteProdEnableRitzPilot false;
discreteProdConvergeFatal    false;              // 不收敛时继续（DIAGNOSTIC MODE）
```

**BFINAL-011 诊断/修复开关（默认全部=基线，生产无需设置）**
```text
solveThermalCouplingFlowAdjoint  true/false   // 标签隔离（诊断用；跳过则不写 Ub/pb）
solvePressureDropFlowAdjoint     true/false   // 标签隔离（诊断用；跳过则不写 Uc/pc）
discreteProdPressurePrecSolver   GAMG;        // GAMG（默认）| PCG（PCG+DIC 备选，非法值 FatalError）
discreteProdPressurePrecScaleCorrection      true;   // GAMG scaleCorrection 覆盖（默认=OF7 对称矩阵默认）
discreteProdPressurePrecNPreSweeps           0;      // GAMG nPreSweeps 覆盖（默认=现状）
discreteProdPressurePrecDirectSolveCoarsest  false;  // GAMG directSolveCoarsest 覆盖（默认=OF7）
```

内层解后 psi 非有限会立即 `FatalError`（含标签、apply 序号；消息串 `non-finite psi`）——这是 BFINAL-011 加的快败仪表，取代「静默 nan 传播 → MMA 门才 abort」。

**pressureGAMG 的 laplacianSchemes 依赖（BFINAL-010 守卫）**：`fvm::laplacian(volScalarField, …)` 的系数面插值词取自 `system/fvSchemes` 中 **`laplacianSchemes`** 条目 ITstream（`Gauss <interp> <snGrad>` 的第二个词；`Gauss` 由 `laplacianScheme::New` 消费，见 OF7 `gaussLaplacianScheme.H`/`laplacianScheme.H:120-135`），**不走 `interpolationSchemes`**。算子 Jᵀ 的 kf 假定 `linear == mesh.weights()`，因此 pressureGAMG 模式下求解器会对 `default` 及任何 `laplacian(prodPressurePrecMobility…` 具名条目强制要求插值词为 `linear`，否则 FatalError（日志标记 `PRODPRECGAMGSCHEME`）。本算例 `default Gauss linear corrected` 满足。数值上矩阵/算子等价性由 `PRODPRECGAMGCHECK` 四组 relL2 硬门最终仲裁。

**oracle 导入/导出（诊断路径）**
```text
discreteExportOnly            true/false   // true=只导出矩阵不迭代求解
discreteUseExplicitSolution   true/false   // 导入外部 SuperLU 精确解
discreteExplicitMatrixFile    "/home/ys/dsH/<case>/explicitJT.mtx";
discreteExplicitRhsFile       "explicitRhs";
discreteExplicitSolutionFile  "/home/ys/b2_case_smoke/explicitSol";  // 真实算例的精确解（只读）
```

### 4.2 fvSolution 关键约定

- 压力 `(p|pa|pb|pc)` 用 `PCG + DIC`，**不是** GAMG（GAMG 聚合路径跨重启不确定，会破坏重复性验收；但生产预条件器内部显式用 GAMG 是另一回事，见 `solveDiscreteFlowAdjointProduction.H`）。
- 速度用 `smoothSolver + GaussSeidel`；温度用 `PBiCGStab + DILU`；过滤/梯度场用 `PCG + DIC`。
- `SIMPLE` 段 `pRefCell 5600`，但**当前 fixedValue 出口 p=0 使 `p.needReference()==false`，pRefCell 实际不生效**（BFINAL-006/007 已确认；唯一参考语义由 `p.needReference()` 门控，见 §5.5）。

---

## 5. src/ 文件内容说明（重点）

所有 `.H` 通过 `#include` 进入 `MTO_HF.C`。按职责分组：

### 5.1 前向主流程
| 文件 | 内容 |
|---|---|
| `MTO_HF.C` | 主程序：编排各 `.H` 块（NS → HeatTransfer → 伴随 → sensitivity → update） |
| `createFields.H` | 场创建与初始化 |
| `opt_initialization.H` | 优化状态与控制 |
| `readTransportProperties.H` / `readThermalProperties.H` | 物性读取 |
| `NS.H` | SIMPLE/RANS 前向流场 + 冻结湍流逻辑（含大量求解诊断输出） |
| `HeatTransfer.H` | 前向温度场 |
| `updateMaterialProperties.H` | Brinkman 阻力 / 热物性插值 |
| `filter_x.H` / `filter_chainrule.H` | Helmholtz 过滤 / Heaviside 投影 与 梯度链式法则 |
| `update.H` / `finalize.H` | 更新与收尾 |

### 5.2 伴随与梯度
| 文件 | 内容 |
|---|---|
| `AdjHeatTransfer.H` | 热伴随（换热目标） |
| `AdjNS_HT.H` / `AdjNS_PD.H` | 换热 / 压降流动伴随块（内部按 `stageB4JacobianProbe` 分流到诊断或生产伴随模块） |
| `AdjNS_FF.H` | 备用伴随入口 |
| `sensitivity.H` | **主梯度装配（含 BFINAL-005 R_x 压力行补丁 +106 行）；已锁定，勿随意改** |
| `rxPressureRowTranspose.H` | 压力/连续性行设计导数贡献（BFINAL-005） |
| `validateDiscreteObjectiveDerivatives.H` | g_w 单元测试（abort 阈值 1e-9） |
| `validateFrozenGradient.H` / `validateGradientChain.H` | 冻结梯度验证编排 |

### 5.3 离散伴随求解器（**当前焦点**）
| 文件 | 内容 |
|---|---|
| `solveDiscreteFlowAdjoint.H` | **大型诊断/oracle 路径**：显式 CSR 导出、matrix-free 测试、dot 测试、BlockDot、ExplicitJToracle、GMRES-ILU、以及 stageB5/6/7/8 探针挂载点 |
| `solveDiscreteFlowAdjointProduction.H` | **生产路径**（优化迭代实际要跑的）：`applyProdFlowJT`、FGMRES/BiCGSTAB、`discreteProdPreconditionerSetup`（pressureGAMG/diagonal/bruteForceL1）、`PRODPRECSETUP`/`PRODPRECGAMGSETUP`/`FGMRES-PROD`/`BCGS-PROD`/`PRODRESID` 日志标记、`discreteProdEnableRitzPilot`。**BFINAL-009/010 焦点文件** |

### 5.4 Stage-B 探针（诊断，默认关，切勿写入生产路径）
| 文件 | 内容 |
|---|---|
| `stageB5BoundaryRelaxOracle.H` | BFINAL-002/003 边界一致 P 行 oracle + dphi/dU 松弛 oracle |
| `stageB6RxDesignOracle.H` | R_x / 切向诊断（RX-A/B/C、加权贡献、prodGsenDPressDrop 复刻） |
| `stageB7PressureBCDiagnostic.H` | 压力边界诊断（BFINAL-007） |
| `stageB8JPPActualResidualFD.H` | 实际冻结-原始压力残差 Jv-FD 探针（BFINAL-008，P1/P3 门） |

### 5.5 压力参考语义（重要，勿回归）
`solveDiscreteFlowAdjoint.H` / `solveDiscreteFlowAdjointProduction.H` 中均有：
```cpp
const bool <前缀>PressureNeedsReference = p.needReference();
```
fixedValue 出口 p=0 时该值为 `false` → **不加任何人工 identity 行、不做 null-space 投影、不 pin 行**；所有物理压力行保持激活。旧代码无条件 `[pRef]=0` 的写法已被静态测试禁止（`test_stage_b_safety_gates.py`）。

### 5.6 验证 / 状态 / 其他
| 文件 | 内容 |
|---|---|
| `validateMmaUnlockGate.H` | MMA 硬门：`mmaUpdateEnabled && !frozenGradientValidated` 时禁止 |
| `evaluateCandidate.H` / `saveAcceptedState.H` / `restoreAcceptedState.H` | 候选评估 / 事务回滚 |
| `saveOptimizerState.H` / `loadOptimizerState.H` | 重启持久化 |
| `validateGate6SSTStrict.H` / `validateSSTDirection.H` / `validateStageB2GradientAmplitude.H` / `validateStageBRepeatability.H` | 各类验收门 |
| `MMA/` | MMA 优化器实现 |
| `tests/test_stage_b_safety_gates.py` | Stage-B 静态回归（5 个测试） |

---

## 6. Python 环境（三套，别混用）

| 环境 | 位置 | 用途 |
|---|---|---|
| WSL `python3` | 3.8.x，无 numpy/scipy | 纯逻辑单元测试（`src/tests/*.py`） |
| **venv_s4** | `/home/ys/dsH/TO-ANISOTROPIC/evidence/agent-group/BFINAL-002/.venv_s4/bin/python`（numpy 2.4.6 + scipy 1.17.1） | **离线稀疏矩阵分析**：splu 直接求解、P2/P8 门、sigma_min 等。**大矩阵离线分析优先用它** |
| Windows Python | `/mnt/c/Users/admin/AppData/Local/Programs/Python/Python313/python.exe`（numpy 2.5.0 + scipy 1.18.0） | 历史 oracle 工作流（`solve_oracle_pd.py` 等），通过 UNC 访问 WSL 文件 |

`explicitJT.mtx` 约 467 MB（134400² 稀疏，nnz≈6.2e6），超过 GitHub 100 MB 单文件限制 → **不要提交**（`.gitignore` 已含 `*.mtx`）。

---

## 7. Git 约定

- 分支：`agent/dsH-stage-b-validation`（当前活动分支）；`agent/p0-prod-preconditioner-onnz` 为本地只读存档分支（勿删）。
- 推送：`git push`（分支 remote 内嵌凭据）；不要 `git push -u origin ...`。
- 提交纪律：**选择性提交**（代码 + 轻量证据；重工件 `*.mtx`、`build/`、`.venv*` 已被 `.gitignore` 排除）。
- 工作区卫生：`git reset --hard` / `git clean -fd` **禁止**（可能破坏他人工作）。
- 大 `.mtx`/日志如需留档：放 `evidence/agent-group/<ROUND>/` 下，仅提交小文件与摘要；或本地保留不提交。

---

## 8. 已知坑与铁律（按踩坑频率排序）

1. **旧二进制陷阱**：`/home/ys/OpenFOAM/ys-7/platforms/.../bin/MTO_HF` 是旧的。运行用绝对路径 `build/bin/MTO_HF` 或 `run_verify_*.sh`。
2. **PATH 污染**：先干净 PATH 再 `source bashrc`，否则静默失败。
3. **`set -e` / tee / SIGPIPE**：不用 `set -e`；日志直接重定向，不用 `tee` 管道（SIGPIPE 会把 wmake 一起带走）。
4. **`FOAM_SIGFPE`**：必须 `unset`，否则伴随 NaN abort。
5. **并行伴随尚未完成**（P2）：当前手写 matrix-free Jᵀ 只遍历本进程内部面，**没有 processor-patch 邻域伴随值交换**。直接删 `if (Pstream::parRun())` 会得到能跑但梯度错误的并行伴随。并行化需要：伴随速度/压力/hA/gradientAdjoint 的 processor halo exchange + processor patch 上 U-U / U-P / P-U / P-P / frozen deviatoric 全部转置项 + 保留全局 reduce 点积，然后在 1/2/4 核上验证（并行 J/Jᵀ 点积 <1e-10、并行 vs 串行伴随解 <1e-8、三方向梯度 <1e-6、最终目标/压降/体积分数梯度通过 FD）。**顺序：先串行生产伴随收敛（BFINAL-011），再做并行。**
6. **BFINAL-009 停止点已由 BFINAL-010 取代**：旧 `PRODPRECGAMGSETUP diagRelL2=0.0408 > 1e-8` 是「裸 fvMatrix 对角（不含边界）」对比「含边界 internalCoeffs 的有效对角」的不相容比较造成的假阳性（`addBoundaryDiag` 语义假设已被证实）。BFINAL-010 已将其替换为四组矩阵/算子等价性分解 + 硬门：`PRODPRECGAMGCHECK (label): interiorDiagRelL2=… boundaryDiagRelL2=… effectiveDiagRelL2=… offDiagRelL2=…`，任一非有限或 effective/offDiag > 1e-8 即 FatalError（阈值不可配置）。实测四值 ≈1e-16（机器精度，两标签一致）。另注意：`fvm::laplacian(volScalarField,…)` 的系数面插值词取自 `laplacianSchemes` 条目 ITstream（不是 interpolationSchemes）——见 §4.1 的 scheme 守卫。b010 观察到的 pressureDrop GAMG 内层 NaN 已由 BFINAL-011 定性为门代码存储物化副作用并修复（见 §8.10），**不要**改插值 scheme 或调 FGMRES 掩盖问题。
7. **MMA 门**：`mmaUpdateEnabled && !frozenGradientValidated` 时 MMA 绝对禁止；不要绕过。
8. **`stageB4JacobianProbe` 与 `stageB2Enabled` 互斥**：同时开会 FatalError。
9. **诊断 vs 生产伴随**：`solveDiscreteFlowAdjoint.H` 是诊断环境（大），`solveDiscreteFlowAdjointProduction.H` 是生产路径（要小、可扩展）。不要把诊断机器复制进生产路径。
10. **`lduMatrix` 非 const 访问会物化对称存储（BFINAL-011 教训）**：`-fvm::laplacian` 组装的矩阵是对称存储（只存 upper，`lower()` 返回 upper 别名）。对非 const 矩阵调用 `.lower()` / `.upper()` 会把 lower 物化为独立副本 → `matrix.symmetric()==false` → PCG 被拒（"Unknown asymmetric matrix solver"）、GAMG 走非对称分支（scaleCorrection 默认 false），足以让某些 rhs（如 pressureDrop 的边界尖峰 P 源）在首个内层 V-cycle 就 NaN。**任何只读检查矩阵的代码必须通过 `const lduMatrix&` 别名访问**（静态测试已加防回归断言）。BFINAL-011 cycle-1 已修复此问题并使两标签生产伴随首次双双收敛（TC 750 iter 9.50e-10 / PD 902 iter 9.93e-10，MTO_RC=0），详见 `evidence/agent-group/BFINAL-011/cycle-1/FINAL_REPORT.md`。

---

## 9. 快速验证命令清单（每次改动后）

```bash
# 1) 静态回归
cd /home/ys/dsH/TO-ANISOTROPIC && python3 src/tests/test_stage_b_safety_gates.py

# 2) 编译（直接重定向，不 tee）
cd src && wclean > /dev/null 2>&1 && wmake > ../wmake_current.log 2>&1
tail -2 ../wmake_current.log          # 期望 WMAKE_EXIT=0

# 3) 生产路径验证运行
bash /home/ys/dsH/run_verify_prod.sh /home/ys/dsH/b8_verify_prod <log>

# 4) 关键日志提取
grep -E "PRODPRECSETUP|PRODPRECGAMGSCHEME|PRODPRECGAMGSETUP|PRODPRECGAMGCHECK|GAMG|FGMRES-PROD|PRODRESID|Production reduced" <log>
```

---

## 10. 当前状态速查（BFINAL-016 cycle-1 第〇阶段完成：量具对齐 + 重要修正，2026-08-20）

- **阶段**：B-final —— 冻结湍流梯度 / 生产伴随闭合。
- **已闭合（求解层）**：J_PU（BFINAL-003）、R_x 压力行（BFINAL-005）、J_PP 实际原始闭包（BFINAL-008）；预条件矩阵/算子四组等价性硬门（BFINAL-010，1e-16）；两标签生产伴随迭代收敛（BFINAL-011）。
- **BFINAL-024（2026-08-20，实现轮）里程碑正式升级声明**：BFINAL-023 证实的第一个生产算子真缺陷 **H7（phiHbyA 压力通道 −interp(rAtU·∇p)·Sf，源于 UEqn 源内一次 −∇p 拷贝）已进入生产 J^T、诊断 J/J^T、b_TC 折叠与显式 CSR（四处一致）**——BFINAL-003/008 的 P 行通量切线语义据此正式升级（原槽位不变，新增 H7 槽）；pressureGAMG 预条件器保持 kf-Laplacian 近似（等价性门四值逐位不变，PRODH7SHARE ~0.70 量化排除份额）。面锚 D2 31.2→22.1%、D1/D3 4.9→4.6-4.8%；gDP 因子 2.145/2.159/2.181 → 2.061/2.007/2.100；J 符号表不变（D1 反）。残余 O(1) 失配方向均匀（~2.0-2.1×），不能归因于 P 行面级切线缺陷 → 上移至 J 链更高层。见 `evidence/agent-group/BFINAL-024/cycle-1/`。
- **BFINAL-012（FAIL）→ 013（定位）→ 014（推导 STOP）→ 015（T2/T3 拆分双缺陷）**。
- **BFINAL-016 cycle-1 第〇阶段（量具对齐，修复未启动）**：(A1) 偏应力块可忽略（3.5e-6）；**(A1b/c) 生产 λ 以 8e-11 满足导出系统——产=导同一算子，BFINAL-015「双算子」系钉扎伪影已更正**；修正后算子量具 λ-direct vs 真值 = **1.769/0.856/2.274**；(B1) **b_TC ≠ T-消除后的 dJ/d(U,p)**（与真流介导差 86%/70%/符号翻转）；**新发现装配矛盾：sensitivity 装配 ≠ λ-direct（D2 147%）**——收缩装配链有独立缺陷，BFINAL-013 单项定位据此修订。A2/B2 留待下一 cycle（顺序：装配审计→算子→源）。见 `evidence/agent-group/BFINAL-016/cycle-1/`。
- **当前阻塞**：三个待修缺陷——算子 P 行语义（产=导）、sensitivity 装配链、b_TC 热消除折叠。
- **下一步任务**：装配审计（gsenshMomentum/rxPressureRowT vs rxc 的 anRUa/assembleWeightedRPa 逐项对比）→ A2 算子修正 → B2 源修正 → BFINAL-012 FD 门全套重验 → C0 串行 MMA smoke。
- **并行伴随（P2）**：未解锁（见 §8.5）。**MMA**：锁定。
