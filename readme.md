# TO-ANISOTROPIC

> **New AI Agent / developer: read [`AI_AGENT_HANDOFF.md`](./AI_AGENT_HANDOFF.md) first.**  
> It is the current authoritative project handoff for goals, scope lock, algorithm path, code architecture, Stage A/B status, BFINAL-009 failure, BFINAL-010 next task, safety gates, and the roadmap to MMA/full-SST validation.  
> The notes below are older implementation notes and may describe historical solver behavior or obsolete commits.

---

源码目录 src/
主入口与主流程编排：

MTO_HF.C — 主程序入口，把所有 .H include 进来按时间步推进
createFields.H、opt_initialization.H、readTransportProperties.H、readThermalProperties.H — 场初始化、优化参数、材料物性
前向（primal）求解链：

NS.H — 流动求解（SIMPLE，含冻结湍流）
HeatTransfer.H — 传热求解
updateThermalFlux.H、updateFrozenTurbulenceFields.H、updateMaterialProperties.H、filter_x.H — 通量/冻结湍流/材料/过滤链
伴随（adjoint）求解链：

AdjNS_HT.H — 热耦合伴随（换热目标）
AdjNS_PD.H — 压降伴随（本次关注的目标）
AdjNS_FF.H — 冷流伴随
sensitivity.H、filter_chainrule.H — 梯度装配 + 过滤链（把伴随场变成对设计变量的梯度 dgdx[0/1]）
离散流动伴随的两个求解器（接手重点）：

solveDiscreteFlowAdjointProduction.H — 生产路径，缩减离散伴随，FGMRES + L1 预条件（会卡在残差 0.157）
solveDiscreteFlowAdjoint.H — B4 诊断/oracle 路径，显式 CSR 组装 + 矩阵自由 deviatoric 转置 + dot-test 校验，可导出 explicitJT.mtx 并导入外部精确解
门控/验证（每个 validate*.H 对应一个验收关卡）：

validateStageB2GradientAmplitude.H — Stage B2/B3 幅值验收（FD vs ADJ 对比，D1/D2/D3 投影）
validateMmaUnlockGate.H — MMA 解锁总闸（含 discreteExplicitSolutionIncludesDeviatoric 等开关）
其余 validate*.H 是各阶段的验收/自检
本次会话改动的关键文件
文件	改动
src/solveDiscreteFlowAdjoint.H	新增 deviatoric 转置显式组装进 CSR（使 explicitJT.mtx 成为完整 J^T 算子）；修复 3 处 bug（dev2 的 2/3 迹系数、边界投影 δ_ij 误乘、neighbor 侧符号颠倒），现已达机器精度 3.9e-16
src/solveDiscreteFlowAdjointProduction.H	上一轮改的：允许伴随非致命退出 + checkpoint 梯度代理
src/validateStageB2GradientAmplitude.H	上一轮改的：冻结前向收敛判据改为窗口均值 + 噪声底容差
本次会话新增的辅助脚本（仓库根目录，未入 git）
构建：build_oracle.sh（干净 PATH + wmake）
运行：run_b4_export.sh（导出完整矩阵）、run_b4_import.sh（导入精确解并算梯度）、run_b2diag.sh、build_minifd.sh
scipy 求解：solve_oracle_pd.py（SuperLU 直接求解 + RHS 重排）
投影/验证：project_pd.py、crosscheck_pd.py（D1/D2/D3 投影 + 体积链交叉验证）
小工具：probe_*.sh、check_*.sh、peek_log.sh、show_cfg.sh、write_C.sh、ps_check.sh
算例目录（在 WSL /home/ys/ 下，仓库外）
/home/ys/b2_case — 原始算例，含旧的 explicitJT.mtx、explicitSol_*.mtx（8 月 9 日的旧 oracle 产物）
/home/ys/b2_case_smoke — 当前测试算例，本次 exact oracle 全在这里跑