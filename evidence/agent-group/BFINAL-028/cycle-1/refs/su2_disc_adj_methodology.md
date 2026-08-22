# SU2 参照：离散伴随方法论（AD 自动微分，无手工 g_f 折叠）

来源：https://github.com/su2code/SU2（shallow clone，/tmp/su2_ref，HEAD bc15466）

## 关键事实

SU2 的 DISCRETE_ADJOINT 走 **算法微分（coDiPack reverse-AD）** 通道：
- CDiscAdjFluidIteration / CDiscAdjHeatIteration（SU2_CFD/src/iteration/）驱动
  反向模式的直接求解器 tape，残差对输入的偏导由 AD 自动生成。
- CHeatSolver.cpp:290-320 `Viscous_Residual`：
  ```cpp
  bool pausePreacc = false;
  if (ReducerStrategy)
    pausePreacc = AD::PausePreaccumulation();
  else
    AD::StartNoSharedReading();
  ... // edge loop
  AD::ResumePreaccumulation(pausePreacc);
  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    if (implicit) Jacobian.SetDiagonalAsColumnSum();
  }
  ```
  (CHeatSolver.cpp:313 注释 "Restore preaccumulation and adjoint evaluation state.")
  这是 preaccumulation（局部 Jacobian 预累积 + SetDiagonalAsColumnSum）策略。

## 方法论差异标注（对照表将引用）

1. **SU2 无 "dR_T/dφ_f 手工转置" 环节**：热残差对通量/速度/压力输入的所有偏导
   由 AD 从算子定义自动获得。因此 SU2 **不提供** 与我方 g0 显式 MINUS
   转置公式逐项可对照的独立符号约定——AD 通道下不存在"手工约定符号出错"的
   空间，也无可引用的符号源。g0/g1 的符号对照必须依赖：连续伴随文献
   （下一参照层）或 AD 自身正确性。
2. **热伴随折叠进流伴随**：CHT 界面（SU2_CFD/src/interfaces/cht/CConjugateHeatInterface.cpp）
   在离散伴随下通过 BGS（Block Gauss-Seidel）在区带间交换伴随温度/伴随焓
   （BGS_ADJ_TEMPERATURE/BGS_ADJ_ENTHALPY，chtMaster.cfg SCREEN_OUTPUT），
   而非在单一 (U,p) 线性系统内折叠 b_TC。即 SU2 的多物理离散伴随 = 分区伴随
   耦合迭代，我方 = 单系统折叠。方法论分支点。
3. 结论：SU2 作为"热折叠项逐项对照"的参照系价值集中在 **泛函形态差异**
   （面积平均 vs 质量加权）与 **方法论对比**（AD 自动 vs 手工转置），
   不提供我方 T1-T8/H7s1/H7s2 槽位符号的逐项独立镜像。

## CHT 边界类型判定（CDriver.cpp:2518 附近）

- 两侧都是 heat：CONJUGATE_HEAT_SS
- 流体侧无能量方程 + 弱耦合热：CONJUGATE_HEAT_WEAKLY_FS / _SF
- nVar=4 的 CConjugateHeatInterface
