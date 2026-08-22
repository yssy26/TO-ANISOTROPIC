# SU2 参照：AVG_TEMPERATURE 泛函（面积平均温度）

来源：https://github.com/su2code/SU2（shallow clone，/tmp/su2_ref，HEAD bc15466）
文件：SU2_CFD/src/solvers/CHeatSolver.cpp（源码路径对应 GitHub 相同路径）

## 关键事实

SU2 的 AVG_TEMPERATURE 泛函是 **AREA 平均**温度，不是质量流量加权温度。

CHeatSolver.cpp:712（HEAT_FLUX 边界上累加）：
```cpp
if ( Boundary == HEAT_FLUX ) {
  AverageT_per_Marker[iMarker] += Twall*config->GetTemperature_Ref()*Area;
}
```

CHeatSolver.cpp:721-731（归一化）：
```cpp
AllBound_AverageT += AverageT_per_Marker[iMarker];
...
SU2_MPI::Allreduce(&Send_Bound_AverageT, &AllBound_AverageT, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
if (Total_HeatFlux_Areas_Monitor != 0.0) {
  Total_AverageT = AllBound_AverageT/Total_HeatFlux_Areas_Monitor;
}
```

CHeatSolver.hpp:353-354（目标函数组合）：
```cpp
case AVG_TEMPERATURE:
  Total_ComboObj = weight * Total_AverageT;
```

## 与我方对照的差异标注

我方泛函（computeObjective.H:72-82, L141-143）是 **质量流量加权出口混合温度**
  Tmix = Σφ_out·T_out/Σφ_out，dJ/dT = −φ_outlet/(M·Tref)，dJ/dφ_f = −(T_f−Tmix)/(M·Tref)。

SU2 无质量加权出口混合温度泛函在此路径；AVG_TEMPERATURE 是面积平均，
其伴随源 d(area-avg)/dT = Area_face/TotalArea（不含通量加权）。
=> g1 出口通量导数项（我方 dJ/dφ_out ≠ 0）在 SU2 面积平均泛函中 **不存在对应物**
   （面积平均的 J 不依赖出口通量 φ）。这是泛函形态差异，不是符号约定差异。

## 测试算例

- TestCases/disc_adj_heat/disc_adj_heat.cfg（2018-11-26，Ole Burghardt，TU Kaiserslautern）：
  WEAKLY_COUPLED_HEAT_EQUATION= YES, MATH_PROBLEM= DISCRETE_ADJOINT,
  OBJECTIVE_FUNCTION= TOTAL_HEATFLUX, INC_ENERGY_EQUATION= NO
- TestCases/coupled_cht/disc_adj_unsteadyCHT_cylinder/：
  fluid.cfg OBJECTIVE_FUNCTION= AVG_TEMPERATURE, MARKER_ANALYZE_AVERAGE= AREA；
  chtMaster.cfg SCREEN_OUTPUT 含 BGS_ADJ_ENTHALPY/BGS_ADJ_TEMPERATURE/SENS_TEMP/SENS_GEO；
  gradient_validation.py（FADO）对 18 个 FFD 设计变量做非定常 CHT 伴随梯度 FD 校验。
