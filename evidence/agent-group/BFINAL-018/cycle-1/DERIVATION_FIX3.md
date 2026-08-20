# BFINAL-018 — 修复③（b_TC 热消除折叠）推导与离线定量

## 1. 问题设定

热目标 `J = -Tmix/Tref`（`Tmix = Σφ_out·T_out/Σφ_out`）依赖（a）出口单元温度 T、
（b）出口面通量 φ_out。冻结验证态下 `J` 对设计的全导数经 T 消除分解为：

```
FD_J = A + B + C
A = J_φ^T·φ'_out            出口通量显式部分（FD_J − thermal_gauge）
B = −Tb^T·(dR_T/dφ)·φ'      T 对通量变化的响应（thermal_gauge − C）
C = −Tb^T·R_T,x             直接热项（dDTDxh 折叠，生产 thermalDiffusionDerivativeDTCell）
```

其中 `thermal_gauge = dJ/dT·T'`（总 T 响应，B16 量具实测）；拉格朗日推导
（本文件 §3）给出 **T-消除后 (U,p) 源向量**：

```
b_TC = M^T·g ,   g_f = dJ/dφ_f − (dR_T/dφ_f)^T·Tb
M = ∂φ/∂(U,p)   （收敛 SIMPLE 通量对 (U,p) 的切线）
```

且完整的 J 设计导数须补两个**不能被 (U,p) 源吸收**的项：

```
dJ/dx = −λ^T·R_flow,x + C + g^T·φ_x ,   λ = J_flow^{-T}·b_TC
```

`φ_x` = 通量对 α 的直接依赖（经 rAU/Sp 的 Brinkman 链，结构与 rxPressureRowT
的 drAU/dHbyA 项相同）。因此量具恒等式的**无歧义形式**为：

```
b_TC^T·w_true == FD_J − C − Gx ,   Gx := g^T·φ_x
```

B16 预注册口径 `b_TC^T·w_true == FD_J − thermal_gauge (=A)` 只可能被
「仅出口项折叠」满足——而该口径在生产伴随链中不可实现（B 无法解析计算，
必须经流行传播）——这是 B16 量具把热介导项「放错边」的根源（B16 审阅要点
「B1 恒等式是否把热介导项放对了边」的实证答案）。

## 2. 缺陷定位：路由（不是 g 公式）

- **g 公式无罪**：内部面 `g_f = −mask·Tb_down·(T_n−T_o)` 是有界迎风残差
  dR_T/dφ 的精确转置（B0.2 逐面 FD 1e-10 验证）；出口 `+dJ/dφ` 项 B0.1 验证；
  入口（fixedValue U）通量非自由变量、跳过正确。
- **路由有罪**：旧代码把 g 路由为 `U 行 += w·Sf·g`（全权重插值）+
  `P(own) −= kf·g; P(nei) += kf·g`——这是**另一个映射**
  `dφ = interp(dU)&Sf + kf(p_n−p_o)` 的转置。而 J 算子 P 行（BFINAL-003
  已验证语义）实现的通量切线是

```
dφ_f = Sf & ( αrel·( w·rAU_o·dH_o + (1−w)·rAU_n·dH_n )
            + (1−αrel)·( w·dU_o + (1−w)·dU_n ) )
     − kf_f·( dp_n − dp_o )
```

  旧路由 (i) 缺 αrel·rAU·dH 的非局部 H^T 块；(ii) 缺 (1−αrel) 因子；
  (iii) **内部 kf 符号与算子自身 J_PP 块相反**（算子转置槽位为
  `out[P(own)] += kf·λPdiff`）。修复 = 逐槽位镜像算子 P 行的 hA/directRel/kf
  三个通道 + deltaH^T（动量 H 转置，含边界 H 对角块）。

## 3. 离线定量（`b18_folding_lab.py`，全部用 b16_states/b15_export 冻结数据）

### 3.1 重建验证

| 项 | 结果 |
|---|---|
| 网格几何（Sf/owner/weights vs b16mesh 导出） | maxdiff 1.4e-19 / 0 |
| 旧 b_TC 重建 vs 导出 explicitRhs_thermalCoupling.mtx | relL2 **2.4e-5** |
| 热算子重建 A_T·T − src − Q | maxabs **2.7e-7**（Q=0） |
| A_T^T·Tb_written − dJdT | relL2 **3.4e-3**（Tb 求解容差级） |
| 我的 C 公式 vs B13 in-pass 探针 thermal（D 配对） | **−1.4367e-4 / −4.7591e-3 / +3.7085e-3 vs −1.442e-4 / −4.759e-3 / 3.698e-3（≤0.4%）** |

### 3.2 C 的真值裁定（直接线性解，不经伴随公式）

`A_T·δT = −(A_x·dx)·T`（纯直接 xh 扰动、流冻结）逐方向求解：
C_true = dJdT·δT = **−0.0303 / −0.0083 / −0.0190**（与生产公式值一致至
0.15–3%；符号按残差约定取负）。生产 thermalDiffusionDerivativeDTCell **无罪**。
（注：b15 写盘 fsenshMeanT 含陈旧导入 λ_TC 的动量项——b15 算例 1/ 目录为
cp -a 残留、写盘 Ub=uniform(0) 与内存态不符——不可用它反演 C；λ_TC15 从未
导出，B17 禁用清单正确。）

### 3.3 量具表（w_true 冻结态，h=1e-3/3e-4 双步长稳定）

| 方向 | FD_J | thermal_g | C | Gx | A=FD−tg | FD−C−Gx（目标） | 旧 b_TC | 新 b_TC(V1) |
|---|---|---|---|---|---|---|---|---|
| D1 | +0.04286 | −0.16687 | −0.03034 | −0.00212 | +0.20972 | **+0.07532** | +0.11251（46%） | +0.05305（−30%） |
| D2 | −0.00347 | +0.08230 | −0.00821 | +0.00244 | −0.08577 | **+0.00230** | −0.02566（+1010%） | −0.00113（符号✗，量级≈0.5×） |
| D3 | −0.00028 | −0.05665 | −0.01979 | +0.00762 | +0.05637 | **+0.01189** | −0.04493（**符号✗**） | +0.03889（+227%） |

结论：
1. 新路由修复 D3 符号翻转、D1 幅度误差减半；D2 残差符号在目标 ≈0.002 的
   近对消区，量级落入 0.5×。
2. 三方向仍不闭合 ≤5%。残余误差与缺陷①（算子 P 行/通量耦合语义，
   产=导同一算子、λ-direct vs 真值 1.768/0.856/2.274）共享同一通量图：
   本轮量具用真值 w_true 检验 M·w_true，而算子的通量图语义属下一轮授权。
   连续性检验 `div(M·w_true + φ_x) = 0` 的信噪比（残差 ~3.2e-4 vs
   |div(φ_x)| ~1.9e-4）不足以在 FD 噪声底之上进一步分解。
3. 分块贡献（新路由，D1/D2/D3）：hA^T (+0.026/−0.005/−0.00003) +
   direct^T (+0.024/+0.003/+0.039) + kf^T (−0.004/−0.001/−0.0001)。

## 4. 代码落点

- `solveDiscreteFlowAdjointProduction.H`：删除旧
  pressure-flux-correction-transpose 辅助与平插值 U 路由；新折叠逐槽位镜像
  applyProdFlowJT 的 P 行通量切线（αrel·prodRAU hA + (1−αrel) direct +
  +kf(P own) −kf(P nei) + 出口 kf_b + deltaH^T + 边界 H 对角块）。
- `solveDiscreteFlowAdjoint.H`：同构镜像（rAUAdj/discreteUpper/Lower）。
- `AdjNS_HT.H`：**未改**（g 公式本身正确；冻结温度/掩码语义保持）。
- `AdjHeatTransfer.H`：**未改**（C 公式经 §3.1/3.2 双重验证无罪）。

## 5. 与生产 J 梯度的关系（下一轮装配审计的输入）

生产 `fsenshMeanT = −dAlphaDxh·(U&Ub)·V + C` 相对完整恒等式仍缺
（a）`−pb^T·R_P,x`（压力行项，rxPressureRowTranspose 目前只对 pc 常驻）
与（b）`Gx = g^T·φ_x` 直接项。二者均属 sensitivity.H 装配层（本轮授权范围
之外）；用冻结态 λ 近似估计其量级为 0.003–0.011（pb 项，λ 未收敛仅供参考）
与 0.002–0.008（Gx 项，已精确定量）。
