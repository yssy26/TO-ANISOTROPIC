# Fira thermalTopO — 连续伴随 CHT-TO 热-动量耦合敏感度（我方 g0 折叠的最近逐项参照）

URL: https://github.com/Fira-Software/thermalTopO
Source file: `src/solvers/thermalAdjointSimple/thermalAdjointSimple.C`
Derivation doc: `docs/derivation.md` (§3.1, §3.2, §6.2)
Clone: `/tmp/fira_thermalTopO` (HEAD captured 2026-08-21)

## 1. 连续伴随体系（对照架构前提）

Lagrangian（derivation.md §3）:

    L = J_P + ∫_Ω (u_a·R_u + p_a R_p + T_a R_T) dΩ

- 温度伴随方程 (§3.1, AT): `-u·∇T_a - ∇·(D_eff ∇T_a) = 0`（边界目标下无体积源；目标从 Γ_pf 通量 BC 进入，见下）
- 动量伴随方程 (§3.2): 新源项 `+ T_a ∇T` (ATC-T)，加在 RHS，与伴随转置对流同侧。
- 温度目标两类：
  - p-norm 目标: `∫_{Γ_pf} P T^{P-1} δT dΓ` → 以 `Γ_pf` 上的通量型 BC `D_eff ∂T_a/∂n = -P T^{P-1}` 进入（目标不直接出现在折叠里）
  - meanTemperature: Γ_pf BC 齐次化，体积源 `-1/|Z|` 进 (AT)
- 我的 D1 目标 J 是出口混合温度（边界泛函）→ Fira 里最接近的是 p-norm 的边界目标路径，但 Fira 用它驱动 T_a BC（固定梯度），不驱动 U/p 源。这是连续伴随/离散伴随架构差异。

## 2. thermalFluxSensitivity() 公式（与 my-side g0 逐项对照的锚点）

实现（L1280-1390，仅支持 Gauss upwind + boundedConvection 开关匹配 fvSchemes）：

内部面:
    AP = CP*Ta[P],  AN = CN*Ta[N]          (C = ρc_p 归一系数，定常属性 C≡1)
    Tf = (phi[facei] >= 0) ? T[P] : T[N]   (迎风原始温度)
    gPhiI[facei] = Tf*(AP - AN) - chi*(AP*T[P] - AN*T[N])
                 = Ta_upwind·(T_own - T_nei)    当 C≡1 且 chi=0(非 bounded)
    chi = boundedConvection_ ? 1 : 0

边界面:
    gPatch[facei] = CP*Ta[celli]*(TPatch[facei] - chi*T[celli])
                 = Ta[celli]*TPatch[facei]      当 C≡1 且 chi=0

结构要点（对照我 g0 = -mask_f·Tb_down·(T_nei - T_own)）:

1. **迎风伴随 × 温差结构相同**：两方都是 `(迎风伴随值) × (own-nei 温差)`。
   - Fira（chi=0）: gPhi = Ta_upwind·(T_own − T_nei)
   - 我方 g0:      g_f  = −Tb_down·(T_nei − T_own) = Tb_down·(T_own − T_nei)
   - **两者完全同形**：Fira 的 Ta_upwind 对应我方 Tb_down（upwind 伴随），温差方向一致，无符号翻转。
2. **差分形式不同**：Fira 用 `Tf*(AP−AN)`（伴随差直接乘迎风 T）；我方用 `Tb_down*(T_N−T_O)`（原始温差乘迎风伴随）。离散上两者是同一连续极限的两种拆分，但离散展开不一样——Fira 的拆法（伴随差×Tf）在逐面意义上对应 `d(phi*T)/dT` 的伴随转置。
3. **bounded 项**：Fira 的 `chi*(AP*T[P] − AN*T[N])` 是 bounded 伴随修正（对 Sp(∇·(Cu), T_a) 的修正），我方离散伴随 g0 无此类 chi 修正项——但这是连续伴随的 bounded-scheme 特有问题，我方采用 exact 离散转置（对 Gauss upwind 的离散算子逐项转置），不适用。
4. **符号对照**：Fira 官方 FD 验证 §5（耦合符号 couplingSign=+1, thermalSensScale=+1），ATC-T 源为 +T_a∇T 进动量 RHS。我方 g0 显式 MINUS 与 Fira 的"迎风伴随×(own−nei)"在能量→动量折叠的同形结果一致，未发现方向性反号。

## 3. projectedFluxMomentumSource()（gPhi 路由进动量源）

（L1578-1800；lambda 投影 + 面-单元路由）

    lambdaSource[P] += -lower[facei]*gPhiI[facei];
    lambdaSource[N] +=  upper[facei]*gPhiI[facei];
    边界: lambdaSource[faceCells[facei]] += intCoeffs[facei]*gPhip[facei];
    solve: fvm::laplacian(rAtU, lambda)  (投影求解，tol 1e-12)
    gFI[facei]  = gPhiI[facei] - (lambda[P] - lambda[N]);
    边界 gFpatch[facei] = gPhip[facei] - lambda[faceCells[facei]];
    路由: S[P] += w[facei]*Gf/V[P];  S[N] += (1-w[facei])*Gf/V[N];
          Gf = gFI[facei]*Sf[facei]
    边界: S[celli] += gFpatch[facei] * cmptMultiply(valueInternalCoeffs, SfPatch)/V[celli]

架构差异：Fira 用 lambda 投影把 gPhi 投影到无散面通量，再以 w*Sf/V 的单元权重路由；我方离散伴随用 M^T 精确转置（αRel·rAU·dH^T + (1-αRel) 直接 + kf P 行槽），无投影、无 1/V。两者是连续伴随（投影）/离散伴随（精确转置）的方法论差，不做逐槽数值对照。

## 4. 结论锚点

- Fira 的 gPhi（chi=0）与我方 g0 在"迎风伴随 × own-nei 温差"上同形同号，官方 FD 门控验证通过（§5 一致收敛）。
- 未见任何参照侧证据支持把 g0 前再加一个负号；BFINAL-026 的 O(1) 反号形态（ADJ=−0.0358 vs FD=+0.0425）不太可能由 g0 自身符号错误引起。
- 候选错项重心转向路由层（T1-T8/H7s1-s2 的 αRel·rAU·dH 槽、kf 槽、H7 压力通道）与参照差异维度（Gx、C 超出授权范围）。

## 5. 补充：docs/atc-t-open-channel.md 的独立佐证

仓库内 `docs/atc-t-open-channel.md` 明确记录：

> "The open-channel example helped expose why **a continuous thermal-to-momentum
> source is not, by itself, a complete discrete transpose of the segregated
> SIMPLE map** in an open-flow topology optimisation setting. The validated
> production path therefore does not rely on the continuous source alone. It
> uses the fixed-point map transpose and the corresponding beta reverse."

即：连续伴随的 `+T_a∇T` 源单独**不是**分离式 SIMPLE 映射的完整离散转置——必须配合 fixed-point map transpose（(I−M_x)^T ψ = J_x, dJ/dβ = J_β + ψ^T M_β）才构成被验证的精确路径。

与我方 BFINAL-018 的结构性结论同构：g0（连续极限意义上的热折叠）本身符号无辜，真正的转置完备性（αRel·rAU·dH^T 等 SIMPLE 映射槽位）在路由层。这为"候选错项在路由层而非 g0"提供参照系侧独立背书。
