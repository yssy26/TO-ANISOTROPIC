# OpenFOAM Othmer 连续伴随 (adjointShapeOptimisationFoam) — 孔隙率(Sp)伴随源耦合符号参照

URL (源码, OpenFOAM-dev 仓库):
https://develop.openfoam.com/Development/OpenFOAM-dev
File: `applications/legacy/incompressible/adjointShapeOptimisationFoam/adjointShapeOptimisationFoam.C`
Clone: `/tmp/of_dev_ref`

## 原始文献
E. Othmer, "A continuous adjoint formulation for the computation of topological and surface sensitivities of ducted flows", IJNMF 58 (2008). (adjointShapeOptimisationFoam 的原始出处)

## 1. 原始(primal) U 方程 — 孔隙率/多孔源

    UEqn == fvModels.source(U)
    其中 UEqn = fvm::div(phi, U) + turbulence->divDevSigma(U) + fvm::Sp(alpha, U)

- `fvm::Sp(alpha, U)` = 多孔(孔隙率)阻力项，对应敏感度通道 α（材料指示函数）。
- 与 SOFO 级 porosity(Brinkman) 源同类：阻力系数正定（Sp 正系数）。

## 2. 伴随 Ua 方程 — 多孔项的伴随转置

    UaEqn = fvm::div(-phi, Ua)
          - adjointTransposeConvection
          + turbulence->divDevSigma(Ua)
          + fvm::Sp(alpha, Ua);
    solve(UaEqn == -fvi::grad(pa));

- 伴随中的多孔项也是 `+ fvm::Sp(alpha, Ua)`：**阻力项 Sp(α,·) 的伴随转置就是它自己**（对角标量项自伴随），符号为 +，无翻转。
- 这一条与 BFINAL-003/018 的结论一致：Brinkman/多孔阻力对伴随系统是自伴随的对角源，不存在跨方程的符号翻转机制。

## 3. 压力投影与出口

    paEqn: fvm::laplacian(rAUa, pa) == fvi::div(phiHbyAa)
    Ua = HbyAa() - rAUa()*fvi::grad(pa)

- 与 SIMPLE 压力修正同构；伴随出口速度 BC 见 `adjointOutletVelocity/adjointOutletVelocityFvPatchVectorField.C`:
  `vectorField::operator=(phiap*patch().Sf()/sqr(patch().magSf()) + UtHat)`
  （伴随出口速度 = phia·Sf/|Sf|² + 切向分量，法向由 phia 决定）

## 4. 设计敏感度 — 目标通过 Ua·U 进入

    alpha += relax*( min(max(alpha + lambda*(Ua & U), zeroAlpha), alphaMax) - alpha )

- 目标敏感度 = λ·(Ua·U)：多孔设计变量敏感度是伴随速度与原始速度的点积（拓扑敏感度标准形），正号。

## 5. 对 BFINAL-028 的参照价值

- 确认"孔隙率/多孔阻力 → 伴随动量"的耦合符号约定为 **+Sp 自伴随**，无隐藏负号。
- 但 Othmer 是**连续伴随**，能量方程不在此框架内（无 dJ/dT→U/p 折叠），且目标为压降类（经 Ua·U 进设计敏感度），**不覆盖** 我方 T-elimination 折叠（dJ/dT → b_TC）的逐项符号。该空位由 Fira thermalTopO 的 gPhi 公式部分填补。
