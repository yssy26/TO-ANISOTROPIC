# BFINAL-018 cycle-1 — 执行摘要（结论先行）

**任务**：合并修复轮——工程缺陷②（伴随容差/停止判据/λ 协议）+ 源缺陷③
（b_TC 热消除折叠），一次编译，验证序列 1–6。

## 一句话结论

**②完全交付（1e-12 双轮双标签 + 泛函稳定判据 + λ 协议，等价性门逐位
1e-16 不变）；③按数学正确口径实现（算子一致通量图转置，实现经块回归
≈1.00 验证），FD 门功能结果：gDP 三方向符号全对、幅度收敛到算子单因子
~2.13×（缺陷①残留，符合预期）；J 符号 2/3 恢复（D3 由错转对、D2 保持、
**D1 仍反**）；gV 与 B12 逐位一致零回归。J(D1) 与 b_TC 量具未闭合部分
如实定位：b_TC^T·w_true 泛函存在近对消结构且对线性化状态极度敏感
（同一代码主循环态/ B2 基线态差 1.8% 使 D1 收缩 0.115→0.018），
残余语义属缺陷①（通量图/算子，本轮禁改）+ sensitivity.H 装配缺口
（pb 压力行项与 Gx 直接项，超出本轮授权文件）。**

## 验证序列结果对照

| 步骤 | 要求 | 实测 | 判定 |
|---|---|---|---|
| 1 静态+编译 | 全过 + WMAKE_EXIT=0 | **11/11 测试 OK**（原 9 + 新 2）；WMAKE_EXIT=0；二进制新于源、含 GRADSTABLE/B18RHSEXPORT | ✅ |
| 2 等价性门 | 四值 1e-16 级逐位不变 | 两标签×两轮 interior/boundary/effective/offDiag = 1.1–1.5e-16 | ✅ |
| 3 双标签 1e-12 | ≤1e-12 或如实报告 | TC 898/894 iter 9.80/9.82e-13；PD 1130/855 iter 9.68/9.36e-13；**GRADSTABLE 4/4**（relChange 5.1e-11…2.8e-12 ≪1e-6）；0 次 Warning | ✅ |
| 3 λ 一致 | 主循环 vs B2 两轮 gradProxy 一致 | 各轮均泛函收敛；但两轮 **持久相差 7.1%**（PD −205194.141 vs −219717.067）——定位为 **B2 模块 full-SST 再基线导致的两轮线性化状态不同**（J=−1.1464 vs −1.16745，b16/b18 逐位复现），非收敛伪影；量具混用路径已由 stage_b2_adjoint_identity.tsv + 协议注释关闭 | ⚠️ 已定位 |
| 4 b_TC 量具 | ≤5% 或如实报告 | 同-pass（第 2 轮）：D1 76%、D2 758%、D3 150% 相对 (FD−C−Gx)；未闭合——见定位分析 | ⚠️ 如实报告 |
| 5 FD 门（功能口径） | gDP 符号全对（幅度属①）；J 符号恢复；gV 不回归 | **gDP 符号 3/3 ✓，ADJ/FD = 2.130/2.113/2.169（单因子）**；**J 符号 0/1/1（D1 仍反，D3 修复）**；gV 符号 3/3、bestRelV 与 B12 **逐位相同** | ⚠️ 部分 |
| 6 诊断模式 | convergeFatal false、记录 MTO_RC | 全程 false；两运行 MTO_RC=0；无 abort | ✅ |

## 修复②内容
1. `discreteFlowAdjointTolerance` 默认 1e-9→**1e-12**（生产+诊断两文件，
   算例同步）；maxIter 保持 4000（实际 855–1130 迭代，充裕）。
2. **GRADSTABLE 泛函稳定判据**（默认开）：FGMRES 残差达标后须最后两个
   cycle 的 gradProxy 相对变化 <1e-6 才准出（`discreteProdGradientStabilityTolerance`
   可调）；否则 Warning 并续迭代至 maxIter。
3. **B2 λ 一致性协议**：B2 模块写 `stage_b2_adjoint_identity.tsv`
   （每轮两标签 gradProxy 指纹）+ 源码协议注释（同-pass 量具、禁止与
   主循环写盘 λ 混用）。
4. 门控 rhs 导出 `stageB18RhsExport`（生产路径物理单位，4 列/单元）。

## 修复③内容
- 旧路由 =「`dφ=interp(dU)·Sf + kf(p_n−p_o)`」的转置（与算子自身 J_PP
  块符号相反、缺 αrel·rAU·H^T 非局部块、缺 (1−αrel) 因子）→
  新路由 = **逐槽位镜像算子 P 行通量切线**（BFINAL-003 语义）：
  hA(αrel·rAU_u)→deltaH^T + direct((1−αrel)·w) + kf(+own/−nei) +
  出口 kf_b。生产与诊断两文件同构修复；g 公式（B0.2 已 1e-10 验证）
  与 AdjHeatTransfer.H 的 C 项（本轮双重裁定无罪）均未改。
- 离线推导（`DERIVATION_FIX3.md` + `b18_folding_lab.py`）：
  - C 真值裁定：A_T·T−src−Q=2.7e-7、A_T^T·Tb−dJdT=3.4e-3 验证热算子
    重建；直接求解 C_true=−0.0303/−0.0083/−0.0190 与生产公式一致
    （B13 in-pass 探针 D 配对 0.4% 互证）。
  - 恒等式澄清：**FD_J = A+B+C；正确源向量目标 = FD−C−Gx**（Gx=通量
    直接 α 项）；B16 预注册口径 FD−thermal_gauge(=A) 只对应「仅出口项
    折叠」，生产不可实现——B16 审阅要点「热介导项放边」的实证答案。
  - 实现保真：单 pass 主循环导出 rhs vs numpy V1 块回归系数
    hA=0.9947/direct=1.0032/kf=1.0010（残差 6.2% = 我方重建输入精度）。

## J(D1) 未恢复的定位（停在此）
1. 生产 J 装配（sensitivity.H，本轮授权外）仍缺 `−pb^T·R_P,x`（压力行项，
   rxPressureRowT 目前只对 pc 常驻）与 `Gx=g^T·φ_x` 直接项——拉格朗日
   推导证明二者不能被任何 (U,p) 源吸收。
2. 缺陷①（算子 P 行/通量耦合语义，本轮禁改）在 TC 源上同样作用——
   gDP 端到端因子 2.13× 即其表现。
3. b_TC 泛函近对消+状态敏感（见上），5% 闭合在①修复前不可达。

## 运行清单
- 编译：`wmake_b18.log`（WMAKE_EXIT=0，全量 ~35 min）。
- 运行：`run_b18fix/`（完整 FD 门 1327 s，MTO_RC=0）、`run_b18main/`
  （单 pass 主循环，实现保真用，MTO_RC=0）。
- 离线：`b18_folding_lab.py/.log/.json`、`b18_gauge.py/.json`、
  `b18_decomp.json`、`DERIVATION_FIX3.md`。
