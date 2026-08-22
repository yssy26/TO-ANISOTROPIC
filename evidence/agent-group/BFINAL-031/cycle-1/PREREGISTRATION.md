# BFINAL-031 预注册（先于对账数据）

> 写就时间：2026-08-22，任何系统性对账数字产出之前。开工时的输入核验意外发现
> （b8/b25 伴随态不一致、rxpr 导出时间戳）已先行记录于 NOTEBOOK §1-§2，作为假设来源
> 而非对账结果。本文件锁定恒等式、预期与判定树。

## 1. 对账恒等式（两标签 × 三方向）

对标签 L ∈ {J(thermalCoupling), gDP(pressureDrop)}，方向 k ∈ {D1,D2,D3}：

```
腿 P（生产 tsv）：  projL_k = Σ_active dgdx·Dk          [stage_b2_summary.tsv，scale=1 已核]
腿 R（离线重算）：  Σ_seg chain(seg_L) 投影 Dk          [seg 由 b25_qgate/1 场离线重算]
腿 S（源收缩）：    b_L^T w_true                        [b18rhs mtx × b26_wstar_h7.npz]

恒等式 A（装配 vs 自身分解，预期机器级闭合）：
   腿P == chain(Σ_seg) 投影   —— 生产中间场逐层自洽
恒等式 B（流介导源恒等式，B26a/B30 口径）：
   momentum_L + pressureRow_L（xh 层投影 z）≈ b_L^T w_true
恒等式 C（全梯度）：腿R(J) == FD_J；腿R(gDP) == FD_gDP（病灶判定的物理基准）
```

## 2. 三口径与 M4 第一道闸（先于真相）

三个仪器口径：
1. 生产 tsv 投影 = <post-chain 设计梯度, Dk>（Dk 为 max=1 原始方向）；
2. rxpr mtx = pre-chain xh 分段场（b8_verify_diag 导出，stageB6RxDesignOracle 开关门）；
3. 源恒等式 = b_L^T w_true。

**预注册判断（M4）**：z（stageB6_rxc_z_analytic）按构造是滤波切线+投影切线的精确
前向链向（oracle L966-998），与生产链式规则互为转置对，故 <pre,z> ≡ <post,D>
是数学恒等式——**"z 投影 vs D 投影"不构成独立口径差**。M4 的两个数量级差
（rxpr_total·z = (−0.056,+0.298,−0.036) vs ADJ_gDP = (−5.230,+1.051,+13.603)，
比值 93/3.5/374 非均匀 ⇒ 排除单一 scale 因子）的候选解释按优先级：

- **H1（主假设，口径差/导出态错位）**：rxpr 导出属于 b8 运行自己的伴随态
  （NOTEBOOK §2：Uc relL2=437、pb=uniform 0、时间戳同秒），而 ADJ_gDP 属于
  b25/B2 门轮伴随态。若逐 cell 复算 rxpr_momentum ≈ −dAlphaDxh·(U&Uc_b8)·V
  （relL2 ≤ 1e-10）而 ≠ b25 态复算 ⇒ H1 成立，M4 判**口径差（导出态错位）**，
  修正口径 = 以 b25 场离线重算分段为准；rxpr mtx 降级为"非 state-B 导出"。
- **H2（真差）**：若 rxpr_momentum ≈ b25 态复算但 <rxpr_total,z> 仍 ≠ ADJ_gDP，
  则生产投影链路（drho/eta/Helmholtz/dfdx 组装）存在真实缺陷 ⇒ 升级为病灶证据。
- **H3（z 失效）**：若 b25 态重算全链（含 Helmholtz 复算）投影 == ADJ_gDP 但
  <重算pre, z_file> ≠ 之，则 z 文件过期（设计态漂移）——仍属口径差，但需重造 z。

判定数据：逐 cell 比较（relL2、maxabs）+ 三口径投影表。**先出此表，再进第二阶段。**

## 3. 逐项对账矩阵（第二阶段）

对每标签每方向输出一行：
`projP | mom | prow | flux(仅J) | C(仅J) | chain后和 | 腿S=b_L^T w_true | FD | 各缺口`

预期（预注册）：
- 恒等式 A 预期机器级闭合（≤1e-10 rel）：生产中间场（fsenshMeanT/gsenshPressureDrop/
  fsensMeanT/gsensPressureDrop/dfdx on-disk）与离线重算逐层互证；若某层破裂 ⇒ 病灶层。
- 恒等式 B 预期闭合（b26a 已证 J/D3 5.7%；gDP 预期 **不闭合**，B30 缺口 0.82/0.37/1.04
  待归因：若 momentum+prow 重算和 ≈ 腿S 但 ≠ 腿P，缺口在链式/投影层；若重算和本身
  ≠ 腿S，缺口在分段装配或源折叠）。
- **D3 28× 预注册病灶位置**：J 链离线四段和（≈FD，b26a 已证）与生产 projJ 的缺口
  ~1.2e-2 必须由矩阵中某一项的异常解释；候选：fluxDirect 生产场 vs 离线 Gx、
  C 生产场 vs 离线 C、链式层放大、dfdx 组装额外项。
- **gDP 缺口预注册**：腿S/腿P = 0.818/0.366/1.036（B30 defcheck）。若 H1 成立后
  b25 态重算 mom+prow 与腿S 闭合而与腿P 有缺口 ⇒ 缺口在链式/投影层（与 M4 同源）；
  若重算 mom+prow 本身 ≠ 腿S ⇒ 分段装配层有独立病灶。

## 4. 裁决映射（判定树，审阅者版）

- 定位到具体项 → 写 B32 修复轮规格（不动生产；含 B13 分解式口径修订声明）；
- 断点在源收缩腿 → 回看 B30 defcheck 口径；
- 断点在生产 tsv 腿 → 仪器问题，回 B2 模块；
- M4 判口径差（H1）→ rxpr 系列从证据链降级，M4 开放差异关闭；
- M4 判真差（H2）→ 升级为装配链病灶主证据。

## 5. 门控（M0，先于一切对账）

复用 b26a/b27 已验收门：mobility <1e-12、V 金字塔 <1e-12。新增：
- 场读取门：b25/1 各伴随场尺寸/有限性检查；
- 导出态门：rxpr_momentum vs b8 态复算、vs b25 态复算的 relL2 双报告。

## 6. 预算与收尾

预算 180 分钟；单步 >15 分钟写 NOTEBOOK；数据齐即收尾。产物 `b31_` 前缀入本目录；
EXECUTOR_SUMMARY.md + EXECUTOR_DONE + 一次选择性提交（不 push）。
