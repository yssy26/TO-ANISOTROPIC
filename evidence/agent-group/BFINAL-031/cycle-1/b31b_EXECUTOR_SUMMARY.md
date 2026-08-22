# BFINAL-031 attempt-2 (b31b) — m2 逐项对账 强制收尾报告

日期：2026-08-22 21:0x（协调人强制收尾令）
仓库：/home/ys/dsH/TO-ANISOTROPIC（分支 agent/dsH-stage-b-validation，HEAD 00a7f27）
范围：完全离线测量工作，零生产修改。本 summary 之外不写任何报告 .md。

## 一行状态
m2 对账矩阵未构建（b31b_m2.py 未写出，因强制收尾先于实现）；全部源码钉扎、代数推导、
数据核验与部分数字 salvage 已完成；压力行项代数与 G0 门根因已在 NOTEBOOK 定稿。

## 已完成（attempt-1 salvage + attempt-2 推导/仪器/核验）
1. 代数推导（NOTEBOOK §6.9-6.11 定稿）：
   - 压力行项 = −T·dAlphaDxh，T=T1+T2 面循环转置（rxPressureRowTranspose.H L223-257）：
     T1[own]+=wf·(Sf&dHbyA[own])·pcDiff，T1[nei]+=(1−wf)·(Sf&dHbyA[nei])·pcDiff；
     T2[own]+=−wf·g0·pcDiff·drAU[own]，nei 对称；边界仅 uAssignable=outlet：
     T1[cell]+=(Sf_b&dHbyA[cell])·λ[cell]，T2[b]=0 精确为零。
   - dflux_w[fi]=g0_int·interp(drAU·w)，边界=0（与 b18 g0 捷径同形）。
   - G0 门根因：生产是 (L_int−V·I)·xp=−V·x；attempt-1 错用 (L−I)xp=−x
     （relL2=2.620e-03，chain fieldrel 5.85/10.1）。修正：diag−V、RHS=−V·x_full、
     chain_apply=lu.solve(−V·mid)·designMask。
2. 数据核验（全部 VERIFIED）：b8==b25 当前态字段全同；pc8=stageB6_lambda.mtx[100800:]
   （L2=19831663.927）、Uc8=lam[:100800].reshape(33600,3)；rxprT=rxprP+rxprM（9.9e-18）；
   b8 门 anchors（Log.verify_diag.txt L2080）：sum(T*w)=−2.696485451749351e-08、
   relErr=2.57680498181e-15、w-support=5040；dAlphaDxh b8==b25。
3. H1 裁决（b31_m0m1_compare.tsv）：rxpr_momentum vs mom_recomp_b8 relL2=2.844e-12（MATCH），
   vs mom_recomp_b25 relL2=1.009（MISMATCH）⇒ rxpr 导出与 b8 运行 lambda 自洽。
4. M4 三口径投影（b31_m4_projections.tsv）：D1 dot_gsensPDdisk_D=−5.2302938463e+00（=ADJ_gDP）、
   dot_dfdx_D=−3.5835040049e-02（=ADJ_J）；D2 dot_gsenshDisk_z=−3.5247766837e+00（反号）；
   D3 dot_gsenshDisk_z=+8.0409993473e+01 vs dot_gsensPDdisk_D=+1.3603422582e+01。
5. attempt-1 跑日志（b31_m2_run.log）：G1 eta5=0.7480051615837531、G2 xh relL2=1.301e-10 通过；
   G0 失败已定位根因；崩溃点 b31_recon.py L430（(134400,)/(33600,) broadcast）。
6. 强制收尾 salvage（b31b_salvage.tsv/.json，纯复用缓存秒级）：
   - b_PD^T w_true = D1 −4.27746124e+00（xADJ_gDP 0.8178）、D2 +3.84462806e-01（0.3660）、
     D3 +1.40963533e+01（1.0362）——与 B30 defcheck 缺口 0.82/0.37/1.04 完全一致。
   - b_TC^T w_true = D1 −4.88401870e-02（xADJ_J 1.3629）、D2 −6.44099895e-03（0.2584）、
     D3 +1.15666175e-02（−0.9636）——J 标签单腿源收缩不能闭合对账（B13 流中介源）。
   - 动量段 pre 投影：mom_gDP preD = −3.128852e+00/+1.409900e+00/+2.193815e+00；
     mom_J preD = +2.044486e-02/−7.990468e-03/−1.206465e-02。
   - D3 28× 复核：FD_J(D3)=−4.4495118952e-04 vs ADJ_J(D3)=−1.2003012269e-02，比值 26.97×。

## 断点位置（如实声明未完成）
- b31b_m2.py 未写出；b31b_matrix.tsv（两标签三方向对账矩阵）未构建；
  b8 pc8 机器精度门未执行；defcheck 缺口与 D3 28× 的定量解释未完成（数值已知，见 salvage）。
- 续作路径完整写在 NOTEBOOK §6.12（8 步实现蓝图，含全部公式与参考行号），
  任何人可从笔记接手。

## 关键文件
- NOTEBOOK.md（§6.9-6.13，350 行）：代数推导、G0 根因、运行计划、当前状态与续作路径。
- b31b_salvage.tsv / b31b_salvage.json：本次 salvage 原始数字。
- b31_recon.py（attempt-1 库；m2() 部分行有 bug 勿复用，§6.12 注明）。
- b31_m0m1_compare.tsv、b31_m4_projections.tsv、b31_stage_cache.npz、b31_m2_run.log。
