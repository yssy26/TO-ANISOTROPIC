# ESI adjointOptimisation 模块不可达记录（参照降级声明，task tier 规则）

按任务书 tier 规则："网络抓取失败 >10 分钟即换源（clone/raw/mirror）或降级到下一个参照"。本节记录对参照系 2（OpenFOAM ESI adjointOptimisation）的抓取失败与降级决定。

## 尝试序列（全部失败，总耗时 < 15 min）

1. `git clone --depth 1 https://github.com/OpenFOAM/Plus.git` — 仓库在 GitHub 上不存在
   (error: "could not read Username for 'https://github.com': terminal prompts disabled")
2. `git clone https://develop.openfoam.com/Development/OpenFOAM-plus.git` — 成功，但为 maintenance-v1812 时代旧版，**不含** adjointOptimisation 模块（仅有 legacy/incompressible/adjointShapeOptimisationFoam）
3. `git clone https://develop.openfoam.com/modules/adjointOptimisation.git` — auth required（仓库非公开）
4. gitlab.com 各候选路径（openfoam/modules/adjointOptimisation 等）— HTTP 403 Forbidden（非公开仓库）

## 结论

- ESI adjointOptimisation 模块（github.com/OpenFOAM/Plus, v2312 中的热/流/孔隙率伴随模块）**无法从本网络环境公开访问**。
- 按 tier 规则降级：
  - tier-2 替代 1：Othmer 连续伴随 duct TO（legacy 随 OpenFOAM-dev 分发）→ `refs/openfoam_othmer_adjoint.md`（已抓）
  - tier-2 替代 2：Fira-Software/thermalTopO（GitHub 公开 CHT-TO 连续伴随，含显式热-动量折叠 gPhi）→ `refs/fira_thermalTopO_flux_sensitivity.md`（已抓）
- 这两者均覆盖了 ESI 模块原本要提供的两个关键面：(a) 孔隙率/Brinkman 伴随源符号（Othmer）；(b) 热目标经伴随温度折入动量源的面敏感度公式（Fira gPhi）。

## 剩余差异（对候选错项排名的影响）

- ESI 模块为离散伴随（连续伴随逆版本），其热折叠若可取得本可提供与我方 T1-T8/H7 槽同类的离散逐槽对照；替代品均为连续伴随，故**离散路由层（αRel·rAU·dH^T、kf、H7 压力通道）无直接逐槽参照**，该层对照仅能依赖内部自洽检查（M1/M2 恒等式）与 BFINAL-018 的槽位镜像论证。
