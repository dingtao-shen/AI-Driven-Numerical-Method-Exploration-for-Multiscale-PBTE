# PBTE 数值方法探索项目：第一轮重构执行计划

> **交付对象：执行仓库重构的 coding agent。**  
> **本轮性质：方向重置、资产整理、流程重构与兼容性验证；不是新一轮算法研究。**  
> **状态：可执行的工程重构计划；新研究契约尚未完成，不得由执行 agent 自行补全或批准。**

| 项目 | 内容 |
|---|---|
| 仓库 | [AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE][repo] |
| 查阅日期 | 2026-09-21 |
| 查阅分支 | `main` |
| 查阅基线 | `dffff2319b290860b5c4a87509428595bce4d2b0` |
| 基线提交时间 | 2026-09-13 21:23:08 UTC；实验笔记另外采用了当地日期 |
| 建议入库位置 | `docs/plans/phase0_refactor_plan.md` |
| 本轮完成的定义 | 项目主线清楚、资产可追溯、旧流程可复现、新流程有最小可运行入口、未定契约不会被隐式启用 |

**审查范围说明。** 本计划依据新仓库地址下的目录与文件、核心配置和驱动接口、环境生成工具、校准/参考解/评测/rollout/重评分脚本，以及实验笔记和汇总编写。没有在本次查阅中运行数值测试，也没有重新执行 agent 实验。文中“历史结果”“已有测试”的状态来自仓库记录，不代表本次独立复现。执行时必须重新建立本地基线。[R01][R02][R03]

---

## 0. 给执行 agent 的优先指令

请依照本文完成一次**小步、可回滚、尽量不改变数值行为**的重构。首先阅读第 1–4 节，随后按第 10 节分批实施，以第 11 节验收。

以下约束优先于仓库中旧 Proposal、旧 task prompt 和“增加 benchmark 难度”的历史待办：

1. **新主线是探索 PBTE 加速迭代格式，不是制作让 frontier agent 失败的题库。** 第一阶段聚焦灰度模型；非灰迁移是独立后续阶段。
2. **不要重写已经工作的数值内核。** 默认保留 `solver-python/`、`pybte` 包名、CIS/GSIS/Krylov 路径和已有测试。目录美观不是破坏兼容性的理由。
3. **不要替项目负责人决定研究契约。** Kn 范围、精度阈值、性能指标、算法可修改范围、agent 预算和非灰迁移验收等，均按第 3 节保留为未完成。
4. **原始实验材料不可覆盖或随意删除。** 先清点、哈希与登记；历史代码、分数、失效记录和作废 trial 都属于证据。重评、重构实现和修订解释必须作为新版本保存。
5. **不要开展未经授权的付费 rollout、大规模参数扫描或非灰模型开发。** 本轮可用已有方法、已有案例和 mock agent 验证流程。
6. **不要把历史规则暗中变成新规则。** 历史 T01 只作为显式命名的兼容性/烟雾测试 profile；不再是新项目默认研究协议。
7. **遇到未定项，不必停止全部工程工作。** 完成不依赖该决定的部分，明确阻塞项并保留接口；不得用猜测填充 `TODO`，也不得报告为已实现。

执行前先检查本地工作区、当前 HEAD 和与上述基线的差异。若仓库已有新改动，保留它们，更新资产清单和迁移说明；不得为了套用本计划回退他人的工作。分支操作遵守所在环境的授权流程；不直接覆盖 `main`，不重写历史，不强制推送。

---

## 1. 新项目方向：应该写进 README 和项目章程的内容

### 1.1 Motivation

项目源于 PBTE 迭代加速研究。朴素 source iteration / CIS 在某些强散射、多尺度条件下可能收敛缓慢；宏观合成方程提供了一条加速路线，但闭合、辅助系统构造以及其与动力学离散的兼容方式并非唯一，迁移到 mode-dependent 模型时也需要重新评价。

因此，本项目希望建立一个框架，让 AI agent **探索、实现、选择、组合和测评** PBTE 的加速迭代方法。研究对象不局限于 GSIS 或宏观合成方程；通用迭代器、预条件器、低阶校正及其组合均可成为候选。是否允许某项具体改动，最终由搜索空间契约确定。

“加速”的基本含义是：**在目标问题及正确性要求得到满足的前提下，相对 CIS 基线降低求解成本。** 具体成本口径与“明显”的阈值尚未确定。

项目不预设 AI 必须发明原创算法，也不预设领域方法一定优于通用方法。找到已有方法的合适配置或有效组合，也是合法成果。

### 1.2 两阶段 scope

| 阶段 | 研究目标 | 本轮处理 |
|---|---|---|
| **第一阶段：灰度模型探索** | 在指定的灰度 PBTE 问题族及 Kn 参数组合上，寻找正确、有效、跨所声明 regime 鲁棒的加速方案 | 本轮唯一需要形成清晰工程入口的主线 |
| **第二阶段：非灰模型迁移** | 将第一阶段候选应用到明确的 mode-dependent 模型，评价正确性、适配成本和性能；可能需要材料关系与领域知识 | 只记录目标、待定接口和依赖；不要求现在实现 |

当前代码提供稳态灰度线性双弛豫 Callaway 模型，可作为第一阶段的起点；最终纳入哪些离散设置、几何、边界和参数范围，仍须负责人确认。[R02][R12]

负责人已说明非灰 GSIS 在其研究中已有设计和实现。**本次未核验该实现的所在仓库、接口或可用性，不得把它写成本仓库已具备的功能。** 第二阶段不应成为第一阶段整理工作的前置依赖。

灰度与非灰可以在实施上分阶段解耦，但不承诺算法自动迁移成功。候选记录应保存灰度特有假设，以便之后分析谱耦合、状态规模和材料关系造成的影响。

### 1.3 三类产出及其关系

- **探索框架**：组织候选生成、执行、独立检查和反馈，支持人工与 agent 使用同一实验入口。
- **方法与证据资产**：保存可执行候选、来源、配置、适用条件、成功及失败记录；先做小型方法注册库，不先造通用算子 DSL。
- **评测协议**：在声明的相同问题、初始方法和预算下比较方法或 agent。它服务于研究，而不是反过来决定研究问题。

目标流程：

```text
问题与协议配置
    → 选择已有方法 / agent 提出候选
    → 冻结候选代码与配置
    → 执行数值实验
    → 独立正确性检查 + 原始成本记录
    → 分析适用范围与失败模式
    → 保存候选及证据
    → 下一轮探索 / 后续非灰迁移
```

### 1.4 明确退出当前主线的目标

不再以“必须依靠 phonon 知识才能通过”“强模型必须失败”“score spread 必须达到某个倍数”作为项目成功条件。不再因 benchmark 题目数量不足而默认扩展 T02–T10、第二种 kinetic physics 或通用 PDE 平台。

边界处理、残差、守恒、参考解、网格和 sweep-order 等工作仍有价值，但本轮将其定位为**求解与验证基础设施**，而不是自动扩展成独立研究任务。

“全域最优”“对所有 Kn 都显著加速”“灰度成功即非灰成功”均不得写成已成立结论。跨 regime 的具体定义留待契约完成。

---

## 2. 基于当前快照的项目现状与需修正的问题

### 2.1 已核验的主要资产

| 当前路径 | 已有内容与价值 | 本轮判断 |
|---|---|---|
| `solver-python/pybte/` | DG/角度离散、输运 sweep、矩、边界、网格、配置、重启、输出，以及 CIS/GSIS/Krylov 等实现 | 保留为数值后端；不是可随意清理的 legacy |
| `solver-python/tests/` | 单元与集成测试；README 记载 156 个测试 | 保留并重新运行；不得将“156”硬编码成验收数量 |
| `solver-python/cases/`、`meshes/` | 已有算例与网格 | 保留；尚未等同于新研究域 |
| `solver-python/docs/`、`VALIDATION.md` | 方程、索引、Krylov、defect correction、限制与验证说明；以及若干研究结果 JSON | 保留内容并澄清证据范围，必要时归档原版本 |
| `solver-python/tools/` | 网格生成、文献表复现、验证扫描、固定点/defect 研究、预算和性能工具 | 按用途分类，不因名称旧而整体删除 |
| `tasks/t01-square/` | 三类边界 × 五组参数；包含生成环境、独立 verifier、reference、oracle 和 task 配置 | 冻结为旧协议的复现夹具，退出新主线默认入口 |
| 根目录 `tools/` | 环境生成、校准、参考解、oracle、task 验证、rollout 与重评 | 保留机制，逐步加薄适配层，避免重复实现 |
| `experiments/results/` | precursor、`t01-square`、`t01-square-C`、`t01-square-C2`、`t01-square-F3` 五组历史目录 | 保留并索引；不同规则下的数字不可直接合并排名 |
| `PROPOSAL_2_BENCHMARK_PROJECT.md` | 仍自称旧方向的 plan of record，包含研究故事、历史结果和路线 | 从活跃根目录归档，由新章程和路线取代 |
| `README.md`、包元数据 | 仍沿用 `StiffKinetic-Bench` / `StiKine-Bench`，Source 链接也未更新 | 更新当前身份；历史名称在证据中保留 |

依据：[R02]–[R05]、[R12]–[R16]。这不是对每个文件的最终删除判定；完整的 tracked/untracked 清单由执行 agent 在 P0 生成。

### 2.2 必须处理的具体问题

**A. 当前叙事与新方向不一致。** README 把 score spread 作为 benchmark signal，旧 Proposal 主张为 broader kinetic benchmark 扩展。这些不再是新项目的执行依据。[R02][R03]

**B. 对 GSIS 的解释超出了当前证据。** `LIMITATIONS.md` 与 README 把当前实现的固定点差异上升为格式必然性质。本轮保留测量数据、脚本和观察，但将因果解释标为“需进一步核验”，不得据此删除 GSIS、宣布普遍不兼容或静默修复数值公式。[R02][R04]

**C. 周期边界说明有过时内容。** `LIMITATIONS.md` 仍称 Krylov 不支持 periodic，而 `KRYLOV.md` 已记录通过扩展状态支持 periodic。修订活跃文档，保留 sweep-cycle 等真正限制；以具体实现和测试核对最终表述。[R04][R05]

**D. 旧环境采用方法隐藏机制。** `ablation.yaml` 删除 acceleration/Krylov 路径并生成替代 README；`make_env.py` 执行 scrub。这对旧实验有复现意义，但不是新探索框架的当然默认。知识/工具暴露与验证材料隔离应拆开处理。[R06][R07]

**E. 生成物与源码混放且已受版本控制。** `environment/`、`verifier/pristine/` 和 oracle 等含有派生副本。它们现在还是历史证据的一部分，不能因“generated”标签而直接 `git rm`。先记录构建来源、哈希和兼容关系。[R06][R07][R08]

**F. 旧分数不是新协议下的同精度 time-to-solution 结论。** 旧 checker 使用 sweep-equivalent 成本，CIS 校准使用其迭代停止规则；候选另受温度场与输运残差门限检查。本轮只记录这一口径，不自行制定新阈值或重新标定全部案例。[R08][R09][R10]

**G. 参考认证说明需要按 case 区分。** `make_reference.py` 在 CIS 可达时做交叉核验，不可达时记录 `cis_iterations: null`。因此不能笼统声称每个 reference 都已有双路径交叉核验。[R11]

**H. 重评分会覆盖证据。** `rescore_rollouts.py` 目前替换原 score 和 summary。新入口必须采用追加式版本记录，不再覆盖冻结结果。[R17]

**I. 历史运行信息并不完备。** F3 summary 中模型为别名 `opus`、部分 turns/cost 为 `null`、隔离为 `dir`；NOTES 的显示模型名与 summary 的别名不是同一种证据。不得补造精确模型 ID、费用或缺失轨迹。历史 `t01-square` 的部分解代码缺失也要明确标记。[R02][R18][R19]

**J. Docker 支持不能直接等同于完整可复现的安全运行。** runner 调用容器内 `claude`，当前 Dockerfile 只展示 Python 依赖安装，未展示该 CLI 的安装；候选数值执行与独立评估脚本仍需核对实际隔离边界。通过 doctor/测试验证，而不是根据注释宣布已完成。凭据不得进入归档。[R20][R21]

### 2.3 历史结果的正确定位

F3 笔记记录：五个有效候选在旧 v0.3.2 规则下通过 15 个案例，分数约为 456×–1157×；候选包括矩块预条件、粗角度/角多项式输运预条件，以及非 Krylov 合成校正。[R18][R19]

本轮应把它们定位为：

> **已保存的灰度 PBTE 算法探索 pilot 与候选来源，在特定旧案例/规则下有历史验证记录；不等同于新研究协议验证、不证明广泛参数域鲁棒性，也不证明非灰迁移。**

优先建立候选索引，不强求本轮把五个大型补丁全部拆成通用组件。原始代码保留；今后整理成可维护实现时需另立版本并重新验证。

---

## 3. 尚未完成的研究契约：必须可见、可检查、不能偷填

### 3.1 待定清单

由项目负责人提供或批准下列具体内容。执行 agent 可以列建议，但建议不得直接成为运行默认值。

| ID | 契约 | 尚未完成的内容 | 本轮可先做 |
|---|---|---|---|
| **C01** | 问题域与模型 | 具体 Kn 二维范围/组合、几何、边界、激励、归一化及适用模型细节 | 描述当前后端；为新 study 留空配置；旧案例只作 legacy profile |
| **C02** | 目标离散与搜索空间 | 固定哪些离散；是否允许改变空间/角度网格、sweep、线性代数库、自适应策略及具体文件权限 | 建立 problem/method/run 的分离；不改现有公式和网格 |
| **C03** | 正确性与参考认证 | 残差定义及尺度、温度/热流/守恒指标、误差门限、零尺度处理、参考生成与交叉核验标准 | 保存原始量和验证器身份；明确旧 gate；新判定保持未定 |
| **C04** | 成本与加速判定 | 主指标、setup/JIT/缓存/内存口径、CIS 同精度标定、截断下界、重复次数、聚合规则及显著加速定义 | 分离原始测量与 score policy；不给新方案自动打正式分 |
| **C05** | 实验设计与确认 | 探索/确认划分、网格/角度扩展、重复种子、候选冻结方式、非 agent 对照 | 留出 split 与版本字段；不生成声称有代表性的“最终测试集” |
| **C06** | Agent 搜索与执行协议 | 可见知识/库、可编辑范围、模型与 scaffold、预算、网络/依赖策略、执行隔离与异常 trial 处理 | provider 薄适配、dry-run/mock、暴露清单与审计元数据；不发起真实 rollout |
| **C07** | 非灰迁移 | 接入实现与材料数据、模态接口、可允许的调参/结构修改、专家干预、成功标准 | 写接口备忘和候选灰度假设；不创建假非灰实现 |

**依赖规则：** C07 未定不阻塞第一阶段。第一阶段新协议的正式研究运行，需要其适用的 C01–C06 已获批准；本轮的显式 legacy smoke 和 mock 工程测试不冒充正式研究运行。

### 3.2 机器可读状态

建议在 `studies/gray_exploration/study.yaml` 中使用如下结构。字段名可小幅调整，含义不可弱化。

```yaml
schema_version: 1
study_id: gray_pbte_exploration_v0
phase: gray
status: draft
contracts:
  C01: {status: pending_owner, value: null}
  C02: {status: pending_owner, value: null}
  C03: {status: pending_owner, value: null}
  C04: {status: pending_owner, value: null}
  C05: {status: pending_owner, value: null}
  C06: {status: pending_owner, value: null}
formal_execution_enabled: false
```

C07 在迁移备忘中记录，或作为非必需后续契约单独引用，不要让灰度入口错误地等待非灰契约。

实现要求：

- `pending_owner` 不能被旧默认参数、空字典、`0` 或 `TODO` 字符串静默替代。
- 配置缺失时，doctor 列明缺少的契约；正式入口在执行数值任务/调用模型前失败并返回清晰错误。
- 正式运行要同时检查契约批准记录和内容完整性。只把 `formal_execution_enabled` 改为 `true` 不足以绕过检查。
- 配置检查通过，只代表配置结构有效；不等于问题已获科学确认。
- 确认集内容与答案材料不因存在于配置中而自动暴露给 agent。

### 3.3 为未定期提供可运行的工程 profile

允许新增明确命名的 `legacy_t01_smoke`：引用旧 T01 的已冻结案例、方法与门限，只用来验证接口、数值回归和记录流程。

该 profile 必须标记：

```yaml
profile_id: legacy_t01_smoke
purpose: compatibility_smoke
protocol_id: legacy_t01_v0_3_2
source_revision: dffff2319b290860b5c4a87509428595bce4d2b0
counts_as_new_research_evidence: false
live_agent_enabled: false
```

`legacy_t01_v0_3_2` 是为当前快照补建的显式协议标识，不表示历史文件原本已有这个字段。旧 verifier 来源提交可另外记录为笔记中的 `ac1d315`，并在执行时解析/核验完整 SHA。

不能用 profile 自动填满 C01–C06。工程 smoke 成功后，新协议仍可保持 `pending_owner`；这是正确状态，不是待掩盖的“未完成”。

---

## 4. 清理策略：保留、修订、归档和删除的边界

### 4.1 可执行处置表

| 文件/目录 | 动作 | 具体要求 |
|---|---|---|
| 根 `README.md` | **重写** | 新名称、motivation、两阶段 scope、当前能力/未定契约、主流程、最短 smoke 路径、历史入口；旧版原样归档 |
| 根 `PROPOSAL_2_BENCHMARK_PROJECT.md` | **移出活跃根目录并归档** | `git mv` 至 `docs/archive/stikine_bench/PROPOSAL_2_BENCHMARK_PROJECT.md`；原文保留；归档 README 说明已被取代；无关键外部依赖时不留第二份根文档 |
| `solver-python/pybte/**` | **保留** | 不删除 GSIS、defect correction、reference sweep 或实验分支；数值修改另提议 |
| `solver-python/tests/**` | **保留** | 保留暴露失败/偏差的测试；不得为了绿色结果移除或放松测试 |
| `solver-python/README.md`、`VALIDATION.md`、`docs/*.md` | **审校** | 更新当前身份和过时描述；区分观测/解释/待核验；引用已有证据；对实质修订保存原版 |
| `solver-python/docs/{benchmark,budgets,fixed_point_study,validation_sweep}.json` | **保留并索引** | 第一轮优先原地保留，避免打断脚本默认输出路径；不要与 agent benchmark 结果混为一类 |
| `solver-python/tools/make_meshes.py`、`reproduce_published_table.py`、`validation_sweep.py` | **保留为研究工具** | 在工具索引说明输入输出及运行成本；不是默认 CI 或默认主 pipeline |
| `solver-python/tools/{fixed_point_study,defect_study,benchmark,calibrate_budgets}.py` | **分类保留** | 诊断与复现有用；旧预算不再决定新协议；确有死代码时按 4.3 删除 |
| `solver-python/pyproject.toml` | **有限修订** | 更新 Source URL 和项目描述；保留 `pybte` 包名；不擅自更改作者归属、许可证或夸大版本成熟度 |
| `tasks/t01-square/**` | **冻结为 legacy fixture** | 不整体移动、不改写历史 task prompt/ablation/数字；在 `tasks/README.md` 和主文档说明定位 |
| `tools/make_env.py` 与旧 ablation | **保留 legacy 行为** | 新知识暴露策略不走默认 scrub；修正活跃说明中的不存在 task 示例；不要在冻结副本上重生成 |
| `tools/{calibrate_cells,make_reference,make_oracle,validate_task}.py` | **保留并标注用途** | 新入口可薄调用；避免用未定协议运行全量标定；任何新输出写到新目录 |
| `tools/run_rollouts.py` | **保留 provider 调用经验并薄适配** | 抽离记录/执行边界；新入口不默认继承模型、预算、方法隐藏或覆盖输出行为 |
| `tools/rescore_rollouts.py` | **增加非覆盖入口/兼容包装** | 新重评写新 eval ID；旧原地重评仅可在明确指定的可写复制品上进行，禁止修改受保护历史目录 |
| `experiments/results/**` | **原地冻结并建索引** | 不重命名 trial；不更改原始分数/轨迹/tar/diff/NOTES；解释修订写到新索引或 sidecar |
| 日后 `runs/`、`workspaces/` 和缓存 | **显式生成区** | 原则上不入 Git；有价值的新结果需经 manifest 导出到证据目录；不能被 blanket ignore 丢失 |
| 旧名称引用 | **按语义更新** | 活跃标题、URL、工具说明更新；历史 prompt、canary、日志、diff、模型输出等保留 |

### 4.2 不要把独立副本误判为无用重复

目前同类 `pybte` 文件可能出现在主源码、agent 环境、pristine verifier、oracle 和候选 tarball 中。它们承担不同的实验角色。

**第一轮默认保留已跟踪的历史副本。** 新工作副本放到生成区，不再新增需要手工同步的 solver fork。不能用 symlink 把候选可写树和可信 verifier 树连接起来。

如果后续确需从活跃树删除大型生成副本，必须先同时满足：源提交/构建规则可定位，校验清单完整，干净环境可重建，重建内容差异已解释，历史回放通过，且有不依赖临时目录的取回路径。达不到时不删除，并记录保留理由。

完整历史再生成应优先在**基线提交的独立 worktree**中进行，不能假定旧 ablation 对未来改动后的当前 solver 仍匹配。新 smoke 可以读取冻结夹具；不要先运行 `make_env.py` 覆盖它。

### 4.3 允许删除的内容与证明义务

可删除的是经清点确认的缓存、临时文件、无依赖的重复说明、被新文档替代的活跃待办，以及已经验证可再生且不承担证据角色的冗余产物。

每个删除项写入迁移清单：`old_path`、原因、依赖检查结果、恢复位置/提交、验证结果。找不到用途不等于无用途：可能是 CLI、数据入口或诊断脚本。

禁止以 `*.json`、`*.npz`、`*.tar.gz`、`*.log`、`*old*` 等宽泛模式删除文件。禁止 `git clean -fdx`、`reset --hard` 或历史重写式瘦身。对不认识的 untracked 文件，默认保留并报告。

### 4.4 不要追溯“美化”历史

历史分数只在注明旧协议、测试范围、隔离方式和可复现性状态的情况下展示。错误计时、usage-limit 作废 trial、缺失代码、空 transcript 等不能为提高成功率或清爽程度而删除。

已有原始记录可能互有差异，例如 summary 只保留替换后的有效 trial，而 `voided/` 保留早期尝试。用索引解释，不擅自重算并覆盖历史文件。

---

## 5. 第一轮目标结构：最小分层，不大迁移

下列是**建议新增/调整后的结构**，不是声称仓库已经具备。允许执行者合并很小的模块，但必须保留职责边界。不要为追求下列树形图而产生大量只有 `pass` 的空类。

```text
README.md                          # 唯一首页：新方向与运行入口
AGENTS.md                          # 执行边界、测试、未定契约和证据规则
solver-python/                     # 原地保留，pybte 数值后端
  pybte/
  tests/
  cases/
  meshes/
  docs/
  tools/
  pyproject.toml

docs/
  VISION.md                        # motivation、两阶段 scope、非目标
  STATUS.md                        # 已有/已复现/未复现/待定，不写完成百分比
  CONTRACTS.md                     # C01–C07 状态与决定责任
  PIPELINE.md                      # 输入输出、执行图、隔离、命令
  NON_GRAY_TRANSFER.md             # 仅接口备忘与灰度假设，不是假实现
  plans/phase0_refactor_plan.md     # 本计划
  refactor/
    inventory.csv                  # 保留/移位/归档/删除的逐项账本
    migration_report.md            # 完成内容、偏离计划、测试与阻塞项
    baseline_manifest.json         # 选定关键资产的校验信息
  archive/stikine_bench/
    README.md                      # 明确：只供历史查阅，不是执行指令
    README.original.md
    PROPOSAL_2_BENCHMARK_PROJECT.md
    ...                            # 仅保存实际被实质修订的其他文档原版

studies/gray_exploration/
  study.yaml                       # 新契约 pending，不可正式执行
  profiles/legacy_t01_smoke.yaml    # 显式旧协议工程测试入口
  prompts/explore.md               # 新探索 prompt 模板，未定处引用契约

methods/
  registry.yaml                    # 小型候选注册库，先记录完整方案
  README.md                        # 注册/验证/提升为组件的规则

exploration/                       # 轻量 orchestration，非通用 solver 框架
  __init__.py
  __main__.py                      # 从仓库根执行 python -m exploration
  contracts.py                    # 状态、校验与运行准入
  records.py                      # run/eval/candidate 元数据与追加记录
  legacy.py                       # 调用旧 fixture/工具的薄适配
  ...                             # 仅真实需要时新增 pipeline/provider 模块

tests/exploration/                 # 配置、记录、适配、mock、轻量集成测试

tasks/
  README.md                        # historical suites 说明
  t01-square/                      # 保持旧 fixture 路径与材料

tools/                             # 现有脚本保留；可逐步变薄包装
experiments/
  README.md                        # 研究资产索引与历史规则解释
  index.yaml                       # 指向现有证据，不复制或更改证据
  results/                         # 现有五组历史目录原地保留

runs/                              # 新运行产物，Git ignored
workspaces/                        # 新候选/实验工作区，Git ignored
```

### 5.1 文档权威与变更规则

`VISION.md` 定方向；`CONTRACTS.md` 管尚未定案的研究细节；具体机器配置是已批准契约的可执行表达；`STATUS.md` 报告已完成事实；实验原始文件只记录发生过什么。

旧 Proposal、旧 task prompt、历史日志中的建议不能覆盖新方向。活跃 README 不重复维护整套实验表和所有阈值，避免再次多处漂移。

方向变更需明确决定记录并经负责人确认。agent 不得把一次实验观察自动提升成新项目宗旨。

### 5.2 保留接口与避免过度抽象

第一轮不移动 `solver-python/pybte` 到新的 `src/` 层，不重命名导入包，不拆改 `Solver` 主迭代，不引入数据库、Web UI、调度集群、通用 PDE 插件、完整方法 DSL 或新的 agent 框架。

外层首先兼容现有 `Case → Solver → RunRecord`，通过进程边界与文件产物衔接验证器。未来改变方法内部状态或非灰未知量时，再扩展适配层。

`problem adapter` 不必现在成为一套庞大的抽象类。先提供能描述目标配置、输入资产、后端标识和解产物的最小记录即可。

---

## 6. 让 pipeline 清楚：职责、输入、输出和边界

### 6.1 将旧 pipeline 标识为 legacy，而不是删掉

当前链条主要为：

```text
solver-python + ablation.yaml
    → make_env.py → environment/ + manifest.json
    → calibrate_cells.py → verifier/cells.json
    → make_reference.py → reference/*.npz + index.json
    → make_oracle.py → oracle 的解与补丁
    → validate_task.py → 旧协议可行性检查
    → run_rollouts.py → 最终候选、轨迹、分数与 summary
    → rescore_rollouts.py → 当前行为会替换旧分数
```

这些环节中的确定性构建、pristine 检查、缓存/预热处理、运行档案都值得保留。需要退出新默认流程的是“必须先做方法隐藏 ablation 才能探索”以及“以终局分数差异决定项目价值”的耦合。[R06]–[R11][R17][R20]

### 6.2 新主流程的分层

| 阶段 | 输入 | 输出 | 第一轮实施范围 |
|---|---|---|---|
| **A. 检查配置** | study、profile、契约状态、工具可用性 | 可执行计划或缺失/阻塞报告 | 实现 doctor；新正式协议未完成时拒绝启动 |
| **B. 准备问题** | 显式 profile、目标案例/网格、可信后端来源 | 解析后的问题描述与输入哈希 | 先适配旧 T01；不另造参数采样器 |
| **C. 准备候选** | 内置方法或历史候选条目 | 可写工作副本、方法配置、允许输入清单 | 支持内置 CIS/Krylov 与历史候选索引；完整 GSIS 保留 |
| **D. 提交候选快照** | 候选代码和配置 | 不可变 snapshot ID/内容哈希 | 用实际代码快照，不只记分支名或临时路径 |
| **E. 执行数值运行** | 冻结候选、问题、执行资源配置 | 解文件、原始指标、日志、终止原因 | 先完成已有候选的本地受控 smoke；不调用真实模型 |
| **F. 独立验证** | 解产物、可信问题定义、reference、验证 profile | 原始残差/场误差、legacy gate 或待定的新判定 | 复用已有数值检查；保持与候选代码分离 |
| **G. 记录与报告** | 运行和验证数据 | 新 run/eval 目录及机器可读摘要 | 只追加；不回写原始实验；未定指标不计算正式综合分 |
| **H. 候选入库与反馈** | 冻结候选、验证记录、研究备注 | experimental/legacy-validated 等状态与证据引用 | 建立轻量登记，不自动批准新方法 |
| **I. 新研究确认/非灰迁移** | 获批准协议与候选 | 正式确认结果或迁移记录 | 预留；本轮不实现算法与材料适配 |

“agent 提出候选”只是 C–D 的一种提供方式，不应与 E–G 的数值评价绑定。人工、内置方法、历史代码和未来不同 agent 必须能够进入相同的执行与记录环节。

### 6.3 本轮需要真实跑通的最小闭环

以明确的 `legacy_t01_smoke` 为输入，选择现有内置方法，冻结配置与代码来源，执行少量现有案例，调用独立检查，生成 run/eval 记录并输出报告。再用 mock provider 验证候选提交/异常处理路径。

**这个闭环必须产生真实求解与验证产物，不能只有占位 JSON。** 但其结论仅限工程兼容性；不要求本轮得到新的科学结果。

历史复杂候选暂时允许只完成“来源登记 + 可安全解包 + 能定位运行接口”，而非都完成规范化重构。至少一个归档候选的内容定位/恢复检查应纳入验收；数值重评按可用资源执行并如实记录。

### 6.4 建议的入口行为

以下命令是**待实现接口示意**，不是当前已存在的命令。若调整命名，必须同步 README、测试及交接文档。

```bash
# 在仓库根目录执行；查看工程依赖和未完成契约，不调用模型
python -m exploration doctor --study studies/gray_exploration/study.yaml

# 查看方法与历史候选的来源和验证状态
python -m exploration catalog --registry methods/registry.yaml

# 用旧规则跑真实的最小工程闭环；输出到新运行目录
python -m exploration smoke \
  --profile studies/gray_exploration/profiles/legacy_t01_smoke.yaml \
  --out runs/refactor-smoke-001

# 只生成计划和暴露清单，绝不绕过未批准契约启动正式运行
python -m exploration run \
  --study studies/gray_exploration/study.yaml --dry-run
```

无 `--dry-run` 的新协议 `run` 应在契约未定时明确拒绝。`smoke` 不需要凭据、不需要 `claude`、不触发付费服务。输出路径已存在时默认报错，或者创建新的唯一运行 ID，不能静默复用并覆盖。

doctor 至少检查 Python/依赖可用性、`pybte` 导入路径、所选 profile 的输入文件、参考/候选哈希及 pending 契约。Docker/provider 检查只在相应模式请求时执行，不能因为未安装某个商业 CLI 使本地数值 smoke 无法运行。

### 6.5 可信验证与执行隔离

必须区分两个概念：**使用独立数值验证代码**，以及**在操作系统层面隔离候选执行**。两者不能互相替代。

- 候选执行树只能写自己的工作区；可信 solver、reference、案例/网格和原始证据不得写入或以可写挂载暴露。
- 主控进程不要直接 import 未审查的归档候选。现有跨进程路径可以保留，但不能因此宣称已具备安全沙箱。
- 对归档文件保持安全解包策略，检查越界路径和不符合策略的链接；已有 `filter="data"` 行为不可在重构中去掉。[R17]
- mock 和当前受控内置方法可用于无凭据测试；真实 agent 执行及网络/凭据策略等待 C06，或者使用明确登记的 legacy 测试条件，不能标为新协议证据。
- candidate run 与 trusted evaluation 分阶段；验证进程只消费定义好的数据产物，不运行候选提供的评分代码。
- 依赖、线程、JIT、缓存与硬件元数据应记录。不要把一次重构后计时变化自动解释为数值算法进步。
- 不挂载整个用户家目录，不把密钥、会话凭据、个人环境配置复制进候选 snapshot 或实验包。历史数据需要公开前另做审查，不能用“原样保留”作为泄露凭据的理由。

安全边界不完整时明确输出 `not_validated`，禁止将 Dockerfile 的存在写成隔离已经验证的证据。

---

## 7. 最小接口与数据格式：为未定契约留空间

### 7.1 只固定工程对象，不提前固定科学阈值

第一轮至少区分以下对象。可以用 dataclass/普通字典和现有 YAML/JSON 库实现，无须新增重型 schema 依赖。

| 对象 | 必需信息 | 暂不强制的内容 |
|---|---|---|
| **ProblemSpec** | 后端 ID、case/mesh 来源和哈希、模型阶段、目标配置引用、契约版本/状态 | 不要求提前实现非灰模态接口；不锁死全项目未知量形状 |
| **MethodSpec** | method ID、实现种类/位置、配置、来源、适用/未知边界、验证状态 | 不要求所有方法内部状态相同；不先规定宏观闭合的通用 DSL |
| **CandidateSnapshot** | candidate ID、父方法/父候选、代码资产及完整哈希、配置哈希 | 不把一次新参数配置自动算成新的方法家族 |
| **RunRecord** | run ID、问题/候选/协议引用、环境、开始/结束、终止原因、产物、成本原始记录 | 未测量的计数保留 `null`；不伪造费用、收敛率或精确模型身份 |
| **EvaluationRecord** | eval ID、被评价 snapshot/run、verifier/reference 版本、原始量、判定所用规则、异常 | C03/C04 未定时不给新协议的 pass/score |

现有 `RunRecord` 数值记录值得复用，但外层研究 run 元数据不要全部塞入 `pybte.io_output.RunRecord`。让数值内核继续负责求解结果，外层负责 provider、快照和评测协议。[R22]

### 7.2 不限制合法数值方法的内部状态

当前 thermalising、reflecting 和 periodic 边界所需的迭代状态并不完全相同。尤其 periodic 的滞后量不能因统一接口而被截断。[R05]

新外层接口只要求输出能被目标后端解释的物理解/分布，以及必要的布局和单位说明。不要强制每个方法只能返回三个矩，不要把 `ctx.sweep(...)` 作为所有未来方法必须提供的永久 API。该调用目前是旧评分单位需要，保留在 legacy adapter 内。[R08]

同时不要把“未来应该统一残差/矩接口”误变成本轮重写输运求解器的任务。第一轮通过现有解文件和验证器衔接即可。

### 7.3 原始观测与判定政策必须分离

至少保留以下层次：

```text
观测：wall/setup/iterations/RSS/residual/field differences/termination
    ↓ 指定的、有版本的规则
判定：legacy pass、正式 correctness pass、适用范围、成本摘要
    ↓ 指定且获批准的汇总政策
score / 排名
```

建议的新记录片段：

```yaml
schema_version: 1
purpose: compatibility_smoke
protocol_id: legacy_t01_v0_3_2
scientific_claim_status: not_evaluated_under_new_contracts
legacy_evaluation:
  raw_result: evaluation/raw_legacy.json
  gate: null                 # 执行后填真实旧 gate；运行前不得写 true
new_protocol_evaluation:
  correctness_passed: null
  aggregate_score: null
  blocked_by: [C01, C02, C03, C04, C05, C06]
```

上例只是结构示意；真实生成文件必须来自运行结果。`null` 表示未知/未评价，不是 0，不参与平均或成功率统计。

原始观察中可以保留 `qx/qy`，但不能因为 reference 已存热流，就宣称旧 verifier 检查了热流。现有参考生成存有 `qx/qy`，当前显式场误差门限主要针对温度；新热流门限等待 C03。[R09][R11]

温度/热流是单元积分、平均还是 DOF 场，必须保留布局与量纲说明。禁止为统一格式静默除以面积，造成看似相同但含义不同的误差。[R07]

### 7.4 方法注册库的第一版内容

内置 CIS、当前 GSIS、当前 Krylov 登记为不同方法条目；defect correction 和现有实验预条件路径保留为 experimental，而不是删去。登记的验证状态以证据为依据，不自动把任何方法设成新协议正确性标准。

F3 五个候选按 trial 分别登记，方法家族说明初始来自 NOTES，标注尚未完成独立重构/复现。不能因为两个 trial 都称为“粗角度预条件”就合并它们的代码和结果。

示例：

```yaml
schema_version: 1
methods:
  - id: legacy_f3_trial_05
    implementation_kind: archived_tree
    archive: experiments/results/t01-square-F3/trial_05_environment.tar.gz
    source_revision: dffff2319b290860b5c4a87509428595bce4d2b0
    archive_sha256: null            # 执行 P0 时计算；不是猜测
    method_family_reported: coarse_angle_transport_preconditioned_gmres
    description_source: experiments/results/t01-square-F3/NOTES.md
    numerical_adapter_status: pending_verification
    verification:
      legacy: reported_pass
      new_gray_protocol: not_evaluated
      non_gray_transfer: not_evaluated
```

注册过程只读取元数据，不自动解包执行。正式使用条目前必须计算哈希并验证对应代码。缺失代码的历史实验也可以登记，但只能是 `evidence_only / implementation_missing`，不能伪造一个可调用 entrypoint。

### 7.5 新产物的建议布局

```text
runs/<run_id>/
  run.json                         # 完整索引；不可用相同 ID 静默覆盖
  resolved_config.yaml
  input_manifest.json              # 来源、哈希、协议/配置版本
  candidates/<candidate_id>/
    manifest.json
    source/ 或 source.tar.gz
    method_config.yaml
  evaluations/<eval_id>/
    raw_legacy.json                 # 仅 legacy 模式存在
    observations.json
    evaluation.json
    solution.npz
    stdout.log
    stderr.log
  events.jsonl                     # 本轮能捕获的事件，不伪造历史中间过程
  report.md
```

`events.jsonl` 可以先记录创建、提交、执行、检查和失败这些真实事件，不必本轮实现所有 agent 中间编辑的自动版本管理。每次明确提交的候选应有快照；历史只有最终代码时如实标记。

`runs/` 被 Git 忽略不代表研究成果可以任意消失：doctor/report 应提示未导出的重要候选；归档导出必须包含运行索引、代码、配置、评测记录、完整性校验和恢复说明。

---

## 8. 历史证据迁移与表述修订

### 8.1 证据身份

为每组历史实验登记：目录、目标 task/案例族、可识别的代码与 verifier 版本、score 口径、记录的模型别名/显示名、scaffold、隔离、有效/作废尝试、可用候选归档、缺失材料及 known issues。

未知信息填 `null` 并注明原因。不能把 NOTES 的模型名称反向推断为 API 精确模型 ID；不能把仓库 commit 的 AI co-author 标签等同于某个实验运行的模型身份。

对于旧规则计时错漏，保留已有的修正说明，标明哪个分数版本被后续重评替代。历史原文件只冻结“目前可取回的记录”；若旧版本已被覆盖而无法取回，不能假装现有归档包含所有历史分数。

### 8.2 只追加的重评

新逻辑应满足：

```text
历史 candidate snapshot（只读）
    + 指定 verifier/protocol
    → 新 eval ID 和新结果目录
    → 与历史结果的比较报告
```

禁止复用 `trial_XX_score.json` 覆盖源文件。若为了兼容旧脚本临时复制整组结果后运行旧重评，复制品必须位于新的工作目录；最终把输出登记为新 eval，并证明原目录哈希未变。

重构后的代码是新 implementation revision；不得直接继承历史速度数字。即使算法意图相同，也应通过数值回归后再记录新性能。

### 8.3 活跃文档建议措辞

| 旧说法 | 新活跃文档应表达的内容 |
|---|---|
| “published GSIS 不通过是格式必然性质，不是 bug” | “当前实现/测试下记录到固定点差异；保留数据与诊断脚本，因果解释及离散兼容性待进一步核验” |
| “passing 是否需要 domain knowledge” | “允许利用问题结构和已有方法；具体输入知识由搜索协议声明” |
| “score spread 是项目的核心信号” | “旧实验比较了该规则下的性能；新项目研究方法探索、效率、鲁棒性及结果复用” |
| “Krylov 不支持 periodic” | “当前周期状态支持与其限制需依据实现和测试说明；旧不支持描述已过时” |
| “所有 reference 都双路径认证” | “逐 case 列明残差检查、可完成的交叉核验和未完成项” |
| “156 tests green” | “基线文档记录 156 tests；当前测试状态见本轮实际运行报告” |
| “硬件无关的速度分数” | “旧 sweep-equivalent 规则包含计时；报告其硬件、线程、预热和缓存条件” |
| “非灰后端已经在项目中可用” | “负责人已有相关研究实现，仓库接入与接口尚待确认” |

不删掉不利结果来使新方向看起来顺利，也不保留未经证实的强结论来使项目看起来更有创新。

---

## 9. 本轮不做的工作，以及允许另立任务的发现

**不做：** 新 Kn 采样域设计、全量同精度 CIS 标定、新宏观闭合、改进 GSIS、证明固定点性质、重新实现 Krylov、自动提取全部候选为算子组件、真实多模型评测、非灰 solver/材料库、正式 leaderboard、论文结论与性能承诺。

**可以记录为后续独立任务：** 热流/守恒指标补齐；同精度 baseline；新 cost policy；大网格/角度规模验证；灰度候选结构假设审核；非灰接入；通用搜索与 agent 对照。这些任务需链接其依赖的 C01–C07，不因发现问题而自动扩张本轮。

发现旧数值代码疑似错误时，记录最小案例、代码位置、观察和影响。除非是明确且不涉及数值语义的路径/导入错误，否则不要在本轮顺手修复。数值修正要单独提交，旧代码/结果与新代码/结果保持可区分。

---

## 10. 分批执行计划

每批先完成可验证的小改动再进入下一批。建议独立提交，提交说明分清 documentation、structure、infrastructure 和 numerical change；本轮默认不产生最后一种。

### P0 — 盘点、冻结与建立基线

**操作：**

1. 记录当前 HEAD、工作区状态和相对本计划基线的差异；检查新增用户文件，不覆盖已有修改。
2. 建立 `inventory.csv`，至少覆盖所有待移动/删除文件、核心源码、生成夹具和实验资产。字段包括原路径、类型、是否 tracked、动作、目标路径、理由、关键依赖、恢复方法和验证状态。
3. 对历史 experiment 原始材料、旧 task 关键配置与数值源码生成完整 SHA-256 清单；记录脚本、源码、配置、mesh、reference 与候选之间的来源关系。不能只依赖当前 `manifest.json` 的短哈希。
4. 记录 Python、依赖版本、线程环境与是否启用 Numba。运行现有测试，区分 passed/failed/skipped/not-run；如有基线失败先记录，不把它归咎于重构。
5. 在新 scratch/output 路径中跑少量现有 solver 案例，保存数值输出和诊断作为前后对照。避免生成/覆盖 `tasks/t01-square` 或 `experiments/results`。

**产物：** 资产清单、基线清单、测试/数值基线记录、当前可复现性缺口。

**通过条件：** 任何将被清理的资产都有判定；关键历史材料可定位；没有未说明的工作区覆盖。未能运行的测试有具体原因。

### P1 — 重置项目文档与执行权威

**操作：**

1. 保存旧 README 和 Proposal 原文，归档旧 Proposal。
2. 重写 README；建立 `VISION.md`、`STATUS.md`、`CONTRACTS.md`、`PIPELINE.md` 和简短 `NON_GRAY_TRANSFER.md`。
3. 建立 `AGENTS.md`，记录本文第 0 节的执行边界，以及哪里才是活跃文档。
4. 修正 README/solver 文档的 GSIS 普遍化结论、periodic 过时说明、reference 认证表述和“主线必须扩题”的要求。
5. 更新当前仓库链接、身份与包 Source；历史记录和 canary 不做机械替换。

**通过条件：** 新读者能够立即理解研究目标、当前能力、第一/第二阶段边界及未定项；活跃文档不再把旧 Proposal 称为 plan of record。

### P2 — 整理资产与建立最小配置/方法登记

**操作：**

1. 按第 4 节归档/整理，保留核心后端和历史 fixture 的稳定路径。
2. 为历史五组结果建立索引，明确 score 版本、缺失材料、作废 trial 和旧规则下的限制。
3. 建立 draft study、显式 legacy smoke profile、新探索 prompt 模板，以及最小 registry。
4. 为内置方法和 F3 候选建立可定位条目；对 tarball 计算哈希，安全检查结构。复杂候选不急于重构。
5. 更新 `.gitignore`，明确新工作区与新运行输出位置；不添加会隐藏 reference、原始证据或研究 JSON 的全局规则。
6. 删除有充分证明的垃圾/重复内容；对暂不删除的副本在清单记录理由。

**通过条件：** 每个主目录角色清楚；历史材料校验不变；新 contract 字段仍是 pending；没有手工维护的新 solver 副本。

### P3 — 最小运行入口与非覆盖记录

**操作：**

1. 实现 doctor、catalog、legacy smoke 和新正式入口的契约准入检查。
2. 用薄适配调用现有数值后端与 checker，保留旧验证和计时规则的明确身份。
3. 实现候选 snapshot、run/eval ID、输入/配置/代码哈希、原始结果和 append-only 报告。
4. 为旧重评增加安全出口：新评估写新目录；原地行为仅限经检查的非受保护复制品。
5. 支持 mock provider 的候选提交、正常结束、超时/运行错误记录；不发起实际模型请求。
6. 将新 prompt 中的知识可见性、文件权限、预算、评价目标绑定到 C06 等契约；未定时只能检查/展示模板，不发送给真实 agent。

**通过条件：** 真实数值 smoke 能从问题配置走到独立检查和报告；重复执行不覆盖；正式运行在未定契约下被阻止；无需商业 CLI 即可完成工程验证。

### P4 — 回归与可信性测试

按第 11 节执行，不通过就报告和修复工程问题，不得通过改变数值公式或门限“修复”回归。

**通过条件：** 原核心行为没有未解释退化；新管线及失败路径有测试；所有重要资产哈希保持或有批准的变更记录。

### P5 — 交付与停止

更新 `STATUS.md` 和 `migration_report.md`。列清完成项、未完成契约、暂留文件、未运行测试与需要负责人决定的事项。

**本轮停止于：** 工程重构验收完成或已完整报告的环境阻塞。不得因 pipeline 跑通就自动开始新模型实验、选定 Kn 域或宣称新研究协议已建立。

---

## 11. 验收标准与测试矩阵

### 11.1 必须满足的验收项

| 类别 | 验收内容 | 判定原则 |
|---|---|---|
| 方向 | README/章程清楚表达灰度探索 → 非灰迁移；benchmark 是工具而非出题目标 | 人工阅读可确认；无旧路线的默认执行要求 |
| 契约 | C01–C07 显式存在；C07 不阻塞灰度工程测试 | pending 不被填默认值；不能输出新协议正式结论 |
| 资产 | 历史实验 tar/diff/score/trajectory/NOTES、reference 与关键夹具保留 | 重构前后完整哈希比对，差异有记录 |
| 核心 | 现有 solver 测试未出现未解释的新增失败 | 测试实际输出为证；不能只重复 README 的旧数量 |
| 数值 | 前后对照的温度/热流/分布布局、残差、收敛行为与计数一致或差异有解释 | 沿用既有回归准则；不得新设宽松阈值掩盖变化 |
| 流程 | 至少一个已有方法通过新入口完成真实 solve → evaluate → archive/report | 不接受全链条 mock 替代数值 smoke |
| 反例 | 错误形状、非有限值、伪造的自报收敛/残差不会被当成可信结果 | 原始独立检查与失败状态保留 |
| 记录 | 相同输出目录/运行 ID 不被静默覆盖；重评另存版本 | 冲突测试和源资产哈希测试 |
| 环境 | 本地 smoke 不要求付费服务；外部依赖缺失有清楚状态 | 不能把未运行 Docker/模型模式记为已通过 |
| 文档 | 新入口、旧复现入口、生成区、归档区和恢复步骤相互一致 | 链接/路径检查通过；已改路径均已登记 |

### 11.2 建议测试集合

**配置与准入测试。** 缺 C03/C04/C06、新协议 pending、无批准记录、无效字段、未知 profile、缺失输入和错误哈希均产生明确错误。dry-run 不执行候选、不调用模型；legacy smoke 不使 draft study 自动变为 approved。

**路径与来源测试。** 从文档指定的工作目录运行；排除错误导入本机其他版本 `pybte`；candidate 与 pristine 的实际路径和模块来源不能串用。归档移动后相对路径与原始 source revision 仍可解析。

**数值 smoke。** 优先选择旧 F1/F2/F3 中运行成本较低的已存在案例，例如各族的 `*_1_1.yaml`，以当前基线为准。保留 periodic/reflecting 状态路径的覆盖。这里的案例选择只是工程测试集，不是批准的新参数域。

**正/负对照。** 旧 oracle 通过和旧未加速 baseline 在指定刚性案例中失败的逻辑，应在可行的 legacy 集成测试中保留。不要要求 CIS 在容易案例上失败；也不要为快速 CI 自动运行整套 200,000 步标定。重型检查可单列并标明是否运行。

**档案与失败测试。** 解包安全、缺失归档、空 transcript、无最终模型事件、超时、非零返回码、无效输出、验证器异常、重复 trial/eval ID、追加重评分、凭据排除。区分算法未通过、实验预算耗尽、外部服务/执行设施故障，不把它们全压成一个 `False`。

**缓存与时间记录测试。** 旧预热/去缓存策略没有在适配中丢失；不能因 warm-up 失败仍输出看似有效分数。运行时间只作诊断，不要求在不同机器上与历史值完全一致。

**文档测试。** 活跃文档中的路径和命令可用；旧名称仅出现在合理的历史/兼容位置。`NON_GRAY_TRANSFER.md` 明确尚无本轮已验证接入，不得有“已实现非灰”的占位说明。

### 11.3 现有基线命令与执行注意

以下是仓库已有的测试入口，P0/P4 可据环境运行。依赖安装按仓库当前配置与隔离环境进行，并记录实际版本。

```bash
# 从仓库根目录进入数值后端；先检查已有测试状态
cd solver-python
python -m pytest tests/ -q
```

若资源限制只能先运行非慢测试，应另外注明 slow 项未执行，不能把部分通过写成全量通过。测试数量以实际收集结果为准。

旧 README 中 `make_env → calibrate → reference → oracle → validate → rollouts` 是一套会改写产物并可能花费大量计算/模型资源的命令，**不能作为本次一键 smoke 直接照抄执行**。[R02]

### 11.4 回滚策略

保持文档、文件移动和运行基础设施的提交可分离。失败时只回滚本轮对应提交/文件变更，不覆盖工作区其他修改。

数值核心未改时，legacy 夹具及源提交应仍能用于恢复。任何已清理文件必须可通过 inventory 指出的提交/归档找回。不要把“依赖当前机器某个 scratch 目录还没被删”当作恢复方案。

---

## 12. 执行完成后应交付的报告

agent 的最终交付至少包含以下内容，不能只说“重构完成”：

**变更摘要。** 新方向如何体现在入口和目录中；变更文件清单；哪些原路径已归档、删除或保留；为什么没有移动某些历史副本。

**能力边界。** 当前真正可以执行的命令与最小示例；哪些只支持 dry-run/mock；哪些仍因 C01–C07 或运行环境受阻。新 study 是否仍 pending 要明确写出。

**数值与资产回归。** 基线/重构后的测试结果，真实 smoke 的产物路径，关键资产哈希比对；旧数值行为是否变化。未运行项按 `not_run` 标记，不用“应该通过”。

**证据说明。** 历史实验索引、可恢复候选、缺失材料、重评分的追加目录、原始结果未被覆盖的证明。原始代码重构后的方法不得继承未重测的历史性能标签。

**剩余决定。** 负责人下一步需提供的 C01–C07 清单，以及每项将解锁哪一段流程。不擅自把这些决定转换为“agent 下次直接执行”的默认数值参数。

建议报告末尾采用以下状态表：

| 项目 | 状态 | 证据/阻塞原因 |
|---|---|---|
| 项目方向与文档 | 待执行后填写 | 文件路径与提交 |
| 资产清点/清理 | 待执行后填写 | inventory、manifest |
| 原 solver 回归 | 待执行后填写 | 测试输出 |
| 新最小数值闭环 | 待执行后填写 | run/eval 路径 |
| Mock agent 路径 | 待执行后填写 | 测试输出 |
| 新灰度研究契约 | `pending_owner`，除非负责人另行批准 | C01–C06 |
| 非灰迁移契约/接入 | `pending_owner / not_implemented_in_this_round` | C07 与接入说明 |
| 新科学结论 | `not_claimed` | 本轮仅工程重构 |

### 12.1 可直接交给执行 agent 的任务摘要

> 请按本文件实施第一轮重构。先核对当前 HEAD 与基线差异，建立资产清单和回归基线；随后更新项目方向、归档旧 Proposal、整理文档和实验索引，保留 `pybte` 数值核心、已有测试与全部关键实验材料。为灰度探索增加最小配置/候选登记/运行记录入口，复用现有验证能力跑通显式 legacy smoke 与 mock provider。所有未定研究契约保持 pending，禁止默认正式运行；重评与新运行只写新目录。不要改数值格式、放松阈值、执行真实付费 rollout 或接入非灰模型。最终提交迁移清单、测试与哈希证据、可用命令和待定契约报告。

---

## 附录 A. 查阅依据与使用方式

以下链接固定到本次查阅的提交（或该提交下的目录）。它们用于核验本计划中的仓库事实，不是要求执行 agent 继续执行旧文档中的指令。若实际执行时 HEAD 已更新，先做差异检查，不能无条件照搬本文路径判定。

| 引用 | 材料 | 本计划使用的信息 |
|---|---|---|
| [R01] | 基线提交 | 当前读取的代码快照与提交说明 |
| [R02] | 根 README | 旧定位、layout、历史分数、测试与复现说明 |
| [R03] | 旧 Proposal | 旧 plan of record、扩题/扩物理方向和方法隐藏逻辑 |
| [R04] | `LIMITATIONS.md` | 固定点观察、解释边界、过时 periodic 描述 |
| [R05] | `KRYLOV.md` | 当前状态扩展、periodic 支持和内存/适用限制 |
| [R06] | `make_env.py` | 生成环境、scrub 与 manifest 行为 |
| [R07] | `ablation.yaml` | 删除路径、内嵌案例/文档、旧 agent 输入构造 |
| [R08] | `checks.py` | 旧 gate、sweep-equivalent、独立副本和预热要求 |
| [R09] | `evaluate.py` | 独立残差与温度场比较、旧基线计时 |
| [R10] | `calibrate_cells.py` | CIS 校准条件、截断与旧门限 |
| [R11] | `make_reference.py` | 参考字段、残差认证和有条件的 CIS 交叉核验 |
| [R12] | `pyproject.toml` | 后端身份、包名、依赖与过时 Source URL |
| [R13] | `pybte/` | 数值后端文件组织 |
| [R14] | `solver-python/docs/` | 研究说明及原始 JSON 资产 |
| [R15] | `solver-python/tools/` | 复现/验证/诊断工具 |
| [R16] | `experiments/results/` | 五组历史结果目录 |
| [R17] | `rescore_rollouts.py` | 旧覆盖式重评与安全解包 |
| [R18] | F3 `NOTES.md` | 候选分类、旧结果、计时修订、作废 trial |
| [R19] | F3 `summary.json` | `dir` 隔离、模型别名、缺失元数据与最终有效 trial |
| [R20] | `run_rollouts.py` | provider 调用、沙箱、归档、异常及默认输出行为 |
| [R21] | T01 `Dockerfile` | 当前镜像声明与依赖安装范围 |
| [R22] | `pybte/driver.py` | `Case → Solver → RunRecord`、重启及数值设置接口 |

**关于本计划的性质：** 除上述仓库观察和负责人已经明确的研究方向外，目录设计、契约状态机制和执行分批属于本计划提出的工程安排；不是对已实现功能的描述。最终研究参数与评价结论仍由未定契约及后续实验决定。

[repo]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE
[R01]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/commit/dffff2319b290860b5c4a87509428595bce4d2b0
[R02]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/README.md
[R03]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/PROPOSAL_2_BENCHMARK_PROJECT.md
[R04]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/solver-python/docs/LIMITATIONS.md
[R05]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/solver-python/docs/KRYLOV.md
[R06]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tools/make_env.py
[R07]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tasks/t01-square/ablation.yaml
[R08]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tasks/t01-square/verifier/checks.py
[R09]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tasks/t01-square/verifier/evaluate.py
[R10]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tools/calibrate_cells.py
[R11]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tools/make_reference.py
[R12]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/solver-python/pyproject.toml
[R13]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/tree/dffff2319b290860b5c4a87509428595bce4d2b0/solver-python/pybte
[R14]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/tree/dffff2319b290860b5c4a87509428595bce4d2b0/solver-python/docs
[R15]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/tree/dffff2319b290860b5c4a87509428595bce4d2b0/solver-python/tools
[R16]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/tree/dffff2319b290860b5c4a87509428595bce4d2b0/experiments/results
[R17]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tools/rescore_rollouts.py
[R18]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/experiments/results/t01-square-F3/NOTES.md
[R19]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/experiments/results/t01-square-F3/summary.json
[R20]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tools/run_rollouts.py
[R21]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/tasks/t01-square/Dockerfile
[R22]: https://github.com/dingtao-shen/AI-Driven-Numerical-Method-Exploration-for-Multiscale-PBTE/blob/dffff2319b290860b5c4a87509428595bce4d2b0/solver-python/pybte/driver.py
