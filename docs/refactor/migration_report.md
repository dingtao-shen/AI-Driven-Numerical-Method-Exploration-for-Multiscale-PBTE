# 第一轮迁移报告

本轮已建立新的灰度方法探索工程入口，保留数值后端及历史证据。
执行分支为 `refactor/phase0-exploration`；初始 HEAD 与计划基线一致：
`dffff2319b290860b5c4a87509428595bce4d2b0`。未推送、未合并 main。
所有研究契约仍待负责人决定，正式研究执行器本轮未实现。

## 变更与恢复

- 首页改为方法探索方向，新增 VISION、CONTRACTS、STATUS、PIPELINE、NON_GRAY_TRANSFER、AGENTS。
  新 study 不继承 T01 的参数域、门限、预算和方法隐藏。
- 原 README 原样保存，Proposal 用 git mv 移入 docs/archive/stikine_bench。
  五份实质修订的 solver 文档保留原版；共七份原文哈希均一致。
  GSIS 数值观察和表格保留，普遍化因果结论改为待核验；periodic 支持按代码与测试修正。
  reference 认证逐 case 说明，三个深扩散 reference 的 CIS 交叉核验仍为 null。
- 包名 pybte、作者、许可证、版本及数值公式不变；包描述与 Source URL 更新。
  旧 make_env 示例改为实际存在的 T01，新增工具用途和成本说明。
- 五组历史结果原地冻结并建立 experiments/index.yaml；F3 五个 trial 各自登记。
  t01-square 丢失的候选树、C 的中断/失效、C2 空轨迹、F3 voided 原尝试与未知模型/usage 均保留。
- 新增 exploration 的 doctor/catalog/run 准入、CIS/Krylov smoke、mock、恢复与导出。
  数值 gate 和成本计算复用旧 checks.py，原 verifier、runner 和模型公式不修改。
  rescore 默认必须写新 eval 目录；受保护目录拒绝原地模式。显式复制品模式保留旧分数副本。
  旧 run_rollouts 保留 provider 机制，但要求新的显式输出目录，不再默认写 experiments/results。
- 新 runs/、workspaces/ 为忽略的生成区；有价值运行包另存本报告 evidence/。
  没有删除实验材料、未知文件、生成 fixture 或 numerical branch。旧独立副本各有角色，
  本轮未验证完整重建，故不移动、不合并、不清理。用户根计划文件保持原样，
  docs/plans/phase0_refactor_plan.md 是字节一致的入库副本。

逐项原路径、动作、依赖与恢复见 [inventory.csv](inventory.csv)；原 404 个 tracked 文件
加用户计划共 405 个文件的完整 SHA-256 在 [baseline_manifest.json](baseline_manifest.json)。
新增文件与工程修改列表见 [changes.json](evidence/changes.json)。
回滚应只撤销本轮文件，原数值资产仍在稳定路径；原文从 archive 或上述基线提交取回。
不要 reset --hard、git clean 或依赖临时目录恢复历史。

## 可执行能力

从仓库根目录执行（已用现有 dev 环境，版本见 environment.json）：

```bash
python -m exploration doctor
python -m exploration catalog
python -m exploration smoke --method krylov --out runs/review-smoke-001
python -m exploration mock --scenario success --out runs/review-mock-001
python -m exploration run --study studies/gray_exploration/study.yaml --dry-run
python -m exploration recover --method legacy_f3_trial_05 --out workspaces/review-f3-001
```

真实输出：runs/phase0-smoke-cis、runs/phase0-smoke-krylov、runs/phase0-mock-success。
各 run.json 索引三个 eval、实际候选 source、方法配置、case/reference/verifier 哈希、
进程日志、温度/热流/分布和原始成本。run/eval ID 冲突报错，无覆盖行为。
新 correctness_passed、aggregate_score 保持 null。

正式 run 的拒绝、重复 run ID、mock timeout/error 的实际 CLI 检查见
[cli_failure_checks.json](evidence/cli_failure_checks.json)。成功 mock 提交实际候选快照，
再执行同一真实数值链；失败场景记录原因，无伪造数值结果。
未知/未批准契约、缺少内容/批准记录、哈希不符、无效 profile 都在执行前拒绝。
即使未来批准全部契约，本轮正式执行器仍会报告 not_implemented_in_phase0。

## 数值与测试证据

原 solver 未改源码或测试。实际收集 161 项，非慢子集前后分别：
121 passed / 40 deselected（29.18 秒、23.26 秒）。全量首次运行停留在首个
nonthermalising CIS 慢测试，消耗约 2589 秒，受控 SIGINT 后报告 no tests ran。
慢测试状态为 not_completed / not_run，而非通过；没有放松测试或改变公式解决耗时。

CIS/Krylov × F1/F2/F3 的 *_1_1 六组前后对照，vdf、temp、qx、qy、residual_history、
mass_history 全部逐位一致，NPZ 也字节一致。CIS 计数 33/75/70，Krylov 21/23/27；
收敛状态和同一 pristine evaluator 重算的残差完全一致。
原直接读取 solver 的 reflecting-wall 残差可能带滞后壁面通量；独立 evaluator 会从 vdf
重建壁面通量。对基线和新结果都这样重算后残差一致，未修改数值格式。

额外使用旧 F1_1e-2 做单例正/负对照：冻结 oracle 228 步通过，残差约 2.77e-15；
CIS 达到原 6000 步 cap，未收敛，残差约 9.32e-6、温度误差约 1.59e-2，旧 gate 拒绝。
这不是新参数域、完整 15-case 验证或同精度 CIS 重新标定。

新增测试覆盖三边界真实 solve/check、伪造收敛和自报残差、错误形状、NaN/Inf、
预热失败、进程超时/返回码/无效 JSON、空 transcript/缺最终事件、快照和来源分离、
危险 tar 链接与越界、缺失档案、重复 ID、原历史哈希、追加重评和导出恢复。
追加重评测试使用 mock grader，验证两版输出及源哈希不变；没有把 mock 分数当成数值结果。
最终 44 项全部通过（75.37 秒），输出见
[exploration_acceptance.log](evidence/exploration_acceptance.log)。55 个活跃文档链接无断链，
Python 编译检查及 git diff --check 通过。

378 个受保护数值/测试/案例/网格、fixture、历史实验和研究 JSON 文件 SHA-256 全部保持，
[asset_audit.json](evidence/asset_audit.json) 记录允许的工程文档/工具修改及恢复来源。
原始分数和轨迹未覆盖。F3 五个归档全部通过结构检查，第五个恢复成功但未数值执行。

## 持久证据与可信边界

[evidence/smoke_cis.tar.gz](evidence/smoke_cis.tar.gz)、
[evidence/smoke_krylov.tar.gz](evidence/smoke_krylov.tar.gz)、
[evidence/mock_success.tar.gz](evidence/mock_success.tar.gz) 各自带相邻 manifest，
保存源码快照、配置、可信 fixture、运行/eval 记录与 NPZ。
已在三个新目录安全解包并验证全部 300/300/301 个文件哈希。
基线 NPZ 与这些包中的对应解字节一致，映射和 SHA 在 numerical_comparison.json。
因此本轮数值证据不依赖 runs/ 或 /tmp 的存续。

独立验证与 OS 隔离是两件事。当前只运行受控内置方法/mock，OS isolation 明确
not_validated。Docker 模式只检查工具与缺失状态，容器/商业 CLI/真实 provider not_run。
未复制 HOME、凭据或个人配置；候选和可信 fixture 使用基线文件白名单，缓存不入包。
真实 agent 的隔离、网络、依赖与异常策略等待 C06。

## 剩余决定与未执行项

| 契约 | 负责人需决定 | 解锁 |
|---|---|---|
| C01 | 模型、Kn 域、几何/边界、激励与归一化 | 正式问题集 |
| C02 | 固定离散、算法/网格/库/文件修改范围 | 正式搜索空间 |
| C03 | 残差/温度/热流/守恒、精度与参考认证 | 新正确性评价 |
| C04 | 成本口径、同精度 CIS、重复/截断/聚合/显著性 | 新加速结论 |
| C05 | 探索/确认划分、分辨率、种子、冻结/对照 | 确认实验 |
| C06 | 知识、模型/scaffold、预算、网络/依赖、隔离、异常 | 真实 agent 执行 |
| C07 | 非灰实现/材料/模态、适配与专家干预、成功标准 | 非灰迁移 |

未运行：40 个原慢测试的完整回归；全部 15-cell 旧验证及重复；大规模校准/研究扫描；
完整旧环境重建；F3 五候选数值重评；Docker 执行、付费模型或非灰接入。
前述耗时限制已记录，不解释成数值缺陷。未因流程可运行就启动新研究。

| 项目 | 状态 | 证据/阻塞 |
|---|---|---|
| 项目方向与文档 | complete | README、VISION、CONTRACTS、PIPELINE |
| 资产清点/整理 | complete | inventory、baseline_manifest、asset_audit |
| 原 solver 回归 | partial_validation | 121 非慢项前后通过；40 慢项未完成 |
| 新最小数值闭环 | complete | 两种内置方法 × 三族；解字节一致；单刚性正/负对照 |
| Mock agent 路径 | complete | 真实成功链与异常记录测试 |
| 新灰度研究契约 | pending_owner | C01–C06；正式研究执行器本轮未实现 |
| 非灰迁移契约/接入 | pending_owner / not_implemented_in_this_round | C07 |
| 新科学结论 | not_claimed | 本轮为工程重构 |
