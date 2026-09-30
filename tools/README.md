# 工具用途与运行成本

新工程入口是仓库根的 `python -m exploration`。下列旧工具保留其复现机制，
未隐式升级为新研究协议。完整旧链条会覆盖产物或消耗大量资源，不能当作 smoke。

| 工具 | 输入 → 输出 | 成本和本轮用途 |
|---|---|---|
| make_env.py | solver-python + ablation → environment/pristine/manifest | 确定性旧环境构建；会写夹具，只在基线独立 checkout 重建 |
| calibrate_cells.py | 旧案例、CIS → cells.json | 最多 200,000 步校准，重型；本轮不运行 |
| make_reference.py | 案例、Krylov/CIS → reference NPZ/index | 每 case 认证，CIS 可达时交叉核验；本轮读取已有 reference |
| make_oracle.py | 旧环境、Krylov → oracle solution/patch | 旧 oracle 构建，不能覆盖冻结版本 |
| validate_task.py | checker、oracle、baseline → 验收统计 | 全族集成耗时分钟级或更长；本轮另记是否运行 |
| run_rollouts.py | 旧 prompt、provider → 候选/轨迹/分数 | 真实 CLI/可能付费，仅 legacy，需显式新目录；phase0 不运行 |
| rescore_rollouts.py | 归档候选、明确 verifier → 新 eval 目录 | 不调用模型，但执行候选，可能数分钟；源证据不覆盖 |

`solver-python/tools/` 的工具从 solver-python 目录按各自 --help/说明运行：

| 工具 | 输入/输出与定位 |
|---|---|
| make_meshes.py | 网格族生成；已有 meshes 保留，不是默认 CI |
| reproduce_published_table.py | 文献测试点与 iteration 表的复现，含长 CIS；本轮不重新做研究扫描 |
| validation_sweep.py | 参数/方法变化 → VALIDATION.md 与 docs/validation_sweep.json；重型诊断 |
| fixed_point_study.py | 离散/参数变化 → --out（默认 fixed_point_study.json）；历史 docs/fixed_point_study.json 原地保留，因果解释待核验 |
| defect_study.py | omega/every/Anderson 等实验设置 → 诊断输出；不改变默认格式 |
| benchmark.py | 求解器性能测量 → 可选 --json；历史 docs/benchmark.json 与 agent 分数分开 |
| calibrate_budgets.py | 旧预算校准 → docs/budgets.json；不决定新 C04/C06 |

恢复旧链条时使用基线 dffff2319b290860b5c4a87509428595bce4d2b0 的独立 checkout，
依次核对其 make_env → calibrate → reference → oracle → validate 的输出哈希。
新入口 smoke 直接读取冻结夹具，不需要重建，且输出只落新目录。
付费 rollout、全量扫描不属于这次重构授权。
