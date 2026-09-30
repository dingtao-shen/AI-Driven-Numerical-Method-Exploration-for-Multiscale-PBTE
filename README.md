# AI-Driven Numerical Method Exploration for Multiscale PBTE

本项目让 AI agent 与研究者使用同一实验入口，探索、实现、选择和组合 PBTE
加速迭代方法。在满足目标问题的正确性要求后，评价相对 CIS 的求解成本。
方法可以是已有算法、合适的配置或组合，不预设原创性，也不预设领域方法一定更优。

第一阶段聚焦灰度模型；第二阶段研究向 mode-dependent 非灰模型的迁移。
当前后端是稳态灰度线性双弛豫 Callaway 模型，保留 DG/角度离散、CIS、GSIS、
Krylov 和实验校正路径。非灰实现尚未在本仓库接入或验证。
Kn 范围、精度、性能、搜索预算等 [C01–C07 契约](docs/CONTRACTS.md)仍为
`pending_owner`，正式研究运行禁止启动。

从仓库根目录执行（Python ≥ 3.10，依赖见 `solver-python/pyproject.toml`；
Numba 用于加速，pytest 用于测试；本轮已用现有 `dev` 环境验证）：

```bash
python -m exploration doctor
python -m exploration catalog
python -m exploration smoke --method krylov --out runs/my-first-smoke
python -m exploration mock --scenario success --out runs/my-first-mock
python -m exploration run --study studies/gray_exploration/study.yaml --dry-run
```

smoke 显式使用 `legacy_t01_smoke`：旧 T01 三类边界各一个已有 `*_1_1` 案例，
真实求解 → 独立验证 → 快照与报告。它验证工程兼容性，不代表新研究协议或新科学结论。
无需模型服务或凭据；输出目录必须不存在。正式 `run` 会在求解或调用模型前拒绝未批准契约。
本轮正式研究执行器本身也尚未实现，批准契约不会自动触发实验。

主流程与边界见 [PIPELINE](docs/PIPELINE.md)，已验证能力见 [STATUS](docs/STATUS.md)，
本轮交付见 [迁移报告](docs/refactor/migration_report.md)。

| 位置 | 作用 |
|---|---|
| [VISION](docs/VISION.md)、[CONTRACTS](docs/CONTRACTS.md) | 研究方向与待负责人决定的契约 |
| [solver-python](solver-python/README.md) | 原地保留的数值后端、测试、案例和研究工具 |
| [studies/gray_exploration](studies/gray_exploration/study.yaml) | draft study、显式兼容 profile、探索 prompt |
| [methods](methods/README.md) | 内置方法与五个 F3 历史候选的来源登记 |
| [exploration](exploration/__main__.py) | 轻量配置、快照、执行与独立评估入口 |
| [experiments](experiments/README.md)、[tasks](tasks/README.md) | 冻结历史证据与旧协议复现夹具 |
| [tools](tools/README.md) | 保留的构建、校准、参考认证、复现与重评工具 |
| `runs/`、`workspaces/` | 新生成区；不入 Git，有价值的结果需显式导出 |
| [旧文档归档](docs/archive/stikine_bench/README.md) | 历史阅读与恢复，已退出执行主线 |

```bash
python -m pytest tests/exploration -q
cd solver-python
python -m pytest tests/ -q
```

旧实验分数只适用于对应的旧案例、规则和计时条件，不能合并成新排行榜。
GSIS 当前实现的固定点差异保留为观察，因果解释和离散兼容性仍需核验。
基准评测服务于方法探索，不再以让强模型失败或扩大题库为成功条件。
