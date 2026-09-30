# 执行入口与记录

在仓库根目录使用配置 → 已有方法/候选提交 → 代码与配置快照 → 求解 →
独立验证 → 原始观测、带版本的判定与报告。外层不改写 `Case → Solver → RunRecord`。

```bash
python -m exploration doctor
python -m exploration catalog
python -m exploration run --study studies/gray_exploration/study.yaml --dry-run
python -m exploration smoke --method krylov --out runs/smoke-001
python -m exploration mock --scenario success --out runs/mock-001
python -m exploration recover --method legacy_f3_trial_05 --out workspaces/recovery-001
python -m exploration export --run runs/smoke-001 --out runs/export/smoke-001.tar.gz
```

默认 smoke profile 是显式命名的 legacy_t01_smoke，包含 F1/F2/F3 的 `*_1_1`。
可选受控内置方法 cis/krylov，默认 cis；当前 GSIS 和实验路径只登记，仍可从原后端运行。
profile 必须通过版本、输入与完整哈希检查，不能用于填充 C01–C06。
这些容易案例是工程覆盖，CIS 通过是预期行为，不是加速收益证据。
正式 run 无 dry-run 时会先拒绝 pending 契约；正式执行器本轮未实现。

运行产物包含 run.json、resolved_config.yaml、input_manifest.json、events.jsonl、
report.md、candidates/submission_001/ 的实际代码和配置，以及 evaluations/<eval_id>/。
每个评估有 solution 对应的 `<case>.npz`、raw_legacy.json、observations.json、
evaluation.json、进程日志、返回码、超时状态、实际 pybte 导入位置。
快照 ID 由代码和方法配置哈希决定；run/eval 目录独占创建，不静默覆盖。
provider 异常会保存原始生命周期，不生成假数值结果。

smoke 通过薄适配调用冻结 checks.py 的 score_cell；候选和 evaluate.py 分进程运行。
独立验证器重算输运残差与温度场，不采信候选自报的残差或收敛标记。
vdf 的轴为 direction、element、nodal DOF；温度及 qx/qy 是单元积分，未除面积。
旧门限检查温度，不意味着旧 gate 检查了热流；热流只作为原始字段保存。
新协议 correctness_passed 和 aggregate_score 始终为 null。

每个 case 使用去缓存的候选副本，旧 runner 保留同族粗网格预热、Case+Solver+run
计时和 fine sweep 单位；候选与 pristine JIT 缓存分开。预热失败不计分。
机器、依赖、线程固定值、JIT/缓存条件和 RSS 随结果记录。
reference 的 CIS 交叉核验逐 case 见冻结 reference/index.json，null 表示未完成。

当前执行仅允许受控内置方法与 mock。原仓库夹具不由候选进程导入；可信副本与候选快照
分别存储，执行前后校验内容。没有可验证的 OS 安全沙箱，状态为 not_validated；
不能据此运行不可信代码。C06 决定真实 agent 的网络、凭据、文件权限和隔离策略。
doctor --mode docker/provider 只检查相关可见工具；容器执行均标记 not_run，
Dockerfile 尚未展示容器内 claude CLI 安装。默认本地 smoke 不需要这些工具。

mock 场景：success、timeout、error、empty_transcript、missing_final、invalid_output。
成功场景提交内置方法快照并走同一真实数值链；其余场景记录异常并停止求解。
dry-run 只显示计划和暴露清单，不执行候选、不调用模型、不发送 prompt。

追加重评示例（会执行归档候选，可能较重，本轮只测试机制，未发起候选全量重评）：

```bash
python tools/rescore_rollouts.py experiments/results/t01-square-F3 \
  --task tasks/t01-square --sandbox-root /tmp/pbte-rescore \
  --trials 5 --out runs/f3-eval-001
```

每次重评使用新的目录和 eval ID，保存源档案/验证器/参考哈希、旧分数比较与新结果。
新入口拒绝受保护历史目录上的原地重评。`--in-place-copy` 只接受明确准备的
workspaces/ 或仓库外普通复制品，拒绝符号/硬链接，并保存被替换文件的副本。
历史 provider 仍在 tools/run_rollouts.py，必须显式指定新的输出目录；本轮不运行它。

export 只打包已知运行产物，去除缓存/执行工作树，附完整 SHA-256 manifest 和恢复说明。
`runs/`、`workspaces/` 被忽略；有价值的导出需明确保存到证据区并纳入版本控制，
不能只依赖临时目录。原始历史材料的恢复来源见 inventory.csv。
导出不是隐私审查；历史轨迹未来公开前仍须单独审核。
