# 当前能力与验证状态

本轮为工程重构，研究方向是灰度 PBTE 加速探索 → 后续非灰迁移。
新 study 仍是 draft，C01–C07 均 pending_owner；不声称新科学结论。

| 项目 | 已观察的状态 | 证据 |
|---|---|---|
| 数值后端与夹具 | 原数值源码、测试、案例/网格、历史结果与研究 JSON 未改 | [378 文件哈希审计](refactor/evidence/asset_audit.json) |
| 原 solver 回归 | 非慢项前后各 121 passed，40 deselected；全量未完成 | [基线](refactor/evidence/solver_baseline_non_slow.log)、[重构后](refactor/evidence/solver_after_non_slow.log) |
| 全量测试 | 首个慢项持续 2589 秒仍未完成，受控中断；不能声称全量通过 | [中断日志](refactor/evidence/solver_baseline.log) |
| 真实数值闭环 | CIS/Krylov × 三类边界均通过旧 gate，六组 NPZ 字节一致，计数和独立残差一致 | [数值对照](refactor/evidence/numerical_comparison.json) |
| 旧刚性正/负对照 | F1_1e-2：oracle 228 步通过，CIS 6000 步未收敛并被拒绝 | [对照](refactor/evidence/legacy_controls.json) |
| 新管线 | 44 passed：配置、快照、独立检查、异常、mock、追加重评与归档恢复 | [验收日志](refactor/evidence/exploration_acceptance.log) |
| 历史候选 | 五个 F3 归档结构与哈希已检查，第五个已安全恢复；数值重评 not_run | [归档检查](refactor/evidence/archive_checks.json)、[恢复](refactor/evidence/recovery.json) |
| 证据导出 | 三个完整运行包已导出并恢复验证，不依赖 scratch | [恢复记录](refactor/evidence/export_recovery.json) |
| 正式研究/非灰 | 正式执行器本轮未实现；C01–C07 pending，非灰未接入 | [契约](CONTRACTS.md)、[迁移备忘](NON_GRAY_TRANSFER.md) |
| Docker/真实 provider | not_run；OS 隔离 not_validated；未发起模型调用 | [doctor](refactor/evidence/doctor_docker.json) |

入口与输出在 [PIPELINE](PIPELINE.md)，变更、原文恢复及未执行项在
[迁移报告](refactor/migration_report.md)。历史研究扫描、完整旧夹具重生成、全 15-case
重复验证和 F3 候选性能重评均未执行；不能把这些状态推断为通过。
