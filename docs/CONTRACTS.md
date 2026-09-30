# 研究契约与准入

负责人尚未批准 C01–C07。机器状态在
[study.yaml](../studies/gray_exploration/study.yaml)，所有 value 和 approval 为 null。
工程 smoke 引用独立旧协议，不会填充或批准新契约。

| ID | 待提供内容 | 解锁流程 |
|---|---|---|
| C01 | 模型、Kn 二维范围/组合、几何、边界、激励、归一化 | 正式问题准备 |
| C02 | 目标离散、允许修改的算法/网格/库、自适应、文件权限 | 正式候选搜索空间 |
| C03 | 残差尺度、温度/热流/守恒、门限、零尺度、参考认证 | 新正确性判定 |
| C04 | 主成本、setup/JIT/缓存/内存、同精度 CIS、截断、重复、聚合、显著性 | 新性能判定 |
| C05 | 探索/确认划分、分辨率扩展、种子、冻结与非 agent 对照 | 正式确认实验 |
| C06 | 可见知识、可编辑范围、模型/scaffold、预算、网络/依赖、隔离、异常 trial | 真实 agent 执行 |
| C07 | 非灰实现/材料、模态接口、允许适配、专家干预和成功标准 | 非灰迁移 |

灰度正式入口要求 C01–C06 全部批准，study status 为 approved，且
formal_execution_enabled 为 true。C07 不参与灰度准入。
只翻转布尔值、空字典、0 或 TODO 不能代替完整契约。

结构字段定义在 [contracts.py](../exploration/contracts.py) 的 `REQUIRED`；这是工程
完整性清单，未预设任何科学数值。例如 C03 value 必须有 residual、field_metrics、
thresholds、zero_scale、reference_certification。契约 value 必须是包含所有必填项的
非空 mapping，不接受空值/占位字符串。具体字段内容由负责人决定。
每个 approval 需记录 approved_by、带时区的 approved_at、decision_ref 和 value_sha256。
哈希为 `digest_value(value)` 的规范 JSON SHA-256；修改内容须重新批准。
结构检查不会证明批准者身份或科学有效性，决定记录仍须人工审核。

本轮只实现准入检查，正式研究执行器等待研究契约后另立任务。
即使结构检查通过，也会报告 not_implemented_in_phase0。
确认集和 reference 不因写入配置就自动成为 agent 可见输入。
