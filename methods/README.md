# 方法登记

[registry.yaml](registry.yaml) 登记 CIS、GSIS、Krylov、实验 defect correction /
宏观预条件路径及 F3 五个独立 trial。内置配置是已有实现示例，不是研究默认值。
来源、适用假设、原始归档 SHA-256 与验证状态独立记录。

F3 方法家族来自 NOTES，状态是 reported legacy pass；未在本轮重构算法或数值重评。
相似方法家族不会合并 trial，也不继承到未来实现的性能标签。

```bash
python -m exploration catalog
python -m exploration recover --method legacy_f3_trial_05 --out workspaces/f3-recovery
```

recover 检查完整哈希、拒绝链接/越界/特殊文件，保留 `filter="data"`，定位 driver；
不导入或执行归档代码。新 smoke 只接受受控的 CIS/Krylov 内置方法。
未来维护版本须另建 implementation revision 和证据，再考虑提取可复用组件。
代码缺失的旧实验只进入 experiments 索引，不伪造可调用入口。
