# Historical suites

`t01-square/` 是冻结的 legacy fixture，包括旧 prompt、ablation、environment、
pristine verifier、reference 和 oracle。路径、历史数字、canary 与独立副本均保留。
它们承担不同实验角色，不是可随意删除或合并的重复文件。

新工程入口只通过显式 `legacy_t01_smoke` 使用已有案例；T01 不再是新研究默认协议。
不要在当前夹具上运行 make_env/calibrate/reference/oracle 覆盖生成。
完整历史重建应在基线提交的独立 checkout/worktree 中执行，并对照哈希和旧 ablation。
执行步骤和成本见 [tools](../tools/README.md)。
