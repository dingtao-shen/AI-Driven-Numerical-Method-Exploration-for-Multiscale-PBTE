# 历史实验与研究证据

[index.yaml](index.yaml) 索引五组历史目录：precursor、t01-square、C、C2、F3。
旧原始文件原地冻结，完整文件哈希在 [baseline_manifest](../docs/refactor/baseline_manifest.json)。
不同 score 版本不能合并排名，也不能视为新协议下同精度 time-to-solution 结论。

F3 为灰度 pilot 与候选来源，NOTES 报告五个有效候选在旧 v0.3.2 的 15 案例通过。
它不证明新参数域鲁棒性或非灰迁移。模型别名 `opus`、NOTES 显示名、缺失 usage
是不同证据；精确模型 ID、缺失费用和轨迹保持 null。
F3 voided/ 的两个原始 quota 失败 trial 保留；当前 summary 只记录替换后的有效 trial。
t01-square 的候选树已丢失；C 中断/失效与 C2 的空 timeout 轨迹均在索引说明。

solver-python/docs/*.json 是求解器诊断资产，原地保留，不与 agent 结果混排。
reference 认证逐 case 读取 index.json；cis_iterations 为 null 的深扩散案例未完成
CIS 交叉核验，不能宣称所有 reference 均双路径认证。

新重评用 `tools/rescore_rollouts.py ... --out runs/<new-eval-id>`。
源目录保持只读语义；新结果附输入/输出哈希和原分数比较，原地模式仅用于明确的非受保护复制品。
详细命令见 [PIPELINE](../docs/PIPELINE.md)。
