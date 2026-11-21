# 量子纠错的两端：Quafu 真机实测与面向容错的 AI 协同设计

第一届 “QUAFU 杯” 量子计算真机挑战赛 · 量子纠错与错误缓解赛道 · 队伍 **Qubit购买力**

本仓库是项目第一部分（QUAFU 真机实验）的代码、汇总数据、真机任务 ID 清单和技术报告。

## 一眼看结果

全部数字均为 QUAFU 真机实测（Quafu-SQC，2026-07-08 至 07-28，20 个运行日，Baihua 156q / Dongling 84q / Shenglian 84q，共 **64,098** 个真机任务）。统计口径为跨日均值与 bootstrap 95% 置信区间。

| 结果 | 数值 | 芯片 | 报告章节 |
| --- | --- | --- | --- |
| 当日自适应选址对冻结选址，d = 7 逻辑错误率改善 | 6.7× [4.9, 8.7] | Baihua | §3.2 |
| d = 7 编码存储 vs 同批最优物理比特 | 0.081 vs 0.149 | Baihua | §4.1 |
| [[4,2,2]] 逻辑 Bell 态 vs 物理 Bell 态（良好态） | 0.940 vs 0.894，接受率–保真度相关 r = 0.95 | Baihua | §4.3 |
| 无中途测量下的两轮 syndrome（SWAP 换出） | 时序特征命中 842 / 850 | Baihua | §7 |
| repetition 码码距扫描 d = 3…13 | 逻辑错误率在 d ≈ 9 触底 | Baihua | §5.1 |

全部实验与任务的对应关系见 [真机任务清单与主要结果.md](真机任务清单与主要结果.md)。

## 目录

```text
README.md                      本文件
运行说明.md                    依赖环境、运行步骤、如何核对数字
真机任务清单与主要结果.md       任务清单说明、各实验主要结果、任务数口径
docs/半决赛技术报告.pdf         技术报告（各结果的方法、对照与统计细节）
code/
  qfab/                        核心库：平台客户端、线路构造、标定解析、解码器、统计
  experiments/                 全部实验脚本（e00–e26）
  tools/                       统计引擎、绘图脚本、离线线路验证器
  requirements.txt             依赖清单
  qrun                         统一运行入口
data/
  task_ids.csv                 纳入统计的 64,098 个真机任务（实验、日期、提交时间、芯片、任务 ID）
  task_summary.csv             按“实验 × 芯片”汇总的任务数与起止时间
  task_ids_supplement.csv      未纳入统计的 1,850 个自测、探针等任务
  stats.md / stats.json        报告全部统计数字
  csv/                         报告主图的同名数据
```

## 快速开始

```bash
cd code
./qrun tools/verify_circuits_offline.py               # 离线验证线路，不需要 Quafu 账户
./qrun tools/freeze_stats.py --cutoff 2026-07-28      # 从原始记录重算全部统计（需原始数据）
```

环境配置和完整步骤见 [运行说明.md](运行说明.md)。每个任务的原始 counts、提交记录和当日标定快照约 1.3 GB，未放入本仓库。

## 第二部分：面向容错的技术储备

项目第二部分（神经解码器、纠错码搜索、综合征线路综合）基于谷歌公开数据与电路级模拟，不在 QUAFU 真机上运行，不在本仓库内。相关工作：

- Rethink the role of neural decoders in quantum error correction（ICML 2026，arXiv:2605.12046）
- Precision-Pyramid: towards real-time neural decoding for fault-tolerant quantum computing（NeurIPS 2026）
- Efficient foundation decoders for fault-tolerant quantum computing（arXiv:2606.27119）
- OmniQEC: discovering practical quantum error-correcting codes by an AI scientist（arXiv:2607.25865）
- Structure-guided synthesis of low-error syndrome extraction circuits（arXiv:2611.25825）
