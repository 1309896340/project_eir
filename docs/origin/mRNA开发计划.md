下面给出一份**16周开发计划**，目标是构建附件所述的 mRNA 设计系统：

> 输入目标蛋白序列 + 物种/细胞/修饰/递送等条件，输出优化后的完整 mRNA 候选（5′UTR + CDS + 3′UTR + polyA），并具备表达、稳定性、免疫风险、可制造性预测，带约束解码、多目标优化和主动学习闭环。

假设团队 4–6 人：算法 2 人、数据/生信 1–2 人、ML 工程 1 人、实验/安全联络 1 人。每周 5 个工作日，周末作为缓冲、文档和风险复盘。

---

# 阶段总览

| 阶段 | 周次 | 核心目标 |
|---|---|---|
| A 需求与数据 | W1–W5 | 明确目标、schema、数据源、特征和 tokenizer |
| B 预训练与生成 | W6–W8 | Transformer backbone、CDS 同义密码子生成、约束解码 |
| C 预测与结构 | W9–W11 | 多任务预测器、结构感知模块、UTR 与条件控制 |
| D 优化与闭环 | W12–W16 | 多目标优化、系统集成、主动学习、实验验证、交付 |

---

# 第1周：项目启动与需求定义

**Milestone**：完成 PRD、输入输出 schema、评价指标、安全策略和系统架构冻结。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 需求工作坊：明确设计目标、应用场景、成功指标 | 附件、业务需求、文献 | PRD v0、系统边界 |
| 2 | 定义输入输出 schema：蛋白、物种、细胞、修饰、递送、目标权重 | PRD v0 | JSON schema、示例样本 |
| 3 | 定义评价指标：表达、半衰期、TE、免疫、可制造性、约束违反率 | 附件第12节 | 指标字典、验收标准 |
| 4 | 安全与合规：敏感序列筛查、受控病原体/毒素过滤 | 生物安全规范 | 安全策略、审查流程 |
| 5 | 架构评审：生成模型 + 多任务预测 + 约束优化 + 主动学习 | 附件第3–11节 | 架构图、技术选型、风险清单 |

---

# 第2周：数据源盘点与数据规范

**Milestone**：确定公开数据源、数据库 schema、数据版本化方案。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 盘点公开数据：GENCODE、RefSeq、Ensembl、CoCoPUTs、Optimus 5-Prime、OpenVaccine、Ribonanza、RMDB、GEO/SRA | 附件数据章节 | 数据源清单、许可记录 |
| 2 | 定义数据库表：sequence、annotation、experiment、label、computed_feature | 附件第8节 | DB schema、ER 图 |
| 3 | 数据许可与伦理审查 | 各数据源条款 | 合规记录、可用数据列表 |
| 4 | 制定预处理标准：U/T、方向、CDS 检查、翻译一致性 | 附件第7节 | QC 流程文档 |
| 5 | 建立数据版本化：DVC/MLflow、目录结构 | QC 流程 | 数据仓库骨架 |

---

# 第3周：天然序列下载与清洗

**Milestone**：预训练序列集 v1，完成按基因/蛋白划分。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 下载 GENCODE human/mouse、RefSeq curated mRNA | 数据源清单 | 原始 FASTA/GTF |
| 2 | 提取 5′UTR、CDS、3′UTR、蛋白序列 | 原始注释 | 结构化序列表 |
| 3 | 清洗：去 N、长度过滤、起始/终止密码子、翻译比对 | 结构化表 | clean transcripts |
| 4 | 划分 train/valid/test：按基因/蛋白划分，防泄漏 | clean transcripts | split 文件 |
| 5 | EDA：GC、长度、密码子使用、UTR 分布 | split 文件 | EDA 报告、图表 |

---

# 第4周：功能标签数据整合

**Milestone**：监督数据 v1，包括 5′UTR、结构、稳定性、Ribo-seq 标签。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 整合 Optimus 5-Prime / 5′UTR reporter 数据 | 公开数据 | 5′UTR-TE 表 |
| 2 | 整合 OpenVaccine、Ribonanza、RMDB | 公开数据 | 结构/降解/reactivity 表 |
| 3 | 整合 GEO/SRA half-life、Ribo-seq 数据 | GEO/SRA | 稳定性/TE 表 |
| 4 | 标签标准化：z-score、rank、相对对照、批次校正 | 各标签表 | normalized labels |
| 5 | 合并 metadata：细胞、剂量、时间、修饰、递送、批次 | 所有表 | 多任务数据集 v1 |

---

# 第5周：特征计算与 Tokenizer

**Milestone**：特征库 + 多粒度 tokenizer + PyTorch DataLoader。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 组成特征：GC、CAI、tAI、CpG、UpA、k-mer、同聚物 | 序列数据 | 组成特征表 |
| 2 | 结构特征：MFE、局部 MFE、配对概率、可及性 | RNA folding 工具 | 结构特征表 |
| 3 | UTR 特征：uAUG/uORF、Kozak、miRNA、ARE、RBP | UTRdb、TargetScan、POSTAR | UTR 注释表 |
| 4 | Tokenizer：核苷酸、密码子、氨基酸、region、条件 token | 附件第1节 | tokenizer 文件 |
| 5 | 构建 Dataset/DataLoader，多任务标签对齐 | 特征+标签 | 训练数据加载器 |

---

# 第6周：Transformer Backbone 预训练

**Milestone**：预训练 backbone v0，完成 masked nucleotide/codon 任务。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 实现混合粒度 embedding：nt + codon + aa + region + position | 附件 Transformer 章节 | 模型代码 |
| 2 | 实现相对位置编码、区域 embedding、锚点相对位置 | tokenizer | 单元测试 |
| 3 | 自监督任务：masked nucleotide、masked codon、region prediction | 预训练集 | 训练脚本 |
| 4 | 启动预训练，监控 loss、学习率、显存 | 训练脚本 | checkpoint |
| 5 | 评估：masked accuracy、region F1、嵌入可视化 | checkpoint | 预训练报告 |

---

# 第7周：CDS 条件生成模型

**Milestone**：Protein → synonymous CDS 生成器 v0，保证氨基酸序列不变。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 实现 Protein Encoder + Codon Decoder | 附件第4节 | 生成模型代码 |
| 2 | 同义密码子 mask 输出层 | 密码子表 | 约束输出层、测试 |
| 3 | Teacher forcing 训练 | protein-CDS pairs | 生成模型 checkpoint |
| 4 | 采样生成候选，翻译验证 | checkpoint | 候选 CDS |
| 5 | 评估：氨基酸保持率、GC、CAI、密码子对偏好 | 候选 CDS | 生成评估报告 |

---

# 第8周：生成模型微调与约束解码

**Milestone**：生成器 v1 + ranking/reward + constrained beam search。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 构建 pairwise ranking 数据 | 同义变体、实验标签 | 偏好对 |
| 2 | 训练 reward/ranking 模型 | 偏好对 | reward 模型 |
| 3 | Reward-guided fine-tuning | reward 模型 | 微调生成器 |
| 4 | Constrained beam search：GC、motif、同聚物、内部 stop | 附件第8节 | 约束解码器 |
| 5 | 生成 top 候选并人工检查 | 约束解码器 | 候选库 v1 |

---

# 第9周：多任务预测器

**Milestone**：表达、稳定性、翻译效率、免疫、可制造性预测器 v0。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 实现共享 encoder + 多任务头 | 附件第6节 | 多任务模型 |
| 2 | 训练表达、稳定性头 | 监督数据 | 指标、checkpoint |
| 3 | 训练翻译效率、免疫风险头 | Ribo-seq、免疫数据 | 指标、checkpoint |
| 4 | 训练可制造性头：规则特征 + 小数据 | IVT/生产数据 | 模型 |
| 5 | 多任务评估、不确定性估计 | checkpoint | 预测报告 |

---

# 第10周：结构感知模块

**Milestone**：结构感知 attention / GNN 集成，提升预测和生成。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 实现 attention bias 加入配对概率 | 结构特征 | 代码、单元测试 |
| 2 | 实现 GNN 分支或结构门控融合 | 结构图 | 融合模型 |
| 3 | 结构辅助预训练/微调 | 结构数据 | 结构感知 checkpoint |
| 4 | 评估：paired/unpaired、reactivity、降解 | 结构标签 | 结构报告 |
| 5 | 集成到预测器，评估增益 | 预测器 | 集成报告 |

---

# 第11周：UTR 与条件控制

**Milestone**：5′UTR/3′UTR 生成/选择模块 + 条件控制 + 完整 mRNA 组装。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 5′UTR 设计：Kozak、uORF、局部结构 | 5′UTR 数据 | 候选 5′UTR |
| 2 | 3′UTR 设计：ARE、miRNA、RBP、polyA 信号 | 3′UTR 注释 | 候选 3′UTR |
| 3 | 条件 token：物种、细胞、修饰、递送 | 附件第8节 | 条件编码 |
| 4 | 训练条件生成/预测模型 | 条件数据 | 条件模型 |
| 5 | 完整 mRNA 组装：5′UTR+CDS+3′UTR+polyA | 各模块 | 完整候选 |

---

# 第12周：多目标优化与 Pareto

**Milestone**：多目标优化器 v1，输出 top-N Pareto 候选。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 定义 Score 函数与权重 | 附件第7节 | 评分模块 |
| 2 | 实现 Pareto 排序/NSGA-II | 评分模块 | 优化器 |
| 3 | 约束过滤：GC、motif、同聚物、酶切位点、剪接信号 | 附件第8节 | 过滤器 |
| 4 | 生成 top-N 候选并评分 | 优化器+过滤器 | 候选排名 |
| 5 | 敏感性分析：权重变化、场景切换 | 候选排名 | 优化报告 |

---

# 第13周：系统集成与 API

**Milestone**：端到端 mRNA 设计平台 v0，输入蛋白输出候选。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 封装推理 pipeline | 各模块 | CLI/API |
| 2 | 构建 Web/批处理接口 | API | 服务 |
| 3 | 集成安全筛查 | 安全策略 | 安全模块 |
| 4 | 集成实验记录与版本追踪 | MLflow/DVC | 跟踪系统 |
| 5 | 端到端测试：蛋白→候选 | 测试用例 | 测试报告 |

---

# 第14周：主动学习闭环与实验设计

**Milestone**：主动学习闭环 v0 + 实验批次设计。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 定义不确定性/多样性采样 | 预测器 | 采样策略 |
| 2 | 选择 100–1000 候选 | 候选库 | 合成订单 |
| 3 | 设计实验：表达、稳定性、免疫、IVT | 实验方案 | 实验 SOP |
| 4 | 建立数据回传格式 | LIMS 模板 | 回传模板 |
| 5 | 模拟闭环迭代 | 历史数据 | 闭环报告 |

---

# 第15周：实验验证与模型微调

**Milestone**：实验数据 v1 + 模型微调 v1。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 接收实验数据，QC | 实验数据 | clean labels |
| 2 | 微调预测器 | clean labels | 预测器 v1 |
| 3 | 微调生成器/reward | clean labels | 生成器 v2 |
| 4 | 重新生成候选并排序 | 生成器 v2 | 候选 v2 |
| 5 | 对比 v0/v1/v2 | 各版本 | 迭代报告 |

---

# 第16周：评估、安全审查与交付

**Milestone**：最终交付包：代码、权重、API、文档、合规报告。

| Day | 任务 | 输入 | 输出 |
|---|---|---|---|
| 1 | 离线评估：氨基酸保持、GC、MFE、motif、预测指标 | 最终模型 | 评估报告 |
| 2 | 安全合规审查：生物安全、序列筛查、数据许可 | 安全策略 | 合规报告 |
| 3 | 文档：模型卡、数据卡、使用说明 | 全部产出 | 文档 |
| 4 | 打包：代码、权重、Docker、API | 文档 | 发布包 |
| 5 | 评审与路线图 | 发布包 | 最终报告、后续计划 |

---

# 关键验收标准

1. **蛋白约束**：生成 CDS 翻译后与目标蛋白 100% 一致。  
2. **表达/稳定性**：在自建或公开测试集上，预测器显著优于基线。  
3. **约束满足**：GC、同聚物、motif、酶切位点等硬约束违反率接近 0。  
4. **多目标**：能输出 Pareto 候选，覆盖高表达、高稳定、低免疫等不同偏好。  
5. **闭环**：实验数据可回传并触发模型微调，形成迭代。  
6. **安全合规**：所有候选通过序列安全筛查。  

---

# 主要风险与调整

- **公开数据不足**：免疫、可制造性数据少，需尽早启动自有实验。  
- **标签批次效应**：必须保留 metadata，使用相对值/rank 归一化。  
- **长序列计算**：尽早采用 local/sparse attention 或 hierarchical Transformer。  
- **安全风险**：从 W1 就内置筛查，避免后期返工。  
- **实验周期**：W14–W15 可并行多批实验，若实验延迟，先用公开数据完成计算闭环。