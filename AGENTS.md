# AGENTS.md

本文件为 AI 编码代理提供本项目的上下文。修改项目内容或新增文档/代码时，请同步更新本文件。

## 项目概述

RNADesign 是一个 **mRNA 序列设计深度学习模型**项目。核心目标不是生成任意 RNA，而是在给定目标蛋白质（或功能目标）的前提下，设计出**稳定、可高效翻译、低免疫原性、易生产且满足安全约束**的 mRNA 序列（5'UTR + CDS + 3'UTR + poly(A)）。

总体技术路线已在设计文档中确定：

```text
条件 Transformer 生成同义 CDS
+ RNA 结构感知模块
+ 多任务表达/稳定性/免疫预测器
+ 约束解码与多目标优化
+ 主动学习实验闭环
```

## 当前状态

- 项目处于**方案设计阶段**，尚无模型实现，未选定训练技术栈/框架。
- `mRNA模型设计.md`：约 3700 行中文设计文档，是当前所有结论的来源。
- `docs/生物学入门.md`：面向零基础生物背景工程师的入门文档（中心法则、遗传密码与同义密码子、mRNA 结构、二级结构、四大优化性质、生物↔ML 概念对照、术语速查表）。
- `data/`：CoCoPUTs（HIVE-CUTs）人类 CDS 组成统计四份 TSV（`o537-Human_CDS.tsv` 密码子计数 77 列、`o537-Human_CDS_Dinuc.tsv` 二核苷酸计数 29 列、`o537-Human_CDS_Junc_Dinuc.tsv` 密码子交界二核苷酸计数 29 列、`o537-Human_CDS_Bicod.tsv` 密码子对计数 4,110 列、约 944 MB，读取需按需选列）。四表同记录、同行序：119,196 条 GRCh38.p13 RefSeq CDS（20,239 个基因；含 13 条线粒体、58,360 条 XP_ 预测、6,358 行 primary/alternate scaffold 重复，使用前需过滤去重；DNA 字母表注意 T↔U，Bicod 表列名为小写 6-mer）。各表字段含义与已验证的计数口径见对应的 `data/*.字段说明.md`。
- `scripts/deduplicate_human_cds.py`：使用标准库处理 `o537-Human_CDS.tsv`。仅保留 `Organelle == genomic`、`Accession` 为 `NC_*` 的行；对每个 `(Gene ID, Protein ID)` 按源文件顺序保留首条。默认生成 `data/o537-Human_CDS.primary_genomic_deduplicated.tsv`（112,454 行数据）；测试位于 `tests/test_deduplicate_human_cds.py`。
- `scripts/download_gencode.cmd`：Windows 批处理脚本，断点续传下载 GENCODE human release_50 数据到 `data/gencode/human/`。默认下载核心三件（pc_transcripts / pc_translations / annotation.gtf + MD5SUMS）；参数 `all` 额外下载基因组 FASTA、RefSeq ID 映射、lncRNA 负样本、polyAs 位点。优先用 PATH 中的 `wget -c`，找不到则降级用 Windows 自带 `curl -C -`。**注意：该文件必须保持纯 ASCII（英文注释）**——UTF-8 中文在中文 Windows 默认 GBK 代码页下会破坏批处理解析。数据源与文件选择依据见 `reference/datasources/` 下 GENCODE 相关文档（注意 `Eir D06 数据源.md` 文件名与内容不符，内容实为 GENCODE 流程）。
- `reference/models/开源模型盘点.md`：开放权重、可直接推理部署的开源模型盘点（2026-09）。密码子优化首选 CodonTransformer（HF `adibvafa/CodonTransformer`，pip 即用）与 LinearDesign（非学习型 baseline）；5'UTR 翻译效率 UTR-LM / Optimus 5-Prime / Framepool；结构感知 RibonanzaNet / RNA-FM；骨干 ESM-2 / Evo-2；分析工具 cubar（R 包，原生读取 CoCoPUTs 表并内置 CodonTransformer 接口）。结论：免疫原性与 IVT 可制造性方向**无**现成开放权重模型，维持规则近似 + 自建数据路线。
- `pyproject.toml` + `.python-version` + `uv.lock`：**Python 环境已用 uv 管理**（Python 3.12——codontransformer 依赖 `numpy<2.0`，而 numpy 1.26.x 的 wheel 只发布到 cp312，故上限 `<3.13`）。**torch 为 CPU 版、设备无关**（默认源 PyPI 轮子，任何机器 `uv sync` 即装）；如需恢复 GPU 版，按 pyproject 内注释配置 `[[tool.uv.index]]` explicit 源——直连 `download.pytorch.org` 实测仅 ~200 B/s，须用上海交大镜像 `mirror.sjtu.edu.cn/pytorch-wheels/cu126`（扁平目录：列表页在 `/cu126/torch/`，文件在 `/cu126/` 下），且 **torch 必须显式写在 `dependencies` 里 `[tool.uv.sources]` 才生效**（uv 缓存中已留有 cu126 轮子）。默认包源：阿里云 PyPI 镜像。本机 GPU：RTX 3080 Laptop 16GB（Ampere, sm_86），驱动 596.21 / CUDA 13.2。已安装：codontransformer 1.6.7（导入名 `CodonTransformer`，注意大小写）、torch 2.12.1+cpu、transformers 4.57.6。
- `.gitignore`：跟踪源码、测试、Markdown 数据说明、`MD5SUMS` 与共享 VS Code 配置；忽略本地环境、缓存、原始/生成数据、模型权重、checkpoint 和实验输出。`CodonTransformer/` 的应用代码与文档纳入 Git，仅 `CodonTransformer/data/` 本地权重及 tokenizer 文件被忽略；未来的 `/models/` 源码目录不再被整体排除。
- `CodonTransformer/data/`：CodonTransformer 预训练权重本地副本（model.safetensors 358MB + config/tokenizer 共 6 个文件，已通过 safetensors 解析校验）。注意 `huggingface_hub` 新版在 hf-mirror.com 上会因 HEAD 请求被 308 重定向回 huggingface.co 而失败（GET 正常），故权重走手动下载，代码从本地路径加载。
- `scripts/codontransformer_demo.py`：CodonTransformer 推理示例（蛋白→宿主为 Homo sapiens 的优化 CDS，`match_protein=True` 即同义约束；含翻译一致性校验、GC 含量与同义变体采样）。运行：`uv run python scripts/codontransformer_demo.py`（自动优先加载本地 `CodonTransformer/data/` 权重）。已验证 CPU 推理通过。
- `CodonTransformer/codontransformer_app.py`：gradio 调试界面（端口 7860），支持 164 物种下拉、确定性/采样解码、temperature/top_p/num_sequences 调节、翻译一致性校验；内置三份测试样例（INS 胰岛素原 110aa / HBB 血红蛋白β 147aa / EPO 促红细胞生成素 193aa，序列取自 `data/gencode/human/gencode.v50.pc_translations.fa.gz` 真实人类蛋白，EPO 样例演示大肠杆菌宿主对比）。运行：`uv run python CodonTransformer/codontransformer_app.py`，或在 VS Code 中 F5（`.vscode/launch.json` 已配置 debugpy 启动项，解释器指向 `.venv`，服务就绪后自动打开浏览器）。使用手册：`CodonTransformer/docs/演示界面使用手册.md`（操作、样例与 FAQ）与 `CodonTransformer/docs/密码子优化核心概念说明.md`（宿主条件、解码模式、同义约束与采样参数）。后者明确当前界面主结果始终为确定性输出，采样仅生成变体表。应用代码和文档受 Git 跟踪，模型权重目录 `CodonTransformer/data/` 保持本地忽略。app 与 demo 两个入口均注册了 `logging.Filter`，按消息过滤 transformers≥4.50 对 `BigBirdForMaskedLM` 的 "has generative capabilities" 误报警告（该提示针对 `.generate()` 能力，而本项目解码走 `predict_dna_sequence` 的单次 MLM forward，不依赖 `.generate()`）。
- 项目仓库根目录为 `RNADesign/`，沿用原 `git/project_eir/.git` 迁移后的 Git 仓库（`main` 分支）；更早的 RNADesign 根目录 `.git` 已移除。
- 原 `git/project_eir/` 的设计与数据源文档已归档至 `docs/origin/`，用于保留迁移前资料。
- 后续可能的演进方向：数据管道搭建（公开数据下载与清洗）→ 模型实现 → 训练评估。

## 设计文档结构（mRNA模型设计.md）

文档由四大部分组成（各部分内部有独立的章节编号）：

1. **模型总体设计**（§1–13）：设计目标、输入输出、三大核心模块（条件生成模型、多任务预测模型、约束优化模块）、目标函数、约束解码、训练策略（预训练→监督微调→生成模型训练→主动学习闭环）、评价指标、关键设计原则。
2. **Transformer 架构调整**（#1–13）：面向 mRNA 场景的定制改造，包括 tokenization、位置编码、attention 机制、编解码结构、输出层、结构信息融合、长序列建模、条件控制、解码策略、预训练任务、block 内部改造（adapter/MoE）。
3. **训练数据体系设计**（#1–15）：五类数据（天然序列、同义密码子设计数据、功能标签数据、结构与生物物理特征、实验闭环数据）、样本 JSON 格式、数据库表结构、数据构建流程、划分策略、QC、MVP 数据集。
4. **公开数据源盘点**（#1–14）：按用途分类的可用公开数据源、推荐组合、常见坑。结论：**没有覆盖完整 mRNA 设计目标的单一公开数据集**，需要多源整合，免疫与可制造性标签基本需自建。

## 核心设计决策（已定，避免重新讨论）

- **蛋白约束优先**：CDS 生成必须保证翻译产物与目标蛋白完全一致，通过同义密码子约束输出层实现（`logits[非同义密码子] = -∞`）。
- **模型架构**：Protein Transformer Encoder + Codon Transformer Decoder（cross-attention），支持自回归与非自回归两种变体，以及 masked refinement 迭代优化。
- **tokenization**：混合粒度——UTR 用单核苷酸 token，CDS 用密码子 token；CDS 位置叠加 codon/amino-acid/position embedding。
- **结构感知**：将 RNA base-pairing probability 作为 attention bias，或采用 Transformer + GNN 双分支融合。
- **多目标评分**：`Score = w1×Translation + w2×Stability − w3×Immunogenicity − w4×StructurePenalty − w5×ManufacturingPenalty − w6×ConstraintViolation`，权重按场景（疫苗/蛋白替代疗法/细胞治疗）调整。
- **硬约束**（约束解码阶段强制满足）：蛋白序列一致、GC 范围、避免长同聚物/指定 motif/限制性位点/内部终止密码子/强 5' 端发夹/隐性剪接与 polyA 信号、安全筛查。
- **训练流程**：自监督预训练（masked nucleotide/codon）→ 监督微调（实验标签）→ 生成模型训练（含 reward-guided）→ 主动学习实验闭环。
- **场景条件化**：宿主/细胞类型/修饰核苷酸/设计目标作为条件 token（如 `[HOST=human] [MOD=m1Ψ]`），配合 adapter/LoRA 或 MoE 做场景适配。

## 关键领域规则（实现时必须遵守）

- 序列统一用 RNA 表示（A/U/G/C），保存时可同时保留 DNA 版本（T ↔ U）；方向统一 5'→3'。
- CDS 长度必须是 3 的倍数，不得含内部终止密码子，翻译结果须与蛋白序列比对一致。
- 不同实验来源的标签尺度不可直接混合，需 within-assay 归一化（z-score / rank / 相对阳性对照），并同时保留 raw_label、normalized_label、batch_id 等 metadata。
- 数据划分至少做按基因/蛋白划分（同蛋白的同义变体不得跨 train/test），避免泄漏；另评估按批次/细胞类型/物种划分。
- 必须保留负样本（低表达/不稳定/高免疫/生产失败变体），否则模型只会生成"像天然"而非"表现更好"的序列。
- GTEx/HPA 等天然表达量不可直接当 mRNA 设计标签（天然表达主要受转录调控影响）。
- OpenVaccine/Ribonanza 是短 RNA 片段的结构/降解数据，只作辅助任务，不代表完整 mRNA 表达。

## 公开数据源速查

| 用途 | 首选数据源 |
|---|---|
| mRNA/CDS/UTR 预训练 | GENCODE (human+mouse)、RefSeq、MANE Select、Ensembl、RNAcentral |
| 密码子/密码子对特征 | CoCoPUTs、HIVE-CUTs、Kazusa（建议自行从高质量 CDS 统计） |
| 5'UTR 翻译效率 | Optimus 5-Prime / Sample et al. reporter 数据集 |
| 结构/reactivity/降解 | OpenVaccine (Kaggle)、Ribonanza、RMDB、bpRNA、ArchiveII |
| 稳定性 (half-life) | GEO/SRA 的 4sU-seq、SLAM-seq、BRIC-seq |
| 翻译效率 (Ribo-seq) | GWIPS-viz、RPFdb、Trips-Viz |
| 3'UTR 调控特征 | TargetScan、miRTarBase、POSTAR/ENCODE eCLIP、UTRdb |
| 免疫刺激 / IVT 可制造性 | 公开数据不足，规则特征近似 + 自建实验数据 |

## 工作约定

- 交流与文档语言：**中文**；序列/模型术语可保留英文（CAI、MFE、UTR 等）。
- 文档为 Markdown，代码块用于流程图、公式和数据格式示例。
- 设计文档中的结论以 `mRNA模型设计.md` 为准；如实现与文档冲突，先在文档中补充决策记录再动代码。
- 新增代码时建议的结构（供参考，建立时再定）：`data/`（下载与清洗管道）、`models/`（生成模型与预测器）、`training/`、`evaluation/`、`configs/`。
- 涉及合成/设计具体序列时注意安全合规：文档明确要求内置序列安全筛查（敏感病原体、毒素、受管制序列过滤），实现时不得省略该模块。
