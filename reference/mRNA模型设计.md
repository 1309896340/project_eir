用于 mRNA 设计的深度神经网络模型，核心目标不是“生成任意 RNA”，而是在给定蛋白质或功能目标的前提下，设计出 **稳定、可高效翻译、低免疫原性、易生产且满足安全约束** 的 mRNA 序列。一个合理的模型通常应设计成“生成模型 + 多任务预测模型 + 约束优化”的组合。

---

## 1. 明确设计目标

mRNA 设计通常要同时优化多个性质：

1. **编码正确蛋白质**
   - ORF 必须翻译成目标氨基酸序列。(ORF open reading frame，开放阅读框 指的是核酸序列中按三联体密码子连续阅读时，从起始密码子ATG/AUG开始、到同框终止密码子TAA,TAG,TGA/UAA,UAG,UGA结束，中间没有终止密码子打断的一段读框。它理论上可以翻译成一条多肽链。)
   - 同义密码子可以变化，但氨基酸序列不能改变。

2. **提高翻译效率**
   - 合理的密码子使用偏好。
   - 合适的起始区结构。
   - 避免过强的 5' 端二级结构影响核糖体扫描。

3. **提高 mRNA 稳定性**
   - 控制 GC 含量。(G：鸟嘌呤Guanine C：胞嘧啶Cytosine;含量越高二级结构越稳定但降低蛋白表达)
   - 优化局部和全局 RNA 二级结构。
   - 设计合适的 5' UTR、3' UTR 和 poly(A) 相关特征。（5′UTR 是起始密码子前的非翻译区，主要调控翻译起始；3′UTR 是终止密码子后的非翻译区，主要调控 mRNA 稳定性、定位和降解。两者都不编码主蛋白，但决定 mRNA 的“寿命”和“翻译效率”。

）

4. **降低不良免疫刺激**
   - 避免某些富集的免疫刺激性序列模式。
   - 控制双链 RNA 倾向。
   - 可考虑修饰核苷酸条件下的设计策略。

5. **提升可制造性**
   - 避免长同聚物、极端 GC 区域、重复序列。
   - 避免不希望出现的限制性酶切位点、剪接样信号、早停信号等。

6. **满足安全与合规约束**
   - 对敏感病原体、毒素或受管制序列设置过滤与审查机制。
   - 设计系统应包含序列安全筛查模块。

---

## 2. 输入与输出设计

### 输入

模型输入可以包括：

```text
目标蛋白氨基酸序列
+ 物种/细胞类型
+ 给药场景或组织类型
+ 是否使用修饰核苷酸
+ 期望优化目标权重
+ 设计约束
```

例如：

```text
Protein sequence: MKT...
Host: human
Cell type: dendritic cell / hepatocyte / general
Objectives: high translation, high stability, low innate immune activation
Constraints: GC 40–60%, avoid motifs, preserve amino acid sequence
```

### 输出

输出可以是完整 mRNA 设计：

```text
5' UTR + CDS + 3' UTR + poly(A) 设计参数
```

其中最关键的是：

- CDS 同义密码子序列；
- 5' UTR；
- 3' UTR；
- 序列性质评分；
- 不确定性评估；
- 是否违反约束。

---

## 3. 模型总体架构

推荐采用模块化架构：

```text
目标蛋白序列
      ↓
蛋白/密码子编码器
      ↓
条件生成模型
      ↓
候选 mRNA 序列
      ↓
多任务性质预测器
      ↓
约束过滤与多目标优化
      ↓
最终候选序列
```

可以拆成三个核心模块：

1. **条件生成模型**
2. **多任务预测模型**
3. **约束与优化模块**

---

## 4. 条件生成模型设计

### 4.1 生成对象

对于 CDS（coding sequence，编码序列。也就是 mRNA 上的蛋白质编码区。它从起始密码子开始，到终止密码子前结束；终止密码子不编码氨基酸。）设计，模型不是自由生成核苷酸，而是在每个氨基酸位置选择一个同义密码子。

例如亮氨酸 L 可由多个密码子编码：

```text
UUA, UUG, CUU, CUC, CUA, CUG
```

模型的任务是：

```text
给定氨基酸序列 A1, A2, ..., An
生成密码子序列 C1, C2, ..., Cn
其中 translate(Ci) = Ai
```

这样可以保证蛋白序列不变。

---

### 4.2 推荐架构：Transformer 编码器-解码器

一个常见设计是：

```text
Protein Encoder: 编码氨基酸上下文
Codon Decoder: 自回归生成同义密码子
```

结构类似：

```text
Amino acid sequence
      ↓
Protein Transformer Encoder
      ↓
Context representation
      ↓
Codon Transformer Decoder
      ↓
Codon distribution at each position
```

在第 i 个位置，模型输出：

```text
P(codon_i | amino_acid_i, sequence_context, previous_codons, constraints)
```

但输出空间只允许该氨基酸对应的同义密码子。

这种设计的优点是：

- 保证编码蛋白不变；
- 能利用长程上下文；
- 可学习密码子对、局部结构、GC 分布等隐含规律；
- 适合与约束解码结合。

---

## 5. 加入 RNA 结构信息

mRNA 的二级结构会显著影响稳定性和翻译效率，因此模型应显式考虑结构。

可以使用两种方式：

### 方式一：结构特征作为输入

对候选序列计算或预测：

- 最小自由能 MFE；
- pairing probability；
- 5' 端局部结构强度；
- 局部窗口 GC 含量；
- unpaired probability；
- accessibility。

然后将这些特征输入预测器。

### 方式二：结构图神经网络

将 RNA 看成图：

- 节点：核苷酸；
- 边：
  - 相邻碱基连接；
  - 碱基配对边；
  - 长程相互作用边。

然后用 GNN 或 Graph Transformer 学习结构表示：

```text
RNA sequence + predicted base-pair graph
          ↓
Graph Neural Network
          ↓
Structure-aware embedding
```

结构模块可以帮助预测：

- 翻译效率；
- 半衰期；
- 免疫刺激风险；
- 降解热点。

---

## 6. 多任务预测模型

生成模型产生候选序列后，需要一个打分模型评估它们。

多任务预测器可以同时预测：

```text
Translation efficiency
mRNA half-life
Protein expression level
Ribosome loading
Innate immune activation risk
GC content quality
Secondary structure penalty
Manufacturability score
```

模型结构可以是：

```text
RNA sequence
   ↓
Nucleotide Transformer / CNN-Transformer hybrid
   ↓
Shared representation
   ↓
多个任务预测头
```

例如：

```text
shared encoder → expression head
               → stability head
               → immunogenicity head
               → structure head
               → manufacturability head
```

每个任务使用不同损失函数：

- 回归任务：MSE、Huber loss；
- 分类任务：cross-entropy；
- 排序任务：pairwise ranking loss；
- 多目标优化：加权 reward 或 Pareto ranking。

---

## 7. 目标函数设计

可以将总体目标定义为：

```text
Score = w1 × Translation
      + w2 × Stability
      - w3 × Immunogenicity
      - w4 × StructurePenalty
      - w5 × ManufacturingPenalty
      - w6 × ConstraintViolation
```

其中：

- `Translation`：预测翻译效率；
- `Stability`：预测 mRNA 半衰期；
- `Immunogenicity`：免疫刺激风险；
- `StructurePenalty`：不利二级结构；
- `ManufacturingPenalty`：生产难度；
- `ConstraintViolation`：违反硬性约束的惩罚。

不同应用场景权重不同。

例如：

### 疫苗场景

可能更重视：

```text
表达水平 + 稳定性 + 适度免疫激活控制
```

### 蛋白替代疗法

可能更重视：

```text
高表达 + 长半衰期 + 低免疫刺激
```

### 细胞治疗场景

可能更重视：

```text
短期表达 + 安全性 + 可控降解
```

---

## 8. 约束解码

单纯生成模型可能产生不符合要求的序列，因此需要约束解码。

常见硬约束包括：

1. 翻译后蛋白序列必须完全一致；
2. 避免内部终止密码子；
3. GC 含量在指定范围；
4. 避免长同聚物，例如：

```text
AAAAAA, GGGGGG
```

5. 避免指定 motif；
6. 避免不希望出现的限制性酶切位点；
7. 避免强发夹结构，尤其是 5' 端附近；
8. 避免隐性剪接位点、polyadenylation-like 信号；
9. 避免高重复性区域；
10. 通过安全序列筛查。

可以采用：

- constrained beam search；
- top-k / nucleus sampling + filter；
- integer programming；
- reinforcement learning；
- genetic algorithm；
- simulated annealing；
- diffusion-based refinement。

推荐流程：

```text
先生成大量候选
→ 快速规则过滤
→ 多任务模型评分
→ 多目标排序
→ 结构与安全复核
→ 选出 top candidates
```

---

## 9. 训练数据设计

模型质量高度依赖数据。可使用的数据包括：

### 9.1 CDS 与表达数据

- 不同 CDS 变体的蛋白表达量；
- ribosome profiling 数据；
- reporter assay 数据；
- codon usage 数据；
- mRNA-seq 与蛋白质组学数据。

### 9.2 稳定性数据

- mRNA half-life 测量；
- degradation rate 数据；
- UTR library 数据；
- RNA-seq time course 数据。

### 9.3 翻译效率数据

- ribosome occupancy；
- polysome profiling；
- reporter gene expression；
- massively parallel reporter assays, MPRA。

### 9.4 结构数据

- RNA 二级结构预测；
- SHAPE-seq、DMS-seq 等实验结构数据；
- base-pairing probability。

### 9.5 免疫相关数据

- TLR/RIG-I/MDA5 相关刺激数据；
- cytokine induction 数据；
- 修饰核苷酸条件下表达数据。

---

## 10. 训练策略

一个实用训练流程可以是：

### 第一步：自监督预训练

用大规模天然 mRNA/CDS 序列做语言模型预训练。

任务包括：

- masked nucleotide prediction；
- masked codon prediction；
- next codon prediction；
- span corruption；
- contrastive learning。

目的：

```text
让模型学习天然 mRNA 的统计规律、密码子上下文和序列模式。
```

---

### 第二步：监督微调

用带实验标签的数据微调：

```text
sequence → expression/stability/translation labels
```

损失函数：

```text
L = L_expression + L_stability + L_translation + L_immunogenicity + ...
```

---

### 第三步：生成模型训练

训练条件生成器：

```text
protein sequence + constraints → optimized codon sequence
```

可以使用：

- teacher forcing；
- sequence-to-sequence training；
- reward-guided fine-tuning；
- reinforcement learning from experimental feedback。

---

### 第四步：主动学习闭环

设计若干候选 mRNA：

```text
模型设计 → 合成测试 → 获得实验结果 → 更新模型
```

这个闭环很重要，因为公共数据往往与实际递送系统、细胞类型和修饰条件不完全匹配。

---

## 11. 推荐的系统框架

一个较完整的 mRNA 设计系统可以这样：

```text
1. 输入目标蛋白序列和设计约束
2. 生成多个 CDS 同义密码子版本
3. 生成或选择 5' UTR / 3' UTR
4. 预测二级结构与可及性
5. 多任务模型预测表达量、稳定性和免疫风险
6. 规则过滤不合格序列
7. 多目标 Pareto 排序
8. 输出 top-N 候选
9. 进行实验验证
10. 将实验结果反馈给模型
```

---

## 12. 评价指标

模型评估不能只看语言模型 loss，还应看生物学指标。

### 离线指标

- amino acid preservation rate；
- GC content distribution；
- codon adaptation index；
- codon pair bias；
- MFE；
- 5' 端结构强度；
- motif violation rate；
- predicted expression score；
- predicted stability score。

### 实验指标

- 蛋白表达量；
- mRNA 半衰期；
- 翻译效率；
- 细胞毒性；
- 免疫因子诱导；
- 批次生产质量；
- 不同细胞类型泛化能力。

---

## 13. 关键设计原则

总结来说，mRNA 设计模型应遵循以下原则：

1. **蛋白约束优先**
   - CDS 生成必须保证翻译产物正确。

2. **多目标优化**
   - 不能只优化密码子适应性，还要考虑结构、稳定性、免疫和生产。

3. **结构感知**
   - RNA 二级结构应成为模型的重要输入或预测目标。

4. **细胞类型特异性**
   - 不同细胞、组织、递送系统下最优 mRNA 可能不同。

5. **生成与筛选结合**
   - 生成模型负责提出候选，预测模型和规则系统负责筛选。

6. **实验闭环**
   - 仅靠计算预测不够，应使用实验数据持续更新模型。

7. **安全约束内置**
   - 从输入、生成、筛选到输出都应进行安全审查。

---

一个较优的深度学习 mRNA 设计模型可以概括为：

```text
条件 Transformer 生成同义 CDS
+ RNA 结构感知模块
+ 多任务表达/稳定性/免疫预测器
+ 约束解码与多目标优化
+ 主动学习实验闭环
```

这样的设计既能保证目标蛋白正确，又能系统性优化 mRNA 的表达、稳定性、安全性和可制造性。

如果选用 **Transformer** 作为 mRNA 设计模型的基础架构，可以从以下几个层面做结构调整，使其更适合 mRNA 序列设计，而不是直接套用普通 NLP Transformer。

---

# 1. 输入表示层调整

## 1.1 从“单核苷酸 token”扩展为多粒度 token

普通 Transformer 通常以字符或词为 token。mRNA 设计中可以采用多种 token 粒度：

### 单核苷酸级别

```text
A, U, G, C
```

优点：

- 分辨率最高；
- 适合学习局部 motif、二级结构、UTR 特征。

缺点：

- 序列较长；
- 对 CDS 的同义密码子约束表达不直接。

---

### 密码子级别

```text
AUG, GCU, UUU, ...
```

优点：

- 特别适合 CDS 设计；
- 每个 token 对应一个氨基酸编码单元；
- 方便施加“同义密码子约束”。

缺点：

- 不适合 UTR，因为 UTR 不按密码子翻译；
- 难以捕捉单碱基级 motif。

---

### 混合粒度表示

推荐使用混合表示：

```text
5' UTR：单核苷酸 token
CDS：密码子 token
3' UTR：单核苷酸 token
```

或者：

```text
CDS 同时使用 nucleotide embedding + codon embedding
```

例如对 CDS 中第 i 个密码子：

```text
embedding_i = codon_embedding_i
            + amino_acid_embedding_i
            + nucleotide_composition_embedding_i
            + position_embedding_i
```

这样模型既知道当前密码子是什么，也知道它编码哪个氨基酸。

---

## 1.2 加入区域类型 embedding

mRNA 不同区域功能不同：

```text
5' cap proximal region
5' UTR
Kozak region
start codon
CDS
stop codon
3' UTR
poly(A) proximal region
```

可以仿照 BERT 的 segment embedding，增加 region embedding：

```text
x_i = token_embedding_i
    + position_embedding_i
    + region_embedding_i
```

例如：

```text
region = {5UTR, CDS, 3UTR, polyA, linker, motif}
```

这能帮助模型区分同样的碱基序列在不同区域中的不同意义。

---

## 1.3 加入氨基酸条件 embedding

对于 CDS 设计，核心约束是：

```text
生成的密码子必须编码给定氨基酸
```

所以每个 CDS 位置可以加入对应的氨基酸 embedding：

```text
x_i = codon_embedding_i
    + amino_acid_embedding_i
    + position_embedding_i
```

对于亮氨酸 L，模型只允许输出：

```text
UUA, UUG, CUU, CUC, CUA, CUG
```

这样 Transformer 会学习：

```text
P(codon_i | amino_acid_i, context)
```

而不是无限制地生成任意三联体。

---

# 2. 位置编码调整

普通 Transformer 的绝对位置编码未必适合 mRNA，因为 mRNA 中很多效应与相对位置、区域边界有关。

## 2.1 使用相对位置编码

推荐引入相对位置编码，例如：

```text
attention_score(i, j) += relative_position_bias(i - j)
```

mRNA 中局部上下文很重要，例如：

- 起始密码子附近；
- 5' 端前几十个碱基；
- stop codon 附近；
- 局部 GC-rich 区域；
- 局部发夹结构。

相对位置编码比单纯绝对位置编码更适合建模这些局部依赖。

---

## 2.2 区域内位置 + 全局位置

可以同时编码：

```text
global_position_embedding
+ local_region_position_embedding
```

例如：

```text
第 20 个 CDS 密码子
```

既有：

```text
global position: 在整条 mRNA 中的位置
region position: 在 CDS 内的位置
```

这对区分以下模式有帮助：

```text
5' UTR 第 20 nt
CDS 第 20 nt
3' UTR 第 20 nt
```

---

## 2.3 关键锚点相对位置编码

mRNA 中有几个关键锚点：

```text
5' end
start codon
stop codon
poly(A) start
```

可以加入相对这些锚点的位置特征：

```text
distance_to_5_end
distance_to_start_codon
distance_to_stop_codon
distance_to_polyA
```

例如：

```text
x_i = token_embedding_i
    + region_embedding_i
    + relative_to_start_embedding_i
    + relative_to_stop_embedding_i
```

这有助于模型学习起始区结构和终止区上下文的影响。

---

# 3. Attention 机制调整

## 3.1 局部 attention + 全局 attention

mRNA 序列较长，完整 self-attention 的复杂度是：

```text
O(n²)
```

对于长 mRNA 不够高效。可以采用：

```text
local window attention + selected global tokens
```

例如：

- 每个位置关注附近 ±k 个 token；
- start codon、stop codon、UTR 边界作为 global token；
- 区域 summary token 参与全局交互。

结构类似 Longformer / BigBird：

```text
局部序列依赖：local attention
长程调控关系：global attention
```

适合处理：

- CDS 内局部密码子偏好；
- 5' UTR 对翻译起始的影响；
- 3' UTR 对稳定性的影响；
- 长程 RNA 二级结构。

---

## 3.2 结构偏置 attention

RNA 的二级结构不是简单线性序列，而是存在碱基配对：

```text
A-U
G-C
G-U wobble
```

可以把预测或实验得到的 base-pairing probability 加入 attention bias：

```text
attention_score(i, j) = Q_i K_j^T / sqrt(d)
                      + position_bias(i, j)
                      + structure_bias(i, j)
```

其中：

```text
structure_bias(i, j) = f(P_pair(i, j))
```

如果 i 和 j 可能形成碱基配对，则提高或调整它们之间的注意力权重。

这样模型不仅关注线性相邻位置，也关注空间结构相邻位置。

---

## 3.3 区域感知 attention mask

不同区域之间的交互强度可以不同。

例如：

- 5' UTR 与起始密码子区域强相关；
- CDS 内相邻密码子强相关；
- 3' UTR 与稳定性相关；
- poly(A) 附近与 3' UTR 相关。

可以设计 region-aware attention bias：

```text
attention_score(i, j) += B[region_i, region_j]
```

其中 `B` 是可学习参数矩阵，例如：

```text
B[5UTR, start_region]
B[CDS, CDS]
B[3UTR, 3UTR]
```

使模型自动学习不同功能区域之间的依赖关系。

---

# 4. 编码器-解码器结构调整

## 4.1 用 Encoder 编码目标蛋白，用 Decoder 生成密码子

对于 CDS 同义优化，推荐使用 encoder-decoder Transformer：

```text
Amino acid sequence
        ↓
Protein Encoder
        ↓
Codon Decoder
        ↓
mRNA CDS sequence
```

模型形式：

```text
P(C | A) = ∏ P(codon_i | codon_<i, amino_acid_sequence)
```

其中：

- `A` 是目标氨基酸序列；
- `C` 是密码子序列。

优点：

- 天然适合“由蛋白生成 mRNA”；
- 解码阶段可以加入同义密码子限制；
- 可以结合 beam search、采样和多目标打分。

---

## 4.2 非自回归 Transformer

因为每个氨基酸位置的密码子可以并行选择，也可以使用非自回归结构：

```text
amino acid sequence
        ↓
Transformer Encoder
        ↓
每个位置输出同义密码子分布
```

即：

```text
P(C | A) ≈ ∏ P(codon_i | A, context)
```

优点：

- 推理速度快；
- 适合大规模候选生成；
- 方便与后处理优化结合。

缺点：

- 对密码子之间的依赖建模弱一些；
- 可通过迭代 refinement 弥补。

---

## 4.3 Masked refinement 结构

可以采用类似 BERT 的迭代式设计：

```text
初始密码子序列
        ↓
mask 部分位置
        ↓
Transformer 预测更优同义密码子
        ↓
重复迭代
```

流程：

```text
1. 根据普通密码子偏好生成初始 CDS
2. 随机或按低置信度 mask 某些密码子
3. Transformer 重新填充
4. 用预测器评分
5. 迭代优化
```

这种结构适合在保持蛋白不变的条件下逐步优化 mRNA。

---

# 5. 输出层调整

## 5.1 同义密码子约束输出层

普通语言模型会输出完整词表：

```text
P(token_i | context)
```

但 CDS 中第 i 个位置只能输出编码对应氨基酸的同义密码子。因此需要 mask logits：

```text
logits[codon not in SynonymousCodons(amino_acid_i)] = -∞
```

例如：

```text
if amino_acid_i = F:
    allowed = {UUU, UUC}
```

则只在这两个密码子中归一化。

这样保证：

```text
translate(generated_CDS) == target_protein
```

---

## 5.2 多头输出

除了生成密码子，还可以让模型同时预测一些辅助属性：

```text
codon choice head
GC content head
local structure head
translation efficiency head
stability head
motif risk head
```

例如：

```text
Transformer hidden states
      ├── codon generation head
      ├── local MFE prediction head
      ├── expression prediction head
      ├── half-life prediction head
      └── immunogenicity risk head
```

这样可以通过多任务学习使表示更有生物学意义。

---

# 6. 结构信息融合方式调整

## 6.1 Transformer + GNN 混合结构

RNA 二级结构可以表示成图：

```text
节点：核苷酸
边：相邻边 + 碱基配对边
```

可采用：

```text
Sequence Transformer
        ↓
Graph Neural Network / Graph Transformer
        ↓
Structure-aware representation
```

或者并行：

```text
序列分支：Transformer
结构分支：GNN
        ↓
融合层
        ↓
预测/生成头
```

融合方式：

```text
h_i = concat(h_i_sequence, h_i_structure)
```

或：

```text
h_i = h_i_sequence + gate_i × h_i_structure
```

其中 gate 可以学习当前结构信息的重要程度。

---

## 6.2 Attention 中直接加入配对概率矩阵

如果有 RNA folding 工具给出的 pairing probability matrix：

```text
P_pair ∈ R^{n×n}
```

可以作为 attention bias：

```text
A = softmax(QK^T / sqrt(d) + B_pos + αP_pair)
```

其中 α 是可学习参数。

这比简单把 MFE 作为全局特征更细粒度。

---

# 7. 长序列建模调整

完整 mRNA 可能从几百到几千甚至上万 nt，标准 Transformer 不一定合适。

可以考虑以下结构。

## 7.1 Sparse Transformer

降低 attention 复杂度：

```text
full attention: O(n²)
sparse attention: O(n√n) 或 O(n log n)
```

适合长 mRNA。

---

## 7.2 Longformer / BigBird 风格结构

设计：

```text
local sliding window attention
+ global attention tokens
+ random sparse attention
```

全局 token 可包括：

```text
[CLS]
[5UTR_SUMMARY]
[START]
[CDS_SUMMARY]
[STOP]
[3UTR_SUMMARY]
```

---

## 7.3 Hierarchical Transformer

将 mRNA 分块：

```text
若干 nucleotide/codon blocks
        ↓
局部 Transformer 编码每个 block
        ↓
block-level Transformer 建模长程关系
```

例如：

```text
第一层：每 30 nt 或 10 codons 一个 block
第二层：block embedding 之间做 attention
```

优点：

- 适合超长 mRNA；
- 兼顾局部 motif 和全局组成；
- 计算量更低。

---

# 8. 条件控制结构调整

mRNA 设计通常依赖使用场景：

```text
宿主物种
细胞类型
递送系统
是否使用修饰核苷酸
目标表达水平
稳定性偏好
免疫风险限制
```

可以把这些条件作为控制 token 输入模型：

```text
[HOST=human] [CELL=hepatocyte] [MOD=pseudouridine] [GOAL=high_expression]
```

结构上类似 conditional Transformer：

```text
condition embedding + sequence embedding → Transformer
```

或者采用 cross-attention：

```text
sequence hidden states attend to condition embeddings
```

这样模型可以根据不同应用场景生成不同风格的 mRNA。

---

# 9. 解码策略结构调整

## 9.1 Constrained beam search

在 beam search 中加入硬约束：

- 同义密码子约束；
- GC 范围；
- 避免指定 motif；
- 避免长同聚物；
- 避免内部 stop codon；
- 控制局部结构强度。

解码时对非法候选直接剪枝。

---

## 9.2 Reward-guided decoding

每生成一部分序列，可以用辅助模型打分：

```text
score = log P_model
      + λ1 expression_score
      + λ2 stability_score
      - λ3 immunogenicity_score
      - λ4 structure_penalty
```

然后选择综合分更高的候选。

---

## 9.3 Pareto 解码

因为 mRNA 设计是多目标问题，不一定存在单一最优解。

可以输出 Pareto front：

```text
候选 A：最高表达
候选 B：最高稳定性
候选 C：最低免疫风险
候选 D：表达/稳定性平衡
```

模型结构上可以输出多个候选，再由多目标排序模块筛选。

---

# 10. 预训练任务调整

Transformer 的能力很大程度来自预训练。mRNA 场景中可以设计专门任务。

## 10.1 Masked nucleotide modeling

随机 mask 单个碱基：

```text
AUGG[MASK]UAC
```

预测被 mask 的碱基。

适合学习局部 motif。

---

## 10.2 Masked codon modeling

对 CDS 随机 mask 一个或多个密码子：

```text
AUG [MASK] GCU UAC
```

预测同义或真实密码子。

更适合 CDS 设计。

---

## 10.3 Amino acid to codon denoising

输入氨基酸序列和被扰动的密码子序列，让模型恢复更自然或更优的密码子选择。

```text
protein + noisy codons → optimized codons
```

---

## 10.4 Structure-aware pretraining

让模型预测：

- 某个位置是否 paired；
- base-pairing probability；
- 局部 MFE；
- unpaired probability；
- accessibility。

例如：

```text
Transformer hidden state → paired/unpaired prediction
```

---

## 10.5 Contrastive learning

构造正负样本：

```text
高表达 mRNA vs 低表达 mRNA
稳定 mRNA vs 不稳定 mRNA
天然 CDS vs 随机同义 CDS
```

让模型学习生物学上有意义的表示。

---

# 11. Transformer Block 内部可改造点

## 11.1 Feed-forward network 加入生物特征门控

普通 Transformer block：

```text
Attention → FFN
```

可以在 FFN 中加入生物特征 gate：

```text
h' = FFN(h)
g = sigmoid(W × bio_features)
output = h + g × h'
```

bio_features 可以包括：

- 局部 GC；
- codon rarity；
- pairing probability；
- region type；
- position relative to start codon。

---

## 11.2 Adapter / LoRA 用于不同细胞类型

不同应用场景可能需要不同模型，但从头训练成本高。

可以采用：

```text
共享主干 Transformer
+ cell-type-specific adapter
+ species-specific adapter
+ delivery-specific adapter
```

例如：

```text
human adapter
mouse adapter
hepatocyte adapter
dendritic-cell adapter
LNP-delivery adapter
```

这样可以在少量数据上快速微调。

---

## 11.3 Mixture-of-Experts

可以设计多个 expert：

```text
expression expert
stability expert
structure expert
immune-risk expert
manufacturability expert
```

通过 gating network 动态选择：

```text
h = Σ gate_k(x) × expert_k(x)
```

适合多目标 mRNA 设计。

---

# 12. 推荐的 Transformer 架构组合

一个比较实用的架构可以是：

```text
Protein Encoder:
    Transformer encoder
    输入：氨基酸序列

Codon Decoder:
    Transformer decoder
    cross-attention 到 protein encoder
    输出：同义密码子分布
    使用 synonymous codon mask

Structure-aware Module:
    attention bias 或 GNN
    融合 RNA 二级结构预测信息

Property Predictor:
    共享 Transformer 表示
    多任务预测：
        expression
        stability
        local structure
        immune risk
        manufacturability

Constrained Decoder:
    beam search / sampling
    加入 GC、motif、结构等约束
```

简化结构如下：

```text
目标蛋白序列
      ↓
Amino Acid Transformer Encoder
      ↓
Codon Transformer Decoder
      ↓
同义密码子约束输出层
      ↓
候选 CDS
      ↓
Structure-aware Transformer / GNN
      ↓
多任务性质预测器
      ↓
约束过滤与多目标排序
      ↓
最终 mRNA 候选
```

---

# 13. 最关键的结构调整总结

如果只挑最重要的几点，建议优先做：

1. **密码子级 tokenization**
   - CDS 用 codon token，而不是普通字符 token。

2. **同义密码子 mask 输出层**
   - 保证生成序列翻译成目标蛋白。

3. **region embedding**
   - 区分 5' UTR、CDS、3' UTR 等功能区域。

4. **相对位置编码**
   - 尤其是相对 start codon、stop codon、5' end 的位置。

5. **结构感知 attention**
   - 将 RNA base-pairing probability 或结构图加入 attention。

6. **长序列 sparse/local attention**
   - 适配几千 nt 的完整 mRNA。

7. **多任务预测头**
   - 同时预测表达、稳定性、免疫风险和可制造性。

8. **条件控制 token**
   - 支持细胞类型、物种、修饰核苷酸、设计目标等条件。

9. **constrained decoding**
   - 在生成阶段硬性满足生物与生产约束。

10. **adapter / MoE**
   - 用于不同组织、细胞类型和递送场景的快速适配。

总体来说，mRNA 设计中的 Transformer 不应只是普通语言模型，而应改造成：

```text
条件控制的、结构感知的、长序列友好的、带同义密码子约束和多任务优化头的 Transformer。
```

训练 mRNA 设计 Transformer 模型的数据，最好不要只是“很多 mRNA 序列”，而应构建成一个 **多层次、多模态、多任务、带实验标签的数据体系**。核心思想是：

```text
天然序列数据用于预训练
实验功能数据用于监督微调
结构/免疫/生产特征用于多任务学习
模型设计-实验验证数据用于主动学习闭环
```

可以把数据分成 5 类：

1. **基础序列数据**
2. **同义密码子设计数据**
3. **功能标签数据**
4. **结构与生物物理特征数据**
5. **实验闭环数据**

---

# 1. 数据总体形式

一个理想样本可以长这样：

```json
{
  "sample_id": "xxx",
  "species": "human",
  "cell_type": "HEK293T",
  "delivery": "LNP / electroporation / transfection",
  "modification": "unmodified / m1Ψ / Ψ",
  "mrna_sequence": "AUGGCU...",
  "utr5": "...",
  "cds": "...",
  "utr3": "...",
  "polyA_length": 120,
  "protein_sequence": "MA...",
  "codon_sequence": ["AUG", "GCU", "..."],
  "region_annotation": ["5UTR", "CDS", "CDS", "3UTR"],
  "labels": {
    "protein_expression": 1.27,
    "mrna_half_life": 6.8,
    "translation_efficiency": 0.92,
    "immune_activation": 0.14,
    "cell_viability": 0.96
  },
  "computed_features": {
    "gc_content": 0.54,
    "cai": 0.78,
    "mfe": -230.5,
    "five_prime_mfe": -11.2,
    "homopolymer_max_len": 4,
    "repeat_score": 0.08
  }
}
```

其中有些字段可以缺失，但最好保留 metadata，方便后续做条件建模。

---

# 2. 第一类：天然 mRNA / CDS 序列数据

这是 Transformer 预训练的基础。

## 2.1 需要的数据

包括：

```text
完整 mRNA 序列
CDS 序列
5' UTR 序列
3' UTR 序列
蛋白质序列
物种信息
基因 ID
转录本 ID
区域注释
```

可来自：

- RefSeq
- Ensembl
- GENCODE
- UniProt 对应蛋白序列
- NCBI nucleotide/protein 数据库
- 各物种基因组注释

对于 mRNA 设计，建议优先收集：

```text
human
mouse
rat
hamster/CHO
monkey
常用表达系统相关物种
```

如果模型用于人类治疗相关设计，人类和哺乳动物数据权重应更高。

---

## 2.2 这些数据的用途

天然序列数据主要用于自监督预训练：

### nucleotide-level 任务

```text
AUGGCUACU → mask 部分碱基 → 预测碱基
```

### codon-level 任务

```text
AUG GCU ACU → mask 一个密码子 → 预测密码子
```

### region prediction

```text
输入序列 → 预测该 token 属于 5' UTR / CDS / 3' UTR
```

### next token / denoising

```text
扰动后的 mRNA → 还原天然 mRNA
```

这些任务能让模型学习：

- 密码子使用偏好；
- UTR 序列模式；
- Kozak 区域规律；
- stop codon 附近上下文；
- GC 分布；
- 天然 mRNA 的局部 motif 和长程依赖。

---

# 3. 第二类：同义密码子设计数据

如果模型目标是：

```text
给定蛋白质序列 → 生成优化 CDS
```

就需要构建如下配对数据：

```text
protein_sequence → native_CDS
```

或者：

```text
amino_acid_sequence → codon_sequence
```

例如：

```json
{
  "protein_sequence": "MALWMRLLPLL...",
  "codon_sequence": ["AUG", "GCU", "CUG", "..."],
  "species": "human",
  "gene_expression_level": "high",
  "tissue": "liver"
}
```

---

## 3.1 可构建多版本同义样本

天然 CDS 只有一个版本，但同一个蛋白可以有大量同义 CDS。可以构建：

```text
天然 CDS
随机同义 CDS
按密码子频率采样的 CDS
高 CAI CDS
低 CAI CDS
GC 优化 CDS
结构优化 CDS
```

这样可以形成对比数据：

```text
same protein, different synonymous CDS, different scores
```

例如：

```json
{
  "protein_sequence": "MKT...",
  "cds_variants": [
    {
      "cds": "AUGAAAACC...",
      "type": "native",
      "score": 1.0
    },
    {
      "cds": "AUGAAGACG...",
      "type": "random_synonymous",
      "score": 0.4
    },
    {
      "cds": "AUGAAGACC...",
      "type": "human_codon_optimized",
      "score": 0.8
    }
  ]
}
```

这些数据可以用于：

- 训练生成模型；
- 训练排序模型；
- 训练对比学习模型；
- 训练 reward model。

---

# 4. 第三类：功能标签数据

仅有天然序列不够，因为设计目标通常是：

```text
高表达
高稳定性
低免疫刺激
好生产
```

所以需要带实验标签的数据。

---

## 4.1 蛋白表达数据

这是最关键标签之一。

样本形式：

```text
mRNA sequence / CDS variant / UTR variant
→ protein expression level
```

常见标签：

```text
荧光强度
luciferase activity
ELISA 蛋白浓度
western blot 定量
flow cytometry MFI
mass spec protein abundance
```

需要记录：

```text
细胞类型
转染方式
剂量
时间点
mRNA 修饰类型
UTR 类型
poly(A) 长度
实验批次
归一化方法
```

推荐标签格式：

```json
{
  "expression_6h": 0.82,
  "expression_24h": 1.35,
  "expression_48h": 0.77,
  "expression_auc": 2.94
}
```

相比单一时间点，**时间序列表达曲线**更有价值。

---

## 4.2 mRNA 稳定性数据

标签可以包括：

```text
mRNA half-life
decay rate
remaining RNA fraction at different time points
```

形式：

```json
{
  "rna_remaining": {
    "0h": 1.0,
    "2h": 0.72,
    "4h": 0.51,
    "8h": 0.28
  },
  "half_life": 3.9
}
```

这种数据用于训练：

```text
sequence → stability
```

模型可以学习：

- UTR 对稳定性的影响；
- CDS 结构对降解的影响；
- AU-rich / GC-rich 模式；
- miRNA 结合位点附近效应；
- poly(A) 邻近区域影响。

---

## 4.3 翻译效率数据

表达量受 mRNA 稳定性和翻译效率共同影响，因此最好单独收集翻译相关数据：

```text
ribosome profiling
polysome profiling
ribosome loading
translation efficiency = ribosome occupancy / mRNA abundance
```

标签形式：

```json
{
  "ribosome_density": 1.18,
  "translation_efficiency": 0.76
}
```

这类数据有助于区分：

```text
高表达是因为 mRNA 更稳定
还是因为单位 mRNA 翻译效率更高
```

---

## 4.4 免疫刺激与细胞毒性数据

对于治疗性 mRNA，免疫刺激风险非常重要。

标签可以包括：

```text
IFN-α
IFN-β
IL-6
TNF-α
ISG expression
TLR activation reporter
RIG-I/MDA5 activation
cell viability
cytotoxicity
```

样本形式：

```json
{
  "immune_activation": {
    "IFNB1_fold_change": 2.1,
    "IL6_pg_ml": 35.2,
    "TNFA_pg_ml": 12.4,
    "ISG_score": 0.31
  },
  "cell_viability": 0.94
}
```

注意需要记录：

```text
是否使用修饰核苷酸
纯化方式
dsRNA 杂质水平
递送系统
细胞类型
剂量
时间点
```

否则模型容易把实验工艺差异误认为序列效应。

---

## 4.5 可制造性数据

可制造性可作为规则特征或监督标签。

可记录：

```text
体外转录产率
完整性
截短产物比例
HPLC 纯化表现
同聚物长度
重复序列
极端 GC 区域
限制性位点
异常二级结构
```

标签示例：

```json
{
  "ivt_yield": 1.43,
  "full_length_fraction": 0.91,
  "truncated_fraction": 0.06,
  "manufacturability_score": 0.84
}
```

---

# 5. 第四类：RNA 结构数据与计算特征

mRNA 结构对表达、稳定性和免疫识别都有影响，所以建议为每条序列计算结构特征。

## 5.1 全局结构特征

```text
MFE
ensemble free energy
partition function
GC content
AU content
base-pairing ratio
average unpaired probability
```

---

## 5.2 局部窗口特征

对序列滑动窗口计算：

```text
局部 GC
局部 MFE
局部 accessibility
局部 pairing probability
局部 homopolymer
局部重复性
```

特别重要的是：

```text
5' UTR
start codon 上下游区域
CDS 前 30–100 nt
stop codon 附近
3' UTR 关键区域
```

样本可包含：

```json
{
  "local_features": [
    {
      "start": 0,
      "end": 50,
      "gc": 0.48,
      "mfe": -5.2,
      "unpaired_prob": 0.71
    },
    {
      "start": 50,
      "end": 100,
      "gc": 0.55,
      "mfe": -8.7,
      "unpaired_prob": 0.63
    }
  ]
}
```

---

## 5.3 base-pairing probability matrix

如果模型使用结构感知 attention，可以保存：

```text
P_pair ∈ R^{L×L}
```

不过完整矩阵很大，可以稀疏保存：

```json
{
  "base_pairs": [
    {"i": 12, "j": 86, "prob": 0.34},
    {"i": 13, "j": 85, "prob": 0.41}
  ]
}
```

---

## 5.4 实验结构数据

如果能获得 SHAPE-seq、DMS-seq 等结构探测数据，可以作为高质量结构标签：

```text
nucleotide → reactivity
paired/unpaired probability
```

标签形式：

```json
{
  "shape_reactivity": [0.12, 0.77, 0.34, ...]
}
```

这类数据可以训练模型更真实地理解 RNA 结构，而不是完全依赖计算预测。

---

# 6. 第五类：UTR 数据

完整 mRNA 设计不仅是 CDS，同样需要 UTR 数据。

## 6.1 5' UTR 数据

5' UTR 影响翻译起始，需记录：

```text
5' UTR 序列
长度
GC content
upstream AUG
uORF
Kozak context
start codon 邻近结构
translation initiation 标签
```

样本：

```json
{
  "utr5": "GGGAGAC...",
  "cds_start_context": "GCCACCAUGG",
  "uorf_count": 0,
  "upstream_aug_count": 0,
  "translation_initiation_score": 0.91
}
```

---

## 6.2 3' UTR 数据

3' UTR 影响稳定性和定位，需要记录：

```text
3' UTR 序列
miRNA target sites
AU-rich elements
CPE elements
polyadenylation signal
RNA binding protein motifs
half-life 标签
```

样本：

```json
{
  "utr3": "UGUAUA...",
  "are_count": 2,
  "mirna_site_count": 5,
  "stability_score": 0.62
}
```

---

# 7. 数据构建流程

推荐按以下流程构建。

---

## Step 1：定义目标任务

先明确模型要做什么：

### 任务 A：只优化 CDS

```text
输入：目标蛋白序列
输出：同义密码子优化 CDS
```

需要重点收集：

```text
protein-CDS pair
codon usage
CDS variant expression data
structure features
```

---

### 任务 B：设计完整 mRNA

```text
输入：目标蛋白 + 设计目标
输出：5' UTR + CDS + 3' UTR + poly(A)
```

需要额外收集：

```text
UTR library data
translation initiation data
stability data
poly(A) 数据
```

---

### 任务 C：性质预测器

```text
输入：mRNA sequence
输出：expression / stability / immune risk / manufacturability
```

需要重点收集：

```text
实验标签数据
metadata
负样本和低表现样本
```

---

## Step 2：收集天然序列并标准化

标准化内容包括：

```text
统一 U/T 表示
统一方向为 5'→3'
去除含 N 或未知碱基过多的序列
检查 CDS 长度是否为 3 的倍数
检查 start codon / stop codon
翻译 CDS 并与蛋白序列比对
切分 5' UTR / CDS / 3' UTR
```

建议保存两套表示：

```text
RNA alphabet: A U G C
DNA alphabet: A T G C
```

模型内部可以统一用 RNA 表示。

---

## Step 3：构建 tokenization

可以同时保存多种 token：

```text
nucleotide tokens
codon tokens
amino acid tokens
region tokens
position tokens
```

例如：

```json
{
  "nt_tokens": ["A", "U", "G", "G", "C", "U"],
  "codon_tokens": ["AUG", "GCU"],
  "aa_tokens": ["M", "A"],
  "region_tokens": ["CDS", "CDS"],
  "codon_index": [0, 1],
  "nt_position_in_codon": [0, 1, 2, 0, 1, 2]
}
```

---

## Step 4：构建监督标签

不同实验来源的标签尺度不一致，需要统一处理。

常见处理：

```text
log transform
z-score within experiment
normalize to positive control
batch correction
rank normalization
quantile normalization
```

例如不同实验中的表达值不能直接混合，可以使用：

```text
relative expression = sample expression / control expression
```

或：

```text
within-assay percentile rank
```

推荐同时保存：

```text
raw_label
normalized_label
normalization_method
batch_id
control_id
```

---

## Step 5：计算特征

为每条序列计算：

### 序列组成特征

```text
GC content
局部 GC
CAI
tAI
codon pair bias
CpG frequency
UpA frequency
k-mer frequency
homopolymer length
repeat score
```

### CDS 特征

```text
rare codon count
codon adaptation index
codon pair score
first 30–50 codons 的结构和稀有密码子分布
```

### UTR 特征

```text
uAUG count
uORF count
Kozak score
miRNA site
ARE motif
RBP motif
polyadenylation signal
```

### 结构特征

```text
MFE
局部 MFE
pairing probability
accessibility
start region ΔG
```

### 安全与过滤特征

```text
限制性位点
长重复序列
长同聚物
潜在剪接位点
不期望 motif
序列相似性筛查结果
```

---

## Step 6：构建正负样本和对比样本

只用“成功序列”训练会有偏差。应构建：

```text
高表达 vs 低表达
稳定 vs 不稳定
强免疫刺激 vs 弱免疫刺激
天然 CDS vs 随机同义 CDS
优化 CDS vs 未优化 CDS
```

特别推荐构建 pairwise 数据：

```json
{
  "protein_sequence": "MKT...",
  "variant_a": "AUGAAA...",
  "variant_b": "AUGAAG...",
  "label": "a_better_than_b",
  "criterion": "expression"
}
```

用于训练 ranking model：

```text
同一个蛋白的多个 mRNA 变体之间排序
```

这比直接预测绝对表达值更稳健。

---

## Step 7：数据划分

普通随机划分会造成信息泄漏。建议多种划分方式同时评估。

### 7.1 按序列随机划分

最简单，但容易高估性能。

```text
train / valid / test = 80 / 10 / 10
```

---

### 7.2 按基因或蛋白划分

同一蛋白的不同同义 CDS 不能同时出现在 train 和 test。

```text
train proteins
valid proteins
test proteins
```

这可以测试模型对新蛋白的泛化能力。

---

### 7.3 按实验批次划分

测试模型跨实验批次泛化。

```text
train assays
test unseen assay batch
```

---

### 7.4 按细胞类型划分

测试模型跨细胞类型泛化。

```text
train: HEK293T, HeLa
test: primary T cells
```

---

### 7.5 按物种划分

如果需要跨物种设计，可做：

```text
train: human/mouse
test: hamster/monkey
```

---

# 8. 数据库表结构建议

可以设计成关系型数据结构。

## 8.1 sequence 表

```text
sequence_id
full_mrna
utr5
cds
utr3
polyA_length
protein_sequence
species
gene_id
transcript_id
source
```

---

## 8.2 annotation 表

```text
sequence_id
region_start
region_end
region_type
feature_type
feature_value
```

例如：

```text
region_type = CDS / 5UTR / 3UTR
feature_type = uORF / miRNA_site / RBP_motif
```

---

## 8.3 experiment 表

```text
experiment_id
cell_type
tissue
delivery_method
dose
time_point
modification
batch_id
assay_type
protocol_summary
```

---

## 8.4 label 表

```text
sequence_id
experiment_id
label_type
raw_value
normalized_value
unit
time_point
```

---

## 8.5 computed_feature 表

```text
sequence_id
feature_name
feature_value
window_start
window_end
method
version
```

---

# 9. 针对不同训练阶段的数据构建

## 9.1 预训练数据

目标：

```text
学习 RNA / CDS / UTR 的统计规律
```

数据：

```text
大规模天然 mRNA、CDS、UTR
多物种转录本
高可信注释
```

训练任务：

```text
masked nucleotide prediction
masked codon prediction
region prediction
denoising
structure prediction
```

样本量：

```text
越大越好，百万级转录本/片段较理想
```

---

## 9.2 生成模型微调数据

目标：

```text
protein → synonymous CDS
```

数据：

```text
protein-CDS pairs
同一蛋白多个同义变体
优化变体与实验表现
```

训练方式：

```text
teacher forcing
sequence-to-sequence
masked codon refinement
ranking loss
reward-guided fine-tuning
```

---

## 9.3 性质预测器数据

目标：

```text
sequence → expression / stability / immune / manufacturability
```

数据：

```text
带实验标签的 mRNA variant library
MPRA 数据
reporter assay 数据
Ribo-seq
RNA half-life
immune activation assay
IVT 生产数据
```

注意：

```text
metadata 比序列本身同样重要
```

否则模型无法区分：

```text
序列导致的表达差异
细胞类型导致的差异
递送方式导致的差异
实验批次导致的差异
```

---

## 9.4 Reward model 数据

目标：

```text
给候选 mRNA 打综合分
```

可用数据：

```text
pairwise preference:
mRNA A 比 mRNA B 表达更高
mRNA A 比 mRNA B 免疫更低
mRNA A 综合优于 mRNA B
```

训练形式：

```text
P(A better than B) = sigmoid(score(A) - score(B))
```

---

# 10. 主动学习数据闭环

最终模型最好通过实验闭环提升。

流程：

```text
1. 模型生成一批候选 mRNA
2. 按不确定性和多样性选择待测序列
3. 合成并实验测试
4. 获得表达、稳定性、免疫等标签
5. 加入训练集
6. 更新生成模型和预测模型
```

选择待测序列时不要只选模型认为最好的，也要选：

```text
高不确定性样本
多样性样本
边界样本
负样本
不同 GC / 结构 / 密码子分布样本
```

这样可以避免模型陷入局部最优。

---

# 11. 数据质量控制

## 11.1 序列 QC

检查：

```text
是否只包含 A/U/G/C
CDS 长度是否为 3 的倍数
是否存在内部 stop codon
翻译是否匹配目标蛋白
UTR/CDS 边界是否正确
是否有过多未知碱基
是否方向正确
```

---

## 11.2 标签 QC

检查：

```text
实验重复一致性
异常值
批次效应
阴性/阳性对照是否正常
时间点是否一致
单位是否统一
归一化是否合理
```

---

## 11.3 数据泄漏检查

避免：

```text
同一 mRNA 变体同时出现在 train 和 test
同一蛋白的高度相似变体跨集合泄漏
同一实验批次产生的重复样本泄漏
```

---

# 12. 最小可行数据集设计

如果从零开始，可以先构建一个 MVP 数据集。

## 阶段 1：预训练集

```text
人类 + 小鼠高质量 mRNA/CDS/UTR
字段：
    sequence
    CDS boundary
    protein sequence
    species
    gene_id
```

用于训练：

```text
masked codon model
masked nucleotide model
```

---

## 阶段 2：监督预测集

收集或构建：

```text
同一 reporter protein 的大量 CDS/UTR 变体
每个变体有表达量、mRNA abundance、半衰期
```

字段：

```text
sequence
variant_type
expression
rna_abundance
translation_efficiency
cell_type
time_point
modification
batch
```

---

## 阶段 3：设计闭环集

模型设计：

```text
100–1000 条同义 CDS 候选
```

实验测试：

```text
表达量
mRNA abundance
细胞毒性
必要时免疫 marker
```

再加入训练集微调。

---

# 13. 推荐的数据样本格式

对于 Transformer 训练，最终可以转成如下格式。

## 13.1 CDS 生成样本

```json
{
  "input": {
    "protein_tokens": ["M", "A", "L", "W"],
    "species": "human",
    "cell_type": "HEK293T",
    "objective": "high_expression"
  },
  "target": {
    "codon_tokens": ["AUG", "GCU", "CUG", "UGG"]
  }
}
```

---

## 13.2 完整 mRNA 预测样本

```json
{
  "input": {
    "nt_tokens": ["G", "G", "A", "A", "U", "G", "..."],
    "region_tokens": ["5UTR", "5UTR", "5UTR", "CDS", "..."],
    "condition_tokens": {
      "species": "human",
      "cell_type": "HEK293T",
      "modification": "m1Ψ",
      "delivery": "LNP"
    }
  },
  "labels": {
    "expression_24h": 1.42,
    "half_life": 5.6,
    "immune_score": 0.18
  }
}
```

---

## 13.3 Pairwise ranking 样本

```json
{
  "protein_sequence": "MKT...",
  "sequence_a": "AUGAAAACC...",
  "sequence_b": "AUGAAGACG...",
  "condition": {
    "cell_type": "HEK293T",
    "modification": "unmodified"
  },
  "label": {
    "preferred": "a",
    "criterion": "protein_expression"
  }
}
```

---

# 14. 需要特别注意的问题

## 14.1 不要只学习天然偏好

天然 CDS 不一定是工程设计最优。天然序列受进化、调控、组织特异性和历史约束影响。

所以数据中应包含：

```text
天然序列
人工优化序列
随机同义序列
失败序列
实验筛选序列
```

否则模型可能只会生成“像天然”的序列，而不是“表达更好”的序列。

---

## 14.2 必须保留负样本

很多数据集只发表高性能序列，但模型需要知道什么是不好的。

负样本包括：

```text
低表达变体
不稳定变体
高免疫刺激变体
结构过强变体
生产失败变体
```

---

## 14.3 metadata 很重要

同一条 mRNA 在不同条件下表现不同：

```text
HEK293T vs primary immune cells
未修饰 mRNA vs m1Ψ 修饰 mRNA
LNP vs electroporation
6h vs 48h
```

所以标签不能脱离实验条件。

---

## 14.4 标签最好用相对值或排序值

不同实验平台的绝对表达值不可直接比较。更稳健的是：

```text
相对阳性对照表达
within-batch rank
pairwise preference
z-score normalized label
```

---

# 15. 总结

训练 mRNA 设计 Transformer 的数据应由以下部分构成：

```text
1. 大规模天然 mRNA/CDS/UTR 序列
   用于自监督预训练

2. protein-CDS 配对数据
   用于学习氨基酸到同义密码子的映射

3. 同一蛋白的多种同义 CDS/UTR 变体
   用于学习序列变化对功能的影响

4. 表达、稳定性、翻译效率、免疫刺激、生产表现等实验标签
   用于监督微调和性质预测

5. RNA 二级结构、GC、CAI、motif、repeat 等计算特征
   用于结构感知和约束优化

6. 主动学习产生的新实验数据
   用于持续提升模型
```

最推荐的数据构建策略是：

```text
天然序列预训练
+ 变体实验数据监督微调
+ 多任务标签建模
+ pairwise/ranking 数据增强
+ 主动学习实验闭环
```

这样训练出的模型才不仅能生成“合法 mRNA”，还能逐步学会生成在具体细胞类型和应用条件下更可能表现良好的 mRNA。

有，但需要注意：**目前基本没有一个“开箱即用、覆盖完整 mRNA 设计所有目标”的公开数据集**。公开数据通常是分散的：有的提供天然 mRNA/CDS 序列，有的提供 5'UTR 翻译效率，有的提供 RNA 结构探测，有的提供 mRNA 半衰期或 Ribo-seq。实际建模时通常要把多个公开数据源整合起来。

可以按用途分成几类。

---

# 1. 天然 mRNA / CDS / UTR 序列数据

这类数据最适合用于 **Transformer 自监督预训练**，例如 masked nucleotide/codon modeling。

| 数据源 | 内容 | 用途 |
|---|---|---|
| **GENCODE** | 人和小鼠高质量基因注释、转录本、CDS、UTR | mRNA/CDS/UTR 预训练 |
| **Ensembl** | 多物种基因组、转录本、蛋白注释 | 多物种模型训练 |
| **NCBI RefSeq** | curated mRNA、CDS、蛋白序列 | 高质量序列库 |
| **MANE Select** | 人类代表性转录本集合 | 减少 isoform 混乱 |
| **CCDS** | 高置信 CDS 注释 | CDS 训练 |
| **UCSC Genome Browser tables** | 基因、转录本、UTR 注释 | 辅助注释 |
| **RNAcentral** | 各类 RNA 序列整合库 | RNA 序列预训练 |

推荐优先使用：

```text
GENCODE human + GENCODE mouse
+ RefSeq curated mRNA/CDS
+ MANE Select human transcripts
```

对于 mRNA 设计，通常要从这些数据库提取：

```text
5' UTR
CDS
3' UTR
protein sequence
transcript ID
gene ID
species
```

---

# 2. 密码子使用与同义优化相关数据

这些数据可用于构建：

```text
protein sequence → native codon sequence
```

以及计算 CAI、codon usage、codon pair bias 等特征。

| 数据源 | 内容 | 用途 |
|---|---|---|
| **Kazusa Codon Usage Database** | 多物种密码子使用频率 | 基础密码子偏好 |
| **HIVE-CUTs** | 大规模物种密码子使用表 | 多物种 codon bias |
| **CoCoPUTs** | codon、codon pair、dinucleotide usage | 密码子对偏好、二核苷酸偏好 |
| **Codon Usage Tables from Ensembl/NCBI CDS** | 可自行统计 | 更可控、更可追溯 |

建议不要只依赖静态 codon usage table，而是自己从高质量 CDS 中统计：

```text
species-specific codon usage
tissue/high-expression gene codon usage
codon pair usage
GC3 distribution
CpG/UpA frequency
```

---

# 3. 5' UTR 翻译效率数据

这类数据对 mRNA 设计非常有价值，尤其用于训练：

```text
5' UTR sequence → translation efficiency
```

## 3.1 Optimus 5-Prime / Sample et al. dataset

这是最常被使用的公开 5'UTR 设计数据之一。

大致特点：

```text
大量合成 5' UTR 序列
reporter assay
标签通常是 mean ribosome load / translation efficiency 相关指标
```

用途：

```text
训练 5' UTR 翻译起始预测模型
学习 uAUG、uORF、GC、局部结构对翻译的影响
```

很多 5'UTR 深度学习模型，例如 Optimus 5-Prime 类模型，都基于类似的大规模 reporter 数据。

---

## 3.2 其他 5'UTR MPRA / reporter 数据

还有不少论文公开过：

```text
5' UTR library
Kozak context variant library
uORF variant library
translation initiation reporter data
```

这些数据通常不一定整理成统一数据库，而是放在：

```text
GEO
SRA
ArrayExpress
论文 supplement
GitHub repository
```

检索关键词可以用：

```text
5' UTR MPRA
5' UTR reporter assay
mean ribosome load
translation initiation library
uORF reporter
Kozak sequence library
```

---

# 4. 3' UTR、miRNA、RBP 调控数据

3'UTR 对 mRNA 稳定性、定位和翻译调控影响很大。

| 数据源 | 内容 | 用途 |
|---|---|---|
| **UTRdb / UTRsite** | UTR 序列与调控元件注释 | UTR motif 特征 |
| **TargetScan** | miRNA 靶位点预测 | 3'UTR 稳定性特征 |
| **miRTarBase** | 实验验证 miRNA-target | miRNA 调控标签 |
| **POSTAR / POSTAR3** | RBP binding sites | RBP motif/CLIP 特征 |
| **ENCODE eCLIP** | RBP 结合实验数据 | RBP 结合建模 |
| **ATtRACT / RBPDB** | RBP binding motifs | motif 特征 |
| **oRNAment** | RBP motif annotation | 3'UTR 调控特征 |

这类数据可以用于给 3'UTR 添加特征：

```text
miRNA site count
RBP motif count
AU-rich element
CPE motif
polyadenylation signal
RBP binding density
```

不过要注意，3'UTR 调控具有很强的细胞类型特异性。

---

# 5. RNA 二级结构与结构探测数据

如果你的 Transformer 要做结构感知 attention 或结构辅助任务，这类数据很重要。

## 5.1 计算/注释结构数据库

| 数据源 | 内容 | 用途 |
|---|---|---|
| **Rfam** | RNA families 与结构注释 | 结构预训练 |
| **RNA STRAND** | 已知 RNA 二级结构 | 结构监督 |
| **ArchiveII** | RNA 二级结构 benchmark | 结构预测任务 |
| **bpRNA / bpRNA-1m** | 大规模 RNA 二级结构注释 | paired/unpaired 预测 |
| **PDB / NDB** | RNA 三维结构 | 高质量结构样本 |

这类数据多为非编码 RNA，但仍可用于让模型学习：

```text
碱基配对
stem-loop
paired/unpaired pattern
局部结构 motif
```

---

## 5.2 实验结构探测数据

| 数据源 | 内容 | 用途 |
|---|---|---|
| **RMDB, RNA Mapping Database** | SHAPE、DMS 等结构探测数据 | 训练 reactivity / accessibility 预测 |
| **Stanford Ribonanza RNA Folding dataset** | 大规模 RNA 化学探测数据 | 结构感知预训练 |
| **Kaggle OpenVaccine dataset** | RNA 降解与结构相关标签 | degradation / reactivity 预测 |

### OpenVaccine 数据集

Kaggle 上的 **Stanford COVID Vaccine Degradation Prediction / OpenVaccine** 数据集比较常用。

它包含：

```text
RNA sequence
predicted structure
loop type
reactivity
deg_Mg_pH10
deg_pH10
deg_Mg_50C
deg_50C
```

适合训练：

```text
sequence → degradation propensity
sequence → reactivity
sequence → local stability
```

但要注意：

```text
它主要是短 RNA 片段，不是完整治疗性 mRNA。
```

因此可作为结构/降解辅助任务数据，而不是完整 mRNA 表达预测数据。

---

# 6. mRNA 稳定性 / half-life 数据

mRNA 半衰期数据公开不少，但常常分散在 GEO/SRA 或论文补充材料中。

常见来源包括：

```text
RNA-seq time course
4sU-seq
SLAM-seq
BRIC-seq
transcription shutoff assay
metabolic labeling
```

代表性方向：

| 数据类型 | 用途 |
|---|---|
| human/mouse mRNA half-life datasets | 训练稳定性预测 |
| yeast mRNA decay datasets | 预训练或跨物种分析 |
| RNA-seq time-course decay data | 拟合 decay rate |
| metabolic labeling 数据 | 更直接估计 half-life |

可检索关键词：

```text
mRNA half-life RNA-seq
mRNA decay rate human
4sU-seq mRNA stability
SLAM-seq mRNA half-life
BRIC-seq mRNA decay
```

常见公共平台：

```text
GEO
SRA
ArrayExpress
Expression Atlas
```

这些数据适合构建：

```text
transcript sequence + cell type → mRNA half-life
```

不过要注意：

```text
天然转录本的 half-life 受到 promoter、RNA binding protein、细胞状态、UTR、poly(A)、亚细胞定位等影响，
不能完全等同于体外合成 mRNA 的稳定性。
```

---

# 7. 翻译效率与 Ribo-seq 数据

翻译效率可以通过：

```text
ribosome profiling / Ribo-seq
polysome profiling
ribosome loading
```

估计。

| 数据源 | 内容 | 用途 |
|---|---|---|
| **GWIPS-viz** | Ribo-seq 数据浏览与整合 | ribosome occupancy |
| **RPFdb** | ribosome profiling 数据库 | 翻译效率 |
| **Trips-Viz** | Ribo-seq 可视化和分析 | ORF translation |
| **GEO/SRA Ribo-seq studies** | 原始 Ribo-seq 数据 | 自行构建 TE 标签 |

可构建标签：

```text
translation_efficiency = Ribo-seq RPKM / RNA-seq RPKM
```

用于训练：

```text
5'UTR + CDS + 3'UTR + condition → translation efficiency
```

注意事项：

```text
Ribo-seq 数据批次效应大
不同实验处理差异大
需要严格标准化
```

---

# 8. 表达量和蛋白丰度数据

这些数据更适合提供背景信息，而不是直接作为 mRNA 设计标签。

| 数据源 | 内容 | 用途 |
|---|---|---|
| **GTEx** | 人类组织 RNA expression | tissue-specific expression context |
| **Human Protein Atlas** | RNA + protein expression | 组织表达背景 |
| **Expression Atlas** | 多物种、多条件表达 | 表达上下文 |
| **recount3** | 大规模 RNA-seq 统一处理数据 | transcript abundance |
| **PaxDb** | 蛋白丰度 | protein abundance |
| **ProteomicsDB** | 蛋白质组数据 | protein expression context |

可用于：

```text
构建 tissue-specific codon usage
分析高表达基因的 CDS/UTR 特征
辅助预训练或条件建模
```

但不建议直接把天然基因表达量当作 mRNA 设计标签，因为天然表达量主要由 promoter、染色质、转录调控等决定，而体外递送 mRNA 不受这些因素控制。

---

# 9. 免疫刺激相关公开数据

这一类相对缺乏，尤其缺少：

```text
完整 mRNA sequence → innate immune activation
```

的标准公开大数据集。

可用资源包括：

| 数据源 | 内容 | 用途 |
|---|---|---|
| **ImmPort** | 免疫学实验数据 | 免疫 response 背景 |
| **GEO/SRA immune stimulation RNA-seq** | TLR/RIG-I/MDA5 刺激转录组 | 免疫通路标签 |
| **IEDB** | 免疫表位数据库 | 主要是蛋白/肽表位，不是 mRNA 先天免疫 |
| **ENCODE / CLIP/RBP 数据** | dsRNA/RBP 相关间接特征 | 辅助特征 |

实际 mRNA 免疫刺激高度依赖：

```text
修饰核苷酸
dsRNA 杂质
纯化方式
递送系统
细胞类型
剂量
时间点
```

所以公开数据只能部分辅助。若模型目标包含低免疫刺激，最好后续补充自有实验数据。

---

# 10. 可制造性 / IVT 产率数据

这一类公开数据最少。

理想标签包括：

```text
in vitro transcription yield
full-length fraction
truncated product fraction
dsRNA impurity
HPLC purification profile
capping efficiency
poly(A) integrity
```

但公开的大规模 sequence-to-IVT-yield 数据非常少。

实际做法通常是：

```text
先用规则特征近似：
    GC content
    local GC
    homopolymer
    repeat
    strong hairpin
    restriction sites
    long U-rich/A-rich regions
    cryptic splice/polyA signals

再通过自有实验积累 IVT 标签。
```

---

# 11. 比较推荐的公开数据组合

如果你要从公开数据起步，我建议按下面组合构建。

## 11.1 预训练数据

```text
GENCODE human/mouse
+ RefSeq mRNA/CDS
+ Ensembl 多物种 CDS
+ RNAcentral/Rfam 可选
```

训练任务：

```text
masked nucleotide modeling
masked codon modeling
region prediction
protein-to-CDS reconstruction
```

---

## 11.2 CDS / 密码子优化数据

```text
GENCODE/RefSeq protein-CDS pairs
+ CoCoPUTs/HIVE-CUTs/Kazusa codon usage
+ 高表达基因子集，例如 GTEx/HPA/recount3 辅助筛选
```

构建：

```text
amino acid sequence → native codon sequence
same protein → multiple synonymous variants
```

---

## 11.3 5'UTR 翻译数据

```text
Optimus 5-Prime / Sample et al. 5'UTR reporter dataset
+ 其他 5'UTR MPRA/reporter datasets from GEO/SRA
```

训练：

```text
5'UTR → translation efficiency
```

---

## 11.4 RNA 结构/降解数据

```text
OpenVaccine Kaggle dataset
+ Ribonanza RNA Folding dataset
+ RMDB
+ bpRNA / RNA STRAND / ArchiveII
```

训练：

```text
sequence → reactivity
sequence → degradation score
sequence → paired/unpaired
sequence → base-pairing/accessibility
```

---

## 11.5 稳定性和翻译效率数据

```text
GEO/SRA mRNA half-life datasets
+ GEO/SRA Ribo-seq datasets
+ GWIPS-viz / RPFdb / Trips-Viz
```

构建：

```text
transcript sequence + cell type → half-life
transcript sequence + condition → translation efficiency
```

---

# 12. 一个实际可行的数据构建方案

可以按三层来做。

## 第一层：大规模无标签预训练集

来源：

```text
GENCODE
RefSeq
Ensembl
RNAcentral
```

样本格式：

```json
{
  "species": "human",
  "transcript_id": "...",
  "utr5": "...",
  "cds": "...",
  "utr3": "...",
  "protein": "...",
  "region_annotation": "..."
}
```

用途：

```text
训练 RNA/CDS Transformer backbone
```

---

## 第二层：公开功能标签微调集

来源：

```text
Optimus 5-Prime
OpenVaccine
Ribonanza
RMDB
Ribo-seq datasets
mRNA half-life datasets
```

样本格式：

```json
{
  "sequence": "...",
  "region": "5UTR/CDS/full_transcript/RNA_fragment",
  "condition": {
    "cell_type": "...",
    "assay": "...",
    "time_point": "..."
  },
  "labels": {
    "translation_efficiency": 0.83,
    "half_life": 4.2,
    "reactivity": [...],
    "degradation": [...]
  }
}
```

用途：

```text
训练多任务预测头
```

---

## 第三层：自有实验闭环数据

公开数据很难覆盖治疗性 mRNA 的真实设计空间，所以最终最好补充：

```text
同一目标蛋白的多种同义 CDS
不同 5'UTR / 3'UTR 组合
不同修饰核苷酸条件
不同细胞类型
表达、稳定性、免疫、可制造性标签
```

用途：

```text
校准模型到自己的应用场景
```

---

# 13. 需要注意的坑

## 13.1 没有统一标签尺度

不同数据集的表达值、稳定性、反应性不能直接混合。

建议做：

```text
within-assay normalization
z-score
rank normalization
relative-to-control normalization
batch correction
```

---

## 13.2 天然表达数据不等于 mRNA 设计标签

GTEx、HPA 这类天然表达数据很有用，但天然基因表达受转录调控影响很大。

对于外源 mRNA 设计，更重要的是：

```text
UTR reporter assay
CDS variant assay
Ribo-seq
mRNA decay assay
synthetic mRNA expression assay
```

---

## 13.3 OpenVaccine/Ribonanza 不是完整 mRNA 表达数据

它们更适合训练：

```text
结构
降解
reactivity
局部稳定性
```

不能直接代表：

```text
完整 mRNA 在细胞内的蛋白表达
```

---

## 13.4 免疫和制造标签公开数据不足

这两类通常需要自建实验数据：

```text
immune activation
IVT yield
dsRNA impurity
full-length fraction
```

公开数据只能提供规则特征或间接特征。

---

# 14. 最简推荐清单

如果只想快速起步，可以先用：

```text
1. GENCODE/RefSeq
   用于 mRNA/CDS/UTR 预训练

2. CoCoPUTs 或 HIVE-CUTs
   用于密码子和密码子对特征

3. Optimus 5-Prime
   用于 5'UTR 翻译效率预测

4. OpenVaccine + Ribonanza + RMDB
   用于结构、reactivity、degradation 相关任务

5. GEO/SRA 中的 mRNA half-life 和 Ribo-seq 数据
   用于稳定性和翻译效率标签

6. TargetScan + POSTAR/ENCODE eCLIP + UTRdb
   用于 3'UTR 调控特征
```

总结一句：

```text
有公开数据，但不是一个完整数据集，而是需要把序列库、UTR reporter 数据、RNA 结构数据、mRNA half-life 数据、Ribo-seq 数据和调控注释库组合起来。
```

对于真正的 mRNA 设计模型，公开数据足够支撑预训练和初步微调，但要达到工程可用，通常还需要针对目标细胞、递送系统和修饰条件建立自有实验闭环数据。
