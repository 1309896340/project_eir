从RefSeq获取数据来训练mRNA设计模型，核心流程是：**按需下载 → 提取区域序列 → 严格QC过滤 → 防泄漏划分**。下面结合公开实践和你的项目目标，详细说明每一步。

---

## 一、从RefSeq获取数据

### 1.1 推荐方式：NCBI Datasets命令行工具

NCBI官方推荐使用 `datasets` 命令行工具下载基因数据包，这是最规范的方式。对于每个蛋白编码基因，可以指定下载包含 **CDS、5′UTR、3′UTR、转录本序列和蛋白序列** 的完整数据包：

```bash
# 安装
conda install -c conda-forge ncbi-datasets-cli

# 下载单个基因的全部区域序列
datasets download gene accession NM_000546.6 \
  --include cds,5p-utr,3p-utr,rna,protein \
  --filename TP53_data.zip

# 批量下载：从accession列表文件读取
datasets download gene accession \
  --inputfile refseq_accessions.txt \
  --include cds,5p-utr,3p-utr,rna,protein \
  --filename batch_data.zip
```

下载后解压，得到 `cds.fna`、`5p_utr.fna`、`3p_utr.fna`、`rna.fna`、`protein.faa` 等文件。

### 1.2 批量下载：FTP + 脚本

如果需要全物种或大规模数据，从NCBI Genomes FTP下载匹配的FASTA和GTF/GFF3文件更高效。RefSeq FTP提供 `*_genomic.gtf.gz` 和 `*_rna.fna.gz` 等文件，序列标识符一致，便于关联。

```bash
# 下载人类RefSeq注释
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz
```

**物种选择建议**：优先下载 **human、mouse、rat**，这三个物种有完整的RefSeq Select注释。如果算力允许，可扩展到115种哺乳动物，覆盖更广的密码子使用和UTR多样性。

### 1.3 使用RefSeq Select减少冗余

一个基因常有多个异构体。RefSeq Select为每个蛋白编码基因提供一个**代表性转录本**，基于临床数据库使用、表达量、编码区保守性、蛋白长度和Swiss-Prot一致性等标准自动选取。人类、小鼠、大鼠的RefSeq Select覆盖所有蛋白编码基因，且人类RefSeq Select与MANE项目整合。

在GFF3/GTF中，可通过 `tag=RefSeq Select` 或 `tag=MANE Select` 标签识别代表性转录本。**使用RefSeq Select可大幅减少冗余，避免同一基因的多个异构体造成数据泄漏。**


## 二、数据清洗

拿到原始数据后，需要经过严格的QC流程。参考mRNA设计模型（如mRNABERT、GEMORNA）的实践，清洗分为以下层次。

### 2.1 格式与方向标准化

- **U/T统一**：所有RNA序列转为 `A/U/G/C`，DNA序列转为 `A/T/G/C`，模型内部统一用RNA字母表。
- **方向标准化**：负链基因必须**反向互补**，统一为5′→3′方向。
- **字符集检查**：仅允许标准碱基，`N` 比例 < 1%，含模糊字符的序列过滤。

### 2.2 CDS结构检查（硬性门槛）

| 检查项 | 标准 | 处理 |
|---|---|---|
| CDS长度 | 是3的倍数 | 不满足则过滤 |
| 起始密码子 | `ATG` / `AUG` | 不满足则过滤或标记 |
| 终止密码子 | `TAA/TAG/TGA` | 不满足则过滤 |
| 内部终止密码子 | 无 | 有则过滤 |
| 移码 | 无 | 有则过滤 |

mRNASyner在预处理中首先**去除重复条目**，然后**排除核苷酸数不是3的倍数的序列**。mRNABERT额外采用 **40% CDS长度过滤**：预测CDS长度低于最大可能长度40%的序列被排除，并去除高同源性冗余序列。

### 2.3 翻译一致性验证（最核心的QC）

**将CDS翻译后，必须与RefSeq提供的蛋白序列完全一致。** 这是区分“真正可翻译的CDS”和“看起来像ORF但实际不编码”的关键。

```python
from Bio.Seq import Seq

def validate_translation(cds_rna, refseq_protein):
    """验证CDS翻译是否与RefSeq蛋白一致"""
    dna = cds_rna.replace("U", "T")
    translated = str(Seq(dna).translate(to_stop=True))
    return translated.rstrip("*") == refseq_protein.rstrip("*")
```

GEMORNA的数据验证标准是：**CDS长度必须恰好是蛋白长度的3倍加3（终止密码子）**，且密码子使用与密码子表一致。不满足的序列被过滤。

### 2.4 序列质量过滤

| 检查项 | 标准 | 依据 |
|---|---|---|
| 全长长度 | 100–15,000 nt | 排除过短/过长 |
| GC含量 | 30%–70% | 排除极端GC |
| 最长同聚物 | ≤ 8 | 避免生产困难 |
| 重复评分 | 低重复 | 避免合成失败 |
| 5′UTR长度 | 通常 < 500 nt | 过长影响翻译 |
| 3′UTR长度 | 通常 < 2,000 nt | 过长影响稳定性 |

GEMORNA的UTR预训练数据来自RefSeq和UTRdb，经过**长度和冗余过滤**，最终得到约800万条5′UTR和200万条3′UTR。

### 2.5 负样本构建

训练编码/非编码判别能力需要负样本：

- **RefSeq中的lncRNA、假基因**：`gene_type != "protein_coding"`
- **随机同义CDS**：对同一蛋白随机选择同义密码子生成的序列
- **框外ORF**：CDS框外的开放阅读框
- **uORF**：5′UTR中的上游ORF

mRNABERT从RefSeq等来源收集了超过**1800万条非冗余mRNA序列**，其中包含未注释数据，使用ORFfinder预测ORF，并做长度过滤以降低假注释风险。


## 三、训练/测试数据集准备

### 3.1 样本结构

每条训练样本应包含完整字段：

```json
{
  "refseq_id": "NM_000546.6",
  "gene_id": "7157",
  "gene_symbol": "TP53",
  "species": "human",
  "utr5": "GCCACC...",
  "cds": "AUGGAGGAG...",
  "utr3": "UGUCUG...",
  "protein_sequence": "MEEPQSD...",
  "full_mrna": "GCCACC...AUG...UGA...UGUCUG",
  "gc_content": 0.52,
  "cds_length": 1182,
  "is_refseq_select": true
}
```

### 3.2 数据划分：防止泄漏是关键

**绝对不能用随机划分。** 同一基因的不同异构体、高度相似的旁系同源基因如果跨train/test，会造成严重的数据泄漏，导致模型评估过于乐观。

推荐划分策略（优先级从高到低）：

| 策略 | 方法 | 适用场景 |
|---|---|---|
| **按基因划分** | 同一基因的所有异构体只出现在一个集合 | 最推荐，防止同源泄漏 |
| **按染色体划分** | train/valid/test 使用不同染色体 | 基础泛化评估 |
| **按物种划分** | train=human，test=mouse | 跨物种泛化 |
| **按实验批次划分** | 不同batch_id分开 | 跨批次泛化 |

mRNABERT将数据划分为 **80/10/10** 的train/validation/test，确保划分基于序列去冗余后的独立集合。ppLM-CO在过滤后使用**固定随机种子42**划分训练和测试集，确保可复现性。

**防泄漏检查清单**：

```python
def check_leakage(train_df, test_df):
    issues = []
    # 基因泄漏
    overlap = set(train_df["gene_id"]) & set(test_df["gene_id"])
    if overlap: issues.append(f"Gene leakage: {len(overlap)}")
    # 蛋白序列泄漏
    overlap = set(train_df["protein_sequence"]) & set(test_df["protein_sequence"])
    if overlap: issues.append(f"Protein leakage: {len(overlap)}")
    # 转录本泄漏
    overlap = set(train_df["refseq_id"]) & set(test_df["refseq_id"])
    if overlap: issues.append(f"Transcript leakage: {len(overlap)}")
    return issues
```

### 3.3 各区域分开准备

由于模型需要同时设计5′UTR、CDS和3′UTR，建议**三个区域分别构建数据集**：

- **CDS数据集**：蛋白序列 → 同义CDS。重点验证翻译一致性，过滤内部终止密码子。
- **5′UTR数据集**：UTR序列 → 翻译效率标签。可结合Optimus 5-Prime等公开reporter数据。
- **3′UTR数据集**：UTR序列 → 稳定性标签。可结合RNA half-life数据。

GEMORNA的微调策略值得参考：预训练用全量RefSeq数据，微调时用PRED-5UTR和PRED-3UTR模型筛选出**高MRL的5′UTR和高稳定性的3′UTR**，随机采样得到约80万条5′UTR和20万条3′UTR微调数据。

### 3.4 条件标签的整合

你的模型需要条件控制（细胞类型、修饰核苷酸、递送方式），因此需要额外收集：

| 条件 | 数据来源 | 处理方式 |
|---|---|---|
| 细胞类型 | GEO/SRA元数据 | 标准化为受控词表 |
| 修饰核苷酸 | 实验记录 | unmodified / m1Ψ / Ψ |
| 递送方式 | 实验记录 | LNP / electroporation |
| 表达量 | reporter assay | 相对对照、z-score |
| 半衰期 | RNA-seq time course | 拟合decay rate |
| 翻译效率 | Ribo-seq | Ribo-seq RPKM / RNA-seq RPKM |


## 四、完整流程总结

```text
RefSeq 数据获取
├── NCBI Datasets CLI → 按accession下载CDS/UTR/蛋白
├── FTP批量下载 → GTF + FASTA + 蛋白FASTA
└── 优先 RefSeq Select / MANE Select 减少冗余

        ↓

数据清洗
├── U/T统一、方向标准化（负链反向互补）
├── CDS结构检查（长度3的倍数、起始/终止/内部终止）
├── 翻译一致性验证（translate(CDS) == RefSeq蛋白）
├── 序列质量过滤（GC、同聚物、重复、N比例）
└── 负样本构建（lncRNA、随机同义、框外ORF）

        ↓

训练/测试集准备
├── 样本结构：utr5 + cds + utr3 + protein + metadata
├── 划分策略：按基因划分（防止同源泄漏）
├── 防泄漏检查：基因/蛋白/转录本/批次
├── 区域分离：CDS、5′UTR、3′UTR 各自成集
└── 条件标签整合：细胞类型、修饰、递送

        ↓

输出：可直接用于Transformer训练的mRNA设计数据集
```


## 五、关键注意事项

1. **RefSeq Select优先**：避免同一基因多个异构体造成冗余和泄漏。
2. **翻译一致性是硬门槛**：不通过则直接过滤，不可妥协。
3. **负样本必不可少**：模型需要学习区分可翻译和不可翻译序列。
4. **按基因划分**：这是防止数据泄漏最有效的方式，优于随机划分。
5. **区域分开处理**：CDS、5′UTR、3′UTR的生成逻辑和标签来源不同，应分别构建数据集。
6. **保留metadata**：细胞类型、修饰、递送条件必须记录，否则模型无法做条件设计。