# 从GENCODE获取数据并准备mRNA设计训练集

## 一、GENCODE概述

GENCODE是国际合作项目，为人类和小鼠基因组提供高质量的结构注释，涵盖蛋白编码基因、非编码RNA、假基因和转录变体。它综合了HAVANA手动注释和Ensembl自动注释，是转录组研究和机器学习模型训练中最常用的注释资源。对于mRNA设计模型，GENCODE主要提供三类数据：

- **基因组序列**（FASTA）：用于提取实际碱基序列；
- **基因注释**（GTF/GFF3）：用于定位5'UTR、CDS、3'UTR等区域；
- **转录本序列**（FASTA）：直接提供成熟mRNA序列，但缺少区域边界注释。

---

## 二、下载所需数据文件

GENCODE数据托管在EBI FTP服务器上。当前人类最新版本为release 49，推荐使用**MANE Select**或**Basic**基因集来减少异构体冗余。

### 2.1 人类数据下载

```bash
# 创建数据目录
mkdir -p data/gencode_human && cd data/gencode_human

# 1. 基因组序列（primary assembly，不含alt/patches）
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# 2. 基因注释（GTF，包含CDS/UTR/exon等feature）
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gunzip gencode.v49.annotation.gtf.gz

# 3. 转录本序列（成熟mRNA，可选，用于快速获取序列）
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
gunzip gencode.v49.transcripts.fa.gz
```

参考。

### 2.2 小鼠数据下载

```bash
mkdir -p data/gencode_mouse && cd data/gencode_mouse

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz
```

### 2.3 文件选择建议

| 文件 | 用途 | 是否必需 |
|---|---|---|
| `GRCh38.primary_assembly.genome.fa.gz` | 提取任意区域序列 | 必需 |
| `gencode.vXX.annotation.gtf.gz` | 定位5'UTR/CDS/3'UTR坐标 | 必需 |
| `gencode.vXX.transcripts.fa.gz` | 快速获取成熟转录本序列 | 可选（GTF+基因组也可提取） |
| `gencode.vXX.pc_translations.fa.gz` | 蛋白序列，用于验证翻译 | 推荐 |

---

## 三、理解GTF格式

GTF文件是制表符分隔的9列格式，GENCODE的格式说明如下：

| 列号 | 内容 | GENCODE中的值 |
|---|---|---|
| 1 | 染色体 | chr1, chr2, ..., chrX, chrY, chrM |
| 2 | 注释来源 | ENSEMBL 或 HAVANA |
| 3 | 特征类型 | gene, transcript, exon, CDS, UTR, start_codon, stop_codon |
| 4 | 起始位置（1-based） | 整数 |
| 5 | 终止位置 | 整数 |
| 6 | 分数 | .（未使用） |
| 7 | 链方向 | + 或 - |
| 8 | 相位（仅CDS） | 0, 1, 2 |
| 9 | 键值对属性 | gene_id, transcript_id, gene_type 等 |

**关键属性字段**：

- `gene_type`：基因生物类型，如 `protein_coding`、`lncRNA`、`miRNA`
- `transcript_type`：转录本类型，如 `protein_coding`、`processed_transcript`
- `gene_id`：ENSG开头的基因ID
- `transcript_id`：ENST开头的转录本ID
- `exon_number`：外显子编号（从5'端计数）
- `level`：注释置信度，1=验证，2=手动注释，3=自动注释

对于mRNA设计模型，**只保留 `gene_type == "protein_coding"` 的转录本**。

---

## 四、从GTF提取区域坐标

mRNA设计模型需要三类区域：5'UTR、CDS、3'UTR。GENCODE GTF中CDS和UTR是分开注释的。

### 4.1 提取各区域BED文件

可以使用现成工具如 `gencode_regions`，从GTF提取3'UTR、5'UTR、CDS等区域并生成BED文件：

```bash
# 使用 gencode_regions 脚本
./create_regions_from_gencode.R gencode.v49.annotation.gtf ./regions_output
# 输出：exons.bed, 3UTR.bed, 5UTR.bed, genes.bed, cds.bed
```

也可以自己用Python解析：

```python
import pandas as pd

def load_gencode_gtf(gtf_path):
    """加载GENCODE GTF，解析第9列属性"""
    records = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = cols
            attr_dict = {}
            for item in attrs.strip().split(";"):
                item = item.strip()
                if not item:
                    continue
                key, _, value = item.partition(" ")
                attr_dict[key] = value.strip('"')
            records.append({
                "chrom": chrom,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "phase": phase,
                "gene_id": attr_dict.get("gene_id", ""),
                "transcript_id": attr_dict.get("transcript_id", ""),
                "gene_type": attr_dict.get("gene_type", ""),
                "transcript_type": attr_dict.get("transcript_type", ""),
                "exon_number": attr_dict.get("exon_number", ""),
            })
    return pd.DataFrame(records)

# 加载并过滤蛋白编码转录本
gtf = load_gencode_gtf("gencode.v49.annotation.gtf")
protein_coding = gtf[gtf["gene_type"] == "protein_coding"]
```

### 4.2 按转录本组装区域边界

对于每个转录本，需要确定：

- **CDS**：所有 `feature == "CDS"` 的外显子拼接
- **5'UTR**：CDS起始位置上游、同一转录本内的外显子区域
- **3'UTR**：CDS终止位置下游、同一转录本内的外显子区域

```python
def get_transcript_regions(gtf_df, transcript_id):
    """获取单个转录本的5'UTR、CDS、3'UTR坐标"""
    tx = gtf_df[gtf_df["transcript_id"] == transcript_id]
    cds = tx[tx["feature"] == "CDS"]
    exons = tx[tx["feature"] == "exon"]

    if cds.empty:
        return None

    strand = tx["strand"].iloc[0]
    cds_start = cds["start"].min()
    cds_end = cds["end"].max()

    # 5'UTR：CDS起始之前的外显子
    # 3'UTR：CDS终止之后的外显子
    if strand == "+":
        utr5_exons = exons[exons["end"] < cds_start]
        utr3_exons = exons[exons["start"] > cds_end]
    else:
        utr5_exons = exons[exons["start"] > cds_end]
        utr3_exons = exons[exons["end"] < cds_start]

    return {
        "utr5": utr5_exons,
        "cds": cds,
        "utr3": utr3_exons,
        "strand": strand
    }
```

---

## 五、提取序列

### 5.1 方法一：使用 gffread（推荐）

`gffread` 是最方便的工具，可以直接从GTF和基因组FASTA提取各区域序列：

```bash
# 提取CDS序列
gffread -x cds.fa -g GRCh38.primary_assembly.genome.fa gencode.v49.annotation.gtf

# 提取成熟转录本序列（拼接外显子）
gffread -w transcripts.fa -g GRCh38.primary_assembly.genome.fa gencode.v49.annotation.gtf

# 提取5'UTR序列
gffread -w utr5.fa -g GRCh38.primary_assembly.genome.fa \
  --region-type 5UTR gencode.v49.annotation.gtf
```



### 5.2 方法二：使用 pyGTF

```python
from pyGTF import GTF
gtf = GTF("gencode.v49.annotation.gtf")
# 提取指定转录本的CDS、UTR序列
seq = gtf.extract_sequence("ENST00000456328.2", "GRCh38.primary_assembly.genome.fa")
```



### 5.3 方法三：从转录本FASTA直接获取

如果只需要成熟mRNA序列，直接使用 `gencode.vXX.transcripts.fa`。FASTA头信息包含转录本ID和基因类型：

```
>ENST00000456328.2|ENSG00000223972.5|...|DDX11L1-202|DDX11L1|1657|processed_transcript|
```

可以提取序列和元数据：

```python
from Bio import SeqIO

records = []
for rec in SeqIO.parse("gencode.v49.transcripts.fa", "fasta"):
    parts = rec.id.split("|")
    tx_id = parts[0]
    gene_id = parts[1]
    gene_name = parts[5]
    tx_type = parts[7]
    length = len(rec.seq)
    records.append({
        "transcript_id": tx_id,
        "gene_id": gene_id,
        "gene_name": gene_name,
        "transcript_type": tx_type,
        "length": length,
        "sequence": str(rec.seq).replace("T", "U")
    })
```

---

## 六、数据清洗与过滤

### 6.1 转录本筛选标准

| 过滤条件 | 标准 | 理由 |
|---|---|---|
| gene_type | `protein_coding` | 只保留编码基因 |
| transcript_type | `protein_coding` | 排除保留内含子等异常转录本 |
| 长度 | 300 nt – 15,000 nt | 排除过短/过长异常 |
| CDS长度 | 是3的倍数 | 保证完整密码子 |
| 起始密码子 | 以ATG开头 | 标准翻译起始 |
| 终止密码子 | 以TAA/TAG/TGA结尾 | 保证正常终止 |
| 内部终止密码子 | 无 | 避免截短蛋白 |
| 未知碱基 | N比例 < 1% | 保证序列质量 |
| transcript_support_level | 1或2 | 提高注释可信度 |



### 6.2 代表性转录本选择

一个基因常有多个异构体。常见策略：

1. **MANE Select**：GENCODE推荐的每个基因一个代表转录本，最简洁
2. **最长CDS**：选择CDS最长的异构体，覆盖最全功能域
3. **最长转录本**：如果CDS长度相同，选最长转录本
4. **最短UTR**：在CDS相同情况下，优先选最短3'UTR和5'UTR的转录本

推荐优先使用MANE Select：

```python
# 从GTF中筛选MANE Select转录本
mane = gtf[gtf["tag"].str.contains("MANE_Select", na=False)]
```

### 6.3 多异构体处理示例

```python
def select_representative(gtf_df, gene_id):
    """为一个基因选择代表性转录本"""
    gene_tx = gtf_df[
        (gtf_df["gene_id"] == gene_id) &
        (gtf_df["feature"] == "CDS")
    ]
    if gene_tx.empty:
        return None

    # 计算每个转录本的CDS总长度
    cds_lengths = gene_tx.groupby("transcript_id").apply(
        lambda x: (x["end"] - x["start"] + 1).sum()
    ).reset_index(name="cds_length")

    # 选CDS最长的；若相同，选转录本最长的
    best = cds_lengths.sort_values(
        ["cds_length"], ascending=False
    ).iloc[0]
    return best["transcript_id"]
```

---

## 七、构建训练样本

### 7.1 样本结构

每条训练样本应包含：

```json
{
  "transcript_id": "ENST00000456328.2",
  "gene_id": "ENSG00000223972.5",
  "gene_name": "DDX11L1",
  "species": "human",
  "chrom": "chr1",
  "strand": "+",
  "utr5": "GGGAAACCC...",
  "cds": "AUGGCU...",
  "utr3": "GGGCCC...",
  "protein_sequence": "MALW...",
  "full_mrna": "GGGAAACCCAUGGCU...GGGCCC",
  "gc_content": 0.52,
  "cds_length": 1200,
  "utr5_length": 150,
  "utr3_length": 800
}
```

### 7.2 序列提取与组装

```python
from Bio.Seq import Seq
from Bio import SeqIO

def build_mrna_sample(tx_id, gtf_df, genome_dict):
    """从GTF和基因组构建完整mRNA样本"""
    tx = gtf_df[gtf_df["transcript_id"] == tx_id]
    strand = tx["strand"].iloc[0]
    chrom = tx["chrom"].iloc[0]

    # 按坐标排序外显子
    exons = tx[tx["feature"] == "exon"].sort_values("start")

    # 拼接外显子序列
    full_seq = ""
    for _, exon in exons.iterrows():
        seq = genome_dict[chrom][exon["start"]-1:exon["end"]]
        full_seq += seq

    if strand == "-":
        full_seq = str(Seq(full_seq).reverse_complement())

    # 确定CDS边界
    cds = tx[tx["feature"] == "CDS"]
    cds_start = cds["start"].min()
    cds_end = cds["end"].max()

    # 从完整序列中定位CDS（需要处理负链坐标转换）
    # ... 具体实现略

    return {
        "transcript_id": tx_id,
        "full_mrna": full_seq.replace("T", "U"),
        "strand": strand,
        "chrom": chrom
    }
```

### 7.3 翻译验证

提取CDS后，必须验证其翻译结果与GENCODE提供的蛋白序列一致：

```python
from Bio.Seq import Seq

def validate_cds(cds_seq, expected_protein):
    """验证CDS翻译是否与预期蛋白一致"""
    # 确保以ATG开头、终止密码子结尾
    if not cds_seq.startswith("ATG"):
        return False
    if cds_seq[-3:] not in ["TAA", "TAG", "TGA"]:
        return False

    # 翻译
    translated = str(Seq(cds_seq).translate(to_stop=True))
    return translated == expected_protein
```

GENCODE提供 `gencode.vXX.pc_translations.fa.gz`，可用于验证。

---

## 八、负样本与非编码数据

### 8.1 为什么要负样本

mRNA设计模型需要区分“可翻译的CDS”和“非编码/不可翻译序列”。GENCODE包含大量lncRNA、假基因等非编码转录本，可作为负样本。

### 8.2 负样本来源

| 来源 | 类型 | 用途 |
|---|---|---|
| GENCODE lncRNA | `gene_type == "lncRNA"` | 训练编码/非编码分类 |
| GENCODE pseudogene | `gene_type == "pseudogene"` | 负样本 |
| 随机同义CDS | 随机生成同义密码子 | 训练ranking模型 |
| 内部ORF | CDS框外ORF | 负样本 |
| uORF | 5'UTR中的ORF | 调控研究 |

### 8.3 构建对比样本

对于同一个蛋白，可以生成多个同义CDS变体：

```python
import random
from Bio.Seq import Seq

# 人类密码子使用表（示例，实际应统计）
codon_table = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    # ... 其他氨基酸
}

def generate_synonymous_variant(protein_seq, codon_table, strategy="random"):
    """生成同义CDS变体"""
    cds = ""
    for aa in protein_seq:
        codons = codon_table.get(aa, [])
        if not codons:
            return None
        if strategy == "random":
            cds += random.choice(codons)
        elif strategy == "most_frequent":
            cds += codons[0]
    return cds
```

---

## 九、特征计算

为每条序列计算以下特征，作为模型的辅助输入或约束：

| 特征 | 说明 | 工具 |
|---|---|---|
| GC含量 | 全局和局部GC | 自写脚本 |
| MFE | 最小自由能 | ViennaRNA / RNAfold |
| 5'端结构强度 | 起始密码子附近ΔG | RNAfold |
| CAI | 密码子适应指数 | CAI工具 |
| 密码子对偏好 | Codon pair bias | CoCoPUTs |
| 同聚物长度 | 最长连续相同碱基 | 自写脚本 |
| 重复序列 | 重复评分 | 自写脚本 |
| motif | AU-rich、剪接样信号 | 正则匹配 |

```bash
# 使用 ViennaRNA 计算 MFE
RNAfold --noPS < sequence.fa > structure.txt
```

---

## 十、数据划分

### 10.1 按染色体划分（推荐）

避免同一基因的变体跨集合泄漏：

```python
train_chroms = ["chr1", "chr5", "chr7", "chr10", "chr13", "chr17", "chr21"]
valid_chroms = ["chr2", "chr9", "chr16"]
test_chroms = ["chr3", "chr8", "chr15"]

train = data[data["chrom"].isin(train_chroms)]
valid = data[data["chrom"].isin(valid_chroms)]
test = data[data["chrom"].isin(test_chroms)]
```

### 10.2 其他划分策略

| 策略 | 用途 | 说明 |
|---|---|---|
| 按基因划分 | 测试新蛋白泛化 | 同一基因不同异构体不跨集合 |
| 按实验批次划分 | 测试跨批次泛化 | 避免批次效应泄漏 |
| 按细胞类型划分 | 测试跨细胞泛化 | 如train=HEK293T，test=T细胞 |
| 按物种划分 | 测试跨物种泛化 | train=human，test=mouse |

---

## 十一、完整数据准备流程

```text
Step 1: 下载 GENCODE 文件
        ├── GRCh38.primary_assembly.genome.fa
        ├── gencode.v49.annotation.gtf
        └── gencode.v49.transcripts.fa

Step 2: 解析 GTF
        ├── 过滤 protein_coding
        ├── 提取 CDS / 5'UTR / 3'UTR 坐标
        └── 选择代表性转录本（MANE Select 或最长CDS）

Step 3: 提取序列
        ├── 使用 gffread 或 pyGTF
        ├── 拼接外显子
        └── 反转录负链

Step 4: 清洗
        ├── 检查起始/终止密码子
        ├── 检查CDS长度是否为3的倍数
        ├── 翻译验证
        └── 过滤低质量序列

Step 5: 构建样本
        ├── 组装 5'UTR + CDS + 3'UTR
        ├── 记录蛋白序列
        └── 添加负样本（lncRNA、随机同义变体）

Step 6: 特征计算
        ├── GC、MFE、CAI、密码子对偏好
        ├── 同聚物、重复、motif
        └── 结构特征

Step 7: 数据划分
        ├── 按染色体划分 train/valid/test
        └── 保存为 JSON/Parquet/FASTA
```

---

## 十二、输出文件结构

```text
data/
├── raw/
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── gencode.v49.annotation.gtf
│   └── gencode.v49.transcripts.fa
├── processed/
│   ├── transcripts_metadata.parquet
│   ├── utr5_sequences.fa
│   ├── cds_sequences.fa
│   ├── utr3_sequences.fa
│   ├── full_mrna_sequences.fa
│   └── protein_sequences.fa
├── features/
│   ├── sequence_features.parquet
│   └── structure_features.parquet
├── splits/
│   ├── train.json
│   ├── valid.json
│   └── test.json
└── negatives/
    ├── lncrna.fa
    └── random_synonymous.fa
```

---

## 十三、注意事项

1. **版本一致性**：基因组FASTA和GTF必须来自同一GENCODE版本，否则坐标不匹配。
2. **链方向**：负链基因需要反转录，CDS边界在负链上的位置与正链相反。
3. **终止密码子**：GENCODE的CDS注释通常**不包括终止密码子**，需要手动添加或从序列中确认。
4. **T/U表示**：训练时统一用RNA字母表（U代替T），输出时保持一致。
5. **数据泄漏**：同一基因的不同异构体不能跨train/test划分。
6. **MANE Select**：如果不需要完整异构体，直接使用MANE Select可大幅简化流程。
7. **验证翻译**：提取CDS后必须与GENCODE蛋白序列比对，确保无误。