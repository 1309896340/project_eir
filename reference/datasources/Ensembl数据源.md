从Ensembl为mRNA设计模型获取数据，核心流程与RefSeq类似，都遵循 **“下载 → 提取 → 清洗 → 划分”** 的范式。但Ensembl在数据组织、访问方式和注释体系上有其特点，需要针对性地调整策略。

---

## 一、从Ensembl获取数据

Ensembl提供三种主要的数据获取方式，各有适用场景：

### 1.1 方式一：BioMart（推荐用于批量提取CDS/UTR序列）

BioMart是Ensembl最灵活的数据检索工具，允许用户通过图形界面或编程接口（R的`biomaRt`包、Python的`pybiomart`）批量下载特定类型的序列。

**核心属性（Attributes）**：
- `cdna`：拼接后的转录本序列，包含UTR（即5′UTR + CDS + 3′UTR）
- `coding`：CDS序列，不含UTR
- `5utr`：5′UTR序列
- `3utr`：3′UTR序列
- `protein`：蛋白序列

**关键提示**：一次查询只能选择一个序列类型，需要分别运行三次查询来获取CDS、5′UTR和3′UTR。建议在属性中同时包含`ensembl_transcript_id`，以便将三条序列正确关联到同一转录本。

**R示例（biomaRt）**：

```r
library(biomaRt)

ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "uswest")

# 分别查询三种序列
cds <- getBM(attributes = c("ensembl_transcript_id", "coding"),
             filters = "ensembl_gene_id",
             values = gene_list, mart = ensembl)

utr5 <- getBM(attributes = c("ensembl_transcript_id", "5utr"),
              filters = "ensembl_gene_id",
              values = gene_list, mart = ensembl)

utr3 <- getBM(attributes = c("ensembl_transcript_id", "3utr"),
              filters = "ensembl_gene_id",
              values = gene_list, mart = ensembl)
```

### 1.2 方式二：REST API（推荐用于程序化按需获取）

Ensembl REST API的`POST /sequence/id`端点支持**批量**获取序列，一次请求最多可提交50个ID。通过`type`参数指定序列类型：

| type参数 | 含义 |
|---|---|
| `cdna` | 拼接转录本 + UTR |
| `cds` | CDS，不含UTR |
| `protein` | 蛋白序列 |

**Python示例**：

```python
import requests

SERVER = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}

def fetch_sequences(ids, seq_type="cds"):
    """批量获取序列，seq_type: cds / cdna / protein"""
    ext = f"/sequence/id?type={seq_type}"
    r = requests.post(SERVER + ext, headers=HEADERS,
                      json={"ids": ids})
    return r.json()
```

对于5′UTR和3′UTR的**单独**获取，REST API不直接支持按`type=5utr`查询。此时需要用`type=cdna`获取完整转录本，再结合GTF中的CDS坐标来切分UTR区域。

### 1.3 方式三：FTP下载（推荐用于全量数据）

Ensembl FTP服务器（`ftp.ensembl.org/pub`）提供按物种组织的数据文件。对于mRNA设计模型，需要下载：

| 文件类型 | 用途 | 示例路径 |
|---|---|---|
| **GTF/GFF3注释** | 定位CDS、UTR、外显子边界 | `current/gtf/homo_sapiens/` |
| **CDS FASTA** | 所有转录本的CDS序列 | `current/fasta/homo_sapiens/cds/` |
| **cDNA FASTA** | 拼接转录本（含UTR） | `current/fasta/homo_sapiens/cdna/` |
| **蛋白 FASTA** | 翻译验证的参考 | `current/fasta/homo_sapiens/pep/` |

**批量下载建议**：Ensembl FTP支持`rsync`，适合大规模数据同步。对于多物种训练（如GEMORNA使用了115种哺乳动物的Ensembl数据），FTP是最可行的方式。

### 1.4 物种与转录本选择

| 策略 | 说明 |
|---|---|
| **物种范围** | 优先human、mouse、rat；可扩展至哺乳动物（GEMORNA用了115种） |
| **代表性转录本** | Ensembl提供**Canonical transcript**标签，每个基因一个代表转录本 |
| **过滤条件** | 只保留`biotype == "protein_coding"`的转录本 |


## 二、数据清洗

### 2.1 格式与方向标准化

- **U/T统一**：Ensembl输出为DNA字母表（`A/T/G/C`），需转为RNA（`A/U/G/C`）用于模型训练。
- **方向标准化**：Ensembl的cDNA和CDS序列**已经是5′→3′方向**（负链已反向互补），无需额外处理。这与从GTF+基因组手动提取时不同。
- **字符集检查**：Ensembl序列中可能包含`N`（用于保持阅读框相位，当转录本起始有相位偏移时）。`N`比例 > 1%的序列应过滤。

### 2.2 CDS结构检查（硬性门槛）

| 检查项 | 标准 | 处理 |
|---|---|---|
| CDS长度 | 是3的倍数 | 不满足则过滤 |
| 起始密码子 | `ATG` | 不满足则过滤或标记 |
| 终止密码子 | `TAA/TAG/TGA` | 不满足则过滤 |
| 内部终止密码子 | 无 | 有则过滤 |
| 翻译一致性 | translate(CDS) == Ensembl蛋白 | 不一致则过滤 |

Ensembl的`ensembldb`包提供了`cds_ok`字段，自动检查CDS长度是否与氨基酸序列长度匹配。在自定义流程中，需要手动实现这一验证。

### 2.3 翻译一致性验证（最核心QC）

**将CDS翻译后，必须与Ensembl提供的蛋白序列完全一致。** GEMORNA的验证标准是：**CDS长度必须恰好是蛋白长度的3倍加3（终止密码子）**，不满足的序列被过滤。

```python
from Bio.Seq import Seq

def validate_ensembl_cds(cds_rna, ensembl_protein):
    """验证Ensembl CDS翻译是否与蛋白一致"""
    dna = cds_rna.replace("U", "T")
    # 去除可能的N（Ensembl用N保持相位）
    dna_clean = dna.replace("N", "")
    if len(dna_clean) % 3 != 0:
        return False, "CDS长度非3的倍数"
    translated = str(Seq(dna_clean).translate(to_stop=True))
    expected = ensembl_protein.rstrip("*")
    if translated.rstrip("*") == expected:
        return True, "一致"
    return False, f"翻译不一致: {translated[:50]} vs {expected[:50]}"
```

**注意**：Ensembl的某些CDS以`N`开头（相位偏移），翻译时需要正确处理。如果CDS以`N`开头且相位为2，需要跳过前2个`N`再翻译。

### 2.4 序列质量过滤

| 检查项 | 标准 | 依据 |
|---|---|---|
| 全长长度 | 100–15,000 nt | 排除过短/过长 |
| CDS长度 | 90–10,000 nt | 保证完整功能域 |
| GC含量 | 30%–70% | 排除极端GC |
| 最长同聚物 | ≤ 8 | 避免生产困难 |
| N比例 | < 1% | 保证序列质量 |
| 重复评分 | 低重复 | 避免合成失败 |

### 2.5 负样本构建

| 来源 | 类型 | 用途 |
|---|---|---|
| Ensembl `biotype != protein_coding` | lncRNA、假基因 | 编码/非编码判别 |
| 随机同义CDS | 同一蛋白随机密码子 | ranking训练 |
| 框外ORF | CDS框外开放阅读框 | 负样本 |
| uORF | 5′UTR中的上游ORF | 调控研究 |


## 三、训练/测试数据集准备

### 3.1 样本结构

```json
{
  "ensembl_transcript_id": "ENST00000269305",
  "ensembl_gene_id": "ENSG00000141510",
  "gene_symbol": "TP53",
  "species": "human",
  "utr5": "GCCACC...",
  "cds": "AUGGAGGAG...",
  "utr3": "UGUCUG...",
  "protein_sequence": "MEEPQSD...",
  "full_mrna": "GCCACC...AUG...UGA...UGUCUG",
  "gc_content": 0.52,
  "cds_length": 1182,
  "is_canonical": true
}
```

### 3.2 数据划分：防止泄漏是关键

**绝对不能用随机划分。** 同一基因的不同异构体、高度相似的旁系同源基因如果跨train/test，会造成严重的数据泄漏。

| 策略 | 方法 | 适用场景 |
|---|---|---|
| **按基因划分** | 同一基因的所有异构体只出现在一个集合 | 最推荐 |
| **按染色体划分** | train/valid/test使用不同染色体 | 基础泛化评估 |
| **按物种划分** | train=human，test=mouse | 跨物种泛化 |
| **相似性聚类划分** | 使用SpanSeq等工具按序列相似性聚类后划分 | 最严格 |

**SpanSeq**是专门为生物序列设计的数据划分工具，通过相似性聚类避免train/test之间的序列泄漏。GEMORNA的数据划分标准是：将蛋白序列聚类到50% identity，然后按80%/10%/10%划分。

```python
# 使用SpanSeq进行相似性聚类划分
# pip install spanseq
from spanseq import SpanSeq

splitter = SpanSeq(identity_threshold=0.5,  # 50% identity
                   n_folds=10)
splits = splitter.split(sequences, labels)
```

**防泄漏检查清单**：

```python
def check_leakage(train_df, test_df):
    issues = []
    for col in ["ensembl_gene_id", "protein_sequence",
                "ensembl_transcript_id"]:
        overlap = set(train_df[col]) & set(test_df[col])
        if overlap:
            issues.append(f"{col} leakage: {len(overlap)}")
    return issues
```

### 3.3 各区域分开准备

| 数据集 | 输入 | 输出/标签 | 关键过滤 |
|---|---|---|---|
| **CDS数据集** | 蛋白序列 | 同义CDS | 翻译一致性、无内部stop |
| **5′UTR数据集** | UTR序列 | 翻译效率 | 长度 < 500 nt |
| **3′UTR数据集** | UTR序列 | 稳定性 | 长度 < 2,000 nt |

GEMORNA的微调策略值得参考：预训练用全量Ensembl数据（115种哺乳动物），微调时用PRED-5UTR和PRED-3UTR模型筛选高MRL的5′UTR和高稳定性的3′UTR，随机采样得到约80万条5′UTR和20万条3′UTR微调数据。

### 3.4 条件标签的整合

你的模型需要条件控制（细胞类型、修饰核苷酸、递送方式），这些标签需要从实验数据中整合：

| 条件 | 数据来源 | 处理方式 |
|---|---|---|
| 细胞类型 | GEO/SRA元数据 | 标准化为受控词表 |
| 修饰核苷酸 | 实验记录 | unmodified / m1Ψ / Ψ |
| 递送方式 | 实验记录 | LNP / electroporation |
| 表达量 | reporter assay | 相对对照、z-score |
| 半衰期 | RNA-seq time course | 拟合decay rate |
| 翻译效率 | Ribo-seq | Ribo-seq RPKM / RNA-seq RPKM |


## 四、Ensembl与RefSeq的关键差异

| 维度 | Ensembl | RefSeq |
|---|---|---|
| **转录本选择** | Canonical transcript标签 | RefSeq Select / MANE Select |
| **序列方向** | cDNA/CDS已标准化为5′→3′ | 同样已标准化 |
| **N碱基** | CDS中可能含N保持相位 | 较少出现 |
| **BioMart** | 原生支持，属性丰富 | 不支持，需用Datasets CLI |
| **REST API** | 支持批量POST | 支持但端点不同 |
| **多物种覆盖** | 更广（115+哺乳动物） | 以模式生物为主 |
| **数据版本** | 定期release（如release-112） | 与基因组组装版本绑定 |


## 五、完整流程总结

```text
Ensembl 数据获取
├── BioMart → 批量提取CDS/5'UTR/3'UTR（R: biomaRt）
├── REST API → 程序化按需获取（POST /sequence/id，批量50 ID）
└── FTP → 全量下载GTF + cDNA + CDS + 蛋白FASTA

        ↓

数据清洗
├── U/T统一、方向确认（Ensembl已标准化）
├── N碱基处理（CDS中可能含N保持相位）
├── CDS结构检查（长度3的倍数、起始/终止/内部终止）
├── 翻译一致性验证（translate(CDS) == Ensembl蛋白）
├── 序列质量过滤（GC、同聚物、N比例）
└── 负样本构建（lncRNA、随机同义、框外ORF）

        ↓

训练/测试集准备
├── 样本结构：utr5 + cds + utr3 + protein + metadata
├── 划分策略：按基因划分 / SpanSeq相似性聚类划分
├── 防泄漏检查：基因/蛋白/转录本/批次
├── 区域分离：CDS、5′UTR、3′UTR 各自成集
└── 条件标签整合：细胞类型、修饰、递送

        ↓

输出：可直接用于Transformer训练的mRNA设计数据集
```


## 六、关键注意事项

1. **Ensembl序列已标准化方向**：cDNA和CDS都是5′→3′，负链已处理，无需手动反向互补。
2. **N碱基需要特殊处理**：Ensembl在CDS中使用`N`保持翻译相位，翻译前需正确处理。
3. **BioMart分三次查询**：CDS、5′UTR、3′UTR需要分别查询，用transcript_id关联。
4. **翻译一致性是硬门槛**：不通过则直接过滤，不可妥协。
5. **SpanSeq是最严格的划分方案**：按序列相似性聚类后划分，最大程度避免泄漏。
6. **保留metadata**：细胞类型、修饰、递送条件必须记录，否则模型无法做条件设计。
7. **版本一致性**：GTF、cDNA、CDS、蛋白FASTA必须来自同一Ensembl release。