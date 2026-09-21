从GEO/SRA获取数据来训练mRNA设计模型，核心在于将其提供的**原始测序数据**，转化为模型可用的**序列—功能标签**，并明确其在你模型中的定位：**mRNA半衰期（稳定性）和翻译效率（TE）标签的主要来源**。

---

## 一、从GEO/SRA获取数据

GEO（Gene Expression Omnibus）是NCBI的基因表达数据库，SRA（Sequence Read Archive）是其底层的原始测序数据存储库。

### 1.1 核心数据资源

| 数据类型 | 典型GEO编号 | 说明 |
|---|---|---|
| **mRNA半衰期（转录抑制）** | GSE208837 | HEK293T和C4-2细胞的mRNA半衰期，使用转录抑制+RNA-seq时间序列 |
| **mRNA半衰期（代谢标记）** | GSE34204 | 7个人类淋巴母细胞系，使用4sU标记计算半衰期 |
| **mRNA半衰期（RATE-seq）** | GSE212990 | iPSC来源神经元发育过程中的转录速率和半衰期变化 |
| **mRNA半衰期（有丝分裂）** | GSE314021 | HeLa细胞在有丝分裂期和间期的mRNA半衰期对比，使用转录抑制时间序列 |
| **Ribo-seq（翻译效率）** | GSE182100 | 大肠杆菌在12种条件下的翻译效率 |
| **Ribo-seq（翻译效率）** | GSE117299 | HeLa细胞的Ribo-seq数据 |
| **Ribo-seq（翻译效率）** | GSE121952 | HepG2细胞的Ribo-seq数据 |
| **Ribo-seq（翻译效率）** | GSE155447 | Huh7细胞的Ribo-seq数据 |
| **mRNA半衰期（人类）** | GSE126520 | 293T、HeLa、RPE细胞的内源mRNA衰减 |

### 1.2 下载流程

**Step 1：获取SRR编号**

在GEO页面找到目标GSE数据集，点击“SRA Run Selector”，下载`SRR_Acc_List.txt`文件，其中包含该数据集所有样本的SRR编号。

**Step 2：安装SRA Toolkit**

```bash
conda install -c bioconda sra-tools
```

**Step 3：下载并转换FASTQ**

```bash
# 下载SRA文件（注意：默认max-size为20GB，大文件需显式指定）
prefetch SRR12345678 --max-size 200G

# 验证下载完整性
vdb-validate SRR12345678

# 转换为FASTQ（双端测序用--split-files）
fasterq-dump --split-files -e 8 SRR12345678
gzip SRR12345678_*.fastq
```



**Step 4：批量下载**

```bash
while read srr; do
    prefetch "$srr" --max-size 200G
    fasterq-dump --split-files -e 8 "$srr"
    gzip "${srr}"_*.fastq
done < SRR_Acc_List.txt
```


## 二、数据清洗

### 2.1 原始数据质量控制

```bash
# 初步QC报告
fastqc SRR12345678_1.fastq.gz SRR12345678_2.fastq.gz

# 汇总多个样本的QC结果
multiqc .
```

FastQC评估的关键指标包括：**每个碱基的序列质量（Phred分数）、每个序列的质量分数、GC含量、N含量、序列长度分布、重复水平、过表达序列和接头含量**。

### 2.2 接头去除与质量过滤

使用`fastp`或`Trim Galore`进行接头去除和低质量碱基过滤：

```bash
# 使用fastp（推荐，速度快）
fastp -i SRR12345678_1.fastq.gz -I SRR12345678_2.fastq.gz \
      -o SRR12345678_1.clean.fastq.gz -O SRR12345678_2.clean.fastq.gz \
      --qualified_quality_phred 20 \
      --length_required 25 \
      --detect_adapter_for_pe \
      --thread 8

# 或使用Trim Galore
trim_galore --quality 20 --length 25 --paired \
    --output_dir clean/ \
    SRR12345678_1.fastq.gz SRR12345678_2.fastq.gz
```



标准过滤条件包括：去除接头序列、过滤Phred分数<20的低质量reads、去除长度<50bp的reads。

### 2.3 rRNA去除

对于RNA-seq数据，rRNA reads会占据大量测序资源，需通过比对去除：

```bash
# 使用bowtie2将reads比对到rRNA参考序列，去除比对的reads
bowtie2 -x rRNA_index -1 clean_1.fastq.gz -2 clean_2.fastq.gz \
    --un-conc-gz rRNA_removed_%.fastq.gz \
    -S /dev/null
```

也可使用SortMeRNA进行rRNA过滤。

### 2.4 比对与定量

将清洗后的reads比对到参考转录组（如GENCODE注释），然后进行定量：

```bash
# 使用Salmon进行转录本定量
salmon index -t transcripts.fa -i salmon_index

salmon quant -i salmon_index -l A \
    -1 clean_1.fastq.gz -2 clean_2.fastq.gz \
    -p 8 --validateMappings -o quant_output
```

**关键点**：使用与训练数据一致的参考转录组（如GENCODE），确保转录本ID与序列数据可以关联。

### 2.5 标签计算

**mRNA半衰期计算**

转录抑制时间序列的标准计算方法：

1. 将各时间点的mRNA水平归一化为RPKM/FPKM
2. 计算降解速率常数：`kdecay = mean(RPKM(t=0) / RPKM(t))`，其中t为各时间点
3. 半衰期：`t1/2 = ln2 / kdecay`

```python
import numpy as np

def compute_half_life(rpkm_time_series):
    """
    rpkm_time_series: dict, {time_hour: rpkm_value}
    返回: 半衰期（小时）
    """
    t0 = rpkm_time_series[0]
    ratios = [t0 / v for t, v in rpkm_time_series.items() if t > 0]
    kdecay = np.mean(ratios)
    if kdecay <= 0:
        return np.nan
    return np.log(2) / kdecay
```

**翻译效率（TE）计算**

TE需要Ribo-seq和RNA-seq的配对数据：

```python
def compute_te(ribo_fpkm, rna_fpkm):
    """
    TE = FPKM(Ribo-seq) / FPKM(RNA-seq)
    """
    if rna_fpkm == 0 or ribo_fpkm == 0:
        return np.nan
    return ribo_fpkm / rna_fpkm
```

TE定义为Ribo-seq的FPKM值与RNA-seq的FPKM值的比值，衡量单位转录本上结合的核糖体数目。

### 2.6 元数据标准化

GEO/SRA数据的元数据（细胞类型、组织、处理条件等）通常不规范，需要标准化：

| 字段 | 标准化方式 | 参考本体 |
|---|---|---|
| **细胞类型** | 模糊匹配到受控词表 | Cellosaurus、EFO |
| **组织** | 统一命名 | Uberon |
| **修饰条件** | 标准化为unmodified/m1Ψ/Ψ | 自定义受控词表 |
| **递送方式** | 标准化为LNP/electroporation等 | 自定义受控词表 |

可使用基于LLM的工具（如CistromeMeta）自动提取和标准化GEO元数据。


## 三、训练/测试数据集准备

### 3.1 样本结构

```json
{
  "geo_accession": "GSE208837",
  "srr_id": "SRR12345678",
  "gene_id": "ENSG00000141510",
  "transcript_id": "ENST00000269305",
  "species": "human",
  "cell_type": "HEK293T",
  "mrna_half_life": 3.2,
  "decay_rate": 0.217,
  "translation_efficiency": 0.85,
  "condition": {
    "modification": "unmodified",
    "delivery": "n/a",
    "treatment": "actinomycin_D",
    "time_points": [0, 2, 4, 6, 8, 12]
  },
  "metadata": {
    "batch": "GSE208837_batch1",
    "platform": "BGISEQ-500"
  }
}
```

### 3.2 数据划分：防止泄漏是关键

GEO/SRA数据来自不同研究、不同实验室，**绝对不能用随机划分**。

| 策略 | 方法 | 说明 |
|---|---|---|
| **按GSE数据集划分** | 同一GSE的所有样本只出现在一个集合 | 最推荐，防止同一研究内部泄漏 |
| **按细胞类型划分** | train=HEK293T+HeLa, test=CHO | 测试跨细胞类型泛化 |
| **按基因聚类划分** | 使用CD-HIT或MMseqs2对基因聚类，同一簇不跨集合 | 防止旁系同源基因泄漏 |
| **按物种划分** | train=human, test=mouse | 测试跨物种泛化 |

**关键原则**：同一GSE研究的样本必须保持在同一个集合中，否则模型可能学到该研究特有的批次效应而非真实的生物学信号。

```python
def split_by_gse(df, test_gses, val_gses):
    """按GSE数据集划分"""
    train = df[~df["geo_accession"].isin(test_gses + val_gses)]
    val = df[df["geo_accession"].isin(val_gses)]
    test = df[df["geo_accession"].isin(test_gses)]
    return train, val, test
```

### 3.3 防泄漏检查

```python
def check_leakage(train_df, test_df):
    issues = []
    # 同一GSE跨集合
    overlap = set(train_df["geo_accession"]) & set(test_df["geo_accession"])
    if overlap:
        issues.append(f"GSE leakage: {len(overlap)} studies")
    # 同一基因跨集合（若按基因划分则需检查）
    overlap = set(train_df["gene_id"]) & set(test_df["gene_id"])
    if overlap:
        issues.append(f"Gene leakage: {len(overlap)} genes")
    return issues
```

### 3.4 在mRNA设计模型中的角色定位

**角色一：mRNA半衰期预测器的主训练数据**

将计算得到的`mrna_half_life`作为回归目标，训练**转录本级稳定性预测器**：

```python
# 输入：mRNA序列 + 细胞类型条件
# 输出：预测半衰期（小时）
# 损失：MSE
```

这是你的模型中**稳定性预测头**的核心标签来源。与OpenVaccine/Ribonanza的局部降解倾向不同，GEO/SRA的half-life数据直接对应**完整mRNA在生理条件下的降解速率**。

**角色二：翻译效率预测器的训练数据**

将`translation_efficiency`作为回归目标，训练**TE预测器**，作为模型优化翻译效率的奖励信号。

**角色三：条件控制的关键标签来源**

GEO/SRA数据涵盖了**多种细胞类型**（HEK293T、HeLa、HepG2、CHO等）和**多种处理条件**（转录抑制、代谢标记等）。这些条件是模型**条件控制能力**的核心训练信号——模型需要学习“在HEK293T中稳定的序列，在T细胞中是否也稳定”。

### 3.5 关键局限性

- **GEO/SRA数据是天然转录本数据**，不是合成mRNA数据。天然mRNA的稳定性受**promoter、染色质状态、RNA结合蛋白、细胞状态**等影响，不能完全等同于外源mRNA的稳定性。
- **标签噪声较大**：不同实验室的实验条件（转录抑制方法、时间点选择、测序深度）差异大，需要严格的批次校正。
- **需要与合成mRNA数据校准**：建议用GEO/SRA数据训练基础稳定性预测器，再用自有的合成mRNA实验数据微调校准。


## 四、完整流程总结

```text
GEO/SRA 数据获取
├── 在GEO找到目标GSE → 下载SRR_Acc_List.txt
├── prefetch下载SRA文件 → fasterq-dump转FASTQ
└── 批量下载多个GSE数据集

        ↓

数据清洗
├── FastQC/MultiQC质控
├── fastp/Trim Galore接头去除和质量过滤
├── rRNA去除（bowtie2/SortMeRNA）
├── 比对和定量（Salmon/kallisto）
├── 半衰期计算（RPKM归一化 → kdecay → t1/2）
├── 翻译效率计算（FPKM(Ribo-seq)/FPKM(RNA-seq)）
└── 元数据标准化（细胞类型、处理条件）

        ↓

训练/测试集准备
├── 样本：gene_id + half_life + TE + cell_type + condition
├── 划分策略：按GSE数据集划分（最推荐）
├── 防泄漏检查：GSE重叠 + 基因重叠
├── 角色定位：半衰期预测主标签 / TE预测主标签 / 条件控制
└── 局限性声明：天然转录本数据 ≠ 合成mRNA稳定性

        ↓

输出：用于稳定性预测和翻译效率预测的核心训练数据集
```


## 五、关键注意事项

1. **按GSE划分是必须的**：同一GSE研究的样本必须保持在同一个集合中，否则模型会学到批次效应而非生物学信号。
2. **半衰期计算需要时间序列**：转录抑制或代谢标记的时间序列数据是计算半衰期的前提。单个时间点无法拟合降解速率。
3. **TE计算需要Ribo-seq和RNA-seq配对**：只有Ribo-seq数据不足以计算TE，必须同时有RNA-seq数据作为分母。
4. **元数据标准化不可忽视**：GEO数据的元数据质量参差不齐，必须标准化后才能作为条件特征使用。
5. **批次校正**：不同GSE研究之间存在显著的批次效应，需要使用ComBat、RUVg等方法校正。
6. **与合成mRNA数据校准**：GEO/SRA数据是天然转录本数据，其稳定性受细胞内调控网络影响。在mRNA设计模型中，应将其作为**基础稳定性预测器**的训练数据，再用自有的合成mRNA实验数据微调校准。
7. **数据许可**：GEO/SRA数据为公共数据，但部分人类数据可能受隐私保护限制，使用前需确认数据使用协议。