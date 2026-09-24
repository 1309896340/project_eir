从Optimus 5-Prime获取数据来训练mRNA设计模型，核心在于将其提供的**5′UTR序列与翻译效率（MRL）的配对数据**，转化为模型可用的**序列-功能标签**，并严格防止数据泄漏。

---

## 一、从Optimus 5-Prime获取数据

Optimus 5-Prime是Sample等人于2019年发表在Nature Biotechnology上的工作，其数据集来自一项大规模并行报告基因分析（MPRA）：约280,000条随机50 nt的5′UTR序列被置于eGFP报告基因上游，在HEK293T细胞中通过多聚核糖体图谱（polysome profiling）测量每条序列的平均核糖体负载（Mean Ribosome Load, MRL）。

### 1.1 官方数据源

| 资源 | 地址 | 内容 |
|---|---|---|
| **GEO（原始数据）** | GSE114002 | 原始测序数据、处理后的CSV文件 |
| **GitHub（官方代码与数据）** | pjsample/human_5utr_modeling | 训练脚本、模型权重、数据加载代码 |
| **HuggingFace（非官方实现）** | multimolecule/optimus5prime | 模型卡片、使用示例 |

原始数据存储在GEO（GSE114002），处理后的CSV文件可直接从GEO下载。官方GitHub仓库的`human_5utrs/data/`目录下包含训练所需的核心数据文件。

### 1.2 核心数据文件

下载后，数据目录中通常包含以下CSV文件：

| 文件名 | 内容 | 用途 |
|---|---|---|
| `egfp_unmod_1.csv` | 未修饰mRNA的MRL测量值 | 主训练集 |
| `random_1.csv` | 随机50nt 5′UTR库的MRL | 补充训练集 |
| `varying_length_25to100.csv` | 可变长度5′UTR（25–100nt） | 长度泛化测试 |
| `human_5utrs.csv` | 天然人类5′UTR及其MRL | 自然序列验证 |
| `snv_*.csv` | 单核苷酸变体的MRL差异 | 变异效应预测 |

每条记录的核心字段是**5′UTR序列**和对应的**MRL值**。官方仓库同时提供了模型训练代码（`training_MRL_CNN.ipynb`）和序列设计代码（`evolve_for_target_MRL.ipynb`），可直接参考其数据加载和预处理逻辑。

### 1.3 批量获取建议

对于需要大规模数据的训练场景，建议：

1. 从GEO下载所有CSV文件，或使用`GEOquery`（R）直接读取GSE114002。
2. 克隆官方GitHub仓库，使用其`human_5utrs/data/`下的数据文件。
3. 若需要原始测序数据重新处理，从GEO下载SRA文件，用`sra-tools`转换后自行计算MRL。


## 二、数据清洗

### 2.1 序列格式标准化

- **字母表统一**：将序列统一为大写RNA字母表（`A/U/G/C`）。Optimus 5-Prime原始数据可能使用DNA字母表（`T`），需转为`U`。
- **长度过滤**：主训练集为固定50 nt序列，可变长度库为25–100 nt。对于mRNA设计模型，建议保留50 nt固定长度作为核心训练集，可变长度序列作为长度泛化验证集。
- **无效字符过滤**：去除含`N`或其他非标准碱基的序列。

### 2.2 测量质量过滤

MPRA数据存在测量噪声，需要过滤低质量测量：

| 过滤条件 | 标准 | 依据 |
|---|---|---|
| **最低reads数** | 每条序列总reads ≥ 200 | 官方重训练时使用的阈值 |
| **重复测量一致性** | 同一序列的重复测量CV < 阈值 | 排除测量不稳定序列 |
| **MRL值范围** | 去除极端离群值（如>3σ） | 避免异常值主导训练 |

### 2.3 序列冗余处理

随机50 nt库中可能存在重复或高度相似的序列：

```python
import pandas as pd
from collections import Counter

def deduplicate_utrs(df, seq_col="utr", mrl_col="mrl"):
    """去除重复序列，保留MRL均值"""
    # 完全重复序列：取MRL均值
    agg = df.groupby(seq_col).agg({
        mrl_col: "mean",
        "reads": "sum"
    }).reset_index()
    # 过滤低reads
    agg = agg[agg["reads"] >= 200]
    return agg
```

对于高度相似但不完全相同的序列（如仅1–2个碱基差异），如果都出现在训练集中，模型可能过拟合到细微差异。建议对训练集内部也做相似性聚类，避免同一簇的序列跨train/test划分。

### 2.4 天然人类5′UTR的验证

`human_5utrs.csv`中的天然序列可用于验证模型在自然序列上的表现。Sample等人测试了35,000条截短的人类5′UTR和3,577个自然变异，模型解释了超过80%的TE变异。清洗时应：

- 确认序列长度分布（天然5′UTR长度差异较大）。
- 对于模型固定50 nt输入的限制，截短或填充至50 nt，或使用支持可变长度的架构（如Framepool）。
- 保留天然序列的基因ID和变异信息，用于后续分析。


## 三、训练/测试数据集准备

### 3.1 训练样本结构

对于你的mRNA设计模型，Optimus 5-Prime数据主要用于训练**5′UTR翻译效率预测器**。每条样本应包含：

```json
{
  "utr5": "AACUCGCUGUAGUAAUUCCAGCGAGAG...",
  "mrl": 1.23,
  "mrl_std": 0.08,
  "reads": 450,
  "cell_type": "HEK293T",
  "modification": "unmodified",
  "library": "egfp_unmod_1",
  "length": 50
}
```

### 3.2 数据划分：防止泄漏是关键

**绝对不能用随机划分。** MPRA库中的序列可能高度相似，随机划分会导致训练集和测试集中出现近乎相同的序列，严重高估模型性能。

推荐划分策略（优先级从高到低）：

| 策略 | 方法 | 适用场景 |
|---|---|---|
| **按序列相似性聚类划分** | 使用CD-HIT或MMseqs2对5′UTR序列聚类，同一簇只出现在一个集合 | 最严格，防止同源泄漏 |
| **按GC含量分层划分** | 按GC含量分箱后分层抽样 | 保证GC分布一致 |
| **按MRL分箱划分** | 按MRL值分箱后分层抽样 | 保证MRL分布一致 |
| **按库来源划分** | train=egfp_unmod_1, test=random_1 | 测试跨库泛化 |

**最推荐：相似性聚类 + 分层抽样。** 先用CD-HIT以90% identity对5′UTR序列聚类，然后按MRL和GC分层，将同一簇的所有序列分配到同一集合。

```python
from sklearn.model_selection import StratifiedGroupKFold

def split_with_clusters(df, seq_col, mrl_col, n_splits=10):
    """基于相似性聚类的分层划分"""
    # 先用CD-HIT聚类（需先安装）
    # cd-hit-est -i utrs.fa -o clusters -c 0.9
    # 得到每条序列的cluster_id

    # 分层：按MRL分箱
    df["mrl_bin"] = pd.qcut(df[mrl_col], q=10, labels=False)

    # 按cluster分组划分
    skf = StratifiedGroupKFold(n_splits=n_splits)
    for train_idx, test_idx in skf.split(df, df["mrl_bin"], groups=df["cluster_id"]):
        # 得到train/test
        pass
```

### 3.3 防泄漏检查

```python
def check_leakage(train_df, test_df, seq_col="utr5"):
    issues = []
    # 完全重复序列
    overlap = set(train_df[seq_col]) & set(test_df[seq_col])
    if overlap:
        issues.append(f"Exact sequence leakage: {len(overlap)}")
    # 高度相似序列（需先聚类）
    train_clusters = set(train_df["cluster_id"])
    test_clusters = set(test_df["cluster_id"])
    overlap_clusters = train_clusters & test_clusters
    if overlap_clusters:
        issues.append(f"Cluster leakage: {len(overlap_clusters)} clusters")
    return issues
```

### 3.4 与mRNA设计模型的集成

Optimus 5-Prime数据在你的模型中的角色：

**角色一：训练5′UTR翻译效率预测头**

将MRL作为回归目标，训练一个预测器：

```python
# 输入：5′UTR序列 + 条件（细胞类型、修饰）
# 输出：预测MRL
# 损失：MSE
```

Sample等人发现，模型在uridine、pseudouridine和m1Ψ修饰条件下保持了稳健的排序一致性，这意味着该数据也可以用于训练**修饰核苷酸条件**下的翻译效率预测。

**角色二：5′UTR生成模型的奖励信号**

在条件生成模型中，Optimus 5-Prime预测器可作为奖励函数的一部分：

```python
def utr5_reward(utr5_sequence, cell_type="HEK293T", modification="unmodified"):
    """基于Optimus 5-Prime的5′UTR奖励"""
    predicted_mrl = optimus_predictor(utr5_sequence, cell_type, modification)
    return predicted_mrl  # 越高越好
```

**角色三：长度泛化验证**

`varying_length_25to100.csv`用于测试模型在不同5′UTR长度上的泛化能力。如果你的模型需要处理可变长度UTR，应使用Framepool等支持可变长度的架构，或在训练时对序列进行截短/填充。

### 3.5 与其他数据源的整合

Optimus 5-Prime只提供HEK293T细胞的未修饰mRNA数据。要训练条件控制模型，需要整合：

| 数据源 | 补充的条件 |
|---|---|
| **细胞类型** | GEO/SRA中其他细胞的MPRA数据 |
| **修饰核苷酸** | 官方论文中m1Ψ/Pseudouridine的条件数据 |
| **天然5′UTR** | `human_5utrs.csv`用于自然序列验证 |
| **变异效应** | `snv_*.csv`用于变异效应预测 |


## 四、完整流程总结

```text
Optimus 5-Prime 数据获取
├── GEO GSE114002 → 下载处理后的CSV文件
├── GitHub pjsample/human_5utr_modeling → 官方数据与代码
└── HuggingFace multimolecule/optimus5prime → 模型卡片与示例

        ↓

数据清洗
├── 序列标准化（U/T统一、长度过滤、无效字符过滤）
├── 测量质量过滤（reads ≥ 200、离群值去除）
├── 序列冗余处理（完全重复去重、相似性聚类）
└── 天然5′UTR验证集整理

        ↓

训练/测试集准备
├── 样本结构：utr5 + mrl + reads + cell_type + modification
├── 划分策略：相似性聚类 + 分层抽样（按MRL和GC）
├── 防泄漏检查：完全重复序列 + 聚类重叠
├── 角色定位：翻译效率预测头 / 生成奖励 / 长度泛化验证
└── 与其他数据源整合：细胞类型、修饰核苷酸

        ↓

输出：用于5′UTR翻译效率预测与生成的训练数据集
```


## 五、关键注意事项

1. **相似性聚类是必须的**：MPRA库中的序列可能存在大量近重复，随机划分会严重高估性能。使用CD-HIT或MMseqs2以90% identity聚类后划分。
2. **reads阈值不可忽视**：官方重训练时过滤了reads < 200的序列，这是保证测量可靠性的最低门槛。
3. **固定长度限制**：Optimus 5-Prime主模型只接受50 nt输入。如果你的模型需要处理可变长度UTR，使用Framepool架构或对序列做截短/填充。
4. **修饰核苷酸条件**：官方数据包含未修饰和m1Ψ/Pseudouridine的条件，可用于训练修饰条件预测，但修饰条件的样本量可能较少。
5. **天然序列与随机序列的分布差异**：随机50 nt序列的碱基组成与天然5′UTR差异显著。如果用随机序列训练，在天然序列上评估时需要注意分布偏移。
6. **版本记录**：记录所用的GEO accession、数据下载日期、官方仓库commit hash，保证可复现。
7. **只用于5′UTR**：Optimus 5-Prime只提供5′UTR数据，不包含CDS或3′UTR。CDS设计需从GENCODE/RefSeq/Ensembl获取，3′UTR稳定性数据需从其他来源（如RNA half-life数据集）补充。