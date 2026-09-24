从RMDB（RNA Mapping Database）获取数据，核心在于将其提供的**RNA结构映射实验的反应性数据**，转化为模型可用的**结构感知预训练信号**。RMDB与你之前了解的Ribonanza、OpenVaccine同属RNA结构探测数据，但RMDB是一个**经过人工审编的归档数据库**，数据来源更广、实验条件更丰富，适合作为结构感知模块的**高质量监督信号**和**跨条件泛化验证集**。

---

## 一、从RMDB获取数据

RMDB是RNA结构映射实验的归档数据库，收录了SHAPE、DMS、CMCT、1M7、mutate-and-map、M2-seq等技术的反应性剖面数据。截至2026年，RMDB包含**超过4,576,919条RNA序列**的化学映射数据，覆盖**1,360个实验**和**176,597个RNA构建体**。

### 1.1 数据组织与下载方式

RMDB的RDAT文件以**GitHub Release资产**形式分发，按五个子类别组织：

| Release | 内容 | 条目数 |
|---|---|---|
| `data-eterna` | Eterna和OpenKnot库实验 | 181 |
| `data-puzzle` | RNA Puzzles盲预测挑战 | 94 |
| `data-riboswitches` | 核糖开关家族（TPP、SAM等） | 253 |
| `data-rna-structures` | 核糖体/tRNA、病毒RNA（CoV-2、HIV等） | 275 |
| `data-general` | 其他异质数据 | 220 |

**下载方式一：GitHub CLI批量下载**

```bash
# 下载单个子类别
gh release download data-riboswitches \
  --repo DasLab/rmdb.github.io \
  --pattern "*.rdat" \
  --dir rdats/

# 下载全部五个子类别
for r in eterna puzzle riboswitches rna-structures general; do
  gh release download "data-$r" \
    --repo DasLab/rmdb.github.io \
    --pattern "*.rdat" \
    --dir rdats/
done
```

每个资产文件命名为`<RMDB_ID>.rdat`，例如`MDLOOP_STD_0000.rdat`。

**下载方式二：网页逐条下载**

访问任一Release页面，点击单个`.rdat`文件即可下载。每条RMDB条目的详情页都有“Download .rdat”按钮，提供永久直链。

**下载方式三：使用RNAGym预处理数据集（推荐用于模型训练）**

Marks-lab在HuggingFace上发布了RNAGym数据集，其中包含**已编译的RMDB化学映射数据**，存储为Parquet格式（brotli压缩），可直接加载训练。

| 文件 | 内容 |
|---|---|
| `compiled_RMDB.tar.gz` | 全部RMDB编译数据 |
| `RMDB_dataset_extra_clean.parquet` | 正常化学映射数据 |
| `RMDB_dataset_extra_cotrans.parquet` | 共转录折叠数据 |
| `RMDB_dataset_extra_degradation.parquet` | 降解测量数据 |
| `RMDB_dataset_extra_invivo.parquet` | 体内化学映射数据 |
| `RMDB_s40_train.parquet` / `RMDB_s40_test.parquet` | 40% identity聚类的train/test划分 |

**优先推荐使用RNAGym的Parquet格式**，因为RDAT原始格式需要专门的`rdat_kit`工具解析，而Parquet可直接用pandas加载。

### 1.2 RDAT文件格式

RDAT（RNA Data）是RMDB的标准格式，采用**基于注释的层级文本结构**，分为三个主要部分：

```text
RDAT_VERSION 0.32
NAME MedLoop
SEQUENCE GGAACGACGAACCGAAAACCGAAGAAAUGAAGAAAGGUUUUCGGUACCGACCUGAAAACCAAAGAAACAACAACAACAAC
STRUCTURE ..........((((((((((...............))))))))))...................................
OFFSET -10
SEQPOS G-9 G-8 A-7 A-6 C-5 G-4 A-3 C-2 G-1 A0 ...
ANNOTATION experimentType:StandardState chemical:Na-HEPES:50mM(pH8.0) temperature:24C
ANNOTATION_DATA:1 modifier:DMS
ANNOTATION_DATA:2 modifier:CMCT
REACTIVITY:1 161.1038 70.2383 75.5198 88.3231 ...
```

核心字段包括：**SEQUENCE**（RNA序列）、**STRUCTURE**（二级结构点括号表示）、**SEQPOS**（位置编号）、**ANNOTATION**（实验条件）、**REACTIVITY**（反应性值）。

RNAGym的Parquet格式将这些字段展平为列：

| 列名 | 说明 |
|---|---|
| `seqID` | 唯一数据ID |
| `sequence` | RNA序列 |
| `modifer` | 化学修饰剂（主要为DMS和2A3） |
| `SNR` | 信噪比 |
| `reads` | 测序reads总数 |
| `temperature` | 实验温度 |
| `chemical` | 缓冲液、盐、配体等 |
| `reverse_transcriptase` | 逆转录酶类型 |
| `reactivity` | 反应性值（字符串格式的数组） |
| `reactivity_err` | 反应性误差 |

**注意**：RMDB数据使用**公共领域CC0许可**，可自由用于商业和非商业用途。


## 二、数据清洗

### 2.1 信噪比过滤

与Ribonanza类似，**SNR < 1.0的序列应被过滤**。RNAGym的预处理流程明确将信噪比<1.0的序列移除。

```python
import pandas as pd
import numpy as np

def filter_low_snr(df, snr_col="SNR", threshold=1.0):
    """过滤低信噪比序列"""
    return df[df[snr_col] >= threshold].copy()
```

### 2.2 首尾无效区域裁剪

实验方法无法捕获序列5′和3′端的化学映射信号，这些位置的反应性为NaN。RNAGym的预处理流程会**裁剪掉序列首尾无化学映射数据的区域**。

```python
def trim_invalid_ends(sequence, reactivity_str):
    """裁剪首尾无反应性数据的区域"""
    react = np.array(eval(reactivity_str), dtype=np.float32)
    valid_mask = ~np.isnan(react)
    if not valid_mask.any():
        return None
    first_valid = np.argmax(valid_mask)
    last_valid = len(valid_mask) - np.argmax(valid_mask[::-1]) - 1
    return sequence[first_valid:last_valid+1], react[first_valid:last_valid+1]
```

### 2.3 缺失值处理

反应性中的NaN位置应被mask掉，不参与损失计算。

```python
def clean_reactivity(reactivity_str, max_val=None):
    """清洗反应性标签"""
    react = np.array(eval(reactivity_str), dtype=np.float32)
    mask = ~np.isnan(react)
    if max_val is not None:
        react = np.clip(react, 0.0, max_val)
    return react, mask
```

### 2.4 按实验条件分层

RMDB数据的一个独特优势是**实验条件丰富**——不同条目可能使用不同的修饰剂（DMS、2A3、CMCT、1M7等）、温度、缓冲液和折叠状态。RNAGym将数据细分为：

| 数据集 | 说明 |
|---|---|
| `extra_clean` | 正常化学映射数据 |
| `extra_cotrans` | 共转录折叠（未重新折叠） |
| `extra_degradation` | 降解测量（非反应性） |
| `extra_invivo` | 体内化学映射 |

清洗时应保留这些条件作为**模型的条件输入**，而非混合训练。

### 2.5 序列去重与聚类

RMDB数据来源多样，可能存在重复或高度相似序列。应进行**完全重复去重**，并使用**MMseqs2以40%序列一致性聚类**，用于后续防泄漏划分。


## 三、训练/测试数据集准备

### 3.1 样本结构

每条训练样本应包含：

```json
{
  "seqID": "MDLOOP_STD_0000",
  "sequence": "GGAACGACGAACCGAAAACCGAAG...",
  "modifier": "DMS",
  "SNR": 3.45,
  "reads": 1250,
  "temperature": 24.0,
  "chemical": "Na-HEPES:50mM(pH8.0)",
  "reverse_transcriptase": "SSII",
  "reactivity": [0.12, 0.45, NaN, ...],
  "reactivity_mask": [true, true, false, ...],
  "cluster_id": 42,
  "fold_type": "refolded"
}
```

### 3.2 数据划分：RNAGym的S40划分方案

RNAGym提供了**标准的40% identity聚类划分**，这是目前RNA结构预测领域广泛采用的做法。流程如下：

1. 从`RMDB_dataset_<integer>.parquet`（Ribonanza train/test）和`RMDB_dataset_extra_clean.parquet`收集所有序列。
2. 过滤SNR < 1.0的序列。
3. 裁剪首尾无数据的区域。
4. 将PDB三级结构序列纳入聚类。
5. 使用**MMseqs2以40%序列一致性聚类**。
6. 按**20%–80%** 划分train/test。

RNAGym的预划分文件可直接使用：`RMDB_s40_train.parquet`和`RMDB_s40_test.parquet`。

**如果需要自定义划分**，使用GroupShuffleSplit确保同一聚类簇不跨集合：

```python
from sklearn.model_selection import GroupShuffleSplit

def split_with_clusters(df, cluster_col="cluster_id",
                        test_size=0.20, val_size=0.15, random_state=42):
    """基于聚类的分层划分"""
    gss = GroupShuffleSplit(n_splits=1, test_size=test_size,
                            random_state=random_state)
    train_val_idx, test_idx = next(gss.split(df, groups=df[cluster_col]))

    train_val = df.iloc[train_val_idx]
    test = df.iloc[test_idx]

    val_ratio = val_size / (1 - test_size)
    gss2 = GroupShuffleSplit(n_splits=1, test_size=val_ratio,
                             random_state=random_state)
    train_idx, val_idx = next(gss2.split(train_val,
                                         groups=train_val[cluster_col]))

    train = train_val.iloc[train_idx]
    val = train_val.iloc[val_idx]
    return train, val, test
```

### 3.3 防泄漏检查

```python
def check_leakage(train_df, test_df, seq_col="sequence", cluster_col="cluster_id"):
    issues = []
    overlap = set(train_df[seq_col]) & set(test_df[seq_col])
    if overlap:
        issues.append(f"Exact sequence leakage: {len(overlap)}")
    if cluster_col in train_df.columns:
        overlap = set(train_df[cluster_col]) & set(test_df[cluster_col])
        if overlap:
            issues.append(f"Cluster leakage: {len(overlap)} clusters")
    return issues
```

### 3.4 在mRNA设计模型中的角色定位

**角色一：结构感知模块的预训练任务**

将`reactivity`作为回归目标，训练**核苷酸级反应性预测头**。RMDB数据的实验条件多样性使其比Ribonanza更适合训练**条件感知的结构预测器**——模型可以学习不同修饰剂、温度、缓冲液条件下的结构差异。

**角色二：与Ribonanza的互补使用**

| 维度 | Ribonanza | RMDB |
|---|---|---|
| 数据规模 | ~200万条序列 | ~457万条序列 |
| 实验条件 | 主要为DMS和2A3 | DMS、2A3、CMCT、1M7等多种 |
| 数据来源 | Eterna玩家设计为主 | 文献归档，来源广泛 |
| 格式 | CSV/Parquet | RDAT/Parquet |
| 审编程度 | 自动处理为主 | 人工审编 |

**建议**：以Ribonanza作为**大规模预训练数据**，以RMDB作为**精调数据和跨条件泛化验证集**。RNAGym的S40划分中，Ribonanza数据作为train/test的基础，RMDB的extra数据作为补充。

**角色三：多修饰剂条件建模**

RMDB包含DMS、2A3、CMCT、1M7等多种修饰剂的数据。可以将`modifier`作为条件输入，训练模型学习**不同化学探测方法下的结构反应性**。这直接对应你的模型需要的**修饰核苷酸条件控制**能力——虽然化学探测修饰与治疗性mRNA的核苷酸修饰（m1Ψ等）不同，但模型学到的“条件-结构”映射关系具有迁移价值。

**角色四：体内 vs 体外结构差异**

RMDB的`extra_invivo`数据集提供了**体内化学映射数据**，这比体外数据更接近真实细胞环境。可以用于验证模型在生理条件下的结构预测能力。


## 四、完整流程总结

```text
RMDB 数据获取
├── GitHub Release → 下载 .rdat 文件（五个子类别）
├── GitHub CLI → 批量下载
└── RNAGym HuggingFace → 直接使用预处理 Parquet（推荐）

        ↓

数据清洗
├── 信噪比过滤（SNR ≥ 1.0）
├── 首尾无效区域裁剪
├── NaN处理：mask，不参与损失
├── 完全重复序列去重
├── 按实验条件分层（clean/cotrans/degradation/invivo）
└── MMseqs2 40% identity 聚类

        ↓

训练/测试集准备
├── 样本：sequence + reactivity + modifier + SNR + conditions
├── 划分策略：RNAGym S40划分（20-80）或自定义聚类划分
├── 防泄漏检查：序列重复 + 聚类重叠
├── 角色定位：结构感知预训练 / 多修饰剂条件建模 / 体内外验证
└── 与Ribonanza互补：Ribonanza预训练，RMDB精调与验证

        ↓

输出：用于结构感知预训练和多条件泛化验证的辅助数据集
```


## 五、关键注意事项

1. **优先使用RNAGym的Parquet格式**：RDAT格式需要`rdat_kit`专门解析，RNAGym已将RMDB编译为可直接加载的Parquet文件，大幅降低工程成本。

2. **SNR过滤是必须的**：与Ribonanza相同，SNR < 1.0的序列测量质量低，应过滤。

3. **首尾区域需要裁剪**：实验方法无法捕获序列首尾的化学映射信号，应裁剪或mask。

4. **按实验条件分层**：RMDB包含共转录、降解、体内等不同类型的数据，不应混合训练。使用`fold_type`或数据集来源作为条件。

5. **40% identity聚类划分是标准做法**：RNAGym基准采用此方案，可直接使用其预划分文件，也可自行用MMseqs2聚类后划分。

6. **RMDB是结构数据，不是稳定性数据**：它是短RNA片段在体外化学探测条件下的反应性剖面，在mRNA设计模型中应定位为**结构感知模块的预训练信号**，而非完整mRNA半衰期的直接预测器。

7. **数据许可宽松**：RMDB数据为CC0公共领域，可自由使用和再分发，但使用时应引用原始作者和RMDB数据库论文。