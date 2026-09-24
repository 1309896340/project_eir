从Ribonanza获取数据来训练mRNA设计模型，核心在于将其提供的**约200万条RNA序列的化学映射（chemical mapping）数据**，转化为模型可用的**序列-结构反应性预测任务**。与OpenVaccine类似，Ribonanza在你的模型中应定位为**RNA二级结构感知与局部稳定性/降解倾向的预训练信号**，而非完整mRNA半衰期的直接标签。

---

## 一、从Ribonanza获取数据

Ribonanza是Stanford Das Lab与Eterna平台合作，于2023年在Kaggle上发起的“Stanford Ribonanza RNA Folding”竞赛数据集。它通过**DMS和2A3两种化学探测方法**，对约200万条多样化RNA序列进行了高通量化学映射测量，提供了核苷酸分辨率的反应性（reactivity）和降解（degradation）剖面。

### 1.1 核心数据资源

| 资源 | 地址 | 内容 |
|---|---|---|
| **Kaggle竞赛（原始数据）** | kaggle.com/competitions/stanford-ribonanza-rna-folding/data | `train_data.csv`、`test_sequences.csv`、`sample_submission.csv` |
| **官方GitHub（数据准备）** | DasLab/DataPrepRibonanzaKaggle2023 | MATLAB脚本、示例数据、数据准备流程 |
| **RibonanzaNet（模型与数据）** | Shujun-He/RibonanzaNet | 训练代码、预训练权重、数据加载示例 |
| **HuggingFace镜像** | multimolecule/ribonanzanet-* | 部分预处理版本和模型卡片 |

官方RibonanzaNet仓库明确指出，训练所需的核心文件是来自Kaggle竞赛页面的 `train_data.csv`、`test_sequences.csv` 和 `sample_submission.csv`。

### 1.2 核心字段结构

`train_data.csv` 的每条记录包含以下关键字段：

| 字段 | 类型 | 说明 |
|---|---|---|
| `sequence` | 字符串 | RNA序列（A/G/U/C），长度范围50–240 nt |
| `experiment_type` | 字符串 | 化学探测方法：`DMS` 或 `2A3` |
| `SN_filter` | 浮点数 | 信噪比过滤器，>1.0表示高质量测量 |
| `reactivity_0001` … `reactivity_XXXX` | 浮点数组 | 每个核苷酸位置的反应性值 |
| `reactivity_error_*` | 浮点数组 | 反应性测量的误差 |
| `deg_*` | 浮点数组 | 降解相关标签 |

**关键点**：序列的**5′端和3′端部分位置可能没有反应性数据**（返回NaN），因为当前实验方法无法捕获序列首尾的化学映射信号。

### 1.3 数据规模与划分

Ribonanza的训练集包含**806,578条序列（长度50–100 nt）**，测试集包含**1,343,823条序列（长度100–240 nt）**。RNAGym基准对Ribonanza数据进行了进一步处理：过滤信噪比<1.0的序列，去除首尾无化学映射数据的区域，最终得到**901k反应性剖面、583k唯一序列、1亿个核苷酸**的数据集，并使用**40%序列一致性聚类**进行20-80的train-test划分。


## 二、数据清洗

### 2.1 信噪比过滤

Ribonanza数据中，部分序列的化学映射信号质量较低。**信噪比（SN_filter）< 1.0的序列应被过滤**。RNAGym的预处理流程明确将信噪比<1.0的序列移除。

```python
import pandas as pd

def filter_low_snr(df, sn_col="SN_filter", threshold=1.0):
    """过滤低信噪比序列"""
    return df[df[sn_col] >= threshold].copy()
```

### 2.2 首尾无效区域裁剪

由于实验方法的限制，序列的**5′和3′端部分位置没有反应性数据**。需要将这些位置裁剪掉，只保留有有效标签的区域。

```python
import numpy as np

def trim_invalid_ends(sequence, reactivity, reactivity_error):
    """裁剪序列首尾无反应性数据的区域"""
    react = np.array(reactivity, dtype=np.float32)
    valid_mask = ~np.isnan(react)
    if not valid_mask.any():
        return None
    first_valid = np.argmax(valid_mask)
    last_valid = len(valid_mask) - np.argmax(valid_mask[::-1]) - 1
    return (
        sequence[first_valid:last_valid+1],
        react[first_valid:last_valid+1],
        np.array(reactivity_error)[first_valid:last_valid+1]
    )
```

### 2.3 缺失值处理

反应性数据中的NaN表示该位置无有效测量。对于降解/反应性预测任务，**NaN位置应被mask掉，不参与损失计算**。

```python
def clean_reactivity(reactivity, max_val=None):
    """清洗反应性标签：裁剪极端值，返回mask"""
    react = np.array(reactivity, dtype=np.float32)
    mask = ~np.isnan(react)
    if max_val is not None:
        react = np.clip(react, 0.0, max_val)
    return react, mask
```

### 2.4 序列去重与聚类

Ribonanza数据包含大量来自不同来源的序列（Eterna玩家设计、病毒窗口、mutate-and-map设计等），可能存在**完全重复或高度相似**的序列。应：

1. **完全重复序列去重**：保留反应性剖面的均值。
2. **相似性聚类**：使用CD-HIT或MMseqs2以**40%序列一致性**对序列聚类，用于后续防泄漏划分。

### 2.5 按实验类型分层

Ribonanza包含**DMS和2A3两种化学探测方法**。两种方法的反应性谱反映不同的结构信息（DMS主要探测未配对腺嘌呤和胞嘧啶，2A3探测未配对腺嘌呤和鸟嘌呤）。建议：

- 保留`experiment_type`字段作为条件。
- 可以训练**多任务模型**同时预测两种方法的反应性。
- 如果模型只关注一种探测方法，则按`experiment_type`过滤。


## 三、训练/测试数据集准备

### 3.1 样本结构

对于结构反应性预测任务，每条训练样本应包含：

```json
{
  "id": "seq_001",
  "sequence": "GGUUGCAG...",
  "experiment_type": "DMS",
  "SN_filter": 1.45,
  "reactivity": [0.12, 0.45, NaN, ...],
  "reactivity_error": [0.02, 0.05, NaN, ...],
  "reactivity_mask": [true, true, false, ...],
  "length": 107,
  "source": "eterna",
  "cluster_id": 42
}
```

### 3.2 数据划分：防止泄漏是关键

Ribonanza序列来源多样，可能存在序列相似性。**绝对不能用随机划分**。推荐策略：

| 策略 | 方法 | 说明 |
|---|---|---|
| **按序列相似性聚类划分** | 使用CD-HIT以40% identity聚类，同一簇只出现在一个集合 | 最严格，推荐 |
| **按来源分层划分** | 按`source`（Eterna、Rfam、病毒窗口等）分层 | 保证来源多样性 |
| **按实验类型划分** | train=DMS，test=2A3 | 测试跨探测方法泛化 |
| **按长度分层划分** | 按序列长度分箱后分层抽样 | 保证长度分布一致 |

**推荐方案**：先用CD-HIT以40% identity对全部序列聚类，然后按**75%:15%:10%** 的比例划分train/validation/test，确保同一簇的序列全部落入同一集合。

```python
from sklearn.model_selection import GroupShuffleSplit

def split_with_clusters(df, cluster_col="cluster_id",
                        test_size=0.10, val_size=0.15, random_state=42):
    """基于聚类的分层划分"""
    gss = GroupShuffleSplit(n_splits=1, test_size=test_size,
                            random_state=random_state)
    train_val_idx, test_idx = next(gss.split(df, groups=df[cluster_col]))

    train_val = df.iloc[train_val_idx]
    test = df.iloc[test_idx]

    val_ratio = val_size / (1 - test_size)
    gss2 = GroupShuffleSplit(n_splits=1, test_size=val_ratio,
                             random_state=random_state)
    train_idx, val_idx = next(gss2.split(train_val, groups=train_val[cluster_col]))

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

Ribonanza数据在你的模型中的**正确定位**是：

**角色一：结构感知模块的预训练任务**

将`reactivity`作为回归目标，训练一个**核苷酸级反应性预测头**：

```python
# 输入：RNA序列 + experiment_type条件
# 输出：每个位置的反应性值
# 损失：MSE（mask掉NaN位置）
```

RibonanzaNet正是以化学映射数据为训练目标，预训练后的模型在RNA二级结构预测和降解预测上达到了state-of-the-art性能。

**角色二：结构感知Attention的辅助信号**

反应性数据可以转化为**碱基配对概率**的代理信号。反应性低的位置倾向于配对，反应性高的位置倾向于未配对。这可以作为结构感知Attention的bias来源，补充或替代计算预测的配对概率矩阵。

**角色三：局部降解倾向的辅助标签**

Ribonanza数据中包含降解相关标签（`deg_*`字段），可以训练模型预测**局部降解倾向**。虽然这仍然是短片段加速条件下的数据，但可以为完整mRNA的稳定性评估提供序列-结构-降解关系的预训练表示。

**角色四：跨探测方法泛化验证**

DMS和2A3两种探测方法的数据可以用于测试模型对**不同化学探测条件**的泛化能力。如果模型在DMS上训练后能在2A3上保持较好性能，说明它学到了通用的RNA结构表示。

### 3.5 关键局限性

**Ribonanza不是完整mRNA稳定性数据。** 它测量的是**50–240 nt短RNA片段**在**体外化学探测条件**下的反应性和降解剖面，而非完整mRNA在生理条件下的半衰期。因此：

- **不能**直接用Ribonanza训练“mRNA半衰期预测器”。
- **可以**用它训练“RNA结构反应性预测器”和“局部降解倾向预测器”，作为结构感知模块的预训练。
- 完整mRNA的半衰期预测仍需从**RNA-seq time course、4sU-seq、SLAM-seq**等数据中获取标签。


## 四、完整流程总结

```text
Ribonanza 数据获取
├── Kaggle竞赛 → train_data.csv / test_sequences.csv
├── DasLab/DataPrepRibonanzaKaggle2023 → 数据准备脚本
└── Shujun-He/RibonanzaNet → 训练代码与预训练权重

        ↓

数据清洗
├── 信噪比过滤（SN_filter ≥ 1.0）
├── 首尾无效区域裁剪（去除无反应性数据的端部）
├── NaN处理：mask，不参与损失
├── 完全重复序列去重
├── 相似性聚类（40% identity）
└── 按实验类型分层（DMS / 2A3）

        ↓

训练/测试集准备
├── 样本：sequence + reactivity + experiment_type + mask
├── 划分策略：聚类分组划分（75:15:10）
├── 防泄漏检查：序列重复 + 聚类重叠
├── 角色定位：结构预训练 / 反应性预测 / 局部降解辅助
└── 局限性声明：短片段+体外探测 ≠ 完整mRNA半衰期

        ↓

输出：用于结构感知预训练和反应性/降解预测的辅助数据集
```


## 五、关键注意事项

1. **信噪比过滤是必须的**：`SN_filter < 1.0`的序列测量质量低，会引入噪声，应过滤。
2. **首尾区域需要裁剪**：实验方法无法捕获序列首尾的化学映射信号，这些位置的反应性为NaN，应裁剪或mask。
3. **相似性聚类划分是必须的**：Ribonanza序列来源多样，可能存在近重复。使用CD-HIT以40% identity聚类后分组划分，这是RNAGym等基准采用的标准做法。
4. **DMS和2A3应作为条件保留**：两种探测方法反映不同的结构信息，应作为模型的条件输入，而非混合训练。
5. **Ribonanza是结构数据，不是稳定性数据**：它是短RNA片段在体外化学探测条件下的反应性/降解剖面，在mRNA设计模型中应定位为**结构感知模块的预训练信号**和**局部降解倾向的辅助标签**，而非完整mRNA半衰期的直接预测器。
6. **RibonanzaNet可直接复用**：如果不想从头训练，可以直接使用RibonanzaNet的预训练权重作为结构感知模块的初始化，然后在你的mRNA设计任务上微调。