从OpenVaccine获取数据来训练mRNA设计模型，核心在于将其提供的**短RNA片段（107 nt）在核苷酸分辨率上的降解速率标签**，转化为模型可用的**序列—降解倾向预测任务**，并明确其在完整mRNA设计体系中的定位：**局部稳定性/降解风险预测器**，而非全长mRNA半衰期的直接预测器。

---

## 一、从OpenVaccine获取数据

OpenVaccine是Stanford Das Lab于2020年在Kaggle上发起的“COVID-19 mRNA Vaccine Degradation Prediction”竞赛数据集，原始数据来自Eterna平台的**3000+条RNA分子**，每条分子在**加速降解条件**（10 mM MgCl₂, 50 mM Na-CHES, pH 10, ~24°C）下测量了多个时间点的降解情况。

### 1.1 核心数据资源

| 资源 | 地址 | 内容 |
|---|---|---|
| **Kaggle竞赛（原始数据）** | kaggle.com/c/stanford-covid-vaccine/data | train.json、test.json、sample_submission.csv |
| **官方GitHub仓库** | github.com/DasLab/KaggleOpenVaccine | 训练集/测试集（.csv和.json）、模型代码、预测脚本 |
| **降解速率汇编** | github.com/DasLab/RNA-deg-rates | `RNA_deg_rates.csv`，汇总的k_deg值、误差和序列位置 |
| **HuggingFace数据集镜像** | 搜索“OpenVaccine” | 部分研究者上传的预处理版本 |

官方GitHub仓库的 `data/Kaggle_RYOS_data/` 目录包含训练集和测试集的CSV与JSON格式文件，可直接下载使用。

### 1.2 核心字段结构

每条样本包含以下关键字段：

| 字段 | 类型 | 说明 |
|---|---|---|
| `id` | string | 样本唯一标识 |
| `sequence` | 107字符字符串 | RNA序列（A/G/U/C），训练集和公共测试集为107 nt |
| `structure` | 107字符字符串 | 二级结构，`(` `)` 表示配对，`.` 表示未配对 |
| `seq_scored` | int | 参与评分的核苷酸数，训练集为68 |
| `seq_length` | int | 序列长度，训练集为107 |
| `reactivity` | 68维浮点向量 | 前68个碱基的反应性，用于推断二级结构 |
| `deg_pH10` | 68维浮点向量 | 无Mg²⁺、高pH（pH 10）条件下的降解值 |
| `deg_Mg_pH10` | 68维浮点向量 | 有Mg²⁺、高pH条件下的降解值 |
| `deg_50C` | 68维浮点向量 | 50°C条件下的降解值 |
| `deg_Mg_50C` | 68维浮点向量 | 有Mg²⁺、50°C条件下的降解值 |

Kaggle竞赛的官方评分仅针对 **`reactivity`、`deg_Mg_pH10` 和 `deg_Mg_50C`** 三列。对于mRNA设计模型，**`deg_Mg_pH10` 和 `deg_Mg_50C` 是最相关的降解标签**，分别对应加速水解条件和热应激条件下的降解倾向。


## 二、数据清洗

OpenVaccine数据相对干净，清洗重点在于**处理序列截断、缺失值编码和异常值**。

### 2.1 序列格式与长度标准化

- **字母表确认**：序列已经是RNA字母表（`A/U/G/C`），无需U/T转换。
- **长度处理**：训练集序列长度统一为107 nt，但**仅有前68个位置有标签**。这意味着每条样本实际可用的监督信号仅限于前68 nt。
- **截断策略**：对于降解预测任务，建议**将序列截取为前68 nt**，与标签长度对齐。剩余的39 nt可用于结构上下文建模，但不参与降解预测的损失计算。

```python
def truncate_to_scored(sequence, labels, seq_scored=68):
    """将序列截取到评分长度"""
    return sequence[:seq_scored], labels[:seq_scored]
```

### 2.2 缺失值与特殊值处理

- **NaN处理**：部分样本可能在个别位置存在NaN。对于降解预测任务，NaN位置应被mask掉，不参与损失计算。
- **极端值**：降解值理论上应在[0,1]范围内，但实验测量可能产生轻微越界值。建议将超出[0, 1]的值裁剪到边界，或标记为异常样本。

```python
import numpy as np

def clean_labels(labels, max_val=1.0):
    """清洗降解标签：裁剪极端值，返回mask"""
    labels = np.array(labels, dtype=np.float32)
    mask = ~np.isnan(labels)  # NaN位置不参与训练
    labels = np.clip(labels, 0.0, max_val)
    return labels, mask
```

### 2.3 重复序列处理

OpenVaccine的3000+条序列可能包含重复或高度相似的分子。**完全相同的序列应去重**，保留降解值的均值。对于高度相似但不完全相同的序列（如仅少数碱基差异），如果都出现在训练集中，可能导致模型过拟合到细微差异。建议对训练集**内部**也做相似性聚类，确保同一簇的序列不跨train/test划分。

### 2.4 与实验条件的绑定

每条样本的降解值对应特定的**实验条件**（Mg²⁺浓度、pH、温度）。这些条件是模型的条件输入，必须作为样本的一部分保留：

```python
# 条件编码示例
condition = {
    "mg_conc": 10.0,      # mM
    "ph": 10.0,
    "temperature": 24.0,  # °C
    "modification": "unmodified"  # 或 "m1Ψ" / "Ψ"
}
```

注意：OpenVaccine主数据集为**未修饰RNA**。Das Lab后续的PERSIST-seq研究包含了**假尿苷（Ψ）和N1-甲基假尿苷（m1Ψ）**修饰条件下的降解数据，可在`RNA-deg-rates`仓库或相关论文补充材料中获取。如果模型需要支持修饰核苷酸条件，应整合这部分数据。


## 三、训练/测试数据集准备

### 3.1 样本结构

对于降解预测任务，每条训练样本应包含：

```json
{
  "id": "id_001",
  "sequence": "GGUUGCAG...",
  "structure": "..(((...))...",
  "reactivity": [0.12, 0.45, ...],
  "deg_Mg_pH10": [0.08, 0.31, ...],
  "deg_Mg_50C": [0.15, 0.42, ...],
  "seq_scored": 68,
  "conditions": {
    "mg_conc": 10.0,
    "ph": 10.0,
    "temperature": 24.0,
    "modification": "unmodified"
  }
}
```

**预测目标的选择**：
- 如果模型关注**加速水解条件下的固有稳定性**，以`deg_Mg_pH10`为主要目标。
- 如果模型关注**热应激下的稳定性**，以`deg_Mg_50C`为主要目标。
- 也可以**多任务同时预测**多个降解条件，让模型学习不同条件下的降解模式差异。

### 3.2 数据划分：防止泄漏是关键

OpenVaccine竞赛本身提供了固定的公共测试集（629条序列），但用于模型训练时，需要从训练集中重新划分train/validation。

**绝对不能用随机划分。** RNA序列之间可能存在序列相似性或结构相似性，随机划分会导致近重复序列跨集合，严重高估性能。

| 策略 | 方法 | 说明 |
|---|---|---|
| **按序列相似性聚类划分** | 使用CD-HIT或MMseqs2对107 nt序列聚类（90% identity），同一簇只出现在一个集合 | 最严格，推荐 |
| **按结构相似性划分** | 使用RNA结构聚类（如RNAclust），同一结构簇不跨集合 | 适合结构感知模型 |
| **按GC含量分层划分** | 按GC含量分箱后分层抽样 | 保证GC分布一致 |
| **官方公共测试集** | 使用竞赛提供的629条公共测试序列作为独立测试集 | 直接可用，但需确认其与训练集的相似性 |

**推荐方案**：先用CD-HIT以90% identity对全部序列聚类，然后按**65%:15%:20%** 的比例划分train/validation/test，确保同一簇的序列全部落入同一集合。如果使用官方公共测试集，应先检查其与训练集的序列相似性，排除潜在泄漏。

```python
from sklearn.model_selection import GroupShuffleSplit

def split_with_clusters(df, cluster_col="cluster_id", 
                        test_size=0.2, val_size=0.15, random_state=42):
    """基于聚类的分层划分"""
    # 先划分出test
    gss = GroupShuffleSplit(n_splits=1, test_size=test_size, 
                            random_state=random_state)
    train_val_idx, test_idx = next(gss.split(df, groups=df[cluster_col]))
    
    train_val = df.iloc[train_val_idx]
    test = df.iloc[test_idx]
    
    # 再从train_val中划分validation
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
def check_leakage(train_df, test_df, seq_col="sequence"):
    issues = []
    # 完全重复序列
    overlap = set(train_df[seq_col]) & set(test_df[seq_col])
    if overlap:
        issues.append(f"Exact sequence leakage: {len(overlap)}")
    # 聚类泄漏
    if "cluster_id" in train_df.columns:
        overlap = set(train_df["cluster_id"]) & set(test_df["cluster_id"])
        if overlap:
            issues.append(f"Cluster leakage: {len(overlap)} clusters")
    return issues
```

### 3.4 在mRNA设计模型中的角色定位

OpenVaccine数据在你的模型中的**正确定位**是：

**角色一：训练局部降解倾向预测器（辅助任务）**

将`deg_Mg_pH10`或`deg_Mg_50C`作为回归目标，训练一个**核苷酸级降解预测头**：

```python
# 输入：RNA序列（68 nt）+ 条件
# 输出：每个位置的降解倾向（68维向量）
# 损失：MSE（mask掉NaN位置）
```

这个预测器可以作为**可制造性/稳定性模块的辅助信号**。虽然它预测的是短片段在加速条件下的降解，但模型学到的序列—结构—降解关系可以部分迁移到完整mRNA的稳定性评估中。

**角色二：结构感知模块的辅助预训练**

`structure`字段和`reactivity`字段提供了**实验推断的二级结构信息**。可以用这些数据训练模型预测**每个位置的配对状态或反应性**，作为结构感知模块的预训练任务。

**角色三：负样本/低稳定性样本的来源**

降解值高的序列可以作为**低稳定性负样本**，用于训练模型的稳定性预测头。将同一序列在不同条件（pH10 vs 50°C）下的降解差异作为对比信号，帮助模型学习条件依赖的稳定性模式。

### 3.5 关键局限性

**OpenVaccine不是完整mRNA稳定性数据。** 它测量的是**107 nt短RNA片段**在**加速降解条件**下的降解速率，而非完整mRNA在生理条件下的半衰期。因此：

- **不能**直接用OpenVaccine训练“mRNA半衰期预测器”。
- **可以**用它训练“局部降解倾向预测器”，作为稳定性评估的辅助特征。
- 完整mRNA的半衰期预测仍需从**RNA-seq time course、4sU-seq、SLAM-seq**等数据中获取标签。


## 四、完整流程总结

```text
OpenVaccine 数据获取
├── Kaggle竞赛 → train.json / test.json
├── 官方GitHub DasLab/KaggleOpenVaccine → CSV/JSON训练集
└── DasLab/RNA-deg-rates → k_deg汇编数据

        ↓

数据清洗
├── 序列截取到前68 nt（与标签对齐）
├── NaN处理：mask，不参与损失
├── 极端值裁剪到[0,1]
├── 完全重复序列去重
├── 条件编码：Mg²⁺、pH、温度、修饰
└── 相似性聚类（用于划分）

        ↓

训练/测试集准备
├── 样本：sequence(68) + structure + deg_Mg_pH10 + deg_Mg_50C
├── 划分策略：聚类分组划分（65:15:20）
├── 防泄漏检查：序列重复 + 聚类重叠
├── 角色定位：局部降解预测辅助任务 / 结构预训练 / 负样本
└── 局限性声明：短片段+加速条件 ≠ 完整mRNA半衰期

        ↓

输出：用于局部降解倾向预测和结构感知预训练的辅助数据集
```


## 五、关键注意事项

1. **标签仅覆盖前68个位置**：序列长107 nt，但只有前68个位置有降解标签。训练时应截取到68 nt，或将后39 nt作为结构上下文但不计算损失。
2. **OpenVaccine是加速降解数据**：10 mM Mg²⁺ + pH 10 + 24°C的条件下测量，降解速率远高于生理条件。模型学到的“降解倾向”需要校准后才能用于生理条件下的稳定性预测。
3. **相似性聚类划分是必须的**：RNA序列之间可能存在近重复，随机划分会严重高估性能。使用CD-HIT以90% identity聚类后分组划分。
4. **官方公共测试集需检查泄漏**：629条公共测试序列如果与训练集存在近重复，不能作为独立测试集使用。
5. **修饰核苷酸数据有限**：主数据集为未修饰RNA，m1Ψ/Pseudouridine修饰条件下的降解数据需要从Das Lab后续研究中补充。
6. **不要将OpenVaccine直接等同于mRNA半衰期数据**：它是短片段加速降解数据，在mRNA设计模型中应定位为**辅助任务和特征来源**，而非稳定性预测的主标签。完整mRNA半衰期数据需从RNA-seq time course等来源获取。