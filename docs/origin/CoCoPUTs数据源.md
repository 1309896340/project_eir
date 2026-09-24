从CoCoPUTs为mRNA设计模型获取数据，核心在于将其提供的**密码子/密码子对使用频率**，转化为模型可用的**序列设计偏好与约束**。它不直接提供用于训练的序列样本，而是提供构建模型训练目标（如密码子优化评分）和约束（如避免稀有密码子）的“金标准”统计数据。

---

## 一、从CoCoPUTs获取数据

CoCoPUTs（Codon and Codon-Pair Usage Tables）由FDA开发，提供基于GenBank/RefSeq的、按物种和人类组织分层的密码子、密码子对及二核苷酸使用频率数据。

### 1.1 核心数据资源

| 数据库 | 内容 | 用途 |
|---|---|---|
| **CoCoPUTs** | 物种级密码子/密码子对使用频率 | 通用密码子优化 |
| **TissueCoCoPUTs** | 52种人类组织特异性密码子使用 | 组织靶向mRNA设计 |
| **CancerCoCoPUTs** | 32种肿瘤类型密码子使用 | 肿瘤新抗原设计 |
| **Cell Lines CoCoPUTs** | 1,866种细胞系的密码子使用 | 体外表达优化 |

### 1.2 下载方式

**方式一：Web界面下载（推荐用于快速获取）**

访问FDA HIVE门户的CoCoPUTs页面（`https://dnahive.fda.gov/dna.cgi?cmd=cuts_main`）。用户可通过搜索窗口选择物种或细胞系，密码子、密码子对、二核苷酸使用数据以独立标签页展示，所有数据集均可下载为Excel文件。Cell Lines CoCoPUTs的数据也可在`https://dnahive.fda.gov/hivecuts/cell-lines/`下载。

**方式二：直接TSV下载（推荐用于程序化处理）**

CoCoPUTs的底层数据以TSV格式存储，可通过URL直接获取。例如，人类RefSeq密码子使用表的直接下载链接为：

```
https://dnahive.fda.gov/dna.cgi?cmd=objFile&ids=537&filename=Refseq_Bicod.tsv&raw=1
```

人类参考密码子使用数据也可通过以下URL获取：

```
https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs
```

**方式三：代码库集成**

VaxPress项目提供了从CoCoPUTs自动下载并解析密码子使用表的Python脚本，可直接参考其实现。BaseBuddy项目也集成了CoCoPUTs数据库，提供了`cocoput_table.tsv`和`cocoput_index.csv`的本地加载逻辑。


## 二、数据清洗

CoCoPUTs是**聚合统计表**，不包含原始序列，因此“清洗”重点在于**格式规范化和数据一致性**，而非序列层面的QC。

### 2.1 格式规范化

CoCoPUTs的TSV文件存在一个已知问题：**包含过多的尾随制表符（trailing tab）**，读取时需要显式处理。同时，不同来源的数据集列名和格式可能不一致，需要统一为标准列名。

```python
import pandas as pd

def load_cocoputs_tsv(filepath):
    """加载并清洗CoCoPUTs TSV文件"""
    with open(filepath, 'r') as f:
        content = f.read()
    # 去除尾随制表符
    content = content.replace('\t\n', '\n')
    from io import StringIO
    df = pd.read_csv(StringIO(content), sep='\t')
    return df

def normalize_codon_table(df):
    """标准化密码子使用表"""
    # 确保密码子列为大写
    df['codon'] = df['codon'].str.upper()
    # 频率列归一化为比例（总和为1）
    freq_col = [c for c in df.columns if 'freq' in c.lower() or 'usage' in c.lower()]
    if freq_col:
        total = df[freq_col[0]].sum()
        df['frequency_normalized'] = df[freq_col[0]] / total
    return df
```

### 2.2 数据一致性检查

| 检查项 | 标准 | 处理 |
|---|---|---|
| 密码子完整性 | 64个密码子齐全 | 缺失则标记或从其他来源补充 |
| 频率总和 | 每种氨基酸的密码子频率之和 = 1 | 不满足则重新归一化 |
| 终止密码子 | 频率通常为0或极低 | 确认不参与氨基酸编码 |
| 物种匹配 | 与目标模型物种一致 | 不一致则切换数据表 |
| 数据版本 | 与GenBank/RefSeq版本对应 | 记录版本号 |

CoCoPUTs数据基于GenBank和RefSeq，**每三个月随新版本更新一次**。训练时应记录所用数据版本，确保可复现。

### 2.3 与序列数据的交叉验证

如果将CoCoPUTs的密码子频率与从GENCODE/RefSeq提取的实际CDS序列进行对比，可以验证统计一致性。例如，检查实际序列中某密码子的观测频率是否与CoCoPUTs报告的理论频率一致，偏差过大的密码子需要排查原因（如数据版本不一致、物种差异）。


## 三、训练/测试数据集准备

CoCoPUTs在mRNA设计模型训练中扮演**特征来源**和**优化目标**的角色，而非直接的训练样本来源。

### 3.1 将CoCoPUTs转化为模型特征

**特征一：密码子适应指数（CAI）**

CAI是衡量CDS密码子使用偏好的经典指标，直接基于CoCoPUTs的密码子频率计算：

```python
import numpy as np

def compute_cai(cds_sequence, codon_freq_table):
    """
    计算密码子适应指数
    codon_freq_table: dict, {codon: frequency}
    """
    codons = [cds_sequence[i:i+3] for i in range(0, len(cds_sequence)-2, 3)]
    # 获取每个氨基酸的最大密码子频率
    aa_max_freq = {}
    codon_to_aa = {...}  # 标准遗传密码
    for codon, freq in codon_freq_table.items():
        aa = codon_to_aa[codon]
        aa_max_freq[aa] = max(aa_max_freq.get(aa, 0), freq)
    
    # 计算几何平均
    log_sum = 0
    for codon in codons:
        aa = codon_to_aa[codon]
        w = codon_freq_table[codon] / aa_max_freq[aa]
        log_sum += np.log(w)
    return np.exp(log_sum / len(codons))
```

**特征二：密码子对偏好评分**

CoCoPUTs的独特价值在于提供**密码子对使用频率**，这是传统密码子表不具备的。密码子对偏好影响翻译延伸效率，可作为模型的特征或约束。

```python
def compute_codon_pair_score(cds, codon_pair_table):
    """计算CDS的密码子对偏好评分"""
    codons = [cds[i:i+3] for i in range(0, len(cds)-2, 3)]
    scores = []
    for i in range(len(codons) - 1):
        pair = codons[i] + codons[i+1]
        scores.append(codon_pair_table.get(pair, 0))
    return np.mean(scores)
```

**特征三：组织/细胞类型特异的密码子使用**

TissueCoCoPUTs和Cell Lines CoCoPUTs提供了组织或细胞系特异性的密码子使用数据。对于条件控制模型，可以将这些数据作为**条件特征**：

```python
# 为不同细胞类型加载对应的密码子表
cell_type_tables = {
    "HEK293T": load_cocoputs("hek293t_codon_table.tsv"),
    "hepatocyte": load_cocoputs("liver_tissue_codon_table.tsv"),
    "dendritic_cell": load_cocoputs("dendritic_cell_codon_table.tsv"),
}
```

### 3.2 构建训练样本

CoCoPUTs本身不提供“输入-输出”配对样本，但可以结合GENCODE/RefSeq的CDS序列构建以下训练任务：

**任务一：密码子优化评分预测**

```json
{
  "protein_sequence": "MKT...",
  "cds_variant": "AUGAAAACC...",
  "cai_score": 0.78,
  "codon_pair_score": 0.65,
  "organism": "human",
  "cell_type": "HEK293T"
}
```

模型学习预测给定CDS变体的密码子使用质量评分，评分由CoCoPUTs数据计算得出。

**任务二：条件密码子生成**

在条件生成模型中，CoCoPUTs的频率数据可作为**软约束**或**奖励信号**：

```python
def codon_usage_reward(cds, cell_type, codon_freq_tables):
    """基于CoCoPUTs的奖励函数"""
    freq_table = codon_freq_tables[cell_type]
    codons = [cds[i:i+3] for i in range(0, len(cds)-2, 3)]
    rewards = [freq_table.get(c, 0.01) for c in codons]
    return np.mean(np.log(rewards))  # 对数似然
```

### 3.3 数据划分策略

由于CoCoPUTs是统计数据，划分的重点在于**避免统计泄漏**：

| 策略 | 方法 | 说明 |
|---|---|---|
| **按物种划分** | train=human+mouse, test=rat | 测试跨物种泛化 |
| **按组织划分** | train=肝脏+肌肉, test=肿瘤 | 测试跨组织泛化 |
| **按细胞系划分** | train=HEK293T+HeLa, test=CHO | 测试跨细胞系泛化 |
| **按数据版本划分** | train=CoCoPUTs 2024Q1, test=2024Q2 | 测试时序泛化 |

**关键原则**：如果使用特定细胞系的密码子表作为条件，该细胞系的CDS序列不应同时出现在训练集中作为优化目标，否则模型可能学会“作弊”。


## 四、完整流程总结

```text
CoCoPUTs 数据获取
├── Web界面 → 按物种/细胞系下载Excel
├── 直接TSV → 通过URL获取原始统计表
└── 代码集成 → 参考VaxPress/BaseBuddy实现

        ↓

数据清洗
├── 去除尾随制表符，标准化列名
├── 密码子频率归一化（每种氨基酸总和=1）
├── 完整性检查（64密码子齐全）
├── 版本记录（与GenBank/RefSeq版本对应）
└── 与序列数据交叉验证

        ↓

训练/测试集准备
├── 特征计算：CAI、密码子对评分、组织特异性频率
├── 奖励函数：基于CoCoPUTs的对数似然
├── 条件特征：细胞类型/组织/肿瘤特异性密码子表
├── 训练样本：结合GENCODE/RefSeq CDS构建
└── 划分策略：按物种/组织/细胞系/版本划分

        ↓

输出：融合CoCoPUTs密码子使用偏好的mRNA设计训练数据集
```


## 五、关键注意事项

1. **CoCoPUTs是统计表，不是序列库**：它提供的是密码子使用频率，需要与GENCODE/RefSeq的CDS序列结合才能构建训练样本。
2. **密码子对数据是独特价值**：CoCoPUTs的密码子对使用频率数据在其他数据库中难以获取，对翻译延伸效率建模非常重要。
3. **组织/细胞类型特异性**：TissueCoCoPUTs和Cell Lines CoCoPUTs为条件控制模型提供了关键的组织特异性密码子偏好数据。
4. **版本一致性**：CoCoPUTs每季度更新，训练时必须记录所用版本，避免与序列数据版本不匹配。
5. **数据格式需处理**：TSV文件存在尾随制表符等格式问题，读取时需要显式清洗。
6. **避免统计泄漏**：如果使用特定细胞系的密码子表作为条件，该细胞系的序列不应同时用于优化目标训练。