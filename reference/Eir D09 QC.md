# 预处理标准与 QC 流程文档

## mRNA 设计数据预处理规范

参考附件第7节与第11节，定义从原始 GENCODE 数据到模型可用训练样本的完整预处理与质量控制标准。核心目标：

> 确保每条进入训练集的序列在 **U/T 表示、方向、CDS 结构、翻译一致性、区域边界、标签质量、数据泄漏** 七个维度上全部合格。

---

# 一、预处理总体流程

```text
原始数据（GENCODE GTF + 基因组 FASTA + 转录本 FASTA + 蛋白 FASTA）
    ↓
Step 1  格式标准化（U/T 统一、大小写、去空白）
    ↓
Step 2  方向标准化（统一 5'→3'，负链反向互补）
    ↓
Step 3  区域切分（5'UTR / CDS / 3'UTR 边界确定）
    ↓
Step 4  CDS 结构检查（长度、起始、终止、内部终止）
    ↓
Step 5  翻译一致性验证（翻译后 == 目标蛋白）
    ↓
Step 6  序列质量过滤（N、长度、同聚物、重复、极端 GC）
    ↓
Step 7  特征计算与记录
    ↓
Step 8  数据划分（按基因/染色体/批次，防泄漏）
    ↓
Step 9  最终 QC 报告与版本冻结
```

---

# 二、Step 1：U/T 表示标准化

## 2.1 标准

- **统一使用 RNA 字母表 `A/U/G/C`** 作为模型内部表示。
- 原始 GENCODE 数据为 DNA 字母表 `A/T/G/C`，需转换。
- 蛋白序列使用标准单字母氨基酸代码 `ACDEFGHIKLMNPQRSTVWY`。
- 不引入 `N`、`X`、`-` 等模糊字符；若存在，按规则处理或过滤。

## 2.2 转换规则

| 数据类型 | 原始表示 | 标准化后 |
|---|---|---|
| 基因组序列 | A/T/G/C | 提取后转 A/U/G/C |
| CDS | A/T/G/C | A/U/G/C |
| 5′UTR | A/T/G/C | A/U/G/C |
| 3′UTR | A/T/G/C | A/U/G/C |
| 蛋白 | 单字母 | 单字母，不变 |

## 2.3 实现

```python
def to_rna(seq: str) -> str:
    """DNA → RNA，统一大写"""
    return seq.upper().replace("T", "U")

def to_dna(seq: str) -> str:
    """RNA → DNA，用于比对 GENCODE 原始数据"""
    return seq.upper().replace("U", "T")

def is_valid_rna(seq: str) -> bool:
    """检查是否只含 A/U/G/C"""
    return set(seq.upper()) <= {"A", "U", "G", "C"}

def is_valid_protein(seq: str) -> bool:
    """检查是否只含标准氨基酸"""
    return set(seq.upper()) <= set("ACDEFGHIKLMNPQRSTVWY")
```

## 2.4 QC 检查项

| 检查项 | 标准 | 处理 |
|---|---|---|
| 是否只含 A/U/G/C | 是 | 否则过滤或记录 |
| 大小写 | 统一大写 | 自动转换 |
| 空白字符 | 无 | 自动去除 |
| N 比例 | < 1% | 超过则过滤 |
| 蛋白字符 | 标准 20 种 | 含 X/* 则过滤或标记 |

---

# 三、Step 2：方向标准化

## 3.1 标准

- **所有序列统一为 5′→3′ 方向。**
- GENCODE 负链基因的序列在提取后需 **反向互补**。
- 5′UTR 在 5′端，CDS 居中，3′UTR 在 3′端。

## 3.2 负链处理

```python
from Bio.Seq import Seq

def strand_normalize(seq: str, strand: str) -> str:
    """将负链序列反向互补，统一为 5'→3'"""
    if strand == "+":
        return seq.upper()
    elif strand == "-":
        return str(Seq(seq).reverse_complement()).upper()
    else:
        raise ValueError(f"Invalid strand: {strand}")

def strand_normalize_rna(seq: str, strand: str) -> str:
    """RNA 版本，输出 A/U/G/C"""
    dna = to_dna(seq)
    normalized = strand_normalize(dna, strand)
    return to_rna(normalized)
```

## 3.3 负链坐标转换

GENCODE 中负链基因的 `start` 和 `end` 是基因组坐标，提取后需注意：

- 负链上 `start > end` 在转录方向上是正常的。
- 5′UTR 在基因组坐标上位于 CDS 的 **下游**（数值更大）。
- 3′UTR 在基因组坐标上位于 CDS 的 **上游**（数值更小）。

```python
def get_regions_by_strand(cds_start, cds_end, exon_starts, exon_ends, strand):
    """根据链方向确定 5'UTR、CDS、3'UTR 的基因组坐标范围"""
    if strand == "+":
        utr5_range = (None, cds_start - 1)
        cds_range = (cds_start, cds_end)
        utr3_range = (cds_end + 1, None)
    else:
        utr5_range = (cds_end + 1, None)   # 负链：5'UTR 在坐标更大处
        cds_range = (cds_start, cds_end)
        utr3_range = (None, cds_start - 1) # 负链：3'UTR 在坐标更小处
    return utr5_range, cds_range, utr3_range
```

## 3.4 QC 检查项

| 检查项 | 标准 | 处理 |
|---|---|---|
| 方向 | 5′→3′ | 负链反向互补 |
| 5′UTR 位置 | 在 CDS 上游 | 检查边界 |
| 3′UTR 位置 | 在 CDS 下游 | 检查边界 |
| 起始密码子 | 在 CDS 开头 | 检查 |
| 终止密码子 | 在 CDS 末尾 | 检查 |

---

# 四、Step 3：区域切分与边界检查

## 4.1 标准

每条序列必须明确切分为：

```text
5'UTR — CDS — 3'UTR
```

- **CDS**：从起始密码子 `AUG` 到终止密码子 `UAA/UAG/UGA`。
- **5′UTR**：CDS 上游、同一转录本内的外显子区域。
- **3′UTR**：CDS 下游、同一转录本内的外显子区域。
- **poly(A)**：若数据包含，记录长度；否则标记为缺失。

## 4.2 边界确定规则

```python
def split_regions(full_mrna: str, cds_start: int, cds_end: int):
    """
    按 CDS 坐标切分区域
    cds_start: CDS 起始在 full_mrna 中的 0-based 位置
    cds_end:   CDS 终止在 full_mrna 中的 0-based 位置（含终止密码子）
    """
    utr5 = full_mrna[:cds_start]
    cds = full_mrna[cds_start:cds_end + 3]  # 含终止密码子
    utr3 = full_mrna[cds_end + 3:]
    return utr5, cds, utr3
```

## 4.3 终止密码子处理

GENCODE 的 CDS 注释 **通常不包括终止密码子**。需注意：

| 来源 | CDS 是否含终止密码子 | 处理 |
|---|---|---|
| GENCODE GTF CDS feature | 不含 | 从序列中补上终止密码子 |
| GENCODE 转录本 FASTA | 含（成熟 mRNA） | 直接从序列读取 |
| 自定义提取 | 需确认 | 显式添加或验证 |

```python
STOP_CODONS = {"UAA", "UAG", "UGA"}
START_CODON = "AUG"

def ensure_stop_codon(cds: str) -> str:
    """确保 CDS 以终止密码子结尾"""
    if cds[-3:] not in STOP_CODONS:
        # 从下游序列中找第一个同框终止密码子
        for i in range(0, len(cds) - 2, 3):
            if cds[i:i+3] in STOP_CODONS:
                return cds[:i+3]
        raise ValueError("No in-frame stop codon found")
    return cds
```

## 4.4 QC 检查项

| 检查项 | 标准 | 处理 |
|---|---|---|
| 5′UTR 存在性 | 可空，但需记录 | 标记 |
| CDS 存在性 | 必须有 | 否则过滤 |
| 3′UTR 存在性 | 可空，但需记录 | 标记 |
| CDS 起始 | `AUG` | 否则过滤或标记 |
| CDS 终止 | `UAA/UAG/UGA` | 否则过滤 |
| 边界一致性 | 与 GENCODE 注释一致 | 否则记录差异 |
| 区域长度 | 与 GTF 计算一致 | 不一致则排查 |

---

# 五、Step 4：CDS 结构检查

## 5.1 标准

| 检查项 | 标准 | 级别 |
|---|---|---|
| CDS 长度 | 是 3 的倍数 | 硬性 |
| 起始密码子 | `AUG` | 硬性 |
| 终止密码子 | `UAA/UAG/UGA` | 硬性 |
| 内部终止密码子 | 无 | 硬性 |
| 移码 | 无 | 硬性 |
| 长度范围 | 通常 90–15,000 nt | 软性 |
| 起始上下文 | Kozak 合理 | 软性 |

## 5.2 实现

```python
def check_cds_structure(cds: str) -> dict:
    """CDS 结构检查"""
    result = {
        "length": len(cds),
        "length_multiple_of_3": len(cds) % 3 == 0,
        "starts_with_aug": cds[:3] == "AUG",
        "ends_with_stop": cds[-3:] in STOP_CODONS,
        "internal_stop": False,
        "internal_stop_positions": [],
        "pass": True
    }

    # 检查内部终止密码子
    for i in range(3, len(cds) - 3, 3):
        codon = cds[i:i+3]
        if codon in STOP_CODONS:
            result["internal_stop"] = True
            result["internal_stop_positions"].append(i)

    # 综合判定
    result["pass"] = (
        result["length_multiple_of_3"]
        and result["starts_with_aug"]
        and result["ends_with_stop"]
        and not result["internal_stop"]
    )
    return result
```

## 5.3 常见问题与处理

| 问题 | 原因 | 处理 |
|---|---|---|
| 长度非 3 倍数 | 注释错误、提取错误 | 过滤 |
| 不以 AUG 开头 | 非标准起始、注释错误 | 过滤或标记 |
| 无终止密码子 | 注释截断 | 从下游找同框终止子 |
| 内部终止密码子 | 假基因、注释错误 | 过滤 |
| CDS 过短 | 微 ORF、注释碎片 | 按长度阈值过滤 |
| CDS 过长 | 异常注释 | 人工审查 |

---

# 六、Step 5：翻译一致性验证

## 6.1 标准

**最核心的硬性检查**：将 CDS 翻译后，必须与 GENCODE 提供的蛋白序列 **完全一致**。

```text
translate(CDS) == GENCODE_protein_sequence
```

不一致则说明：

- CDS 边界错误；
- 起始/终止密码子错误；
- 存在移码或内部终止；
- 注释版本不匹配；
- 提取过程出错。

## 6.2 实现

```python
from Bio.Seq import Seq

def translate_cds(cds: str) -> str:
    """翻译 CDS，去除终止密码子"""
    dna = to_dna(cds)
    protein = str(Seq(dna).translate(to_stop=True))
    return protein

def check_translation_consistency(cds: str, expected_protein: str) -> dict:
    """翻译一致性检查"""
    translated = translate_cds(cds)

    # 去除可能的起始 Met 差异
    translated_clean = translated.rstrip("*")
    expected_clean = expected_protein.upper().rstrip("*")

    is_match = translated_clean == expected_clean

    # 若不一致，定位第一个差异位置
    diff_pos = None
    if not is_match:
        min_len = min(len(translated_clean), len(expected_clean))
        for i in range(min_len):
            if translated_clean[i] != expected_clean[i]:
                diff_pos = i
                break
        if diff_pos is None:
            diff_pos = min_len  # 长度差异

    return {
        "match": is_match,
        "translated_length": len(translated_clean),
        "expected_length": len(expected_clean),
        "first_diff_position": diff_pos,
        "translated_preview": translated_clean[:50],
        "expected_preview": expected_clean[:50],
    }
```

## 6.3 不一致的常见原因

| 原因 | 表现 | 处理 |
|---|---|---|
| CDS 边界偏移 | 翻译产物整体移码 | 修正边界 |
| 缺少终止密码子 | 翻译到序列末尾 | 补充终止密码子 |
| 内部终止密码子 | 翻译提前截断 | 过滤或排查 |
| 版本不匹配 | 蛋白序列与 CDS 来自不同版本 | 统一版本 |
| 负链未反向互补 | 翻译产物完全不同 | 反向互补 |
| 起始密码子非 AUG | 翻译起始错误 | 过滤或标记 |
| 注释错误 | 局部差异 | 人工审查 |

## 6.4 一致性等级

| 等级 | 定义 | 处理 |
|---|---|---|
| A | 完全一致 | 通过 |
| B | 仅末尾终止密码子差异 | 通过，记录 |
| C | 起始 Met 差异 | 通过，记录 |
| D | 局部 1–2 个氨基酸差异 | 人工审查 |
| E | 移码、大量差异 | 过滤 |
| F | 完全不一致 | 过滤 |

---

# 七、Step 6：序列质量过滤

## 7.1 过滤标准

| 检查项 | 标准 | 级别 |
|---|---|---|
| 字符集 | 仅 A/U/G/C | 硬性 |
| N 比例 | < 1% | 硬性 |
| 全长长度 | 100–15,000 nt | 软性 |
| CDS 长度 | 90–10,000 nt | 软性 |
| GC 含量 | 30%–70% | 软性 |
| 最长同聚物 | ≤ 8（训练集可放宽） | 软性 |
| 重复评分 | < 阈值 | 软性 |
| 未知碱基 | 无 | 硬性 |

## 7.2 实现

```python
def sequence_qc(seq: str) -> dict:
    """序列质量检查"""
    seq = seq.upper()
    n_count = seq.count("N")
    gc = (seq.count("G") + seq.count("C")) / len(seq) if len(seq) > 0 else 0

    # 最长同聚物
    max_homo = 1
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            current += 1
            max_homo = max(max_homo, current)
        else:
            current = 1

    return {
        "length": len(seq),
        "valid_chars": is_valid_rna(seq),
        "n_count": n_count,
        "n_ratio": n_count / len(seq) if len(seq) > 0 else 0,
        "gc_content": gc,
        "max_homopolymer": max_homo,
        "pass": (
            is_valid_rna(seq)
            and n_count == 0
            and 30 <= gc * 100 <= 70
            and max_homo <= 8
        )
    }
```

## 7.3 重复序列检测

```python
def repeat_score(seq: str, k: int = 6) -> float:
    """简单重复评分：k-mer 重复频率"""
    if len(seq) < k:
        return 0.0
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    from collections import Counter
    counts = Counter(kmers)
    repeated = sum(c - 1 for c in counts.values() if c > 1)
    return repeated / len(kmers) if kmers else 0.0
```

---

# 八、Step 7：特征计算与记录

## 8.1 必算特征

| 类别 | 特征 | 方法 |
|---|---|---|
| 组成 | GC 含量、局部 GC | 自写脚本 |
| 组成 | CAI、tAI | CAI 工具 |
| 组成 | 密码子对偏好 | CoCoPUTs |
| 组成 | CpG、UpA 频率 | 自写脚本 |
| 组成 | 同聚物最长长度 | 自写脚本 |
| 组成 | 重复评分 | 自写脚本 |
| 结构 | MFE、局部 MFE | ViennaRNA |
| 结构 | 5′端局部 MFE | RNAfold |
| 结构 | 配对概率、可及性 | RNAfold / LinearFold |
| UTR | uAUG、uORF 数量 | 自写脚本 |
| UTR | Kozak 评分 | 自写脚本 |
| UTR | miRNA 位点 | TargetScan |
| UTR | ARE 数量 | 正则匹配 |
| UTR | RBP motif | POSTAR / ENCODE |
| 安全 | 限制性酶切位点 | Biopython |
| 安全 | 剪接样信号 | 正则匹配 |

## 8.2 特征记录

所有特征写入 `computed_feature` 表，必须包含：

```json
{
  "sequence_id": 1,
  "feature_name": "gc_content",
  "feature_value": 0.52,
  "window_start": null,
  "window_end": null,
  "method": "custom",
  "tool_version": "1.0.0",
  "computed_at": "2026-09-15T10:00:00Z"
}
```

局部特征必须记录窗口：

```json
{
  "sequence_id": 1,
  "feature_name": "mfe_window",
  "feature_value": -8.7,
  "window_start": 50,
  "window_end": 100,
  "method": "ViennaRNA",
  "tool_version": "2.6.4"
}
```

---

# 九、Step 8：数据划分与防泄漏

## 9.1 划分策略

| 策略 | 用途 | 实现 |
|---|---|---|
| 按染色体 | 基础泛化 | train/valid/test 用不同染色体 |
| 按基因 | 新蛋白泛化 | 同一基因不跨集合 |
| 按实验批次 | 跨批次泛化 | 不同 batch_id 分开 |
| 按细胞类型 | 跨细胞泛化 | 不同 cell_type 分开 |
| 按物种 | 跨物种泛化 | human/mouse 分开 |

## 9.2 防泄漏检查

```python
def check_leakage(train_df, test_df):
    """检查 train/test 是否有泄漏"""
    issues = []

    # 同一基因跨集合
    train_genes = set(train_df["gene_id"])
    test_genes = set(test_df["gene_id"])
    overlap_genes = train_genes & test_genes
    if overlap_genes:
        issues.append(f"Gene leakage: {len(overlap_genes)} genes")

    # 同一转录本跨集合
    train_tx = set(train_df["transcript_id"])
    test_tx = set(test_df["transcript_id"])
    overlap_tx = train_tx & test_tx
    if overlap_tx:
        issues.append(f"Transcript leakage: {len(overlap_tx)} transcripts")

    # 同一蛋白跨集合
    train_prot = set(train_df["protein_sequence"])
    test_prot = set(test_df["protein_sequence"])
    overlap_prot = train_prot & test_prot
    if overlap_prot:
        issues.append(f"Protein leakage: {len(overlap_prot)} proteins")

    # 同一实验批次跨集合
    train_batch = set(train_df["batch_id"].dropna())
    test_batch = set(test_df["batch_id"].dropna())
    overlap_batch = train_batch & test_batch
    if overlap_batch:
        issues.append(f"Batch leakage: {len(overlap_batch)} batches")

    return issues
```

## 9.3 划分记录

每次划分必须保存：

```json
{
  "split_id": "split_v1",
  "strategy": "by_chromosome",
  "train_chroms": ["chr1", "chr5", "chr7", "chr10", "chr13", "chr17", "chr21"],
  "valid_chroms": ["chr2", "chr9", "chr16"],
  "test_chroms": ["chr3", "chr8", "chr15"],
  "train_count": 45000,
  "valid_count": 5000,
  "test_count": 5000,
  "leakage_check": "pass",
  "created_at": "2026-09-15T10:00:00Z"
}
```

---

# 十、Step 9：最终 QC 报告

## 10.1 报告结构

```json
{
  "qc_id": "qc_v1",
  "source": "GENCODE v49",
  "species": "human",
  "total_input": 100000,
  "total_pass": 85000,
  "total_fail": 15000,
  "pass_rate": 0.85,
  "checks": {
    "u_t_standardization": {
      "pass": 100000,
      "fail": 0
    },
    "strand_normalization": {
      "pass": 100000,
      "fail": 0
    },
    "region_splitting": {
      "pass": 98000,
      "fail": 2000,
      "reasons": {
        "missing_cds": 1500,
        "boundary_mismatch": 500
      }
    },
    "cds_structure": {
      "pass": 90000,
      "fail": 10000,
      "reasons": {
        "not_multiple_of_3": 3000,
        "no_start_codon": 2000,
        "no_stop_codon": 2500,
        "internal_stop": 2500
      }
    },
    "translation_consistency": {
      "pass": 88000,
      "fail": 2000,
      "reasons": {
        "frameshift": 800,
        "local_mismatch": 700,
        "complete_mismatch": 500
      }
    },
    "sequence_quality": {
      "pass": 85000,
      "fail": 3000,
      "reasons": {
        "n_content": 500,
        "extreme_gc": 1500,
        "long_homopolymer": 1000
      }
    }
  },
  "final_dataset": {
    "train": 68000,
    "valid": 8500,
    "test": 8500
  },
  "version": "1.0.0",
  "created_at": "2026-09-15T10:00:00Z"
}
```

## 10.2 关键指标

| 指标 | 目标 | 说明 |
|---|---|---|
| 总通过率 | > 80% | 从原始到最终 |
| U/T 标准化通过率 | 100% | 必须 |
| 方向标准化通过率 | 100% | 必须 |
| CDS 结构通过率 | > 90% | 蛋白编码转录本 |
| 翻译一致性通过率 | > 98% | 核心指标 |
| 序列质量通过率 | > 95% | 过滤低质量 |
| 数据泄漏 | 0 | 必须 |

---

# 十一、QC 工具链

| 工具 | 用途 | 版本 |
|---|---|---|
| Biopython | 序列处理、翻译 | 1.81+ |
| pyGTF / gffread | GTF 解析、序列提取 | 最新 |
| ViennaRNA | 二级结构、MFE | 2.6+ |
| LinearFold | 长序列折叠 | 最新 |
| CAI 工具 | 密码子适应指数 | 最新 |
| CoCoPUTs | 密码子对偏好 | 最新 |
| TargetScan | miRNA 位点 | 最新 |
| POSTAR / ENCODE | RBP motif | 最新 |
| pandas / polars | 数据处理 | 最新 |
| pydantic | Schema 校验 | 2.x |
| Great Expectations | 数据质量检查 | 最新 |

---

# 十二、QC 流程总结

```text
输入：GENCODE GTF + 基因组 FASTA + 蛋白 FASTA
    ↓
[QC1] U/T 标准化         → 仅 A/U/G/C，大写
    ↓
[QC2] 方向标准化         → 统一 5'→3'，负链反向互补
    ↓
[QC3] 区域切分           → 5'UTR / CDS / 3'UTR 边界
    ↓
[QC4] CDS 结构检查       → 长度、起始、终止、内部终止
    ↓
[QC5] 翻译一致性         → translate(CDS) == 蛋白
    ↓
[QC6] 序列质量过滤       → N、GC、同聚物、重复
    ↓
[QC7] 特征计算           → 组成、结构、UTR、安全
    ↓
[QC8] 数据划分           → 按基因/染色体/批次，防泄漏
    ↓
[QC9] QC 报告            → 通过率、失败原因、版本
    ↓
输出：模型可用训练样本 + QC 报告
```

**核心原则**：

1. **翻译一致性是硬门槛**，不一致直接过滤。
2. **U/T 和方向必须 100% 标准化**，否则模型学到错误表示。
3. **CDS 结构必须完整**，长度、起始、终止、内部终止四项全查。
4. **区域边界必须与 GENCODE 一致**，否则 UTR 特征错误。
5. **数据泄漏必须为零**，同一基因/蛋白不跨集合。
6. **所有 QC 步骤可追溯**，记录版本、参数、通过率。
7. **失败样本保留原因**，便于后续排查和模型负样本构建。