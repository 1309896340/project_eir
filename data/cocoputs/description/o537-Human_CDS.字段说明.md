# o537-Human_CDS.tsv 字段说明

> 人类 RefSeq CDS 的**密码子组成统计表**（来源 CoCoPUTs / HIVE-CUTs）。每行对应一条蛋白编码序列（CDS），给出其基因/基因组元数据、密码子总数、GC 含量，以及 64 种密码子在该 CDS 中的出现次数。
>
> 本文档中的"已验证"结论均基于对全表 119,196 行的实际校验（2026-09-22），非仅凭字段名推测。

## 文件基本信息

| 属性 | 值 |
|---|---|
| 路径 | `data/o537-Human_CDS.tsv` |
| 大小 | 30,957,713 字节（约 29.5 MB） |
| 行数 | 119,197 行 = 1 行表头 + 119,196 行数据（每行一条 CDS） |
| 列数 | 77（数据行行尾多带一个制表符，解析时会多出一个空列，见文末注意） |
| 分隔符 | 制表符 `\t`，无注释行 |
| 物种 | Homo sapiens（Taxid 9606），参考基因组 GCF_000001405.39（GRCh38.p13） |
| 数据源 | CoCoPUTs（HIVE-CUTs 的持续更新版），FDA DNAHIVE：<https://dnahive.fda.gov> |
| 文件名前缀 | `o537` 为 CoCoPUTs 下载库内部的物种文件编号（非 NCBI Taxid；Taxid 见数据内第 6 列）。同前缀的另外三个文件为同物种的二核苷酸 / 密码子对（bicodon）统计表 |

## 字段说明

### 一、元数据字段（第 1–8 列）

| 列号 | 字段 | 类型 | 含义 | 示例 |
|---|---|---|---|---|
| 1 | Gene ID | string | RefSeq 基因符号：有官方 symbol 的用 symbol，否则用 `LOC*` 形式的编号。全表共 20,239 个不同基因 | `OR4F5`、`LOC112268260` |
| 2 | Protein ID | string | RefSeq 蛋白 accession.version，按前缀分质量层级：`NP_` 人工审校（60,249 条）、`XP_` 计算预测的模型转录本（58,360 条）、`YP_` 线粒体基因组编码蛋白（13 条）；另有 574 行为空（免疫球蛋白/TCR 基因段等无蛋白 accession 的记录） | `NP_001005484.1` |
| 3 | Accession | string | CDS 所在基因组 contig 的 accession（**不含版本号**）。`NC_000001`–`NC_000024` 为染色体（`NC_012920` = 线粒体 rCRS），`NT_*` 为 alternate / 未定位 scaffold | `NC_000001` |
| 4 | Division | string | 来源数据库分区，全表恒为 `refseq` | `refseq` |
| 5 | Assembly | string | 参考基因组组装 accession，全表恒为 `GCF_000001405.39`（GRCh38.p13） | `GCF_000001405.39` |
| 6 | Taxid | int | NCBI 分类学 ID，全表恒为 9606 | `9606` |
| 7 | Species | string | 物种名，全表恒为 `Homo sapiens` | `Homo sapiens` |
| 8 | Organelle | string | 序列来源细胞器：`genomic`（核基因组，119,183 条）/ `mitochondrion`（线粒体，13 条） | `genomic` |

### 二、总量与 GC 统计字段（第 9–13 列）

| 列号 | 字段 | 类型 | 含义 | 示例 |
|---|---|---|---|---|
| 9 | # Codons | int | 该 CDS 的密码子总数，**含终止密码子**（已验证：ND1 蛋白 318 aa + 终止 = 319；OR4F5 305 aa + 终止 = 306；且 64 列计数之和与之相等，全表零例外） | `306` |
| 10 | GC% | float | 整条 CDS 的 GC 百分比，≈ (GC1% + GC2% + GC3%) / 3（已抽样验证） | `42.81` |
| 11 | GC1% | float | 密码子第 1 位上的 GC 百分比 | `44.44` |
| 12 | GC2% | float | 密码子第 2 位上的 GC 百分比 | `34.64` |
| 13 | GC3% | float | 密码子第 3 位上的 GC 百分比。第 3 位是同义替换的主要自由度，是密码子优化 / GC 约束解码的关键指标 | `49.35` |

注：933 行的 CDS 无终止密码子（多为不完整的预测注释），361 行含多个终止密码子（多为硒蛋白等内部 TGA 注释特例）。

### 三、密码子计数列（第 14–77 列，共 64 列）

- 每列为对应 **DNA 三联体密码子**在该 CDS 中的出现次数，整数**原始计数**（非频率、非千分比）。
- 密码子用 **DNA 字母表（T）**，本项目统一用 RNA（U），使用时注意 T↔U 转换。
- 列排列规则：按「第二位碱基（T→A→C→G）× 第一位碱基（T→C→A→G）× 第三位碱基（T→C→A→G）」排列，与同前缀二核苷酸文件的 16 列顺序兼容（密码子前两位 = 二核苷酸）。
- 终止密码子 TAA/TAG/TGA 的计数包含在第 32/33/64 列中，且计入 `# Codons`。

| 列号 | DNA 密码子 | RNA 密码子 | 编码氨基酸 | 列号 | DNA 密码子 | RNA 密码子 | 编码氨基酸 |
|---|---|---|---|---|---|---|---|
| 14 | TTT | UUU | Phe (F) | 46 | TCT | UCU | Ser (S) |
| 15 | TTC | UUC | Phe (F) | 47 | TCC | UCC | Ser (S) |
| 16 | TTA | UUA | Leu (L) | 48 | TCA | UCA | Ser (S) |
| 17 | TTG | UUG | Leu (L) | 49 | TCG | UCG | Ser (S) |
| 18 | CTT | CUU | Leu (L) | 50 | CCT | CCU | Pro (P) |
| 19 | CTC | CUC | Leu (L) | 51 | CCC | CCC | Pro (P) |
| 20 | CTA | CUA | Leu (L) | 52 | CCA | CCA | Pro (P) |
| 21 | CTG | CUG | Leu (L) | 53 | CCG | CCG | Pro (P) |
| 22 | ATT | AUU | Ile (I) | 54 | ACT | ACU | Thr (T) |
| 23 | ATC | AUC | Ile (I) | 55 | ACC | ACC | Thr (T) |
| 24 | ATA | AUA | Ile (I) | 56 | ACA | ACA | Thr (T) |
| 25 | ATG | AUG | Met (M)，起始密码子 | 57 | ACG | ACG | Thr (T) |
| 26 | GTT | GUU | Val (V) | 58 | GCT | GCU | Ala (A) |
| 27 | GTC | GUC | Val (V) | 59 | GCC | GCC | Ala (A) |
| 28 | GTA | GUA | Val (V) | 60 | GCA | GCA | Ala (A) |
| 29 | GTG | GUG | Val (V) | 61 | GCG | GCG | Ala (A) |
| 30 | TAT | UAU | Tyr (Y) | 62 | TGT | UGU | Cys (C) |
| 31 | TAC | UAC | Tyr (Y) | 63 | TGC | UGC | Cys (C) |
| 32 | TAA | UAA | **终止** | 64 | TGA | UGA | **终止**（硒蛋白中可编码 Sec/U） |
| 33 | TAG | UAG | **终止** | 65 | TGG | UGG | Trp (W) |
| 34 | CAT | CAU | His (H) | 66 | CGT | CGU | Arg (R) |
| 35 | CAC | CAC | His (H) | 67 | CGC | CGC | Arg (R) |
| 36 | CAA | CAA | Gln (Q) | 68 | CGA | CGA | Arg (R) |
| 37 | CAG | CAG | Gln (Q) | 69 | CGG | CGG | Arg (R) |
| 38 | AAT | AAU | Asn (N) | 70 | AGT | AGU | Ser (S) |
| 39 | AAC | AAC | Asn (N) | 71 | AGC | AGC | Ser (S) |
| 40 | AAA | AAA | Lys (K) | 72 | AGA | AGA | Arg (R) |
| 41 | AAG | AAG | Lys (K) | 73 | AGG | AGG | Arg (R) |
| 42 | GAT | GAU | Asp (D) | 74 | GGT | GGU | Gly (G) |
| 43 | GAC | GAC | Asp (D) | 75 | GGC | GGC | Gly (G) |
| 44 | GAA | GAA | Glu (E) | 76 | GGA | GGA | Gly (G) |
| 45 | GAG | GAG | Glu (E) | 77 | GGG | GGG | Gly (G) |

## 已验证的数据性质（全表校验，2026-09-22）

- 64 列密码子计数之和 = `# Codons`：119,196 行零例外。
- 与 `o537-Human_CDS_Dinuc.tsv`、`o537-Human_CDS_Junc_Dinuc.tsv`、`o537-Human_CDS_Bicod.tsv` **同记录、同行序**（Gene ID + Protein ID 逐行比对零错位），可按行号直接拼接四个表的特征。
- 三表中 GC% / GC1% / GC2% / GC3% 完全一致。
- 密码子计数与两张二核苷酸表逐类型精确互推（闭合校验 1,886,432 项，99.99% 通过；例外行为含内部终止密码子的硒蛋白类注释），口径见另三份字段说明。
- **(Gene ID, Protein ID) 键不唯一**：6,358 行重复，主要因同一 CDS 同时注释在主染色体与 alternate scaffold（例：`ACACA|NP_942133.1` 同时出现在 `NC_000017` 与 `NT_187614`）。聚合前需按 (Gene ID, Protein ID, Accession) 去重。

## 使用注意事项

1. **DNA 字母表**：所有密码子基于 DNA（T）；本项目统一 RNA（U），注意 T↔U 转换。
2. **原始计数而非频率**：跨基因比较需先归一化——除以 `# Codons` 得频率，或按氨基酸家族计算 RSCU（相对同义密码子使用度）/ CAI 权重。
3. **构建"人类密码子偏好参考"时建议过滤**：保留 `Organelle == genomic`（剔除 13 条线粒体 CDS）；按需只取 `NP_`（人工审校）子集；先处理 6,358 行重复；XP_ 为预测模型，质量分层使用。
4. **组成统计 ≠ 表达标签**：本表只有序列组成信息，与表达量/稳定性无直接对应，仅用作先验分布与参考统计（与设计文档中"天然表达数据不可直接当设计标签"的原则一致）。
5. **数据划分提示**（对应设计文档防泄漏规则）：同一 Gene ID 的所有行（含其同物蛋白变体、重复行）必须进入同一数据划分。
6. 解析坑：数据行行尾多一个制表符，会多读出一个空列。

## 读取示例

```python
import pandas as pd

df = pd.read_csv("data/o537-Human_CDS.tsv", sep="\t")
df = df.dropna(axis=1, how="all")      # 去掉行尾多余制表符产生的空列
meta_cols   = df.columns[:8]           # 元数据
stat_cols   = df.columns[8:13]         # # Codons, GC%, GC1-3%
codon_cols  = df.columns[13:]          # 64 个密码子计数（DNA 字母表）
```

## 与另外三份文件的关系

`o537-Human_CDS_Dinuc.tsv`（16 种二核苷酸计数）、`o537-Human_CDS_Junc_Dinuc.tsv`（密码子交界二核苷酸计数）与 `o537-Human_CDS_Bicod.tsv`（4,096 种密码子对计数）与本文件**同记录、同行序**，可按行拼接得到每条 CDS 的完整组成特征。字段说明见各自 md。

## 来源与引用

- CoCoPUTs / TissueCoCoPUTs / CancerCoCoPUTs，FDA DNAHIVE：<https://dnahive.fda.gov>（"Available Files to Download" 页附 readme.txt，为各表官方格式定义）
- Athey J. et al. *A new and updated resource for codon usage tables.* Nucleic Acids Research, 2017（HIVE-CUTs）
- Alexaki A. et al. *Codon and Codon-Pair Usage Tables (CoCoPUTs): facilitating genetic variation analyses and standardization.* Nucleic Acids Research, 2019
- 对应设计文档：`mRNA模型设计.md` 公开数据源盘点中"密码子/密码子对特征：CoCoPUTs、HIVE-CUTs"
