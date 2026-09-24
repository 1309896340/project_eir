# o537-Human_CDS_Junc_Dinuc.tsv 字段说明

> 人类 RefSeq CDS 的**密码子交界二核苷酸统计表**（来源 CoCoPUTs / HIVE-CUTs）。每行对应一条蛋白编码序列（CDS），只统计**跨越密码子边界**的二核苷酸——即上一密码子第 3 位与下一密码子第 1 位组成的有序碱基对。与 `o537-Human_CDS.tsv` **同记录、同行序**（同一批 119,196 条 CDS）。
>
> 本文档中的"已验证"结论均基于对全表 119,196 行的实际校验（2026-09-22），非仅凭字段名推测。

## 文件基本信息

| 属性 | 值 |
|---|---|
| 路径 | `data/o537-Human_CDS_Junc_Dinuc.tsv` |
| 大小 | 18,635,035 字节（约 17.8 MB） |
| 行数 | 119,197 行 = 1 行表头 + 119,196 行数据 |
| 列数 | 29（数据行行尾多带一个制表符，解析时会多出一个空列） |
| 分隔符 | 制表符 `\t` |
| 数据源 | CoCoPUTs（HIVE-CUTs 持续更新版），FDA DNAHIVE：<https://dnahive.fda.gov> |

## 字段说明

### 一、元数据字段（第 1–8 列）

与 `o537-Human_CDS.tsv` 完全相同（逐行比对一致），详见该文件的字段说明文档。摘要：Gene ID（基因符号）、Protein ID（`NP_`/`XP_`/`YP_`，574 行为空）、Accession（contig，不带版本号）、Division（`refseq`）、Assembly（`GCF_000001405.39`，GRCh38.p13）、Taxid（9606）、Species（`Homo sapiens`）、Organelle（`genomic` / `mitochondrion`）。

### 二、总量与 GC 统计字段（第 9–13 列）

| 列号 | 字段 | 类型 | 含义 | 示例 |
|---|---|---|---|---|
| 9 | # Junction Dinucleotides | int | 交界二核苷酸总数。**已验证恒等于 #Codons − 1**（每对相邻密码子产生 1 个交界对；306 个密码子 → 305 个），全表 119,196 行零例外 | `305` |
| 10–13 | GC% / GC1% / GC2% / GC3% | float | 与 `o537-Human_CDS.tsv` 完全一致（三表逐行比对零差异） | `42.81` 等 |

### 三、交界二核苷酸计数列（第 14–29 列，共 16 列）

每列为对应**交界二核苷酸**（DNA 字母表，5'→3' 有序对：`[上一密码子第 3 位][下一密码子第 1 位]`）在该 CDS 中的出现次数（整数原始计数）。

这一统计与**密码子对使用偏向（codon-pair bias, CPB）**直接对应：同义密码子重排 CDS 时，密码子内部二核苷酸不变，但交界二核苷酸组成会改变——这正是 codon-pair 优化（如减慢/加快翻译延伸速率、调节 mRNA 稳定性）的主要作用面。

列顺序与 `o537-Human_CDS_Dinuc.tsv` 相同，按「第一位碱基（T→C→A→G）× 第二位碱基（T→C→A→G）」排列：

| 列号 | 二核苷酸(DNA) | RNA 记法 | 列号 | 二核苷酸(DNA) | RNA 记法 |
|---|---|---|---|---|---|
| 14 | TT | UU | 22 | AT | AU |
| 15 | TC | UC | 23 | AC | AC |
| 16 | TA | UA（UpA，稳定性相关） | 24 | AA | AA |
| 17 | TG | UG | 25 | AG | AG |
| 18 | CT | CU | 26 | GT | GU |
| 19 | CC | CC | 27 | GC | GC（CpG，免疫相关） |
| 20 | CA | CA | 28 | GA | GA |
| 21 | CG | CG（CpG） | 29 | GG | GG |

## 已验证的数据性质（全表校验，2026-09-22）

- 16 列计数之和 = `# Junction Dinucleotides`：119,196 行零例外。
- `# Junction Dinucleotides = #Codons − 1`：全表零例外。
- 与另两表同记录、同行序；GC 四列逐行一致。
- **跨表闭合校验**（1,886,432 项）：`Dinuc 表计数[XY] = 本表[XY] + 密码子内部对[XY]`（密码子内部对由 CDS 表密码子计数精确推导，扣除最后一个密码子的 2 个内部对）在 99.99% 的项上精确成立；极少数例外为含内部终止密码子的硒蛋白类注释行。三张表计数口径完全自洽。

## 使用注意事项

1. **DNA 字母表**（T），本项目统一 RNA（U），注意转换。
2. **原始计数**：归一化时除以 `# Junction Dinucleotides`（或 `# Codons − 1`）。
3. 本表是三表中与 **codon-pair bias / 交界序列工程** 最直接相关的一份，可用作：
   - 密码子对优化（CPB）打分的参考分布（人类天然交界二核苷酸频率；完整密码子对计数见 `o537-Human_CDS_Bicod.tsv`，本表可由其精确推导，已验证）；
   - 约束解码中交界处 UpA/CpG 控制的先验统计；
   - 与 Dinuc 表相减可还原"密码子内部二核苷酸"组成（两种成分的调控意义不同）。
4. 行尾多余制表符会多读出一个空列；过滤与去重建议同 `o537-Human_CDS.tsv`（线粒体、重复行、XP_ 分层）。

## 读取示例

```python
import pandas as pd

cds  = pd.read_csv("data/o537-Human_CDS.tsv", sep="\t").dropna(axis=1, how="all")
junc = pd.read_csv("data/o537-Human_CDS_Junc_Dinuc.tsv", sep="\t").dropna(axis=1, how="all")
# 两表行序一致，可直接按列拼接
assert (cds["Gene ID"] == junc["Gene ID"]).all() and (cds["Protein ID"] == junc["Protein ID"]).all()
junc_cols = junc.columns[13:]     # 16 个交界二核苷酸计数
junc_freq = junc[junc_cols].div(junc["# Junction Dinucleotides"], axis=0)
```

## 来源与引用

- CoCoPUTs，FDA DNAHIVE：<https://dnahive.fda.gov>（"Available Files to Download" 页附 readme.txt）
- Athey J. et al. *A new and updated resource for codon usage tables.* Nucleic Acids Research, 2017（HIVE-CUTs）
- Alexaki A. et al. *Codon and Codon-Pair Usage Tables (CoCoPUTs).* Nucleic Acids Research, 2019
