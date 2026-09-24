# o537-Human_CDS_Dinuc.tsv 字段说明

> 人类 RefSeq CDS 的**二核苷酸组成统计表**（来源 CoCoPUTs / HIVE-CUTs）。每行对应一条蛋白编码序列（CDS），统计 16 种相邻二核苷酸的出现次数。与 `o537-Human_CDS.tsv` **同记录、同行序**（同一批 119,196 条 CDS）。
>
> 本文档中的"已验证"结论均基于对全表 119,196 行的实际校验（2026-09-22），非仅凭字段名推测。

## 文件基本信息

| 属性 | 值 |
|---|---|
| 路径 | `data/o537-Human_CDS_Dinuc.tsv` |
| 大小 | 19,679,475 字节（约 18.8 MB） |
| 行数 | 119,197 行 = 1 行表头 + 119,196 行数据 |
| 列数 | 29（数据行行尾多带一个制表符，解析时会多出一个空列） |
| 分隔符 | 制表符 `\t` |
| 数据源 | CoCoPUTs（HIVE-CUTs 持续更新版），FDA DNAHIVE：<https://dnahive.fda.gov> |

## 字段说明

### 一、元数据字段（第 1–8 列）

与 `o537-Human_CDS.tsv` 完全相同（逐行比对一致），详见该文件的字段说明文档。摘要：

| 列号 | 字段 | 含义 |
|---|---|---|
| 1 | Gene ID | RefSeq 基因符号（symbol 或 `LOC*` 编号） |
| 2 | Protein ID | RefSeq 蛋白 accession.version（`NP_`/`XP_`/`YP_`，574 行为空） |
| 3 | Accession | 所在 contig accession（不带版本号；`NC_*` 染色体，`NT_*` scaffold） |
| 4 | Division | 恒为 `refseq` |
| 5 | Assembly | 恒为 `GCF_000001405.39`（GRCh38.p13） |
| 6 | Taxid | 恒为 9606 |
| 7 | Species | 恒为 `Homo sapiens` |
| 8 | Organelle | `genomic`（核基因组）/ `mitochondrion`（线粒体，13 条） |

### 二、总量与 GC 统计字段（第 9–13 列）

| 列号 | 字段 | 类型 | 含义 | 示例 |
|---|---|---|---|---|
| 9 | # Dinucleotides | int | 二核苷酸计数总数。**已验证恒等于 3 × #Codons − 3**（即计数范围为序列上起始位置 1 … 3N−3 的相邻碱基对，最后一个密码子内部的 2 个相邻对不计入），全表 119,196 行零例外 | `915`（306 密码子的 CDS） |
| 10–13 | GC% / GC1% / GC2% / GC3% | float | 与 `o537-Human_CDS.tsv` 完全一致（三表逐行比对零差异） | `42.81` 等 |

### 三、二核苷酸计数列（第 14–29 列，共 16 列）

每列为对应**二核苷酸**（DNA 字母表，5'→3' 有序对）在该 CDS 中的出现次数（整数原始计数）。计数口径（已通过跨表闭合校验验证）：

```text
Dinuc[XY] = Junc[XY]（密码子交界对，见 Junc_Dinuc 表）
          + 密码子内部对（每个密码子贡献第 1–2 位、第 2–3 位两个对）
          − 最后一个密码子的 2 个内部对（不计）
```

即：**序列上除最后一个密码子内部 2 对以外的全部相邻碱基对**（含跨密码子边界的对）。列顺序按「第一位碱基（T→C→A→G）× 第二位碱基（T→C→A→G）」排列，与密码子表的"前两位"顺序一致。

| 列号 | 二核苷酸(DNA) | RNA 记法 | 备注 | 列号 | 二核苷酸(DNA) | RNA 记法 | 备注 |
|---|---|---|---|---|---|---|---|
| 14 | TT | UU | | 22 | AT | AU | |
| 15 | TC | UC | | 23 | AC | AC | |
| 16 | TA | UA | **UpA**：含量与 RNA 降解速率正相关，mRNA 稳定性设计常用特征 | 24 | AA | AA | |
| 17 | TG | UG | | 25 | AG | AG | |
| 18 | CT | CU | | 26 | GT | GU | |
| 19 | CC | CC | | 27 | GC | GC | **CpG**：与先天免疫刺激（如 TLR9）相关，mRNA 免疫原性设计常用特征 |
| 20 | CA | CA | | 28 | GA | GA | |
| 21 | CG | CG | CpG，见上 | 29 | GG | GG | |

注：RNA 记法中 TA→UA、TG→UG（T 换 U）；"UpA/CpG" 为文献常用称呼。

## 已验证的数据性质（全表校验，2026-09-22）

- 16 列计数之和 = `# Dinucleotides`：119,196 行零例外。
- `# Dinucleotides = 3 × #Codons − 3`：全表零例外。
- 与另三表同记录、同行序；GC 四列逐行一致。
- **跨表闭合校验**（1,886,432 项）：`Dinuc[XY] = Junc[XY] + 密码子内部对[XY]`（密码子内部对由 CDS 表的密码子计数精确推导）在 99.99% 的项上精确成立；极少数例外为含内部终止密码子的硒蛋白类注释行。此校验确认三张表计数口径完全自洽。

## 使用注意事项

1. **DNA 字母表**（T），本项目统一 RNA（U），注意转换。
2. **原始计数**：跨基因比较需除以 `# Dinucleotides`（或 `# Codons × 3`）归一化为频率。
3. mRNA 设计中最常用的是 **TA（UpA）与 CG（CpG）** 两列：UpA 偏低通常利于稳定性，CpG 含量影响免疫刺激特性；设计文档中免疫/稳定性打分可将二核苷酸频率作为输入特征。
4. 行尾多余制表符会多读出一个空列；过滤与去重建议同 `o537-Human_CDS.tsv`（线粒体、重复行、XP_ 分层）。

## 读取示例

```python
import pandas as pd

df = pd.read_csv("data/o537-Human_CDS_Dinuc.tsv", sep="\t").dropna(axis=1, how="all")
dinuc_cols = df.columns[13:]   # 16 个二核苷酸计数，顺序 TT,TC,TA,TG,CT,CC,CA,CG,AT,AC,AA,AG,GT,GC,GA,GG
freq = df[dinuc_cols].div(df["# Dinucleotides"], axis=0)   # 归一化为二核苷酸频率
```

## 来源与引用

- CoCoPUTs，FDA DNAHIVE：<https://dnahive.fda.gov>（"Available Files to Download" 页附 readme.txt）
- Athey J. et al. *A new and updated resource for codon usage tables.* Nucleic Acids Research, 2017（HIVE-CUTs）
- Alexaki A. et al. *Codon and Codon-Pair Usage Tables (CoCoPUTs).* Nucleic Acids Research, 2019
