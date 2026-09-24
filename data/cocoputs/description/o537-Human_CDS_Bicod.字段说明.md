# o537-Human_CDS_Bicod.tsv 字段说明

> 人类 RefSeq CDS 的**密码子对（bicodon）组成统计表**（来源 CoCoPUTs / HIVE-CUTs）。每行对应一条蛋白编码序列（CDS），统计全部 4,096 种相邻密码子对的出现次数。与 `o537-Human_CDS.tsv` **同记录、同行序**（同一批 119,196 条 CDS）。
>
> 本文档中的"已验证"结论均基于对全表 119,196 行的实际校验（2026-09-22），非仅凭字段名推测。

## 文件基本信息

| 属性 | 值 |
|---|---|
| 路径 | `data/o537-Human_CDS_Bicod.tsv` |
| 大小 | 989,716,154 字节（约 944 MB，四表中最大，读取方式见文末） |
| 行数 | 119,197 行 = 1 行表头 + 119,196 行数据 |
| 列数 | 4,110 = 8 元数据 + 1 总数 + 4 GC + **4,096 密码子对计数** + 1 个行尾制表符产生的空列 |
| 分隔符 | 制表符 `\t` |
| 数据源 | CoCoPUTs（HIVE-CUTs 持续更新版），FDA DNAHIVE：<https://dnahive.fda.gov> |

## 字段说明

### 一、元数据字段（第 1–8 列）

与 `o537-Human_CDS.tsv` 完全相同（逐行比对一致），详见该文件的字段说明文档。摘要：Gene ID（基因符号）、Protein ID（`NP_`/`XP_`/`YP_`，574 行为空）、Accession（contig，不带版本号）、Division（`refseq`）、Assembly（`GCF_000001405.39`，GRCh38.p13）、Taxid（9606）、Species（`Homo sapiens`）、Organelle（`genomic` / `mitochondrion`）。

### 二、总量与 GC 统计字段（第 9–13 列）

| 列号 | 字段 | 类型 | 含义 | 示例 |
|---|---|---|---|---|
| 9 | # Codon Pairs | int | 相邻密码子对总数。**已验证恒等于 #Codons − 1**（CDS 内每对框内相邻密码子 (i, i+1) 计 1 对；306 个密码子 → 305 对），全表 119,196 行零例外 | `305` |
| 10–13 | GC% / GC1% / GC2% / GC3% | float | 与 `o537-Human_CDS.tsv` 完全一致（逐行比对零差异） | `42.81` 等 |

### 三、密码子对计数列（第 14–4,109 列，共 4,096 列）

- 列名为**小写六联体 `abcdef`**（注意：本表小写，`o537-Human_CDS.tsv` 密码子列大写），表示框内相邻的两个密码子拼接：`abc` = 第 i 个密码子，`def` = 第 i+1 个密码子（DNA 字母表，5'→3'）。值为该密码子对在该 CDS 中的出现次数（整数**原始计数**）。
- 列排列规则（已程序化验证全部 4,096 列）：第 k 列（k = 1…4,096，对应文件第 13+k 列）= `CODONS[⌊(k−1)/64⌋] + CODONS[(k−1) mod 64]`，其中 `CODONS` 即 `o537-Human_CDS.tsv` 第 14–77 列的 64 密码子顺序（外层 = 第一个密码子，内层 = 第二个密码子，同序循环）。
  - 前 16 列：`tttttt … tttgtg`（第一个密码子固定 ttt，第二个遍历前 16 个）；
  - 最后 4 列：`gggggt、gggggc、ggggga、gggggg`。
- **终止密码子参与配对**：CDS 末位的终止密码子作为第二个密码子计入最后一个密码子对（如 `…|taa`）；硒蛋白类含内部 TGA 的注释行同样按实际序列计数。
- 与交界二核苷酸的精确关系（已验证）：每个 bicodon `abc|def` 的第 3、4 位 `c·d` 即一个交界二核苷酸，因此 `Junc[XY] = Σ（第 3 位=X 且第 4 位=Y 的 bicodon 计数）`——Junc 表可由本表完全推导。

## 已验证的数据性质（全表校验，2026-09-22）

- 4,096 列计数之和 = `# Codon Pairs`：119,196 行零例外。
- `# Codon Pairs = #Codons − 1`（与 `o537-Human_CDS.tsv` 逐行比对）：零例外。
- 与另外三表同记录、同行序（Gene ID + Protein ID 逐行零错位）；GC 四列逐行一致。
- **跨表闭合校验**：抽样 2,384 行 × 16 类型 = 38,144 项，由 bicodon 第 3+4 位推导的交界二核苷酸计数与 `o537-Human_CDS_Junc_Dinuc.tsv` 完全一致，零差异——证实"密码子对"就是框内相邻密码子 (i, i+1)，四张表计数口径完全自洽。
- 列名 6-mer 规则与排列顺序：4,096 列逐一程序化验证通过。

## 使用注意事项

1. **文件近 1 GB，按需读取**：整表载入内存（int64 × 4,096 列 × 119,196 行 ≈ 3.9 GB）通常不必要——计数极稀疏（每行非零值 ≈ #CodonPairs ≪ 4,096），多数分析只需部分列或元数据（见读取示例）。
2. **DNA 字母表（T）**，本项目统一 RNA（U），注意转换；列名为小写，匹配列名时勿与 CDS 表大写混淆。
3. **原始计数**：归一化时除以 `# Codon Pairs`（或 `# Codons − 1`）得密码子对频率；做 codon-pair bias（CPB）打分时通常还要对"期望频率"（由单密码子频率估计）取比值/对数。
4. **终止密码子对**：建模或统计人类密码子对先验时，按需剥离以终止密码子为第二个元素的对（`taa/tag/tga` 结尾的 3×64=192 列中的相关部分）。
5. 过滤与去重建议同 `o537-Human_CDS.tsv`（剔除线粒体 13 条、按 (Gene ID, Protein ID, Accession) 去重 6,358 行重复、XP_ 预测质量分层）。
6. 行尾多余制表符会多读出一个空列。

## 读取示例

```python
import pandas as pd

# 1) 只读元数据 + 总量（不全表载入）
meta = pd.read_csv("data/o537-Human_CDS_Bicod.tsv", sep="\t", usecols=range(13))

# 2) 需要全量计数时：用最小整数类型省内存（uint16 单元格上限 65,535，足够覆盖单条 CDS 的对计数）
df = pd.read_csv("data/o537-Human_CDS_Bicod.tsv", sep="\t",
                 dtype={c: "uint16" for c in pd.read_csv(
                     "data/o537-Human_CDS_Bicod.tsv", sep="\t", nrows=0).columns[13:]})
df = df.dropna(axis=1, how="all")   # 去掉行尾多余制表符产生的空列
pair_cols = df.columns[13:]         # 4,096 个密码子对计数（小写 6-mer，abc|def）
```

## 与其他三表的关系

四表同记录、同行序（按行号可直接拼接）。本表是信息最全的一份：单密码子计数（CDS 表）、交界二核苷酸（Junc_Dinuc 表）均可由本表结合/对照推导（交界二核苷酸已验证可精确推导）。密码子对统计直接对应**密码子对使用偏向（codon-pair bias）**，是同义密码子重排中改变翻译延伸速率与 mRNA 稳定性的主要作用面。

## 来源与引用

- CoCoPUTs，FDA DNAHIVE：<https://dnahive.fda.gov>（"Available Files to Download" 页附 readme.txt，为各表官方格式定义）
- Athey J. et al. *A new and updated resource for codon usage tables.* Nucleic Acids Research, 2017（HIVE-CUTs）
- Alexaki A. et al. *Codon and Codon-Pair Usage Tables (CoCoPUTs): facilitating genetic variation analyses and standardization.* Nucleic Acids Research, 2019
- 对应设计文档：`mRNA模型设计.md` 公开数据源盘点中"密码子/密码子对特征：CoCoPUTs、HIVE-CUTs"
