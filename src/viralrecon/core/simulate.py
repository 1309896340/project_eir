"""模拟数据生成:合成参考基因组 + 样本真值变异 + Illumina 风格 reads。

本模块是演示 app 的"测序数据源":不使用任何真实病原体序列,
参考基因组按冠状病毒样 GC 含量(~38%)随机合成,变异与 reads 均
由种子驱动、可复现。reads 的质量值遵循 Illumina 特征(前端高、
末端线性衰减 + 少量低质量 reads),测序错误概率与 Phred 质量自洽
(p_err = 10^(-Q/10)),因此质控图与变异检出频率相互印证。
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

BASES = "ACGT"
BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
# 冠状病毒样组成:AT 偏多,目标 GC ~38%
_BASE_WEIGHTS = (0.31, 0.19, 0.19, 0.31)

# 变异位点之间的最小间隔(且离基因组两端各留 100bp),
# 避免同一 read 内多个变异叠加导致坐标解释复杂化
_VARIANT_SPACING = 60
_EDGE_MARGIN = 100


@dataclass(frozen=True)
class SimParams:
    """一次模拟运行的全部参数。"""

    genome_len: int = 5000
    coverage: int = 100            # 目标平均测序深度(x)
    read_len: int = 150
    error_rate_scale: float = 1.0  # 测序错误率整体缩放(1.0 即严格按 Q 值)
    n_snps: int = 20
    n_indels: int = 6
    indel_max_len: int = 5
    low_quality_frac: float = 0.08  # 低质量 reads 占比(Q 整体 -12)
    coverage_mode: str = "amplicon"  # uniform | amplicon
    q_start: float = 38.0          # read 首端基线质量
    q_end: float = 32.0            # read 末端基线质量
    seed: int = 42


@dataclass(frozen=True)
class TruthVariant:
    """样本相对参考的"真值"变异(模拟注入,用于评估检出效果)。"""

    pos: int      # 0-based 参考坐标
    ref: str      # 参考等位:snp 单碱基;del 为被删序列;ins 为空串
    alt: str      # 替代等位:snp 单碱基;ins 为插入序列;del 为空串
    kind: str     # "snp" | "ins" | "del"


@dataclass
class SampleBuild:
    """参考基因组 + 注入变异后的样本基因组及坐标映射。"""

    reference: str
    sample_genome: str
    truth: list[TruthVariant]
    # sample_genome 每个碱基对应的参考坐标;插入碱基挂靠前一位,
    # 因此该数组仅在插入处出现相邻相等(其余严格递增)
    ref_of_sample: list[int]
    # 删除接缝:sample 坐标系中 del 之后的首个位置 -> 变异(按坐标排序)
    del_seams: list[tuple[int, TruthVariant]] = field(default_factory=list)


def generate_reference(length: int, rng: random.Random) -> str:
    """按目标碱基组成随机合成参考基因组(非真实病原体序列)。"""
    return "".join(
        rng.choices(BASES, weights=_BASE_WEIGHTS)[0] for _ in range(length)
    )


def _pick_positions(rng: random.Random, span: int, n: int) -> list[int]:
    """在留边区间内取 n 个满足最小间隔的位点(等距骨架 + 有界抖动)。

    等距骨架保证在容量内必然成功,抖动保留随机性且不破坏最小间距。
    """
    if n == 0:
        return []
    lo = _EDGE_MARGIN
    hi = span - _EDGE_MARGIN
    step = (hi - lo) / n
    if step < _VARIANT_SPACING:
        raise ValueError("基因组长度不足以容纳指定数量的间隔变异")
    jitter = max(0.0, (step - _VARIANT_SPACING) / 2)
    picked = [
        min(hi - 1, int(lo + step * k + rng.uniform(0, jitter)))
        for k in range(n)
    ]
    return sorted(set(picked)) if len(set(picked)) == n else picked


def _random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(BASES, weights=_BASE_WEIGHTS)[0]
                   for _ in range(length))


def build_sample(params: SimParams, rng: random.Random) -> SampleBuild:
    """向参考基因组注入 SNP/Indel,构建样本真值基因组与坐标映射。"""
    reference = generate_reference(params.genome_len, rng)
    positions = _pick_positions(rng, params.genome_len,
                                params.n_snps + params.n_indels)
    rng.shuffle(positions)
    snp_pos = sorted(positions[: params.n_snps])
    indel_pos = sorted(positions[params.n_snps:])

    snps: dict[int, TruthVariant] = {}
    for p in snp_pos:
        alt = rng.choice([b for b in BASES if b != reference[p]])
        snps[p] = TruthVariant(p, reference[p], alt, "snp")

    indels: dict[int, TruthVariant] = {}
    for p in indel_pos:
        if rng.random() < 0.5:
            length = rng.randint(1, params.indel_max_len)
            indels[p] = TruthVariant(p, reference[p: p + length], "", "del")
        else:
            length = rng.randint(1, params.indel_max_len)
            indels[p] = TruthVariant(p, "", _random_seq(rng, length), "ins")

    parts: list[str] = []
    ref_of_sample: list[int] = []
    seams: list[tuple[int, TruthVariant]] = []
    truth = [snps[p] for p in snp_pos] + [indels[p] for p in indel_pos]
    truth.sort(key=lambda v: v.pos)

    i = 0
    while i < len(reference):
        if i in snps:
            parts.append(snps[i].alt)
            ref_of_sample.append(i)
            i += 1
            continue
        if i in indels:
            v = indels[i]
            if v.kind == "ins":
                parts.append(reference[i])
                ref_of_sample.append(i)
                for ch in v.alt:
                    parts.append(ch)
                    ref_of_sample.append(i)  # 插入碱基挂靠锚点
                i += 1
                continue
            # del:锚点起的整段参考碱基从样本基因组中移除
            seams.append((len(ref_of_sample), v))
            i += len(v.ref)
            continue
        parts.append(reference[i])
        ref_of_sample.append(i)
        i += 1

    return SampleBuild(reference, "".join(parts), truth, ref_of_sample, seams)


def coverage_weights(params: SimParams, sample_len: int,
                     rng: random.Random) -> list[float]:
    """read 起点的采样权重。

    amplicon 模式模拟引物扩增的波浪形覆盖(周期 ~1200bp,每个扩增子
    随机效率),是 viralrecon 真实报告中最典型的覆盖度形态;
    uniform 模式为均匀覆盖。
    """
    n_pos = max(1, sample_len - params.read_len + 1)
    if params.coverage_mode == "uniform":
        return [1.0] * n_pos
    amp_len = 1200
    phase = rng.uniform(0, 2 * math.pi)
    weights: list[float] = []
    amp_eff = 1.0
    for s in range(n_pos):
        if s % amp_len == 0:
            amp_eff = min(1.6, max(0.4, rng.gauss(1.0, 0.3)))
        wave = 0.30 + 0.70 * 0.5 * (1 + math.sin(
            2 * math.pi * s / amp_len + phase))
        weights.append(wave * amp_eff)
    return weights


def sample_read_starts(params: SimParams, sample_len: int,
                       rng: random.Random) -> list[int]:
    """按覆盖度权重采样全部 read 起点(样本坐标系)。"""
    n_reads = max(1, round(params.genome_len * params.coverage / params.read_len))
    weights = coverage_weights(params, sample_len, rng)
    return rng.choices(range(len(weights)), weights=weights, k=n_reads)


__all__ = [
    "BASES",
    "BASE_IDX",
    "SimParams",
    "TruthVariant",
    "SampleBuild",
    "generate_reference",
    "build_sample",
    "coverage_weights",
    "sample_read_starts",
]
