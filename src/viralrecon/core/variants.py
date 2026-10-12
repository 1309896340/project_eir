"""变异检测:从 pileup 按 ivar 风格阈值调用 SNP/Indel 并过滤。

对应真实流程中的 `ivar variants`(低频变异过滤:depth/qual 阈值)。
过滤规则与 viralrecon 默认参数一致:
    depth < 10   -> LOW_DEPTH
    meanQ < 20   -> LOW_QUAL
    AF < 0.75    -> LOW_FREQ(ivar 共识/变异的默认频率阈值)
其余为 PASS。另与模拟注入的真值集对照,给出 TP/FP/FN。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .pileup import Pileup
from .simulate import TruthVariant

# 只有频率达到该阈值的位点才进入结果表(纯测序噪声远低于此,不展示)
_TABLE_MIN_FREQ = 0.05


@dataclass(frozen=True)
class VariantFilters:
    min_depth: int = 10
    min_qual: float = 20.0
    min_freq: float = 0.75


@dataclass(frozen=True)
class CalledVariant:
    pos: int          # 1-based 展示坐标
    ref: str
    alt: str
    kind: str         # snp | ins | del
    depth: int
    af: float
    qual: float
    filter_status: str


def call_variants(pileup: Pileup, reference: str,
                  filters: VariantFilters) -> list[CalledVariant]:
    """逐位点扫描 pileup,按阈值产出变异表(含未通过过滤的行)。"""
    counts = pileup.counts_array()
    qual_mean = pileup.qual_mean_array()
    qual_sums = np.array(pileup.qual_sums, dtype=np.float64).reshape(
        pileup.ref_len, 4
    )
    dels = np.array(pileup.del_count, dtype=np.int64)
    ins = np.array(pileup.ins_count, dtype=np.int64)
    called: list[CalledVariant] = []

    def position_qual(pos: int) -> float:
        """indel 位点质量:该位碱基计数加权均值;纯 del 位无任何碱基
        计数(reads 全部跨过)时取邻近覆盖位均值。"""
        for p in (pos, pos - 1, pos + 1, pos - 2, pos + 2):
            if 0 <= p < pileup.ref_len:
                n = int(counts[p].sum())
                if n > 0:
                    return float(qual_sums[p].sum() / n)
        return 0.0

    for pos in range(len(reference)):
        depth = int(counts[pos].sum() + dels[pos])
        if depth == 0:
            continue
        ref_idx = {"A": 0, "C": 1, "G": 2, "T": 3}[reference[pos]]

        # SNP/错配:非参考碱基中的最高频者
        alt_counts = counts[pos].copy()
        alt_counts[ref_idx] = 0
        alt_idx = int(np.argmax(alt_counts))
        alt_n = int(alt_counts[alt_idx])
        if alt_n > 0:
            af = alt_n / depth
            if af >= _TABLE_MIN_FREQ:
                called.append(_make_row(
                    pos, reference[pos], "ACGT"[alt_idx], "snp",
                    depth, af, float(qual_mean[pos][alt_idx]), filters,
                ))

        # 缺失:连续 del 段合并为一行(VCF 风格,段首行代表整段)
        if dels[pos] > 0 and dels[pos] / depth >= _TABLE_MIN_FREQ:
            prev_depth = int(counts[pos - 1].sum() + dels[pos - 1]) if pos > 0 else 0
            prev_frac = dels[pos - 1] / prev_depth if prev_depth else 0.0
            if prev_frac < _TABLE_MIN_FREQ:  # 仅在段首输出
                end = pos
                while end < len(reference):
                    d = int(counts[end].sum() + dels[end])
                    if d and dels[end] / d >= _TABLE_MIN_FREQ:
                        end += 1
                    else:
                        break
                called.append(_make_row(
                    pos, reference[pos:end], "", "del",
                    depth, dels[pos] / depth,
                    position_qual(pos), filters,
                ))

        # 插入:挂靠该位点的插入串
        if ins[pos] > 0 and ins[pos] / depth >= _TABLE_MIN_FREQ:
            seq, _ = pileup.top_insertion(pos)
            called.append(_make_row(
                pos, reference[pos], seq, "ins",
                depth, ins[pos] / depth,
                position_qual(pos), filters,
            ))

    return called


def _make_row(pos: int, ref: str, alt: str, kind: str, depth: int,
              af: float, qual: float, f: VariantFilters) -> CalledVariant:
    if depth < f.min_depth:
        status = "LOW_DEPTH"
    elif qual < f.min_qual:
        status = "LOW_QUAL"
    elif af < f.min_freq:
        status = "LOW_FREQ"
    else:
        status = "PASS"
    return CalledVariant(pos + 1, ref, alt, kind, depth, af, qual, status)


def score_against_truth(called: list[CalledVariant],
                        truth: list[TruthVariant]) -> dict[str, int]:
    """PASS 变异与真值对照:(pos, kind) 匹配即 TP。"""
    truth_keys = {(v.pos + 1, v.kind) for v in truth}
    pass_keys = {(c.pos, c.kind) for c in called if c.filter_status == "PASS"}
    tp = len(truth_keys & pass_keys)
    return {
        "truth_total": len(truth_keys),
        "tp": tp,
        "fn": len(truth_keys - pass_keys),
        "fp": len(pass_keys - truth_keys),
    }
