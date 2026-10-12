"""共有序列(consensus)调用:按 ivar consensus 的多数碱基逻辑。

每位置取众数碱基;支持 reads 的碱基比例 >= min_frac 才调用,
否则写 N(真实流程中 ivar 默认 0.75)。高频缺失使该位从共识中
移除,高频插入按众数插入串补回,因此共识长度可不同于参考。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .pileup import Pileup


@dataclass(frozen=True)
class ConsensusResult:
    seq: str
    n_count: int
    diffs_vs_ref: int   # 与参考的差异位点数(含 N 与缺失)
    length: int


def build_consensus(pileup: Pileup, reference: str,
                    min_frac: float = 0.75) -> ConsensusResult:
    counts = pileup.counts_array()
    dels = np.array(pileup.del_count, dtype=np.int64)
    ins = np.array(pileup.ins_count, dtype=np.int64)

    parts: list[str] = []
    n_count = 0
    diffs = 0
    for pos in range(len(reference)):
        depth = int(counts[pos].sum()) + int(dels[pos])
        if depth == 0 or int(dels[pos]) / depth >= min_frac:
            # 无覆盖写 N;高频缺失则从共识中移除该位(ivar 行为)
            if depth == 0:
                parts.append("N")
                n_count += 1
            diffs += 1
            continue
        base_counts = counts[pos]
        top_idx = int(np.argmax(base_counts))
        frac = int(base_counts[top_idx]) / depth
        if frac >= min_frac:
            base = "ACGT"[top_idx]
        else:
            base = "N"  # 有覆盖但碱基分歧,不足以调用
            n_count += 1
        if base != reference[pos]:
            diffs += 1
        parts.append(base)
        # 高频插入:按众数插入串补回
        if int(ins[pos]) / depth >= min_frac:
            seq, _ = pileup.top_insertion(pos)
            parts.append(seq)

    seq = "".join(parts)
    return ConsensusResult(seq, n_count, diffs, len(seq))


def to_fasta(name: str, seq: str, width: int = 60) -> str:
    lines = [f">{name}"]
    for i in range(0, len(seq), width):
        lines.append(seq[i: i + width])
    return "\n".join(lines) + "\n"
