"""pileup 累加器:模拟"比对"后的逐位点碱基/质量/Indel 计数。

演示采用理想化比对(reads 的参考坐标由模拟器直接给出,真实流程
用 bowtie2/minimap2 做序列比对),但 pileup 统计口径与
samtools mpileup 一致:
- depth[pos] = 覆盖该位的 reads 数 = 碱基计数(匹配+错配) + 跨位缺失(*)
- 插入(+A..)挂靠前一位点计数
"""

from __future__ import annotations

from collections import Counter

import numpy as np

FloatArray = np.ndarray  # type: ignore[type-arg,misc]



class Pileup:
    """参考基因组逐位点计数。counts/qual_sums 按 (pos*4+base) 平铺存储。"""

    def __init__(self, ref_len: int) -> None:
        self.ref_len = ref_len
        self.counts = [0] * (ref_len * 4)
        self.qual_sums = [0.0] * (ref_len * 4)
        self.del_count = [0] * ref_len
        self.ins_count = [0] * ref_len
        self.ins_seqs: list[dict[str, int]] = [{} for _ in range(ref_len)]

    def add_base(self, ref_pos: int, base_idx: int, qual: int) -> None:
        k = ref_pos * 4 + base_idx
        self.counts[k] += 1
        self.qual_sums[k] += qual

    def add_del_span(self, ref_pos: int, length: int) -> None:
        for p in range(ref_pos, min(ref_pos + length, self.ref_len)):
            self.del_count[p] += 1

    def add_insertion(self, ref_pos: int, seq: str) -> None:
        self.ins_count[ref_pos] += 1
        bucket = self.ins_seqs[ref_pos]
        bucket[seq] = bucket.get(seq, 0) + 1

    # ---------------- 汇总(numpy) ----------------

    def counts_array(self) -> FloatArray:
        return np.array(self.counts, dtype=np.int64).reshape(self.ref_len, 4)

    def qual_mean_array(self) -> FloatArray:
        counts = np.maximum(self.counts_array(), 1)
        quals = np.array(self.qual_sums, dtype=np.float64).reshape(self.ref_len, 4)
        return quals / counts

    def depth_array(self) -> FloatArray:
        """depth = 碱基计数和 + 跨位缺失计数(与 mpileup 口径一致)。"""
        return self.counts_array().sum(axis=1) + np.array(
            self.del_count, dtype=np.int64
        )

    def top_insertion(self, ref_pos: int) -> tuple[str, int]:
        """该位点最高频的插入串及其次数。"""
        counter = Counter(self.ins_seqs[ref_pos])
        if not counter:
            return "", 0
        seq, n = counter.most_common(1)[0]
        return seq, n
