"""FastQC 风格的质控统计(在模拟 reads 的单遍扫描中流式收集)。

对应真实流程中的 FastQC/MultiQC 报告:per-cycle 质量分布、
每条 read 的 GC 含量分布、碱基总量。
"""

from __future__ import annotations

import numpy as np

# 每 cycle 保留的质量值样本上限(随机抽稀,画分位数足够)
_MAX_CYCLE_SAMPLES = 4000


class QcCollector:
    """单遍流式收集质控指标。"""

    def __init__(self, read_len: int) -> None:
        self.n_reads = 0
        self.total_bases = 0
        self.gc_hist = [0] * 101
        self.cycle_quality: list[list[int]] = [[] for _ in range(read_len)]
        self._cycle_seen = [0] * read_len

    def add_read(self, quals: list[int], gc_count: int, n_bases: int) -> None:
        self.n_reads += 1
        self.total_bases += n_bases
        self.gc_hist[min(100, round(gc_count * 100 / max(1, n_bases)))] += 1
        for cycle, q in enumerate(quals):
            seen = self._cycle_seen[cycle]
            if len(self.cycle_quality[cycle]) < _MAX_CYCLE_SAMPLES:
                self.cycle_quality[cycle].append(q)
            elif seen % 97 == 0:  # 已满后随机替换,保持近似均匀
                slot = (seen * 31 + 7) % _MAX_CYCLE_SAMPLES
                self.cycle_quality[cycle][slot] = q
            self._cycle_seen[cycle] = seen + 1

    def report(self) -> dict[str, object]:
        """输出 JSON 友好的报告结构(per-cycle 质量分位数等)。"""
        q25: list[float] = []
        q50: list[float] = []
        q75: list[float] = []
        for samples in self.cycle_quality:
            if not samples:
                continue
            s = np.array(samples, dtype=np.float64)
            q25.append(float(np.percentile(s, 25)))
            q50.append(float(np.percentile(s, 50)))
            q75.append(float(np.percentile(s, 75)))
        gc_total = sum(self.gc_hist) or 1
        return {
            "n_reads": self.n_reads,
            "total_bases": self.total_bases,
            "read_len": len(self.cycle_quality),
            "mean_gc": float(
                sum(i * n for i, n in enumerate(self.gc_hist))
                / (100 * gc_total)
            ),
            "gc_hist": self.gc_hist,
            "cycle_q25": q25,
            "cycle_q50": q50,
            "cycle_q75": q75,
        }
