"""变异检测与共识调用单元测试:手工构造 pileup 验证阈值逻辑。"""

from __future__ import annotations

import unittest

from viralrecon.core.consensus import build_consensus
from viralrecon.core.pileup import Pileup
from viralrecon.core.variants import VariantFilters, call_variants


def _seed_base_counts(pileup: Pileup, pos: int, base: str, n: int,
                      qual: int = 35) -> None:
    idx = "ACGT".index(base)
    for _ in range(n):
        pileup.add_base(pos, idx, qual)


class VariantCallTests(unittest.TestCase):
    def setUp(self) -> None:
        self.ref = "A" * 30
        self.filters = VariantFilters()

    def _rows(self, pileup: Pileup):
        return call_variants(pileup, self.ref, self.filters)

    def test_high_frequency_snp_pass(self) -> None:
        p = Pileup(30)
        _seed_base_counts(p, 10, "A", 3)
        _seed_base_counts(p, 10, "G", 17)  # AF=17/20=0.85
        rows = [r for r in self._rows(p) if r.kind == "snp"]
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual((r.ref, r.alt, r.filter_status), ("A", "G", "PASS"))
        self.assertEqual(r.pos, 11)  # 1-based
        self.assertAlmostEqual(r.af, 0.85, places=3)

    def test_low_depth_filtered(self) -> None:
        p = Pileup(30)
        _seed_base_counts(p, 5, "A", 1)
        _seed_base_counts(p, 5, "G", 4)  # depth=5 < 10
        rows = [r for r in self._rows(p) if r.kind == "snp"]
        self.assertEqual(rows[0].filter_status, "LOW_DEPTH")

    def test_low_qual_filtered(self) -> None:
        p = Pileup(30)
        _seed_base_counts(p, 5, "A", 2)
        _seed_base_counts(p, 5, "G", 18, qual=15)  # depth=20,qual<20
        rows = [r for r in self._rows(p) if r.kind == "snp"]
        self.assertEqual(rows[0].filter_status, "LOW_QUAL")

    def test_low_freq_filtered_and_noise_hidden(self) -> None:
        p = Pileup(30)
        _seed_base_counts(p, 5, "A", 10)
        _seed_base_counts(p, 5, "G", 10)  # AF=0.5:进表但 LOW_FREQ
        _seed_base_counts(p, 8, "A", 190)
        _seed_base_counts(p, 8, "G", 1)   # AF~0.005:纯噪声不进表
        rows = self._rows(p)
        freq_row = [r for r in rows if r.pos == 6]
        self.assertEqual(freq_row[0].filter_status, "LOW_FREQ")
        self.assertFalse(any(r.pos == 9 for r in rows))

    def test_deletion_called_as_single_segment(self) -> None:
        p = Pileup(30)
        _seed_base_counts(p, 10, "A", 4)
        _seed_base_counts(p, 15, "A", 2)
        for pos in (12, 13, 14):
            _seed_base_counts(p, pos, "A", 2)
        for _ in range(10):  # 10 个 reads 跨过 del 区间 12-14:frac=10/12
            p.add_del_span(12, 3)
        rows = [r for r in self._rows(p) if r.kind == "del"]
        self.assertEqual(len(rows), 1)  # 连续 del 段合并为一行
        r = rows[0]
        self.assertEqual(r.filter_status, "PASS")
        self.assertEqual(r.alt, "")
        self.assertEqual(r.ref, "AAA")  # 被删参考串

    def test_insertion_called_with_top_seq(self) -> None:
        p = Pileup(30)
        _seed_base_counts(p, 10, "A", 10)
        for _ in range(15):
            p.add_insertion(10, "TGG")
        rows = [r for r in self._rows(p) if r.kind == "ins"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].alt, "TGG")
        self.assertEqual(rows[0].filter_status, "PASS")


class ConsensusTests(unittest.TestCase):
    def test_majority_and_n(self) -> None:
        ref = "A" * 20
        p = Pileup(20)
        _seed_base_counts(p, 5, "A", 10)
        _seed_base_counts(p, 5, "G", 10)   # 分歧 50% -> N
        _seed_base_counts(p, 6, "C", 20)   # 100% C -> 调用
        _seed_base_counts(p, 9, "A", 9)    # AF=1.0 但 ref 是 A -> 同参考
        # 位点 0-4/7-8/10-19 无覆盖 -> N
        cons = build_consensus(p, ref)
        self.assertEqual(cons.seq[6], "C")
        self.assertEqual(cons.seq[5], "N")
        self.assertEqual(cons.seq[0], "N")
        self.assertGreater(cons.n_count, 0)

    def test_high_freq_deletion_removes_base(self) -> None:
        ref = "ACGT" * 5
        p = Pileup(20)
        for pos in range(20):
            _seed_base_counts(p, pos, ref[pos], 10)
        for _ in range(50):  # 跨过位点 8-9 的 reads 计缺失
            p.add_del_span(8, 2)
        cons = build_consensus(p, ref)
        # 高频缺失:共识移除位点 8-9(1-based 9-10)
        self.assertEqual(len(cons.seq), 18)
        self.assertNotIn("N", cons.seq)

    def test_high_freq_insertion_restored(self) -> None:
        ref = "A" * 10
        p = Pileup(10)
        for pos in range(10):
            _seed_base_counts(p, pos, "A", 10)
        for _ in range(30):
            p.add_insertion(4, "TT")
        cons = build_consensus(p, ref)
        self.assertEqual(cons.seq, "A" * 5 + "TT" + "A" * 5)


if __name__ == "__main__":
    unittest.main()
