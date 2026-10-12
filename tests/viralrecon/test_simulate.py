"""数据模拟层测试:基因组合成、变异注入、reads 采样。"""

from __future__ import annotations

import unittest

from viralrecon.core.simulate import (
    SimParams,
    build_sample,
    generate_reference,
    sample_read_starts,
)


class ReferenceTests(unittest.TestCase):
    def test_seed_reproducible(self) -> None:
        import random

        r1 = generate_reference(500, random.Random(7))
        r2 = generate_reference(500, random.Random(7))
        self.assertEqual(r1, r2)
        self.assertEqual(len(r1), 500)
        self.assertTrue(set(r1) <= set("ACGT"))

    def test_gc_content_in_target_range(self) -> None:
        import random

        ref = generate_reference(50000, random.Random(1))
        gc = sum(1 for b in ref if b in "GC") / len(ref)
        # 目标 ~38%,大样本下波动很小
        self.assertGreater(gc, 0.33)
        self.assertLess(gc, 0.43)


class SampleBuildTests(unittest.TestCase):
    def setUp(self) -> None:
        import random

        self.params = SimParams(
            genome_len=3000, coverage=40, read_len=100,
            n_snps=8, n_indels=4, seed=99, coverage_mode="uniform",
        )
        self.build = build_sample(self.params, random.Random(self.params.seed))

    def test_truth_counts(self) -> None:
        self.assertEqual(len(self.build.truth), 12)
        self.assertEqual(
            sum(1 for v in self.build.truth if v.kind == "snp"), 8
        )
        self.assertEqual(
            sum(1 for v in self.build.truth if v.kind in ("ins", "del")), 4
        )

    def test_variant_spacing(self) -> None:
        positions = sorted(v.pos for v in self.build.truth)
        for a, b in zip(positions, positions[1:]):
            self.assertGreaterEqual(b - a, 60)
        self.assertGreaterEqual(positions[0], 100)
        self.assertLess(positions[-1], 3000 - 100)

    def test_sample_genome_length_accounting(self) -> None:
        ins_bases = sum(len(v.alt) for v in self.build.truth if v.kind == "ins")
        del_bases = sum(len(v.ref) for v in self.build.truth if v.kind == "del")
        self.assertEqual(
            len(self.build.sample_genome),
            self.params.genome_len + ins_bases - del_bases,
        )
        self.assertEqual(len(self.build.ref_of_sample),
                         len(self.build.sample_genome))

    def test_reproducible(self) -> None:
        import random

        other = build_sample(self.params, random.Random(99))
        self.assertEqual(self.build.sample_genome, other.sample_genome)
        self.assertEqual(self.build.truth, other.truth)


class ReadStartsTests(unittest.TestCase):
    def test_read_count_matches_coverage(self) -> None:
        import random

        params = SimParams(genome_len=2000, coverage=50, read_len=100,
                           n_snps=2, n_indels=1, coverage_mode="uniform", seed=3)
        build = build_sample(params, random.Random(3))
        starts = sample_read_starts(params, len(build.sample_genome),
                                    random.Random(3))
        self.assertEqual(len(starts), round(2000 * 50 / 100))
        self.assertTrue(all(0 <= s < len(build.sample_genome) - 100 + 1
                            for s in starts))


if __name__ == "__main__":
    unittest.main()
