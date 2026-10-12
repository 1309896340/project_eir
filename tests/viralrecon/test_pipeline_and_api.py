"""端到端管线与 FastAPI 接口测试(纯模拟,秒级)。"""

from __future__ import annotations

import unittest
from typing import Any

from fastapi.testclient import TestClient

from viralrecon import app as app_module
from viralrecon.core.pipeline import run_pipeline
from viralrecon.core.simulate import SimParams


def _uniform_params(**overrides: Any) -> SimParams:
    """uniform 覆盖的小参数:深度均匀,真值变异应全部可检出。"""
    values: dict[str, Any] = dict(
        genome_len=2000, coverage=60, read_len=100,
        n_snps=10, n_indels=4, coverage_mode="uniform", seed=7,
    )
    values.update(overrides)
    return SimParams(**values)


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = run_pipeline(_uniform_params())

    def test_report_structure_complete(self) -> None:
        r = self.report
        for key in ("params", "timings", "qc", "coverage", "variants", "consensus"):
            self.assertIn(key, r)
        self.assertIn("fasta", r["consensus"])
        self.assertIn("depth", r["coverage"])

    def test_qc_accounting(self) -> None:
        qc = self.report["qc"]
        self.assertEqual(qc["n_reads"], round(2000 * 60 / 100))
        self.assertEqual(qc["total_bases"], qc["n_reads"] * 100)
        self.assertEqual(len(qc["cycle_q50"]), 100)
        # Illumina 前高后低:首 cycle 中位质量应高于末 cycle
        self.assertGreater(qc["cycle_q50"][0], qc["cycle_q50"][-1])

    def test_coverage_matches_target(self) -> None:
        cov = self.report["coverage"]
        self.assertAlmostEqual(cov["mean"], 60.0, delta=12.0)
        # 基因组两端 read 无法完整覆盖,深度自然衰减(与真实比对一致),
        # 故全域阈值放至 0.98;内部区域应为高覆盖
        self.assertGreater(cov["breadth"], 0.98)
        self.assertGreater(cov["pct_ge_10x"], 0.98)
        self.assertLessEqual(len(cov["depth"]), 4000)

    def test_all_truth_variants_recovered(self) -> None:
        """均匀覆盖 + 充足深度:注入变异应全部 PASS 检出,且无假阳。"""
        v = self.report["variants"]
        self.assertEqual(v["truth_total"], 14)
        self.assertEqual(v["tp"], 14)
        self.assertEqual(v["fn"], 0)
        self.assertEqual(v["fp"], 0)

    def test_consensus_quality(self) -> None:
        cons = self.report["consensus"]
        n_pct = cons["n_count"] / cons["length"]
        self.assertLess(n_pct, 0.05)
        self.assertGreaterEqual(cons["diffs_vs_ref"], 10)  # 至少包含注入 SNP

    def test_seed_reproducible(self) -> None:
        other = run_pipeline(_uniform_params())
        self.assertEqual(
            self.report["consensus"]["fasta"], other["consensus"]["fasta"]
        )
        self.assertEqual(
            self.report["variants"]["rows"], other["variants"]["rows"]
        )


class ApiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.client = TestClient(app_module.app)

    def test_presets_and_pipeline_meta(self) -> None:
        res = self.client.get("/api/presets")
        self.assertEqual(res.status_code, 200)
        keys = {p["key"] for p in res.json()["presets"]}
        self.assertEqual(keys, {"quick", "standard", "coronavirus"})
        steps = self.client.get("/api/pipeline").json()
        self.assertEqual(len(steps), 5)

    def test_run_quick(self) -> None:
        res = self.client.post("/api/run", json={"preset": "quick"})
        self.assertEqual(res.status_code, 200)
        r = res.json()
        self.assertGreater(r["qc"]["n_reads"], 0)
        self.assertTrue(r["consensus"]["fasta"].startswith(">"))
        self.assertGreaterEqual(r["variants"]["truth_total"], 10)

    def test_run_with_overrides_and_seed(self) -> None:
        res = self.client.post("/api/run", json={
            "preset": "quick",
            "params": {"coverage": 30, "coverage_mode": "uniform"},
            "seed": 123,
        })
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.json()["params"]["seed"], 123)
        self.assertEqual(res.json()["params"]["coverage"], 30)

    def test_run_rejects_unknown_param(self) -> None:
        res = self.client.post("/api/run", json={
            "preset": "quick", "params": {"nonexistent": 1},
        })
        self.assertEqual(res.status_code, 400)

    def test_run_rejects_budget_overflow(self) -> None:
        res = self.client.post("/api/run", json={
            "preset": "coronavirus", "params": {"coverage": 200},
        })
        self.assertEqual(res.status_code, 400)
        self.assertIn("预算", res.json()["detail"])

    def test_run_clamps_out_of_range(self) -> None:
        res = self.client.post("/api/run", json={
            "preset": "quick", "params": {"coverage": 99999},
        })
        # clamp 到上限 200:2000x200 = 40 万碱基 < 预算,应正常完成
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.json()["params"]["coverage"], 200)


if __name__ == "__main__":
    unittest.main()
