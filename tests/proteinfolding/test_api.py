"""ProteinFolding FastAPI 后端测试(不启动真实服务,不跑真实模拟)。"""

from __future__ import annotations

import sys
import tempfile
import unittest
from typing import Any
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).parents[2]))

from fastapi.testclient import TestClient  # noqa: E402

from ProteinFolding import app as app_module  # noqa: E402

FIXTURE = Path(__file__).parent / "fixtures" / "ubq20.pdb"
SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRT"


class ApiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.client = TestClient(app_module.app)
        app_module.JOBS.clear()

    def test_validate_ok(self) -> None:
        res = self.client.post("/api/validate", json={"sequence": SEQ})
        self.assertEqual(res.status_code, 200)
        data = res.json()
        self.assertTrue(data["ok"])
        self.assertLess(20, data["length"])

    def test_validate_rejects_short(self) -> None:
        res = self.client.post("/api/validate", json={"sequence": "MALWM"})
        self.assertEqual(res.status_code, 200)
        self.assertFalse(res.json()["ok"])

    def test_samples_and_presets(self) -> None:
        self.assertEqual(len(self.client.get("/api/samples").json()), 3)
        presets = self.client.get("/api/presets").json()
        self.assertEqual(set(presets), {"fast", "standard", "high"})

    def test_submit_rejects_bad_sequence(self) -> None:
        res = self.client.post("/api/jobs", json={"sequence": "MALWM"})
        self.assertEqual(res.status_code, 400)

    def test_submit_rejects_bad_params(self) -> None:
        res = self.client.post(
            "/api/jobs", json={"sequence": SEQ, "params": {"nope": 1}}
        )
        self.assertEqual(res.status_code, 400)

    def test_job_not_found(self) -> None:
        res = self.client.get("/api/jobs/deadbeef")
        self.assertEqual(res.status_code, 404)

    def test_result_file_rejects_traversal(self) -> None:
        res = self.client.get("/api/results/..%2F..%2Fsecret/native.pdb")
        self.assertIn(res.status_code, (400, 404))

    def test_job_lifecycle_with_mock_pipeline(self) -> None:
        """mock 掉 run_pipeline,验证 排队→运行→完成 的状态机与缓存直返。"""
        fake_result = {
            "key": "deadbeefcafebabe",
            "cached": False,
            "dir": "nowhere",
            "meta": {"rmsd_A": [30.0, 4.2], "platform": "CPU"},
            "predictor": "mock",
        }
        with mock.patch.object(
            app_module, "run_pipeline", return_value=fake_result
        ) as mocked:
            res = self.client.post(
                "/api/jobs", json={"sequence": SEQ, "preset": "standard"}
            )
            self.assertEqual(res.status_code, 200)
            job_id = res.json()["job_id"]
            status: dict[str, Any] = {"status": "queued"}
            for _ in range(100):
                status = self.client.get(f"/api/jobs/{job_id}").json()
                if status["status"] in ("done", "error"):
                    break
                import time

                time.sleep(0.05)
            self.assertEqual(status["status"], "done")
            result = status["result"]
            assert result is not None
            self.assertEqual(result["key"], "deadbeefcafebabe")
            mocked.assert_called_once()

    def test_cache_hit_short_circuit(self) -> None:
        """缓存目录齐备时提交直接返回 done,不进队列。"""
        with tempfile.TemporaryDirectory() as td:
            # 将缓存根指到临时目录
            with mock.patch.object(
                app_module, "find_cached", return_value=Path(td) / "fake"
            ):
                res = self.client.post(
                    "/api/jobs",
                    json={"sequence": SEQ, "preset": "standard"},
                )
            self.assertEqual(res.status_code, 200)
            job_id = res.json()["job_id"]
            status = self.client.get(f"/api/jobs/{job_id}").json()
            self.assertEqual(status["status"], "done")
            self.assertTrue(status["result"]["cached"])


if __name__ == "__main__":
    unittest.main()
