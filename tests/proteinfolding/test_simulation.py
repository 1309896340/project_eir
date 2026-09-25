"""ProteinFolding targeted MD 模拟集成测试(小体系,CPU,~1 分钟)。"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2]))

from ProteinFolding.core.config import SimParams, resolve_params  # noqa: E402
from ProteinFolding.core.linear import make_linear_pdb  # noqa: E402
from ProteinFolding.core.simulation import (  # noqa: E402
    choose_platform,
    run_targeted_md,
)

FIXTURE = Path(__file__).parent / "fixtures" / "ubq20.pdb"


class SimulationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.native = FIXTURE.read_text()
        cls.linear, _ = make_linear_pdb(cls.native)

    def test_presets_exist(self) -> None:
        for key in ("fast", "standard", "high"):
            p = resolve_params(key, None)
            self.assertGreater(p.guide_steps, 0)
            self.assertGreater(p.relax_steps, 0)

    def test_override_rejects_unknown_key(self) -> None:
        with self.assertRaises(ValueError):
            resolve_params("standard", {"nonexistent": 1})

    def test_choose_platform_invalid(self) -> None:
        with self.assertRaises(RuntimeError):
            choose_platform("TPU")

    def test_end_to_end_md(self) -> None:
        """小体系全流程:线性链 → GBn2 → 引导+松弛 → DCD 产物。

        步数取标准档的一半量级,足以验证数值稳定与 RMSD 收敛趋势。
        """
        params = SimParams(
            guide_steps=3500,
            relax_steps=1500,
            n_frames=100,
            minimization_iters=500,
            k_max=5000.0,
            seed=42,
        )
        with tempfile.TemporaryDirectory() as td:
            meta = run_targeted_md(
                self.linear, self.native, params, Path(td), platform_pref="CPU"
            )
            self.assertTrue((Path(td) / "topology.pdb").exists())
            self.assertTrue((Path(td) / "trajectory.dcd").exists())
            self.assertTrue((Path(td) / "meta.json").exists())
            self.assertGreater(meta["n_atoms"], 300)
            self.assertGreaterEqual(len(meta["rmsd_A"]), 90)
            # 起点为全延伸链,与天然态距离应远(>25 Å)
            self.assertGreater(meta["rmsd_start_A"], 25.0)
            # 关键帧形变 + 松弛后应收敛到演示级贴合(<6 Å)
            self.assertLess(meta["rmsd_end_A"], 6.0)
            self.assertEqual(meta["platform"], "CPU")


if __name__ == "__main__":
    unittest.main()
