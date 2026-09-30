"""ProteinFolding 线性链构建与缓存层测试。

fixture ubq20.pdb:泛素(1UBQ)前 20 个残基的重原子子集。
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))

from ProteinFolding.core.cache import cache_key, find_cached  # noqa: E402
from ProteinFolding.core.config import SimParams  # noqa: E402
from ProteinFolding.core.linear import FloatArray, make_linear_pdb  # noqa: E402

FIXTURE = Path(__file__).parent / "fixtures" / "ubq20.pdb"


def _atom_coords(pdb_text: str) -> dict[int, dict[str, FloatArray]]:
    coords: dict[int, dict[str, FloatArray]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        res = int(line[22:26])
        name = line[12:16].strip()
        xyz = np.array([float(line[30 + 8 * i : 38 + 8 * i]) for i in range(3)])
        coords.setdefault(res, {})[name] = xyz
    return coords


class LinearChainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.native = FIXTURE.read_text()
        cls.linear, cls.sequence = make_linear_pdb(cls.native)

    def test_residue_count_preserved(self) -> None:
        self.assertEqual(len(self.sequence), 20)

    def test_extended_backbone_span(self) -> None:
        coords = _atom_coords(self.linear)
        bb_n = np.array([coords[r]["N"] for r in sorted(coords)])
        span = float(np.linalg.norm(bb_n[-1] - bb_n[0]))
        # 全延伸时每残基 ~3.3-3.5 Å
        self.assertGreater(span, 19 * 3.0)
        self.assertLess(span, 19 * 3.8)

    def test_chirality_matches_native(self) -> None:
        native_coords = _atom_coords(self.native)
        linear_coords = _atom_coords(self.linear)
        dets = []
        for coords in (native_coords, linear_coords):
            r = coords[5]
            m = np.column_stack([r["N"] - r["CA"], r["C"] - r["CA"], r["CB"] - r["CA"]])
            dets.append(float(np.linalg.det(m)))
        # Kabsch 旋转保手性:行列式符号一致
        self.assertEqual(np.sign(dets[0]), np.sign(dets[1]))

    def test_no_heavy_atom_clashes(self) -> None:
        coords = _atom_coords(self.linear)
        xyz = np.array(
            [v for r in sorted(coords) for v in coords[r].values()]
        )
        d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        self.assertGreater(float(d.min()), 1.15)

    def test_oxt_added_when_missing(self) -> None:
        self.assertIn("OXT", _atom_coords(self.linear)[20])

    def test_sidechain_atom_inventory_preserved(self) -> None:
        native_coords = _atom_coords(self.native)
        linear_coords = _atom_coords(self.linear)
        for res in sorted(native_coords):
            self.assertEqual(
                set(native_coords[res]) - {"O", "OXT"},
                set(linear_coords[res]) - {"O", "OXT"},
                f"残基 {res} 侧链原子清单不一致",
            )


class CacheTests(unittest.TestCase):
    def test_key_is_stable_and_sensitive(self) -> None:
        k1 = cache_key("MALWM", SimParams())
        k2 = cache_key("MALWM", SimParams())
        k3 = cache_key("MALWN", SimParams())
        k4 = cache_key("MALWM", SimParams(guide_steps=99))
        self.assertEqual(k1, k2)
        self.assertNotEqual(k1, k3)
        self.assertNotEqual(k1, k4)
        self.assertLessEqual(len(k1), 16)

    def test_find_cached_missing(self) -> None:
        self.assertIsNone(find_cached("nonexistent0000"))


if __name__ == "__main__":
    unittest.main()
