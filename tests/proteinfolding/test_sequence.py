"""ProteinFolding 序列清洗与校验测试。"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2]))

from ProteinFolding.core.sequence import clean_and_check, strip_fasta  # noqa: E402


class StripFastaTests(unittest.TestCase):
    def test_plain_sequence(self) -> None:
        self.assertEqual(strip_fasta("ACDE\nFGHI"), "ACDEFGHI")

    def test_fasta_single(self) -> None:
        self.assertEqual(strip_fasta(">sp|123|INS\nMALWM\nRLLPL"), "MALWMRLLPL")

    def test_fasta_takes_first(self) -> None:
        text = ">a\nMALWM\n>b\nRLLPL"
        self.assertEqual(strip_fasta(text), "MALWM")

    def test_digits_removed_in_clean(self) -> None:
        check = clean_and_check("MA1LWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALY")
        self.assertTrue(check.ok)
        self.assertEqual(check.sequence[:6], "MALWMR")


class CleanAndCheckTests(unittest.TestCase):
    def test_valid_sequence(self) -> None:
        seq = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT" * 2
        check = clean_and_check(seq)
        self.assertTrue(check.ok)
        self.assertEqual(check.sequence, seq)

    def test_too_short(self) -> None:
        check = clean_and_check("MALWM")
        self.assertFalse(check.ok)
        self.assertIn("超出支持范围", check.errors[0])

    def test_too_long(self) -> None:
        seq = "A" * 401
        check = clean_and_check(seq)
        self.assertFalse(check.ok)

    def test_nonstandard_residue_rejected(self) -> None:
        seq = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKX" * 1 + "A" * 5
        check = clean_and_check(seq + "A" * (100 - len(seq)))
        self.assertFalse(check.ok)
        self.assertIn("非标准氨基酸", check.errors[0])

    def test_lowercase_uppercased(self) -> None:
        seq = ("malwmrllpllallalwgpdpaaafvnqhlcgshlvealylvcgergffytpkt" * 2)[:60]
        check = clean_and_check(seq)
        self.assertTrue(check.ok)
        self.assertEqual(check.sequence, seq.upper())

    def test_long_sequence_warns(self) -> None:
        seq = "A" * 320
        check = clean_and_check(seq)
        self.assertTrue(check.ok)
        self.assertTrue(any("300" in w for w in check.warnings))


if __name__ == "__main__":
    unittest.main()
