"""Tests for the human CDS primary-reference deduplication utility."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).parents[1] / "scripts" / "deduplicate_human_cds.py"
SPEC = importlib.util.spec_from_file_location("deduplicate_human_cds", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
deduplicate_human_cds = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(deduplicate_human_cds)


HEADER = ["Gene ID", "Protein ID", "Accession", "Organelle", "# Codons"]


class DeduplicateFileTests(unittest.TestCase):
    def write_tsv(self, path: Path, rows: list[list[str]], header: list[str] = HEADER) -> None:
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(header)
            writer.writerows(rows)

    def read_tsv(self, path: Path) -> list[list[str]]:
        with path.open(newline="", encoding="utf-8") as stream:
            return list(csv.reader(stream, delimiter="\t"))

    def test_keeps_first_primary_record_per_gene_protein_key(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.tsv"
            output_path = Path(directory) / "output.tsv"
            self.write_tsv(
                input_path,
                [
                    ["GENE1", "NP_1.1", "NC_000001", "genomic", "10"],
                    ["GENE1", "NP_1.1", "NC_000002", "genomic", "11"],
                    ["GENE2", "", "NC_000003", "genomic", "12"],
                    ["GENE2", "", "NC_000004", "genomic", "13"],
                ],
            )

            summary = deduplicate_human_cds.deduplicate_file(input_path, output_path)

            self.assertEqual(
                self.read_tsv(output_path),
                [
                    HEADER,
                    ["GENE1", "NP_1.1", "NC_000001", "genomic", "10"],
                    ["GENE2", "", "NC_000003", "genomic", "12"],
                ],
            )
            self.assertEqual(summary.kept_rows, 2)
            self.assertEqual(summary.duplicate_rows, 2)

    def test_excludes_mitochondrial_and_non_primary_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.tsv"
            output_path = Path(directory) / "output.tsv"
            self.write_tsv(
                input_path,
                [
                    ["MT-ND1", "YP_1.1", "NC_012920", "mitochondrion", "319"],
                    ["GENE1", "NP_1.1", "NT_123456", "genomic", "10"],
                    ["GENE1", "NP_1.1", "NC_000001", "genomic", "10"],
                ],
            )

            summary = deduplicate_human_cds.deduplicate_file(input_path, output_path)

            self.assertEqual(
                self.read_tsv(output_path),
                [HEADER, ["GENE1", "NP_1.1", "NC_000001", "genomic", "10"]],
            )
            self.assertEqual(summary.non_genomic_rows, 1)
            self.assertEqual(summary.non_primary_rows, 1)

    def test_rejects_input_without_required_columns(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.tsv"
            output_path = Path(directory) / "output.tsv"
            self.write_tsv(input_path, [["GENE1", "NP_1.1", "NC_000001"]], HEADER[:3])

            with self.assertRaisesRegex(ValueError, "缺少必需列"):
                deduplicate_human_cds.deduplicate_file(input_path, output_path)


if __name__ == "__main__":
    unittest.main()
