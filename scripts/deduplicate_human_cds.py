"""保留人类主染色体 CDS 的首条代表记录。

规则：仅保留 ``Organelle == genomic`` 且 ``Accession`` 以 ``NC_``
开头的记录；对每个 ``(Gene ID, Protein ID)`` 组合，按源文件顺序保留
第一次出现的记录。
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import NamedTuple, Sequence


REQUIRED_COLUMNS = ("Gene ID", "Protein ID", "Accession", "Organelle")
PROJECT_ROOT = Path(__file__).parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "data" / "o537-Human_CDS.tsv"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "o537-Human_CDS.primary_genomic_deduplicated.tsv"


class DeduplicationSummary(NamedTuple):
    """处理结果计数。"""

    input_rows: int
    non_genomic_rows: int
    non_primary_rows: int
    duplicate_rows: int
    kept_rows: int


def _normalise_header(header: list[str]) -> list[str]:
    """移除 CoCoPUTs 行末制表符产生的空表头列。"""
    return header[:-1] if header and header[-1] == "" else header


def deduplicate_file(input_path: Path, output_path: Path) -> DeduplicationSummary:
    """按项目规则筛选并去重一个 Human CDS TSV 文件。"""
    with input_path.open("r", newline="", encoding="utf-8-sig") as input_stream:
        reader = csv.reader(input_stream, delimiter="\t")
        try:
            header = _normalise_header(next(reader))
        except StopIteration as error:
            raise ValueError("输入文件为空，未找到表头") from error

        missing_columns = [column for column in REQUIRED_COLUMNS if column not in header]
        if missing_columns:
            raise ValueError(f"输入文件缺少必需列：{', '.join(missing_columns)}")

        column_index = {column: header.index(column) for column in REQUIRED_COLUMNS}
        seen_keys: set[tuple[str, str]] = set()
        kept_rows: list[list[str]] = []
        input_rows = non_genomic_rows = non_primary_rows = duplicate_rows = 0

        for line_number, row in enumerate(reader, start=2):
            input_rows += 1
            row = row[: len(header)]
            if len(row) != len(header):
                raise ValueError(f"第 {line_number} 行的列数与表头不一致")

            if row[column_index["Organelle"]] != "genomic":
                non_genomic_rows += 1
                continue
            if not row[column_index["Accession"]].startswith("NC_"):
                non_primary_rows += 1
                continue

            key = (row[column_index["Gene ID"]], row[column_index["Protein ID"]])
            if key in seen_keys:
                duplicate_rows += 1
                continue

            seen_keys.add(key)
            kept_rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_stream:
        writer = csv.writer(output_stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(kept_rows)

    return DeduplicationSummary(
        input_rows=input_rows,
        non_genomic_rows=non_genomic_rows,
        non_primary_rows=non_primary_rows,
        duplicate_rows=duplicate_rows,
        kept_rows=len(kept_rows),
    )


def parse_args(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="筛选并去重人类主染色体 CDS 统计表。")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="输入 TSV 路径")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="输出 TSV 路径")
    return parser.parse_args(arguments)


def main(arguments: Sequence[str] | None = None) -> int:
    args = parse_args(arguments)
    summary = deduplicate_file(args.input, args.output)
    print(f"输入数据行数：{summary.input_rows}")
    print(f"排除非 genomic 行数：{summary.non_genomic_rows}")
    print(f"排除非 NC_* 行数：{summary.non_primary_rows}")
    print(f"排除重复主染色体行数：{summary.duplicate_rows}")
    print(f"最终保留行数：{summary.kept_rows}")
    print(f"输出文件：{args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
