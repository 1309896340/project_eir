"""氨基酸序列的输入清洗与校验。"""

from __future__ import annotations

import re
from dataclasses import dataclass

# 力场 amber14 只认识标准 20 种残基;越界序列一律拒绝而不是静默替换
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

MIN_LENGTH = 20
MAX_LENGTH = 400

_FASTA_HEADER = re.compile(r"^>.*$")


@dataclass
class SequenceCheck:
    sequence: str
    errors: list[str]
    warnings: list[str]

    @property
    def ok(self) -> bool:
        return not self.errors


def strip_fasta(text: str) -> str:
    """接受裸序列或 FASTA 文本;FASTA 时仅取首条记录。"""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if lines and lines[0].startswith(">"):
        body: list[str] = []
        for ln in lines[1:]:
            if ln.startswith(">"):
                break  # 遇到第二条 header 即停止
            body.append(ln)
        return "".join(body)
    return "".join(lines)


def clean_and_check(text: str) -> SequenceCheck:
    seq = strip_fasta(text).upper()
    seq = re.sub(r"[\s0-9\-*\.]", "", seq)

    errors: list[str] = []
    warnings: list[str] = []

    invalid = sorted(set(seq) - STANDARD_AA)
    if invalid:
        errors.append(
            "包含非标准氨基酸字符: {}。本工具的物理模拟基于标准残基力场,"
            "请先手动处理非标准残基(如 U 沿硒代半胱氨酸、修饰残基等)。".format(
                " ".join(invalid)
            )
        )
    if not MIN_LENGTH <= len(seq) <= MAX_LENGTH:
        errors.append(
            f"序列长度 {len(seq)} aa 超出支持范围 [{MIN_LENGTH}, {MAX_LENGTH}] aa。"
            "上限主要受 ESMFold/ESM3 预测质量与模拟耗时约束。"
        )
    elif len(seq) > 300:
        warnings.append(
            "序列超过 300 aa,模拟耗时将显著增加(默认 CPU 下约 10~30 分钟);"
            "建议先用『快速预览』档体验。"
        )

    return SequenceCheck(sequence=seq, errors=errors, warnings=warnings)
