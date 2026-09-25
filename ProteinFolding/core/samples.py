"""内置演示样例,与 CodonTransformer app 共用同一批 GENCODE v50 真实人类蛋白。"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Sample:
    key: str
    label: str
    sequence: str


SAMPLES: list[Sample] = [
    Sample(
        key="ins",
        label="INS 胰岛素原 (110 aa)",
        sequence=(
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGG"
            "GPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        ),
    ),
    Sample(
        key="hbb",
        label="HBB 血红蛋白β (147 aa)",
        sequence=(
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVL"
            "GAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGV"
            "ANALAHKYH"
        ),
    ),
    Sample(
        key="epo",
        label="EPO 促红细胞生成素 (193 aa)",
        sequence=(
            "MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVP"
            "DTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLR"
            "ALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR"
        ),
    ),
]

SAMPLES_BY_KEY = {s.key: s for s in SAMPLES}
