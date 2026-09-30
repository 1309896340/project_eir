"""CodonTransformer 推理示例:对给定蛋白生成宿主=Homo sapiens 的优化 CDS。

运行(需先 `uv sync`):
    uv run python scripts/codontransformer_demo.py

权重优先从本地 `CodonTransformer/data/` 加载(已下载),否则回退 HuggingFace。
输出为 DNA(T),同时按项目约定给出 RNA(U)版本。
"""

import logging
from pathlib import Path
from typing import cast

import torch
from transformers import AutoTokenizer, BigBirdForMaskedLM

from CodonTransformer.CodonData import get_amino_acid_sequence
from CodonTransformer.CodonEvaluation import get_GC_content
from CodonTransformer.CodonPrediction import predict_dna_sequence
from CodonTransformer.CodonUtils import DNASequencePrediction


class _DropGenerativeCapabilityWarning(logging.Filter):
    """丢弃 transformers 的 "has generative capabilities" 误报(详见 codontransformer_app.py)。"""

    def filter(self, record: logging.LogRecord) -> bool:
        return "has generative capabilities" not in record.getMessage()


logging.getLogger("transformers.modeling_utils").addFilter(_DropGenerativeCapabilityWarning())

PROTEIN = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLA"
ORGANISM = "Homo sapiens"
MODEL_DIR = Path(__file__).resolve().parent.parent / "CodonTransformer" / "data"
MODEL_ID = "adibvafa/CodonTransformer"


def main() -> None:
    source = MODEL_DIR if (MODEL_DIR / "model.safetensors").exists() else MODEL_ID
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"model:  {source}")
    print(f"device: {device} ({torch.cuda.get_device_name(0) if device.type == 'cuda' else 'CPU'})")

    tokenizer = AutoTokenizer.from_pretrained(str(source))
    model = cast(
        "torch.nn.Module", BigBirdForMaskedLM.from_pretrained(str(source))
    )
    model.to(device)
    model.eval()

    # 1) 确定性解码:全局最优的单条设计(match_protein=True 强制同义约束,
    #    实现为非同义密码子 logits 置 -inf,与设计文档的硬约束方案一致)
    # deterministic=True + num_sequences=1 时运行时返回单条预测
    output = cast(
        "DNASequencePrediction",
        predict_dna_sequence(
            protein=PROTEIN,
            organism=ORGANISM,
            device=device,
            tokenizer=tokenizer,
            model=model,
            attention_type="original_full",
            deterministic=True,
            match_protein=True,
        ),
    )
    dna = output.predicted_dna

    # 翻译一致性校验(蛋白序列一致为硬约束)
    # 注:return_correct_seq=False 时运行时返回 str,但库注解为宽联合类型
    translated = cast("str", get_amino_acid_sequence(dna, stop_symbol="", codon_table=1))
    assert translated.rstrip("_") == PROTEIN, "翻译产物与输入蛋白不一致!"

    print(f"\norganism:      {output.organism}")
    print(f"protein ({len(PROTEIN)} aa):  {PROTEIN}")
    print(f"DNA ({len(dna)} nt):    {dna}")
    print(f"RNA:           {dna.replace('T', 'U')}")
    print(f"GC content:    {get_GC_content(dna):.1f}%")
    print(f"translate ok:  {translated == PROTEIN}")

    # 2) 非确定性采样:同义变体候选(用于后续多目标打分/主动学习)
    # num_sequences=3 时运行时返回 List[DNASequencePrediction]
    variants = predict_dna_sequence(
        protein=PROTEIN,
        organism=ORGANISM,
        device=device,
        tokenizer=tokenizer,
        model=model,
        attention_type="original_full",
        deterministic=False,
        temperature=0.5,
        top_p=0.95,
        num_sequences=3,
        match_protein=True,
    )
    variant_list = cast("list[DNASequencePrediction]", variants)
    print(f"\nsampled {len(variant_list)} synonymous variants (temperature=0.5):")
    seen = {dna}
    for i, v in enumerate(variant_list, 1):
        tag = "" if v.predicted_dna not in seen else " (dup)"
        seen.add(v.predicted_dna)
        print(f"  variant {i}: GC={get_GC_content(v.predicted_dna):.1f}%{tag}")


if __name__ == "__main__":
    main()
