"""CodonTransformer 推理示例:对给定蛋白生成宿主=Homo sapiens 的优化 CDS。

运行(需先 `uv sync`):
    uv run python scripts/codontransformer_demo.py

权重优先从本地 `CodonTransformer/data/` 加载(已下载),否则回退 HuggingFace。
输出为 DNA(T),同时按项目约定给出 RNA(U)版本。
"""

from pathlib import Path

import torch
from transformers import AutoTokenizer, BigBirdForMaskedLM

from CodonTransformer.CodonData import get_amino_acid_sequence
from CodonTransformer.CodonEvaluation import get_GC_content
from CodonTransformer.CodonPrediction import predict_dna_sequence

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
    model = BigBirdForMaskedLM.from_pretrained(str(source)).to(device).eval()

    # 1) 确定性解码:全局最优的单条设计(match_protein=True 强制同义约束,
    #    实现为非同义密码子 logits 置 -inf,与设计文档的硬约束方案一致)
    output = predict_dna_sequence(
        protein=PROTEIN,
        organism=ORGANISM,
        device=device,
        tokenizer=tokenizer,
        model=model,
        attention_type="original_full",
        deterministic=True,
        match_protein=True,
    )
    dna = output.predicted_dna

    # 翻译一致性校验(蛋白序列一致为硬约束)
    translated = get_amino_acid_sequence(dna, stop_symbol="", codon_table=1)
    assert translated.rstrip("_") == PROTEIN, "翻译产物与输入蛋白不一致!"

    print(f"\norganism:      {output.organism}")
    print(f"protein ({len(PROTEIN)} aa):  {PROTEIN}")
    print(f"DNA ({len(dna)} nt):    {dna}")
    print(f"RNA:           {dna.replace('T', 'U')}")
    print(f"GC content:    {get_GC_content(dna):.1f}%")
    print(f"translate ok:  {translated == PROTEIN}")

    # 2) 非确定性采样:同义变体候选(用于后续多目标打分/主动学习)
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
    print(f"\nsampled {len(variants)} synonymous variants (temperature=0.5):")
    seen = {dna}
    for i, v in enumerate(variants, 1):
        tag = "" if v.predicted_dna not in seen else " (dup)"
        seen.add(v.predicted_dna)
        print(f"  variant {i}: GC={get_GC_content(v.predicted_dna):.1f}%{tag}")


if __name__ == "__main__":
    main()
