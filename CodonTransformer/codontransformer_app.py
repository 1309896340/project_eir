"""CodonTransformer 调试界面:蛋白序列 -> 宿主特异优化 CDS。

运行(需先 `uv sync`,权重从同目录 `data/` 加载):
    uv run python CodonTransformer/codontransformer_app.py
浏览器自动打开 http://127.0.0.1:7860

测试样例取自 GENCODE v50 人类蛋白(INS/HBB/EPO)。
"""

import logging
import time
from pathlib import Path

import gradio as gr
import torch
from transformers import AutoTokenizer, BigBirdForMaskedLM

from CodonTransformer.CodonData import get_amino_acid_sequence
from CodonTransformer.CodonEvaluation import get_GC_content
from CodonTransformer.CodonPrediction import predict_dna_sequence
from CodonTransformer.CodonUtils import ORGANISM2ID, DNASequencePrediction


class _DropGenerativeCapabilityWarning(logging.Filter):
    """丢弃 transformers 的 "has generative capabilities" 误报。

    transformers>=4.50 将 BigBirdForMaskedLM 移出 GenerationMixin(失去 .generate),
    且模型加载期间 can_generate() 被调用多次、每次都重复打印该警告。本应用解码
    走 predict_dna_sequence 的单次 MLM forward,不调 .generate(),与该提示无关,
    故按消息内容精确过滤,不影响 transformers 的其他警告。
    """

    def filter(self, record: logging.LogRecord) -> bool:
        return "has generative capabilities" not in record.getMessage()


logging.getLogger("transformers.modeling_utils").addFilter(_DropGenerativeCapabilityWarning())

MODEL_DIR = Path(__file__).resolve().parent / "data"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# ---- 测试样例:GENCODE v50 真实人类蛋白 ----
EX_INS = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"  # 胰岛素原 110aa
EX_HBB = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"  # 血红蛋白β 147aa
EX_EPO = "MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR"  # 促红细胞生成素 193aa

# ---- 全局加载模型(一次性) ----
tokenizer = AutoTokenizer.from_pretrained(str(MODEL_DIR))
model = BigBirdForMaskedLM.from_pretrained(str(MODEL_DIR))
# transformers 用 functools.wraps 包装了 Module.to,pyright 解析 _Wrapped
# 绑定调用时会把位置实参错配到 self(reportArgumentType,误报),故忽略
model.to(DEVICE)  # pyright: ignore[reportArgumentType]
model.eval()


def optimize(protein: str, organism: str, mode: str, temperature: float,
             top_p: float, num_sequences: int, match_protein: bool):
    protein = "".join(protein.split()).upper()
    if not protein:
        raise gr.Error("请输入蛋白序列(单字母氨基酸)")

    t0 = time.perf_counter()
    best = predict_dna_sequence(
        protein=protein, organism=organism, device=DEVICE,
        tokenizer=tokenizer, model=model,
        attention_type="original_full",
        deterministic=True, match_protein=match_protein,
    )
    assert isinstance(best, DNASequencePrediction)  # deterministic=True 必返回单条
    dna = best.predicted_dna
    translated = get_amino_acid_sequence(dna, stop_symbol="", codon_table=1)
    assert isinstance(translated, str)  # return_correct_seq=False 必返回 str
    ok = translated.rstrip("_") == protein
    elapsed = time.perf_counter() - t0

    summary = (
        f"**宿主** {organism} | **长度** {len(dna)} nt | **GC** {get_GC_content(dna):.1f}%"
        f" | **翻译一致** {'✅' if ok else '❌'} | **耗时** {elapsed:.1f}s"
    )

    variants_info = ""
    variant_rows: list[list] = []
    if mode == "采样(多条变体)":
        t1 = time.perf_counter()
        variants = predict_dna_sequence(
            protein=protein, organism=organism, device=DEVICE,
            tokenizer=tokenizer, model=model,
            attention_type="original_full",
            deterministic=False, temperature=temperature, top_p=top_p,
            num_sequences=int(num_sequences), match_protein=match_protein,
        )
        assert isinstance(variants, list)  # num_sequences>1 必返回列表
        for i, v in enumerate(variants, 1):
            vt = get_amino_acid_sequence(v.predicted_dna, stop_symbol="", codon_table=1)
            assert isinstance(vt, str)
            variant_rows.append([
                i,
                f"{get_GC_content(v.predicted_dna):.1f}%",
                "✅" if vt.rstrip("_") == protein else "❌",
                v.predicted_dna,
                v.predicted_dna.replace("T", "U"),
            ])
        variants_info = (
            f"采样 {len(variant_rows)} 条同义变体"
            f"(temperature={temperature}, top_p={top_p}, {time.perf_counter() - t1:.1f}s)"
        )

    return (
        summary,
        dna,
        dna.replace("T", "U"),
        translated,
        variants_info,
        variant_rows,
    )


def build_ui() -> gr.Blocks:
    with gr.Blocks(title="CodonTransformer 调试台") as demo:
        gr.Markdown(
            "# CodonTransformer 密码子优化调试台\n"
            "输入蛋白序列(单字母氨基酸),选择宿主与解码模式,生成优化 CDS。\n"
            "`match_protein=True` 即同义约束(非同义密码子 logits 置 -inf),"
            "对应设计文档的蛋白序列一致硬约束。"
        )
        with gr.Row():
            with gr.Column():
                protein_in = gr.Textbox(label="蛋白序列(单字母氨基酸,5'→3' 编码方向)",
                                        lines=4, placeholder="以 M 开头…")
                organism_in = gr.Dropdown(label="宿主(164 物种,可输入过滤)",
                                          choices=sorted(ORGANISM2ID.keys()),
                                          value="Homo sapiens")
                mode_in = gr.Radio(["确定性(最优单条)", "采样(多条变体)"],
                                   value="确定性(最优单条)", label="解码模式")
                match_in = gr.Checkbox(value=True, label="match_protein(强制同义约束)")
                with gr.Accordion("采样参数(仅采样模式生效)", open=False):
                    temp_in = gr.Slider(0.1, 0.9, 0.5, step=0.1, label="temperature")
                    topp_in = gr.Slider(0.5, 1.0, 0.95, step=0.05, label="top_p")
                    nseq_in = gr.Slider(1, 5, 3, step=1, label="num_sequences")
                btn = gr.Button("生成优化 CDS", variant="primary")
            with gr.Column():
                summary_out = gr.Markdown(label="概览")
                dna_out = gr.Textbox(label="优化 DNA (5'→3', T)")
                rna_out = gr.Textbox(label="RNA 版本(项目约定, T→U)")
                translated_out = gr.Textbox(label="翻译产物(应为原蛋白+终止符)")
                variants_info_out = gr.Markdown()
                variants_out = gr.Dataframe(
                    headers=["#", "GC", "一致", "DNA", "RNA"],
                    wrap=True,
                    interactive=False,
                    label="同义变体(采样模式)",
                )

        btn.click(
            optimize,
            inputs=[protein_in, organism_in, mode_in, temp_in, topp_in, nseq_in, match_in],
            outputs=[summary_out, dna_out, rna_out, translated_out, variants_info_out, variants_out],
        )

        gr.Examples(
            examples=[
                [EX_INS, "Homo sapiens", "确定性(最优单条)", 0.5, 0.95, 3, True],
                [EX_HBB, "Homo sapiens", "采样(多条变体)", 0.5, 0.95, 3, True],
                [EX_EPO, "Escherichia coli general", "采样(多条变体)", 0.2, 0.95, 3, True],
            ],
            inputs=[protein_in, organism_in, mode_in, temp_in, topp_in, nseq_in, match_in],
            label="测试样例(真实人类蛋白,GENCODE v50):INS 胰岛素原 110aa / HBB 血红蛋白β 147aa / EPO 促红细胞生成素 193aa(演示大肠杆菌宿主对比)",
        )
    return demo


if __name__ == "__main__":
    build_ui().launch(server_name="127.0.0.1", server_port=7860, inbrowser=True)
