# 数据版本化与数据仓库骨架

## 基于 DVC + MLflow 的 mRNA 设计数据仓库

参考 QC 流程文档，建立从原始数据到模型训练样本的完整数据版本化体系。核心目标：

> **数据可追溯、版本可复现、实验可对比、QC 可审计。**

---

# 一、总体设计原则

| 原则 | 说明 |
|---|---|
| 数据与代码分离 | 大文件不进 Git，由 DVC 管理 |
| 版本不可变 | 每次 QC 或处理产生新版本，不覆盖旧版本 |
| 血缘可追溯 | 记录 raw → processed → features → splits 的完整链路 |
| 实验可对比 | MLflow 记录每次训练的数据版本、参数、指标 |
| QC 可审计 | 每次 QC 报告与数据版本绑定 |
| 环境可复现 | Docker + conda/uv 锁定依赖 |
| 远程可共享 | DVC remote 支持 S3/MinIO/SSH |

---

# 二、目录结构

```text
mrna-designer/
│
├── .dvc/                          # DVC 内部配置
│   ├── config
│   └── plots/
│
├── .dvcignore                     # DVC 忽略规则
├── .gitignore                     # Git 忽略规则
├── .mlflow/                       # MLflow 本地配置
│
├── configs/                       # 配置文件
│   ├── config.yaml                # 主配置
│   ├── data/
│   │   ├── gencode_human.yaml     # GENCODE 人类数据配置
│   │   ├── gencode_mouse.yaml     # GENCODE 小鼠数据配置
│   │   └── splits.yaml            # 数据划分配置
│   ├── model/
│   │   ├── generator.yaml         # 生成模型配置
│   │   ├── predictor.yaml         # 预测模型配置
│   │   └── structure.yaml         # 结构模块配置
│   ├── qc/
│   │   ├── sequence_qc.yaml       # 序列 QC 配置
│   │   ├── cds_check.yaml         # CDS 检查配置
│   │   └── translation.yaml       # 翻译一致性配置
│   └── train/
│       ├── pretrain.yaml
│       ├── finetune.yaml
│       └── predictor.yaml
│
├── data/                          # 数据目录（DVC 管理）
│   ├── raw/                       # 原始数据（只读，不修改）
│   │   ├── gencode/
│   │   │   ├── human/
│   │   │   │   ├── GRCh38.primary_assembly.genome.fa
│   │   │   │   ├── gencode.v49.annotation.gtf
│   │   │   │   ├── gencode.v49.transcripts.fa
│   │   │   │   └── gencode.v49.pc_translations.fa
│   │   │   └── mouse/
│   │   │       ├── GRCm39.primary_assembly.genome.fa
│   │   │       ├── gencode.vM36.annotation.gtf
│   │   │       └── gencode.vM36.pc_translations.fa
│   │   ├── external/
│   │   │   ├── cocoputs/
│   │   │   ├── optimus_5prime/
│   │   │   ├── openvaccine/
│   │   │   ├── ribonanza/
│   │   │   ├── rmdb/
│   │   │   └── geo_sra/
│   │   └── README.md
│   │
│   ├── interim/                   # 中间数据（DVC 管理）
│   │   ├── parsed_gtf/
│   │   │   ├── human_cds.parquet
│   │   │   ├── human_utr5.parquet
│   │   │   ├── human_utr3.parquet
│   │   │   └── human_metadata.parquet
│   │   ├── extracted_sequences/
│   │   │   ├── human_cds.fa
│   │   │   ├── human_utr5.fa
│   │   │   ├── human_utr3.fa
│   │   │   └── human_full_mrna.fa
│   │   └── strand_normalized/
│   │       └── human_normalized.parquet
│   │
│   ├── processed/                 # 处理后数据（DVC 管理）
│   │   ├── v1/
│   │   │   ├── sequences.parquet
│   │   │   ├── annotations.parquet
│   │   │   ├── experiments.parquet
│   │   │   ├── labels.parquet
│   │   │   └── features.parquet
│   │   └── v2/
│   │       └── ...
│   │
│   ├── features/                  # 计算特征（DVC 管理）
│   │   ├── composition/
│   │   │   ├── gc_content.parquet
│   │   │   ├── cai.parquet
│   │   │   └── codon_pair.parquet
│   │   ├── structure/
│   │   │   ├── mfe.parquet
│   │   │   ├── pairing_prob.parquet
│   │   │   └── accessibility.parquet
│   │   ├── utr/
│   │   │   ├── uorf.parquet
│   │   │   ├── kozak.parquet
│   │   │   └── mirna_site.parquet
│   │   └── safety/
│   │       ├── restriction_site.parquet
│   │       └── splice_site.parquet
│   │
│   ├── splits/                    # 数据划分（DVC 管理）
│   │   ├── v1/
│   │   │   ├── train.parquet
│   │   │   ├── valid.parquet
│   │   │   ├── test.parquet
│   │   │   └── split_metadata.json
│   │   └── v2/
│   │       └── ...
│   │
│   ├── negatives/                 # 负样本（DVC 管理）
│   │   ├── lncrna.fa
│   │   ├── pseudogene.fa
│   │   └── random_synonymous.fa
│   │
│   └── qc/                        # QC 报告（DVC 管理）
│       ├── v1/
│       │   ├── qc_report.json
│       │   ├── qc_summary.md
│       │   ├── failures/
│       │   │   ├── cds_structure_fail.parquet
│       │   │   ├── translation_fail.parquet
│       │   │   └── sequence_quality_fail.parquet
│       │   └── plots/
│       │       ├── gc_distribution.png
│       │       ├── length_distribution.png
│       │       └── pass_rate.png
│       └── v2/
│           └── ...
│
├── models/                        # 模型（DVC 管理）
│   ├── pretrained/
│   │   ├── backbone_v0.pt
│   │   └── backbone_v1.pt
│   ├── generator/
│   │   ├── cds_generator_v0.pt
│   │   └── cds_generator_v1.pt
│   ├── predictor/
│   │   ├── multitask_v0.pt
│   │   └── multitask_v1.pt
│   └── structure/
│       └── structure_module_v0.pt
│
├── src/                           # 源代码
│   ├── data/
│   │   ├── download.py            # 下载 GENCODE 数据
│   │   ├── parse_gtf.py           # 解析 GTF
│   │   ├── extract_sequence.py    # 提取序列
│   │   ├── strand_normalize.py    # 方向标准化
│   │   ├── split_regions.py       # 区域切分
│   │   ├── check_cds.py           # CDS 结构检查
│   │   ├── check_translation.py   # 翻译一致性
│   │   ├── sequence_qc.py         # 序列质量过滤
│   │   ├── compute_features.py    # 特征计算
│   │   ├── build_splits.py        # 数据划分
│   │   └── build_negatives.py     # 负样本构建
│   ├── models/
│   │   ├── generator.py
│   │   ├── predictor.py
│   │   ├── structure.py
│   │   └── tokenizer.py
│   ├── train/
│   │   ├── pretrain.py
│   │   ├── finetune_generator.py
│   │   ├── train_predictor.py
│   │   └── reward_model.py
│   ├── eval/
│   │   ├── evaluate_generator.py
│   │   ├── evaluate_predictor.py
│   │   └── metrics.py
│   ├── optimize/
│   │   ├── constrained_decode.py
│   │   └── pareto.py
│   ├── safety/
│   │   └── biosecurity_screen.py
│   └── utils/
│       ├── io.py
│       ├── logging.py
│       └── config.py
│
├── pipelines/                     # DVC 流水线定义
│   ├── dvc.yaml                   # 主流水线
│   ├── download.dvc
│   ├── parse_gtf.dvc
│   ├── extract_sequence.dvc
│   ├── qc.dvc
│   ├── features.dvc
│   ├── splits.dvc
│   └── train.dvc
│
├── scripts/                       # 便捷脚本
│   ├── setup_env.sh
│   ├── download_all.sh
│   ├── run_qc.sh
│   ├── run_training.sh
│   └── upload_remote.sh
│
├── notebooks/                     # 探索性分析
│   ├── 01_data_overview.ipynb
│   ├── 02_qc_analysis.ipynb
│   ├── 03_feature_eda.ipynb
│   └── 04_model_eval.ipynb
│
├── tests/                         # 单元测试
│   ├── test_data/
│   ├── test_models/
│   └── test_qc/
│
├── docs/                          # 文档
│   ├── data_dictionary.md
│   ├── qc_process.md
│   ├── db_schema.md
│   └── model_card.md
│
├── mlruns/                        # MLflow 本地运行记录（.gitignore）
├── .gitignore
├── .dvcignore
├── dvc.yaml
├── dvc.lock
├── pyproject.toml
├── requirements.txt
├── environment.yml
├── Dockerfile
├── Makefile
└── README.md
```

---

# 三、DVC 配置

## 3.1 初始化 DVC

```bash
# 在项目根目录初始化
git init
dvc init

# 配置远程存储（以 MinIO 为例）
dvc remote add -d storage s3://mrna-designer-dvc/data
dvc remote modify storage endpointurl http://minio.example.com:9000
dvc remote modify storage access_key_id ${AWS_ACCESS_KEY_ID}
dvc remote modify storage secret_access_key ${AWS_SECRET_ACCESS_KEY}
```

## 3.2 .dvcignore

```text
# 忽略不需要 DVC 管理的文件
.git/
.mlruns/
mlruns/
notebooks/
*.pyc
__pycache__/
.env
```

## 3.3 .gitignore

```text
# 数据由 DVC 管理，Git 只跟踪 .dvc 文件
data/raw/
data/interim/
data/processed/
data/features/
data/splits/
data/negatives/
data/qc/
models/

# MLflow
mlruns/
.mlruns/

# Python
__pycache__/
*.pyc
*.egg-info/
.pytest_cache/

# 环境
.env
.venv/
venv/

# IDE
.idea/
.vscode/
```

## 3.4 DVC 流水线定义

```yaml
# dvc.yaml
stages:
  download_gencode:
    cmd: python src/data/download.py --config configs/data/gencode_human.yaml
    deps:
      - src/data/download.py
      - configs/data/gencode_human.yaml
    outs:
      - data/raw/gencode/human/
    params:
      - configs/data/gencode_human.yaml:
          - version
          - species

  parse_gtf:
    cmd: python src/data/parse_gtf.py --config configs/data/gencode_human.yaml
    deps:
      - src/data/parse_gtf.py
      - data/raw/gencode/human/gencode.v49.annotation.gtf
    outs:
      - data/interim/parsed_gtf/
    params:
      - configs/data/gencode_human.yaml:
          - gene_type
          - transcript_type

  extract_sequence:
    cmd: python src/data/extract_sequence.py --config configs/data/gencode_human.yaml
    deps:
      - src/data/extract_sequence.py
      - data/raw/gencode/human/GRCh38.primary_assembly.genome.fa
      - data/interim/parsed_gtf/
    outs:
      - data/interim/extracted_sequences/

  strand_normalize:
    cmd: python src/data/strand_normalize.py
    deps:
      - src/data/strand_normalize.py
      - data/interim/extracted_sequences/
    outs:
      - data/interim/strand_normalized/

  qc:
    cmd: python src/data/sequence_qc.py --config configs/qc/sequence_qc.yaml
    deps:
      - src/data/sequence_qc.py
      - src/data/check_cds.py
      - src/data/check_translation.py
      - data/interim/strand_normalized/
    outs:
      - data/qc/v1/
    params:
      - configs/qc/sequence_qc.yaml:
          - n_ratio_max
          - gc_min
          - gc_max
          - max_homopolymer

  compute_features:
    cmd: python src/data/compute_features.py
    deps:
      - src/data/compute_features.py
      - data/qc/v1/
    outs:
      - data/features/

  build_splits:
    cmd: python src/data/build_splits.py --config configs/data/splits.yaml
    deps:
      - src/data/build_splits.py
      - data/features/
    outs:
      - data/splits/v1/
    params:
      - configs/data/splits.yaml:
          - strategy
          - train_chroms
          - valid_chroms
          - test_chroms

  build_negatives:
    cmd: python src/data/build_negatives.py
    deps:
      - src/data/build_negatives.py
      - data/raw/gencode/human/
    outs:
      - data/negatives/

  pretrain:
    cmd: python src/train/pretrain.py --config configs/train/pretrain.yaml
    deps:
      - src/train/pretrain.py
      - src/models/generator.py
      - data/splits/v1/
    outs:
      - models/pretrained/backbone_v0.pt
    params:
      - configs/train/pretrain.yaml:
          - epochs
          - lr
          - batch_size
```

## 3.5 DVC 常用命令

```bash
# 添加数据到 DVC
dvc add data/raw/gencode/human/

# 提交 .dvc 文件到 Git
git add data/raw/gencode/human.dvc .gitignore
git commit -m "Add GENCODE human data v49"

# 推送数据到远程
dvc push

# 拉取数据
dvc pull

# 查看数据版本
dvc list . --dvc-only

# 切换数据版本
git checkout <commit>
dvc checkout

# 查看流水线
dvc dag

# 复现流水线
dvc repro

# 比较数据版本
dvc diff

# 查看指标
dvc metrics show
```

---

# 四、MLflow 配置

## 4.1 MLflow 服务启动

```bash
# 本地启动
mlflow server \
  --backend-store-uri postgresql://user:pass@localhost:5432/mlflow \
  --default-artifact-root s3://mrna-designer-mlflow/artifacts \
  --host 0.0.0.0 \
  --port 5000
```

## 4.2 MLflow 集成代码

```python
# src/train/pretrain.py
import mlflow
import mlflow.pytorch
from omegaconf import OmegaConf
import subprocess

def get_dvc_data_version(path: str) -> str:
    """获取 DVC 数据版本"""
    result = subprocess.run(
        ["dvc", "get", "--json", path],
        capture_output=True, text=True
    )
    return result.stdout.strip()

def train(config_path: str):
    config = OmegaConf.load(config_path)

    # 设置 MLflow
    mlflow.set_tracking_uri("http://mlflow.example.com:5000")
    mlflow.set_experiment("mrna-designer-pretrain")

    with mlflow.start_run(run_name=f"pretrain_{config.version}"):
        # 记录配置
        mlflow.log_params(OmegaConf.to_container(config))

        # 记录数据版本
        mlflow.log_param("data_version", get_dvc_data_version("data/splits/v1"))
        mlflow.log_param("dvc_commit", subprocess.check_output(
            ["git", "rev-parse", "HEAD"]
        ).decode().strip())
        mlflow.log_param("gencode_version", config.data.gencode_version)

        # 记录代码版本
        mlflow.log_param("git_commit", subprocess.check_output(
            ["git", "rev-parse", "HEAD"]
        ).decode().strip())

        # 训练
        model = build_model(config)
        for epoch in range(config.epochs):
            train_loss = train_epoch(model, train_loader)
            val_loss = validate(model, val_loader)

            # 记录指标
            mlflow.log_metrics({
                "train_loss": train_loss,
                "val_loss": val_loss,
                "epoch": epoch
            }, step=epoch)

        # 记录模型
        mlflow.pytorch.log_model(model, "model")

        # 记录 DVC 数据版本
        mlflow.log_artifact("data/splits/v1/split_metadata.json")

        # 记录 QC 报告
        mlflow.log_artifact("data/qc/v1/qc_report.json")

        # 记录特征重要性
        mlflow.log_artifact("data/features/")

        return model
```

## 4.3 MLflow 实验组织

```text
Experiments:
├── mrna-designer-pretrain
│   ├── run_001: backbone_v0, GENCODE v49
│   ├── run_002: backbone_v0, GENCODE v49 + mouse
│   └── run_003: backbone_v1, GENCODE v49
│
├── mrna-designer-generator
│   ├── run_001: cds_generator_v0
│   ├── run_002: cds_generator_v1 + reward
│   └── run_003: cds_generator_v2 + constrained decode
│
├── mrna-designer-predictor
│   ├── run_001: multitask_v0
│   ├── run_002: multitask_v1 + structure
│   └── run_003: multitask_v2 + condition
│
└── mrna-designer-reward
    ├── run_001: reward_v0
    └── run_002: reward_v1
```

## 4.4 MLflow 模型注册

```python
# 注册模型到 Model Registry
mlflow.register_model(
    model_uri=f"runs:/{run_id}/model",
    name="mRNA-CDS-Generator"
)

# 阶段转换
client = mlflow.tracking.MlflowClient()
client.transition_model_version_stage(
    name="mRNA-CDS-Generator",
    version=1,
    stage="Staging"
)
client.transition_model_version_stage(
    name="mRNA-CDS-Generator",
    version=1,
    stage="Production"
)
```

---

# 五、数据版本化策略

## 5.1 版本命名规范

| 数据类型 | 版本格式 | 示例 |
|---|---|---|
| 原始数据 | `{source}_{version}` | `gencode_v49` |
| 中间数据 | `{stage}_{date}` | `parsed_20260915` |
| 处理后数据 | `v{major}.{minor}` | `v1.0` |
| 特征 | `feat_v{major}` | `feat_v1` |
| 划分 | `split_v{major}` | `split_v1` |
| QC | `qc_v{major}` | `qc_v1` |
| 模型 | `{type}_v{major}` | `generator_v1` |

## 5.2 版本升级规则

| 变更类型 | 版本升级 | 示例 |
|---|---|---|
| 数据源更新 | major | GENCODE v49 → v50 |
| QC 标准变更 | major | GC 范围 40–60 → 35–65 |
| 特征新增 | minor | 新增 MFE 特征 |
| 划分策略变更 | major | 按染色体 → 按基因 |
| Bug 修复 | patch | 修复负链处理 |
| 模型架构变更 | major | Transformer → GNN |
| 超参数调整 | minor | lr 1e-4 → 5e-5 |

## 5.3 数据血缘追踪

```json
{
  "data_version": "v1.0",
  "created_at": "2026-09-15T10:00:00Z",
  "lineage": {
    "raw": {
      "gencode_human": "v49",
      "gencode_mouse": "vM36",
      "cocoputs": "2024-01",
      "optimus_5prime": "2023-06"
    },
    "interim": {
      "parsed_gtf": "parsed_20260915",
      "extracted_sequences": "extracted_20260915",
      "strand_normalized": "normalized_20260915"
    },
    "processed": {
      "sequences": "v1.0",
      "annotations": "v1.0",
      "labels": "v1.0"
    },
    "features": {
      "composition": "feat_v1",
      "structure": "feat_v1",
      "utr": "feat_v1",
      "safety": "feat_v1"
    },
    "splits": {
      "train": "split_v1",
      "valid": "split_v1",
      "test": "split_v1"
    },
    "qc": {
      "report": "qc_v1",
      "pass_rate": 0.85
    }
  },
  "dvc_commit": "a1b2c3d4e5f6",
  "git_commit": "f6e5d4c3b2a1"
}
```

---

# 六、QC 与数据版本绑定

## 6.1 QC 报告结构

```json
{
  "qc_id": "qc_v1",
  "data_version": "v1.0",
  "created_at": "2026-09-15T10:00:00Z",
  "source": {
    "gencode_human": "v49",
    "gencode_mouse": "vM36"
  },
  "checks": {
    "u_t_standardization": {"pass": 100000, "fail": 0},
    "strand_normalization": {"pass": 100000, "fail": 0},
    "region_splitting": {"pass": 98000, "fail": 2000},
    "cds_structure": {"pass": 90000, "fail": 10000},
    "translation_consistency": {"pass": 88000, "fail": 2000},
    "sequence_quality": {"pass": 85000, "fail": 3000}
  },
  "final_dataset": {
    "train": 68000,
    "valid": 8500,
    "test": 8500
  },
  "dvc_commit": "a1b2c3d4e5f6",
  "git_commit": "f6e5d4c3b2a1"
}
```

## 6.2 QC 与 DVC 集成

```bash
# QC 报告作为 DVC 输出
dvc add data/qc/v1/

# QC 报告版本与数据版本绑定
dvc metrics add data/qc/v1/qc_report.json

# 查看 QC 指标
dvc metrics show data/qc/v1/qc_report.json
```

## 6.3 QC 与 MLflow 集成

```python
# 训练时记录 QC 报告
with mlflow.start_run():
    mlflow.log_artifact("data/qc/v1/qc_report.json")
    mlflow.log_metric("qc_pass_rate", 0.85)
    mlflow.log_param("qc_version", "qc_v1")
```

---

# 七、数据仓库骨架使用流程

## 7.1 初始化

```bash
# 克隆项目
git clone https://github.com/example/mrna-designer.git
cd mrna-designer

# 创建环境
conda env create -f environment.yml
conda activate mrna-designer

# 初始化 DVC
dvc pull  # 拉取所有数据

# 启动 MLflow
mlflow server --host 0.0.0.0 --port 5000
```

## 7.2 数据准备

```bash
# 复现整个数据流水线
dvc repro download_gencode
dvc repro parse_gtf
dvc repro extract_sequence
dvc repro qc
dvc repro compute_features
dvc repro build_splits
dvc repro build_negatives

# 或一键复现
dvc repro
```

## 7.3 训练

```bash
# 预训练
dvc repro pretrain

# 或手动运行并记录 MLflow
python src/train/pretrain.py --config configs/train/pretrain.yaml
```

## 7.4 版本切换

```bash
# 切换到数据版本 v1
git checkout v1.0
dvc checkout

# 切换到数据版本 v2
git checkout v2.0
dvc checkout
```

## 7.5 共享与协作

```bash
# 推送数据到远程
dvc push

# 拉取最新数据
dvc pull

# 查看数据状态
dvc status
```

---

# 八、Makefile 便捷命令

```makefile
.PHONY: setup download qc features splits train eval clean

setup:
	conda env create -f environment.yml
	dvc init --no-scm
	dvc remote add -d storage s3://mrna-designer-dvc/data
	pre-commit install

download:
	dvc repro download_gencode

qc:
	dvc repro parse_gtf
	dvc repro extract_sequence
	dvc repro strand_normalize
	dvc repro qc

features:
	dvc repro compute_features
	dvc repro build_negatives

splits:
	dvc repro build_splits

pretrain:
	dvc repro pretrain

finetune:
	dvc repro finetune_generator

predictor:
	dvc repro train_predictor

train: pretrain finetune predictor

eval:
	python src/eval/evaluate_generator.py
	python src/eval/evaluate_predictor.py

mlflow:
	mlflow server --host 0.0.0.0 --port 5000 \
	  --backend-store-uri postgresql://user:pass@localhost:5432/mlflow \
	  --default-artifact-root s3://mrna-designer-mlflow/artifacts

push:
	dvc push
	git push

pull:
	git pull
	dvc pull

clean:
	dvc remove --all
	rm -rf data/interim data/processed data/features data/splits
```

---

# 九、数据版本化最佳实践

| 实践 | 说明 |
|---|---|
| 原始数据只读 | `data/raw/` 永不修改，只添加新版本 |
| 中间数据可重建 | `data/interim/` 可通过流水线重建 |
| 处理后数据版本化 | `data/processed/` 每次 QC 产生新版本 |
| 特征与序列绑定 | 特征的 `sequence_id` 与数据版本绑定 |
| QC 报告不可变 | 每次 QC 报告与数据版本一一对应 |
| 划分固定 | 一旦确定，不再修改，保证可对比 |
| 模型与数据绑定 | MLflow 记录训练数据版本 |
| 远程备份 | DVC remote 定期备份 |
| 版本标签 | 重要版本打 Git tag |
| 文档同步 | 每次版本升级更新文档 |

---

# 十、总结

本数据仓库骨架以 **DVC + MLflow** 为核心，提供：

1. **完整目录结构**：从 raw 到 splits，覆盖 QC 流程每一步。
2. **DVC 流水线**：`dvc.yaml` 定义可复现的数据处理流程。
3. **MLflow 集成**：记录实验、数据版本、模型、指标。
4. **版本化策略**：命名规范、升级规则、血缘追踪。
5. **QC 绑定**：QC 报告与数据版本一一对应。
6. **便捷命令**：Makefile 一键复现。
7. **协作共享**：DVC remote 支持团队协作。

该骨架可直接支撑从 GENCODE 数据准备到模型训练的完整流程，保证数据可追溯、版本可复现、实验可对比、QC 可审计。