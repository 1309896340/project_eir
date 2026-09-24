下面给出可直接用于 API 契约、数据校验和前后端联调的 **输入/输出 JSON Schema** 与 **示例样本**。设计原则：

- 输入：患者基因样本相关元数据 + 目标蛋白 + 物种/细胞 + 修饰 + 递送 + 目标权重 + 约束。
- 输出：top-N mRNA 候选，包含完整序列、分区、评分、约束报告、安全报告、不确定性和人工审核标记。
- 序列统一使用 RNA 字母表 `A/U/G/C`；蛋白使用单字母氨基酸。
- 权重建议归一化，总和为 1。
- 所有候选必须经过安全筛查和人工审核。

---

# 一、输入 JSON Schema

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://example.com/schemas/mrna-design-input.schema.json",
  "title": "mRNA Design Input",
  "type": "object",
  "required": [
    "request_id",
    "patient",
    "target_protein",
    "delivery",
    "modification",
    "objectives"
  ],
  "properties": {
    "request_id": {
      "type": "string",
      "description": "请求唯一 ID"
    },
    "patient": {
      "type": "object",
      "required": ["patient_id", "species", "cell_type"],
      "properties": {
        "patient_id": {
          "type": "string",
          "description": "脱敏患者 ID"
        },
        "species": {
          "type": "string",
          "enum": ["human", "mouse", "rat", "monkey", "hamster", "cho", "other"]
        },
        "cell_type": {
          "type": "string",
          "examples": [
            "hepatocyte",
            "dendritic_cell",
            "T_cell",
            "HEK293T",
            "CHO",
            "cardiomyocyte",
            "other"
          ]
        },
        "tissue": {
          "type": "string",
          "examples": ["liver", "muscle", "tumor", "blood", "other"]
        },
        "genetic_sample_id": {
          "type": "string",
          "description": "患者基因样本 ID"
        },
        "notes": {
          "type": "string"
        }
      }
    },
    "target_protein": {
      "type": "object",
      "required": ["name", "sequence"],
      "properties": {
        "name": {
          "type": "string"
        },
        "sequence": {
          "type": "string",
          "pattern": "^[ACDEFGHIKLMNPQRSTVWY]+$",
          "description": "目标蛋白氨基酸序列，单字母代码"
        },
        "uniprot_id": {
          "type": "string"
        },
        "source": {
          "type": "string",
          "enum": ["patient_genetic_sample", "reference", "synthetic", "other"]
        },
        "notes": {
          "type": "string"
        }
      }
    },
    "delivery": {
      "type": "object",
      "required": ["method"],
      "properties": {
        "method": {
          "type": "string",
          "enum": [
            "LNP",
            "electroporation",
            "lipofection",
            "naked_mRNA",
            "viral_vector",
            "other"
          ]
        },
        "route": {
          "type": "string",
          "enum": [
            "intramuscular",
            "intravenous",
            "subcutaneous",
            "intratumoral",
            "ex_vivo",
            "other"
          ]
        },
        "notes": {
          "type": "string"
        }
      }
    },
    "modification": {
      "type": "object",
      "required": ["nucleoside"],
      "properties": {
        "nucleoside": {
          "type": "string",
          "enum": ["unmodified", "m1Ψ", "Ψ", "5mC", "m6A", "other"]
        },
        "capping": {
          "type": "string",
          "enum": ["CleanCap", "ARCA", "enzymatic", "none", "other"]
        },
        "polyA_length": {
          "type": "integer",
          "minimum": 0,
          "maximum": 300
        },
        "notes": {
          "type": "string"
        }
      }
    },
    "objectives": {
      "type": "object",
      "required": ["expression", "stability", "immunogenicity"],
      "properties": {
        "expression": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        },
        "stability": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        },
        "translation_efficiency": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        },
        "immunogenicity": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        },
        "manufacturability": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        },
        "safety": {
          "type": "number",
          "minimum": 0,
          "maximum": 1
        }
      },
      "additionalProperties": false,
      "description": "目标权重，建议总和为 1"
    },
    "constraints": {
      "type": "object",
      "properties": {
        "gc_range": {
          "type": "array",
          "items": {
            "type": "number",
            "minimum": 0,
            "maximum": 1
          },
          "minItems": 2,
          "maxItems": 2
        },
        "max_homopolymer": {
          "type": "integer",
          "minimum": 1,
          "maximum": 10
        },
        "avoid_motifs": {
          "type": "array",
          "items": {
            "type": "string",
            "pattern": "^[AUGC]+$"
          }
        },
        "forbidden_enzymes": {
          "type": "array",
          "items": {
            "type": "string"
          }
        },
        "avoid_internal_stop": {
          "type": "boolean"
        },
        "min_length": {
          "type": "integer",
          "minimum": 1
        },
        "max_length": {
          "type": "integer",
          "minimum": 1
        },
        "safety_screen": {
          "type": "boolean"
        }
      }
    },
    "metadata": {
      "type": "object"
    }
  }
}
```

---

# 二、输出 JSON Schema

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://example.com/schemas/mrna-design-output.schema.json",
  "title": "mRNA Design Output",
  "type": "object",
  "required": [
    "request_id",
    "status",
    "candidates",
    "model_version",
    "generated_at"
  ],
  "properties": {
    "request_id": {
      "type": "string"
    },
    "status": {
      "type": "string",
      "enum": ["success", "partial_success", "failure"]
    },
    "candidates": {
      "type": "array",
      "items": {
        "type": "object",
        "required": [
          "candidate_id",
          "rank",
          "mrna_sequence",
          "regions",
          "scores",
          "constraints_report",
          "safety_report",
          "uncertainty",
          "review_required"
        ],
        "properties": {
          "candidate_id": {
            "type": "string"
          },
          "rank": {
            "type": "integer",
            "minimum": 1
          },
          "mrna_sequence": {
            "type": "string",
            "pattern": "^[AUGC]+$"
          },
          "regions": {
            "type": "object",
            "required": ["utr5", "cds", "utr3", "polyA"],
            "properties": {
              "utr5": {
                "type": "string",
                "pattern": "^[AUGC]*$"
              },
              "cds": {
                "type": "string",
                "pattern": "^[AUGC]+$"
              },
              "utr3": {
                "type": "string",
                "pattern": "^[AUGC]*$"
              },
              "polyA": {
                "type": "string",
                "pattern": "^A*$"
              }
            }
          },
          "scores": {
            "type": "object",
            "required": [
              "expression",
              "stability",
              "translation_efficiency",
              "immunogenicity",
              "manufacturability",
              "safety",
              "composite_score"
            ],
            "properties": {
              "expression": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              },
              "stability": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              },
              "translation_efficiency": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              },
              "immunogenicity": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              },
              "manufacturability": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              },
              "safety": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              },
              "composite_score": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              }
            }
          },
          "constraints_report": {
            "type": "object",
            "required": [
              "gc_content",
              "internal_stop",
              "homopolymer_max",
              "motif_violations",
              "enzyme_sites",
              "length",
              "pass"
            ],
            "properties": {
              "gc_content": {
                "type": "number",
                "minimum": 0,
                "maximum": 1
              },
              "internal_stop": {
                "type": "boolean"
              },
              "homopolymer_max": {
                "type": "integer",
                "minimum": 0
              },
              "motif_violations": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "enzyme_sites": {
                "type": "array",
                "items": {
                  "type": "string"
                }
              },
              "length": {
                "type": "integer",
                "minimum": 1
              },
              "pass": {
                "type": "boolean"
              }
            }
          },
          "safety_report": {
            "type": "object",
            "required": [
              "biosecurity_pass",
              "controlled_sequence",
              "patient_privacy_ok",
              "review_required"
            ],
            "properties": {
              "biosecurity_pass": {
                "type": "boolean"
              },
              "controlled_sequence": {
                "type": "boolean"
              },
              "patient_privacy_ok": {
                "type": "boolean"
              },
              "review_required": {
                "type": "boolean"
              },
              "notes": {
                "type": "string"
              }
            }
          },
          "uncertainty": {
            "type": "number",
            "minimum": 0,
            "maximum": 1
          },
          "review_required": {
            "type": "boolean"
          },
          "warnings": {
            "type": "array",
            "items": {
              "type": "string"
            }
          },
          "provenance": {
            "type": "object"
          }
        }
      }
    },
    "errors": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "warnings": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "model_version": {
      "type": "string"
    },
    "generated_at": {
      "type": "string",
      "format": "date-time"
    },
    "schema_version": {
      "type": "string"
    }
  }
}
```

---

# 三、示例输入样本

```json
{
  "request_id": "REQ-2026-001",
  "patient": {
    "patient_id": "P001",
    "species": "human",
    "cell_type": "hepatocyte",
    "tissue": "liver",
    "genetic_sample_id": "GS-001",
    "notes": "罕见遗传病，需补充正常酶"
  },
  "target_protein": {
    "name": "示例酶",
    "sequence": "MKT",
    "uniprot_id": "EXAMPLE",
    "source": "patient_genetic_sample",
    "notes": "由患者基因样本分析得到的目标蛋白"
  },
  "delivery": {
    "method": "LNP",
    "route": "intravenous",
    "notes": "肝靶向 LNP"
  },
  "modification": {
    "nucleoside": "m1Ψ",
    "capping": "CleanCap",
    "polyA_length": 120,
    "notes": "降低免疫原性"
  },
  "objectives": {
    "expression": 0.4,
    "stability": 0.3,
    "immunogenicity": 0.2,
    "manufacturability": 0.1
  },
  "constraints": {
    "gc_range": [0.4, 0.6],
    "max_homopolymer": 4,
    "avoid_motifs": ["AUUUA"],
    "forbidden_enzymes": ["EcoRI", "BamHI"],
    "avoid_internal_stop": true,
    "min_length": 100,
    "max_length": 5000,
    "safety_screen": true
  },
  "metadata": {
    "requested_by": "clinician_001",
    "priority": "high"
  }
}
```

---

# 四、示例输出样本

```json
{
  "request_id": "REQ-2026-001",
  "status": "success",
  "candidates": [
    {
      "candidate_id": "C001",
      "rank": 1,
      "mrna_sequence": "GGGAAACCCAUGAAAACCUAAGGGCCCAAAAAAAAAAAAAAAAAAAA",
      "regions": {
        "utr5": "GGGAAACCC",
        "cds": "AUGAAAACCUAA",
        "utr3": "GGGCCC",
        "polyA": "AAAAAAAAAAAAAAAAAAAA"
      },
      "scores": {
        "expression": 0.91,
        "stability": 0.87,
        "translation_efficiency": 0.89,
        "immunogenicity": 0.12,
        "manufacturability": 0.84,
        "safety": 0.95,
        "composite_score": 0.88
      },
      "constraints_report": {
        "gc_content": 0.52,
        "internal_stop": false,
        "homopolymer_max": 3,
        "motif_violations": [],
        "enzyme_sites": [],
        "length": 47,
        "pass": true
      },
      "safety_report": {
        "biosecurity_pass": true,
        "controlled_sequence": false,
        "patient_privacy_ok": true,
        "review_required": true,
        "notes": "需分子生物学家和临床医生审核"
      },
      "uncertainty": 0.08,
      "review_required": true,
      "warnings": [],
      "provenance": {
        "model_version": "0.1.0",
        "pipeline": "generator-v0 + predictor-v0",
        "data_version": "2026-09-15"
      }
    },
    {
      "candidate_id": "C002",
      "rank": 2,
      "mrna_sequence": "GGGAAACCCAUGAAAACCUAAGGGCCCAAAAAAAAAAAAAAAAAAAA",
      "regions": {
        "utr5": "GGGAAACCC",
        "cds": "AUGAAAACCUAA",
        "utr3": "GGGCCC",
        "polyA": "AAAAAAAAAAAAAAAAAAAA"
      },
      "scores": {
        "expression": 0.85,
        "stability": 0.93,
        "translation_efficiency": 0.82,
        "immunogenicity": 0.10,
        "manufacturability": 0.88,
        "safety": 0.96,
        "composite_score": 0.86
      },
      "constraints_report": {
        "gc_content": 0.52,
        "internal_stop": false,
        "homopolymer_max": 3,
        "motif_violations": [],
        "enzyme_sites": [],
        "length": 47,
        "pass": true
      },
      "safety_report": {
        "biosecurity_pass": true,
        "controlled_sequence": false,
        "patient_privacy_ok": true,
        "review_required": true
      },
      "uncertainty": 0.10,
      "review_required": true,
      "warnings": []
    }
  ],
  "errors": [],
  "warnings": [],
  "model_version": "0.1.0",
  "generated_at": "2026-09-15T10:00:00Z",
  "schema_version": "1.0.0"
}
```

---

# 五、校验与版本说明

- **输入必填**：`request_id`、`patient`、`target_protein`、`delivery`、`modification`、`objectives`。
- **输出必填**：`request_id`、`status`、`candidates`、`model_version`、`generated_at`。
- **权重**：`objectives` 中各权重建议总和为 1；JSON Schema 不强制求和，可在业务层校验。
- **序列**：mRNA 使用 `A/U/G/C`；蛋白使用标准单字母氨基酸。
- **安全**：所有候选必须 `biosecurity_pass = true`，且 `review_required = true` 进入人工审核。
- **版本**：`schema_version`、`model_version`、`data_version` 均需记录，保证可追溯。
- **扩展**：可在 `metadata`、`provenance` 中追加实验批次、递送配方、HLA 类型等字段，不影响核心契约。