# Script templates / 脚本模板

Runnable, minimally-opinionated templates that implement the playbook. Adapt paths and
parameters to your data. Python = `.py`, R = `.R`. Ordered by pipeline stage.

可直接运行的模板，实现 `playbook.md` 的流程。请按你的数据调整路径与参数。

| Stage / 阶段 | Python | R |
|---|---|---|
| 1. Ambient RNA / 环境RNA | `01_ambient_cellbender.sh` | `01_ambient_soupx.R` |
| 2. Doublets + QC / 双细胞+质控 | `02_qc_doublets.py` | `02_qc_doublets.R` |
| 3. Integration / 整合 | `03_integrate_scvi.py`, `03_integrate_harmony.py` | — |
| 4. Annotation / 注释 | `04_annotate.py` | — |
| 4b. Malignant / 恶性判定 | — | `04b_infercnv.R` |
| 5. Pseudobulk DE / 伪bulk差异 | `05_pseudobulk_de.py` | `05_pseudobulk_de.R` |
| 6. Module scores / 模块打分 | `06_module_scores.py` | `06_module_scores.R` |
| 7. Neighborhood / 邻域 | `07_neighborhood_scores.py` | — |
| 8. Cell–cell comm (secondary) / 细胞通讯（次要） | `08_ccc_liana.py` | `08_ccc_cellchat.R` |
| 9. Join clinical / 对接临床 | `09_join_clinical.py` | — |
| gene sets / 基因集 | `gene_sets.py` | `gene_sets.R` |

All scripts set seeds and are safe to read top-to-bottom as documentation.
所有脚本均设定随机种子，可作为文档从上到下阅读。
