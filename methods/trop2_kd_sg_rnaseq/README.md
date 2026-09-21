# TACSTD2 knockdown and sacituzumab RNA-seq

Public GEO matrices only. The question is whether CLDN4, a classical NHEJ core, the cGAS–STING axis, and a compact ISG set move after TACSTD2/Trop2 knockdown or sacituzumab govitecan (IMMU-132).

```bash
pip install -r methods/trop2_kd_sg_rnaseq/requirements.txt
python3 methods/trop2_kd_sg_rnaseq/analyze.py
```

Processed matrices are downloaded from the NCBI GEO FTP into `results/trop2_kd_sg_rnaseq/data/` (gitignored). Tables and figures land in `results/trop2_kd_sg_rnaseq/`.

Gene lists are fixed in `gene_sets.py` before the contrasts are scored. The write-up is `results/trop2_kd_sg_rnaseq/RESULTS.md`.
