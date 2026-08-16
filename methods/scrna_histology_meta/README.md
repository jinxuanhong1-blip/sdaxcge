# methods/scrna_histology_meta

Additive public-only slice: **split-and-pool lung tumor scRNA by histology**.

For every open series with **LUAD and/or LUSC labels** and **TACSTD2 / CLDN4 in epithelium**, compute sample/patient Spearman of **malignant (or epithelial) TACSTD2 or CLDN4 vs T/NK fraction separately in LUAD and in LUSC**, then Fisher-z meta **within histology**. Forest by histology. Honest n / ρ / p.

This is the larger-n single-cell counterpart of the bulk TLS histology interaction (`notes/opus_tls`). It does **not** replace that analysis.

**Open GEO only.** EGA / dbGaP / DAC matrices are catalogued, not downloaded.

| Path | Role |
| --- | --- |
| [playbook.md](playbook.md) | Estimand, gates, what not to claim (EN + 中文) |
| [config/series_registry.csv](config/series_registry.csv) | Series audit |
| [scripts/](scripts/) | Extract / Spearman / meta / forest |
| [harvested/](harvested/) | Sample-level tables (sibling extracts + GSE241934) |
| [results/](results/) | Effects, meta, forests, writeup |

```bash
python3 -m pip install -r methods/scrna_histology_meta/requirements.txt
# optional: download GSE241934 MTX into /tmp/gse241934 then
python3 methods/scrna_histology_meta/scripts/extract_gse241934.py \
  --data-dir /tmp/gse241934 \
  --out methods/scrna_histology_meta/harvested/gse241934_per_sample.tsv
python3 methods/scrna_histology_meta/scripts/01_compute_cohort_effects.py
python3 methods/scrna_histology_meta/scripts/02_meta_and_forest.py
```

Primary result (do not round away the CI): **LUAD TACSTD2 vs T/NK** random-effects ρ = **−0.21** (95% CI −0.44 to +0.05), **p = 0.11**, k = 5, n = 72, I² = 0%. **LUSC has no primary cohort with n ≥ 6.**
