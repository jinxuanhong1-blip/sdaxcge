# Additive GSE148071 REAL Slingshot/PAGA scored by CLDN4

Public advanced NSCLC malignant-like epithelium (Wu et al. 2021, GSE148071). **CLDN4 only.** Not a TACSTD2 redo. No dual-high gate.

PR #394 barrier ρ=+0.590 is given and is not re-audited. This folder runs **real slingshot** (R) plus PAGA and asks where CLDN4, barrier (CLDN4 held out), and IFN sit along pseudotime. Patient is the unit.

- Finding: [`FINDING.md`](FINDING.md)
- Lineage table: [`tables/lineage_level.tsv`](tables/lineage_level.tsv)
- Patient table: [`tables/patient_level.tsv`](tables/patient_level.tsv)

```bash
bash methods/gse148071_slingshot_real_cldn4/scripts/install_tools.sh
bash methods/gse148071_slingshot_real_cldn4/scripts/download.sh /tmp/gse148071_sling_data
python3 methods/gse148071_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse148071_sling_data \
  --out /tmp/gse148071_sling_data/epithelium.h5ad
python3 methods/gse148071_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse148071_sling_data/epithelium.h5ad \
  --outdir methods/gse148071_slingshot_real_cldn4 \
  --finding methods/gse148071_slingshot_real_cldn4/FINDING.md
```

Honest n is written by the analyze script. GEO n=42 is the catalog, not the Spearman n.
