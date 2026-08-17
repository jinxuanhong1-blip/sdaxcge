# Inputs

- `GSE131907_samples.tsv` / `GSE205335_patients.tsv` — PR #320 author-malignant CLDN4 %pos and T/NK fractions (winning pair).
- `GSE131907_malignant_counts.tsv.gz` / `GSE205335_malignant_counts.tsv.gz` — patient UMI-sum of author malignant cells (from the public scRNA GSEA extra).
- `GSE205335_gsm_sample_metadata.csv` — GSM → patient / orig.ident map (A3 GSE205335).
- `a8_sets.json` — Hallmark IFN, custom MHC-I, KEGG/GO tight junction.
- `GSE131907_tnk_counts.tsv.gz` / `GSE205335_tnk_counts.tsv.gz` — author T/NK UMI-sum from the same winning-pair units (`build_tnk_pseudobulk.py`).
