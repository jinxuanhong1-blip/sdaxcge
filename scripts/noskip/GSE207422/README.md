# GSE207422 no-skip: malignant TACSTD2 vs MPR/NMPR and T/NK

Hu et al., Genome Medicine 2023 (PMID 36869384). Neoadjuvant PD-1 + chemo NSCLC scRNA, ~92k cells (BD Rhapsody).

Honest result: claimed per-patient malignant TACSTD2 vs T/NK ρ −0.40 to −0.50 is **not supported** (post-treatment ρ = −0.021, p = 0.95, n = 12). NMPR > MPR is directional only (p = 0.21). TACSTD2 is malignant-restricted vs T/NK (paired p = 6.1e-5).

## Inputs (processed GEO only)

- `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (175 MB; 24,292 genes × 92,330 cells)
- `GSE207422_NSCLC_scRNAseq_metadata.xlsx` (sample-level clinical labels)

No barcode-level author annotation file is deposited on GEO or in paper Additional files 1–4. Labels are reconstructed from the authors' published canonical-marker scheme (Fig. 1B / Methods) and their malignant rule (epithelial, then exclude normal lung programs; CopyKAT-like CNV score as support).

## Run

```bash
python3 scripts/noskip/GSE207422/download.py --outdir data/GSE207422
python3 scripts/noskip/GSE207422/analyze.py \
  --matrix data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --metadata data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  --gene-chr scripts/noskip/GSE207422/gene_chr.tsv \
  --outdir results/noskip/GSE207422
```
