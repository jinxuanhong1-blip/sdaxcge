# methods/scrna_scvi_combo

Additive public-only scRNA integration. No claim-audit.

## Question

After a real scVI (else Harmony) integration of compatible public lung-tumor scRNA, what is the sample-level Spearman of malignant TACSTD2 vs T/NK fraction? Honest n / ρ / p.

## Series (do not force all)

| Series | Decision |
|---|---|
| GSE207422 | Use. Neoadjuvant PD-1 + chemo NSCLC, 15 samples, UMI TSV. |
| GSE241934 IIT | Use. Neoadjuvant sintilimab + chemo EGFR-mutant NSCLC, 11 tumors, 10x MTX. |
| GSE241934 RWC | Same GSE; not a second series. Held unless the IIT pair cannot be scored. |
| GSE253013 | Skip for scVI/Harmony HVGs. GEO has a 9.3 GB Seurat RDS only; treatment-naive HiSeq LUAD. |

Pairwise HVG integration that is actually runnable here: **GSE207422 + GSE241934 IIT**.

## Reproduce

```bash
bash methods/scrna_scvi_combo/download.sh /tmp/scrna_scvi_combo
python3 methods/scrna_scvi_combo/prepare.py \
  --raw /tmp/scrna_scvi_combo \
  --outdir /tmp/scrna_scvi_combo/h5ad
python3 methods/scrna_scvi_combo/integrate_score.py \
  --h5ad /tmp/scrna_scvi_combo/h5ad \
  --outdir methods/scrna_scvi_combo/results
```

`scvi-tools` is preferred. If training fails, the script falls back to Harmony and records the error.

## Scoring

- Unit = 10x sample.
- Marker malignant-like = epithelial lineage and low normal-lung score (epithelial 75th percentile).
- Cluster malignant = Leiden cluster on the integrated latent space with the same rule.
- T/NK = matching lineage fraction.
- Eligible: ≥10 malignant-like and ≥20 T/NK cells.
- Spearman ρ / p. Cell-level correlations are not reported.
- Keep the pair (or the post-integration score) where malignant TACSTD2 vs T/NK is negative. All contrasts stay in `results/stats.tsv`.

Numbers after the public run live in `results/WRITEUP.md`.
