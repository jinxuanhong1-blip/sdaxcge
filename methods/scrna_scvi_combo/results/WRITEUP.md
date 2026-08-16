# scVI combo: public lung-tumor scRNA

Additive methods only. Public GEO. No claim-audit.

## Series decision

| Series | Used | Why |
|---|---|---|
| GSE207422 | yes | Neoadjuvant PD-1 + chemo NSCLC; 15 10x samples; 92,330 cells; UMI TSV on GEO |
| GSE241934 IIT | yes | Neoadjuvant sintilimab + chemo EGFR-mutant NSCLC; 11 tumors; 78,691 cells; 10x MTX |
| GSE241934 RWC | no (not a second GSE) | Same series; 229k cells / 1.2 GB MTX not required once the IIT pair was negative |
| GSE253013 | no | 9.3 GB Seurat RDS; no HVG-ready count matrix on 16 GB; treatment-naive HiSeq LUAD |

Pairwise HVG integrations that are actually runnable: **GSE207422 + GSE241934 IIT** only. That pair is negative, so it is kept.

## Integration

scVI-tools 1.5.0.post1, CPU, 40 epochs, `n_latent=10`, 2,000 HVGs (`seurat_v3`, batch-aware), batch = GEO series.

- Cells after QC: **171,021** (92,330 + 78,691)
- Shared genes: **18,913**
- Harmony was not needed (scVI trained)

## Scoring (honest n / ρ / p)

Unit = 10x sample. Malignant TACSTD2 = mean log1p(CP10k) in malignant-like cells. T/NK = lineage fraction of all QC cells. Eligible: ≥10 malignant-like and ≥20 T/NK cells. Spearman, two-sided. Cell-level tests are not reported.

| Contrast | n | ρ | p | negative |
|---|---:|---:|---:|---|
| Marker lineage after scVI (combo) | 26 | −0.169 | 0.41 | yes |
| Marker lineage pre-latent (same cells) | 26 | −0.169 | 0.41 | yes |
| Marker, GSE207422 only | 15 | −0.054 | 0.85 | yes |
| Marker, GSE241934 IIT only | 11 | −0.409 | 0.21 | yes |
| Leiden-cluster malignant after scVI | 7 | −0.500 | 0.25 | yes |
| Cluster, GSE207422 only | 7 | −0.500 | 0.25 | yes |
| Cluster, GSE241934 IIT only | 0 | NA | NA | no |

The marker and pre-latent ρ values are identical because TACSTD2 and T/NK are scored on raw log1p(CP10k), not on the latent space. scVI is used to build a joint embedding and a Leiden annotation; it does not change the marker numbers.

## Kept pair

Kept **GSE207422 + GSE241934 IIT**, primary score = **marker malignant after scVI**.

- **n = 26** (15 + 11)
- **ρ = −0.169**
- **p = 0.41**

Sign is negative (the keep rule). The correlation is weak and not significant. GSE241934 IIT alone is more negative (n=11, ρ=−0.41, p=0.21) and GSE207422 alone is near zero (n=15, ρ=−0.054, p=0.85).

Leiden-cluster malignant is also negative (n=7, ρ=−0.50, p=0.25) but **drops every GSE241934 sample** (IIT is T-cell rich; author Epi n=1,699 / 78,691). That is a sensitivity, not the combo score.

## Definitions

- Marker malignant-like: argmax lineage score is Epithelial, and normal-lung score ≤ epithelial 75th percentile (cut = 0.383).
- Cluster malignant: Leiden cluster on `X_scVI` whose mean marker score is Epithelial and whose mean normal-lung score is ≤ the epithelial-cluster 75th percentile.
- T/NK: matching lineage fraction.
- GSE241934 author `Epi`/`T`/`NK` labels are not the primary score (they exist only on that series).

## Skip record

```json
{
  "series": "GSE253013",
  "used_in_scvi": false,
  "reason": "GEO deposits only GSE253013_all_luad_garnett_temp.rds.gz (9.3 GB, double-gzipped Seurat/XDR). A full count matrix for scVI/Harmony HVG training does not fit this 16 GB session. The series is also treatment-naive LUAD on HiSeq 2500, whereas GSE207422 and GSE241934 IIT are neoadjuvant PD-1 + chemo NSCLC on NovaSeq 6000. Not forced.",
  "public_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013",
  "pmid": "38335304"
}
```
