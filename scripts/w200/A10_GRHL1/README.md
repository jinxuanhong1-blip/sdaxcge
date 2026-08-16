# A10: GRHL1 vs TACSTD2 / CLDN4 in public lung RNA

Honest public-data test of the claim that **GRHL1** is a shared correlate of **TACSTD2 (TROP2)** and **CLDN4** in lung.

The test is co-expression only (Spearman + Pearson). It can support or refute a *shared-correlate* claim. It does **not** prove direct transcriptional regulation.

Primary cohort is pre-specified: **TCGA NSCLC primary tumors (LUAD + LUSC)**. LUAD, LUSC, and adjacent-normal slices are reported as sensitivity, not as a search for a nicer number.

```bash
python3 scripts/w200/A10_GRHL1/download.py
python3 scripts/w200/A10_GRHL1/analyze.py
```

Outputs land in `results/w200/A10_GRHL1/`.

## Support rule (fixed before looking at the numbers)

A pair supports the claim if Spearman ρ ≥ 0.30, p ≤ 0.05, and the sign is positive.

- both pairs pass → `SUPPORTED_FOR_BOTH`
- only one pair passes → `PARTIAL_ONLY_ONE_TARGET`
- neither pair passes → `NOT_SUPPORTED`

GRHL2 / GRHL3 and CLDN3 / CLDN7 are context only (the closest paralog and neighboring claudins). They are not used to rescue the claim.

Pooling LUAD+LUSC can cancel a within-LUAD GRHL1–CLDN4 correlation. That split is reported as sensitivity. The primary score stays on mixed NSCLC.

## Data

UCSC Xena TCGA hub, HiSeqV2 RSEM gene-level `log2(norm_count+1)`:

- [TCGA.LUAD.sampleMap/HiSeqV2.gz](https://tcga-xena-hub.s3.dualstack.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz)
- [TCGA.LUSC.sampleMap/HiSeqV2.gz](https://tcga-xena-hub.s3.dualstack.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap/HiSeqV2.gz)

Raw matrices stay in `/tmp/xena_tcga_lung/`. The repo stores the small gene extract and derived tables.
