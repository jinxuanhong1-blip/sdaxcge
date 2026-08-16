# TISCH2 NSCLC pool: epithelial TACSTD2/CLDN4 vs T/NK fraction

**Slice:** `methods/tisch_nsclc_pool/` only. Additive. Does not rewrite `notes/fable_tisch/` or other scRNA slices.

**Inputs:** public TISCH2 `expression.h5` + `CellMetainfo_table.tsv` only
(`https://tisch.compbio.cn/static/data/<DS>/<DS>_expression.h5`). No extra git clones, no FASTQ, no dbGaP.

**Question:** at the sample/patient level, does epithelial TACSTD2 or CLDN4 track T/NK fraction? Then pool. Note which objects carry ICI labels.

## Catalog (TISCH2 NSCLC gallery, 17 objects)

User-requested: `NSCLC_GSE151537`, `GSE146100`, `GSE117570`, `GSE127465`, `GSE131907`, `GSE148071`, plus every other TISCH NSCLC object that can be scored for TACSTD2.

| TISCH id | Gallery treatment | Epithelium + T/NK in CellMetainfo | ICI labels in CellMetainfo | Role |
|---|---|---|---|---|
| NSCLC_GSE151537 | Immunotherapy | T-only (0 epi) | none | honest skip (GEO is treatment-naive sorted T cells) |
| NSCLC_GSE146100 | Immunotherapy | yes (1 patient, 3 nodules) | none | sample-level table; n=3 too small for Spearman |
| NSCLC_GSE117570 | None | yes | none | pool candidate |
| NSCLC_GSE127465 | None | yes | none | pool candidate (tumor, not PBMC) |
| NSCLC_GSE131907 | None | yes (no Malignant call) | none | pool candidate (tumor-tissue epithelial proxy) |
| NSCLC_GSE148071 | None | yes | none | pool candidate |
| NSCLC_EMTAB6149 | None | yes | none | skip pool: no Patient/Sample column |
| NSCLC_GSE143423 | None | malignant + few T | none | try; likely underpowered |
| NSCLC_GSE149655 | None | alveolar/club, no Malignant | none | try; n_patients=2 |
| NSCLC_GSE150660 | None | yes | none | try; n_patients=2 |
| NSCLC_GSE153935 | None | yes | none | pool candidate if n eligible ≥5 |
| NSCLC_GSE162498 | None | yes | none | pool candidate if n eligible ≥5 |
| NSCLC_GSE127471 | None | PBMC only | none | skip (no epithelium) |
| NSCLC_GSE139555 | None | T/immune-sorted | none | skip (no epithelium) |
| NSCLC_GSE99254 | None | T-sorted | none | skip (no epithelium) |
| NSCLC_GSE176021_aPD1 | Immunotherapy | T-sorted (~817k) | none in downloadable CellMetainfo | skip (no epithelium; do not download huge h5) |
| NSCLC_GSE179373 | None (`Treatment=Systemic therapy`) | T-only | Treatment only | skip (no epithelium) |

**ICI honesty:** TISCH gallery tags GSE151537 / GSE146100 / GSE176021_aPD1 as Immunotherapy. Downloadable CellMetainfo for those objects has **no Response / RECIST / MPR column**. This slice cannot test ICI benefit. GSE151537 and GSE176021 cannot even score epithelium.

## Methods (short)

1. TISCH major-lineage, labels stripped. Malignant if ≥20 cells in the unit; else epithelial-like (`Epithelial`, `Alveolar`, `Basal`, `Ciliated`, `Club`).
2. T/NK = `CD8T`, `CD8Tex`, `CD4Tconv`, `Treg`, `Tprolif`, `TMKI67`, `NK`, `Tcell`, `NKT`, `ILC`.
3. Tumor-like tissue only (`Tumor` / metastasis / mLN / mBrain / pleural). Drop NAT / Normal / PBMC / nLN.
4. Unit = `Sample` when it is a real sample (not a per-cell id); else `Patient`.
5. Eligible unit: ≥20 scored epithelial cells and ≥20 T/NK.
6. Scores: mean TISCH `log2(TPM/10+1)` and % positive (`>0`) for TACSTD2 and CLDN4 in the scored epithelial compartment. T/NK fraction = n_TNK / n_cells in the unit.
7. Per-dataset Spearman (n, ρ, p). Spearman only if n≥5 eligible units.
8. Pool: Fisher-z inverse-variance of dataset ρ (n≥6); secondary = rank-within-dataset then one Spearman. Report I². Do not tune filters to a target ρ.

## Reproduce

```bash
python3 -m pip install -r methods/tisch_nsclc_pool/requirements.txt
python3 methods/tisch_nsclc_pool/download.py
python3 methods/tisch_nsclc_pool/analyze.py
python3 methods/tisch_nsclc_pool/pool.py
```

Raw h5/tsv stay in `data/tisch_nsclc_pool/` (gitignored).
