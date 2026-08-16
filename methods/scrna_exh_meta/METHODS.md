# Methods — T/NK cytotoxicity / exhaustion vs malignant TACSTD2 and CLDN4

Additive extra analysis. Public GEO objects only. Patient is the unit.

## Question

In public lung ICI or neoadjuvant-ICI scRNA that still has a malignant (or author epithelial) compartment **and** T/NK cells, how do per-patient malignant **TACSTD2** and **CLDN4** relate to T/NK **cytotoxicity** and **exhaustion**?

## Gene sets

| Score | Compartment | Genes |
|---|---|---|
| Cytotoxicity | T/NK | GZMB, PRF1, GNLY, NKG7 |
| Exhaustion | T/NK | PDCD1, HAVCR2, LAG3, TIGIT, TOX |
| TACSTD2 | malignant / author Epi | TACSTD2 |
| CLDN4 | malignant / author Epi | CLDN4 |

Per cell: `log1p(CP10k) = log1p(UMI / library_size × 10⁴)`.  
Per patient: mean of the cell-level gene (or mean of the genes in the set) inside the named compartment.  
A gene missing from a matrix is dropped; the remaining genes in the set are averaged. Present genes are recorded in `results/summary.json`.

## Patient filter

Keep a patient when **≥10 malignant (or author-Epi) cells** and **≥10 T/NK cells**.  
Spearman uses the patients that still have both scores after that filter.  
Meta-analysis uses cohorts with **n ≥ 5**.

## Cohorts

| Cohort | GEO | Setting | Malignant definition | T/NK definition |
|---|---|---|---|---|
| GSE207422 | Hu et al., *Genome Med* 2023 | Neoadjuvant PD-1 + chemo; post-resection | Marker: epithelial lineage **and** tumor-epi module > normal-lung module **and** normal-lung < 0.4. CopyKAT IDs are not on GEO. | Marker lineage max (CD3D/E/G, TRAC, CD2, NKG7, GNLY, KLRD1) |
| GSE205335 | Park/Ahn/Lee lung ICI atlas | Palliative ICI; tumor / LN / effusion | Author `lineage.sub == Malignant cells` | Author `lineage.total == T/NK cells` |
| GSE241934 IIT | NEOTIDE / CTONG2104 | EGFR-mutant neoadjuvant sintilimab + chemo | Author `major.cell.type == Epi` | Author T or NK |
| GSE241934 RWC | same series, real-world arm | EGFR-WT neoadjuvant PD-1 + chemo | Author Epi | Author T or NK |
| GSE291670 | Xia et al., *J Transl Med* 2025 | Neoadjuvant anlotinib + camrelizumab | Same marker rule as GSE207422 | Same marker T/NK |
| GSE233203 | leftover | Pre-ABCP pleural effusion (3 R / 4 NR) | Same marker rule | Same marker T/NK; **enter the meta only if ≥5 patients pass the cell filter** |

GSE241934 IIT and RWC are separate studies (different EGFR status and enrollment). They are not pooled before the meta.

GSE205335 patients are the mean of non-`Normal*` samples. RECIST is carried but is not the primary contrast.

## Lineage modules (marker cohorts)

Epithelial: EPCAM, KRT8/18/19/5/7/17, ELF3, CDH1, MUC1  
T/NK: CD3D, CD3E, CD3G, TRAC, CD2, NKG7, GNLY, KLRD1  
Normal lung: SFTPA1/A2/B/D, AGER, NAPSA, SCGB1A1, SCGB3A2, TPPP3, FOXJ1, CAPS  
Tumor-epi: EPCAM, KRT8/18/19/7, CEACAM5/6, MUC1, ELF3  

QC: library size ≥ 200 UMI.

## Statistics

1. Per cohort: Spearman ρ and two-sided p on the patient table.  
2. Fisher z = artanh(ρ), SE = 1/√(n−3).  
3. Inverse-variance fixed-effect and DerSimonian–Laird random-effect meta.  
4. Report k, N, ρ, 95% CI, p, I². No cell-level p-values (those are pseudo-replicates).

## Leftovers that were not scored in the meta

See `leftover_catalog.tsv`. Short version:

- GSE146100 — both compartments in the paper, **n=1 patient**
- GSE243013 — public object is CD45+ immune-only
- GSE266035, GSE176022, GSE162498 — T/CD3-sorted, no epithelium
- GSE131907, GSE253013, GSE150938, GSE189357 — not ICI / not neoadjuvant ICI, or n=3 unlabeled tumors
- GSE325414 — PEF ablation, not ICI

## How to run

```bash
python3 methods/scrna_exh_meta/download.py --datadir /tmp/scrna_exh_meta
python3 methods/scrna_exh_meta/extract.py  --datadir /tmp/scrna_exh_meta
python3 methods/scrna_exh_meta/analyze.py  --datadir /tmp/scrna_exh_meta \
  --outdir methods/scrna_exh_meta/results
```

Outputs: `methods/scrna_exh_meta/results/` and a copy of the extra figure under `results/scrna_exh_meta/`.
