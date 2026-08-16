# methods/scrna_scenic — ELF3 / GRHL1 / KLF4 regulons vs TACSTD2 / CLDN4

Additive public scRNA GRN. **A10 ELF3–CLDN4 bulk RNA is taken as given** and is not re-tested. This folder asks a different question: do **ELF3 / GRHL1 / KLF4 regulon activities** track **TACSTD2 / CLDN4** in public lung epithelium?

**No ChIP peaks are invented.** There is no public human lung ELF3 ChIP (A10). Motif cisTarget databases were not downloaded. Regulons are public curated TF–target edges and/or co-expression in this matrix.

Outputs: [`results/scrna_scenic/`](../../results/scrna_scenic/).

## Question (pre-specified)

1. Score ELF3, GRHL1, and KLF4 regulons with **AUCell** in epithelial / malignant-like cells.
2. Correlate AUCell with TACSTD2 and CLDN4.
3. Combinatorial split: **all histologies vs LUAD-only**.
4. Report **honest n**: sample-level is primary; cell-level ρ is exploratory (pseudoreplication).

## Dataset

**GSE207422** (Hu et al., *Genome Med* 2023, PMID 36869384). Neoadjuvant PD-1 + chemo NSCLC. Public processed UMI matrix (175 MB gzip, 92,330 cells) + sample metadata with **Adeno vs Squamous**.

| Alternative | Why not primary |
| --- | --- |
| GSE253013 | 9.3 GB RDS; LUAD-only (no LUAD vs all contrast) |
| GSE131907 | LUAD atlas, no ICI; 2.9 GB log2TPM txt / 390 MB raw UMI. Not required once GSE207422 provides the split |

Author **CopyKAT malignant labels are not on GEO**. Malignant-like = epithelial lineage (marker argmax) AND normal-lung score ≤ 75th percentile of epithelial cells. That is a proxy, not CopyKAT.

## GRN method

pySCENIC **imports** in this environment. Full SCENIC (`grn → ctx → aucell`) was **not** run: cisTarget feather rankings are multi-GB and would be used to *imply* motif support at TACSTD2/CLDN4. That is too close to a binding claim we cannot make.

What *was* run:

| Layer | What it is | What it is not |
| --- | --- | --- |
| Public prior | TRRUST v2 + DoRothEA (OmniPath) + CollecTRI (OmniPath) union per TF | Lung ChIP, motif scan, or a binding map |
| Pearson | TF vs genes in GSE207422 epithelial log1p-CP10k; top 50 with r ≥ 0.10 | Causal regulation |
| GRNBoost2 (optional) | arboreto, three TFs only, subsampled cells | cisTarget-pruned SCENIC regulon |
| AUCell | Aibar 2017 recovery-curve score on those gene sets | The R `AUCell` binary / pySCENIC CLI |

**TACSTD2 and CLDN4 are held out of every regulon** before AUCell so the correlation is not circular. Public priors **do not list** TACSTD2 or CLDN4 as ELF3/GRHL1/KLF4 targets (see `resources/tf_targets_public.tsv`). GRHL1 has only 2–3 prior genes; its prior AUCell is under-specified and the Pearson regulon is the informative GRHL1 test.

## Units and verdict rules (locked before results)

- **Primary:** Spearman ρ of *sample means* (malignant-like cells; sample kept if ≥20 such cells).
- **Exploratory:** cell-level Spearman (report n; do not treat p as a patient-level claim).
- **SUPPORT:** ρ ≥ 0.30 and p < 0.05 for **both** TACSTD2 and CLDN4, same regulon, sample-level.
- **PARTIAL:** one target or one split (all vs LUAD) meets SUPPORT.
- **NOT_SUPPORTED:** neither target meets SUPPORT.
- **UNDERPOWERED:** n_samples < 6.

## Run

```bash
python3 methods/scrna_scenic/scripts/download_gse207422.py --outdir data/scrna_scenic/GSE207422
python3 methods/scrna_scenic/scripts/analyze.py
# optional: SKIP_GRNBOOST=1 python3 methods/scrna_scenic/scripts/analyze.py
```

Refresh priors (not required; committed snapshot is enough):

```bash
python3 methods/scrna_scenic/scripts/download_priors.py
```

## Results

**Verdict: `PARTIAL`.** Full table: [`results/scrna_scenic/README.md`](../../results/scrna_scenic/README.md).

| split | n samples | call | one-line |
| --- | ---: | --- | --- |
| all (Adeno+Squamous) | 12 | SUPPORT | Pearson regulons and KLF4/combinatorial priors track both genes; **ELF3_prior does not** (ρ=0.08 / 0.01) |
| LUAD-only | 6 | PARTIAL | ELF3 RNA still tracks (A10); regulon AUCell mostly fails p<0.05 |

Honest n: sample-level is the claim. Cell-level malignant-like n=9782 (all) / 2360 (LUAD) inflates p (ELF3_prior vs CLDN4: cell ρ=0.11 p=1e-28 vs sample ρ=0.007 n=12).

GRNBoost2 did not run (arboreto/dask `diagnostics_port` error on this stack). Intersection ELF3 prior ∩ Pearson = 0 genes.
