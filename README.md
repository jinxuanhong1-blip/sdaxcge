# TACSTD2 surface-gene co-expression in TCGA-BRCA (B1 analog)

Honest ranking of how well a candidate surface marker co-expresses with an
anchor gene, across the entire cell-surface proteome, in TCGA breast cancer.

**Driving question:** In TCGA-BRCA, is `CLDN4` the *top* co-expression partner
of `TACSTD2` (TROP2) among surface genes?

**Answer: No.** Using Spearman correlation across 1,097 primary-tumour samples
and the 2,618 surfaceome genes present in the matrix, `CLDN4` ranks **#4**
(ρ = 0.348, bootstrap 95% CI 0.29–0.40, FDR q ≈ 1e-29). The genes ahead of it
are `EFNA1` (ρ = 0.429), `PVRL4` / NECTIN4 (0.351) and `EFNA4` (0.350).

That #4 rank is stable: one-vial-per-patient is a no-op on this extract,
ABSOLUTE-purity partial correlation barely moves the pair (0.346 → 0.343),
and CLDN4 does not become #1 inside any PAM50 stratum (closest: Basal, #2
behind `SLC39A2`). Among *claudins* CLDN4 *is* #1; among all expressed genes
it is #30 / 17,659. Adjacent-normal ρ (0.80) is *higher* than the tumour pair,
so the coupling looks like a normal-epithelium program diluted in tumours.

Full write-up and tables: [`results/w200/B1_BRCA/`](results/w200/B1_BRCA/).

## What "honest ranking" means here

Every surface gene is scored against the anchor and the focus gene's true
position is reported. No hand-picked neighbourhood, both Spearman (primary)
and Pearson, the actual #1 gene named next to the focus rank, plus the
robustness checks used in the PAAD analog (bootstrap CI, purity-partial,
subtype strata, adjacent-normal, transcriptome-wide rank).

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/w200/B1_BRCA/download_data.py   # → data/ (git-ignored)
python3 scripts/w200/B1_BRCA/run_analysis.py    # → results/w200/B1_BRCA/
python3 scripts/w200/B1_BRCA/test_ranking.py    # synthetic-data unit checks
```

## Data sources

- **Expression:** UCSC Xena `TCGA.BRCA.sampleMap/HiSeqV2` (log2 norm_count+1, HGNC symbols).
- **Subtypes:** Xena `BRCA_clinicalMatrix` (`PAM50Call_RNAseq`).
- **Purity:** PanCanAtlas ABSOLUTE (GDC open file `4f277128-f793-4354-a13d-30cc7fe9f6b5`).
- **Surfaceome:** Bausch-Fluck et al., *PNAS* 2018 table S3 (GitHub mirror).
