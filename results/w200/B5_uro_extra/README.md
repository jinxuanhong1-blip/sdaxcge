# B5 extra: Mariathasan / Snyder urothelial ICI RNA

Public pretreatment RNA for **TACSTD2** (TROP2) and **CLDN4** versus
atezolizumab response in the two named urothelial cohorts. Written to be
conservative. Nulls are nulls.

## Bottom line

IMvigor210 (298 RECIST-evaluable of 348 RNA samples): TACSTD2 AUC 0.556
(p=0.16), CLDN4 AUC 0.553 (p=0.19). Median-split response ORs are 1.47 and
1.59 — high expression, if anything, slightly *better* response, not worse.
Snyder 2017 (21 evaluable of 25) is uninformative (AUC 0.53 and 0.45; wide
CIs). Immune-gene controls (CXCL9, CD8A) work on the same IMvigor210 matrix.

See `RESULTS.md`.

## Cohorts

1. **IMvigor210** — Mariathasan et al., *Nature* 2018
   ([doi:10.1038/nature25501](https://doi.org/10.1038/nature25501)).
   Atezolizumab in metastatic urothelial carcinoma
   (NCT02108652 / NCT02951767). Public object: Genentech
   `IMvigor210CoreBiologies` 1.0.1
   ([research-pub.gene.com](http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/)).
   Counts and RECIST / OS / subtype annotations come from the package `cds`
   CountDataSet (EGA EGAS00001002556 is the raw archive; this analysis uses
   the public processed package, not a new EGA download).

2. **Snyder 2017** — Snyder et al., *PLoS Med* 2017
   ([doi:10.1371/journal.pmed.1002309](https://doi.org/10.1371/journal.pmed.1002309)).
   Atezolizumab in metastatic urothelial carcinoma. Public kallisto gene
   counts and clinical table from
   [hammerlab/multi-omic-urothelial-anti-pdl1](https://github.com/hammerlab/multi-omic-urothelial-anti-pdl1).

## Methods (prespecified)

- Expression: `log2(CPM + 1)` from raw counts (IMvigor210) or kallisto
  `est_counts` (Snyder) and the sample library size. One gene at a time;
  library-size CPM is enough for a between-sample rank test.
- Primary endpoint: investigator/IRF RECIST CR/PR vs SD/PD. IMvigor210 `NE`
  (n=50) and Snyder “only scanned at/baseline” (n=4) are excluded.
- Tests: two-sided Mann–Whitney; ROC AUC with 2000-bootstrap 95% CI;
  median-split 2×2 Fisher exact + Woolf-style OR CI via `Table2x2`.
- Secondary, labeled as such: CR/PR vs PD only; DESeq size-factor scale;
  Cox OS per log2CPM; Snyder paper `is_benefit`; TCGA subtype / immune
  phenotype strata.
- Positive controls on the same matrix: CD8A, CXCL9, mean(GZMA, PRF1).

No threshold was optimized. No model was trained. Exploratory strata are
not discovery.

## Reproduce

```sh
python3 -m pip install -r requirements.txt
python3 analyze.py
```

This reads the small tables in `processed/` (gene-level extracts + clinical).
To rebuild those extracts from the public sources (needs `Rscript` and
network):

```sh
python3 analyze.py --rebuild-processed
```

Source checksums are in `processed/provenance.tsv`. The full 70 MB
IMvigor210 package is not stored in git.

## Outputs

- `RESULTS.md` — interpretation
- `tables/summary.json` — headline numbers
- `tables/tests.csv` — every test that was run
- `tables/sample_imvigor210.csv`, `tables/sample_snyder2017.csv` — auditable
  per-sample values
- `tables/imvigor210_subtype_context.csv` — ORR and gene medians by TCGA
  subtype
- `figures/response_boxplots.png` — primary boxplots
- `figures/imvigor210_os_km.png` — OS by median split
- `figures/median_split_or.png` — response ORs including controls

## References

- Mariathasan S, Turley SJ, Nickles D, et al. TGFβ attenuates tumour
  response to PD-L1 blockade by contributing to exclusion of T cells.
  *Nature*. 2018;554:544–548.
- Snyder A, Nathanson T, Funt SA, et al. Contribution of systemic and
  somatic factors to clinical response and resistance to PD-L1 blockade in
  urothelial cancer: An exploratory multi-omic analysis. *PLoS Med*.
  2017;14:e1002309.
