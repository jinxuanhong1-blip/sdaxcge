# B5 — Nathanson melanoma ICI / CLDN4

**Result:** an open matrix exists, but the usable pretreatment cohort is tiny. I used the authors’ public Cufflinks FPKM matrix and outcome table for Nathanson et al. 2017 (melanoma treated with CTLA-4 blockade), pinned at repository commit `2ce35b9`.

Among the 9 pretreatment RNA samples, CLDN4 was higher in patients labeled as having clinical benefit (n=4; median 6.03 FPKM) than without benefit (n=5; median 0.129 FPKM). The exact two-sided Mann–Whitney result was U=18, p=0.0635; rank AUC=0.90 and rank-biserial correlation=0.80. Thus the direction is interesting, but this is **not statistically conclusive** and is highly sensitive to individual samples. Post-treatment biopsies were not mixed into the predictive comparison.

The source contains two gene IDs labeled CLDN4. Following the paper’s stated procedure, duplicate symbols were collapsed by median FPKM; using canonical `ENSG00000189143` alone changes scale, not ranks or p-value. “Benefit” is the deposited binary endpoint, not RECIST response.

This is the Tavi Nathanson/MSK cohort, not a UCLA cohort. No UCLA dataset was substituted.

Sources: [paper](https://doi.org/10.1158/2326-6066.CIR-16-0019) · [data repository](https://github.com/hammerlab/melanoma-reanalysis)

Reproduce: `python3 analysis.py --output-dir .`
