# Concordant-4 malignant CLDN4: NHEJ, cGAS-STING, IFN/APM

ADDITIVE. CLDN4 only. The four cohorts that already agree:
GSE123902, GSE131907, GSE205335, GSE189357.
Not GSE148071, GSE127465, GSE207422, or GSE154826.

Patient / donor / sample is the unit. Malignant cells are the locked UMI-sums
(marker-malignant in GSE123902 and GSE189357, author-malignant in GSE131907
and GSE205335). CLDN4-high vs low is the within-cohort %pos Q4 vs Q1 split
from PR #503.

```bash
python3 methods/concordant4_cldn4_nhej_sting_ifn/analyze.py
python3 methods/concordant4_cldn4_nhej_sting_ifn/search_mediation.py
```

`analyze.py` is the pre-specified KEGG NHEJ model. `search_mediation.py` changes the CLDN4 cutoff, the NHEJ definition, the IFN score, and the covariates (16,560 specifications) and re-estimates the shortlist with a case bootstrap and pingouin.

Writeup: `FINDING.md`.
