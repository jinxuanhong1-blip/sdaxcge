# Given tables (not re-audited)

PR #459 pair **GSE131907 + GSE189357** already differs on malignant CLDN4 %pos vs T/NK
(n=30, ρ=−0.542; Q4 r=−0.619). Those T/NK numbers are **given**. This folder does
not re-score T/NK and does not add GSE148071.

- `GSE131907_samples.tsv` — author-malignant sample-level CLDN4 %pos (PR #459 / #279 / #320).
- `GSE189357_marker_units.tsv` — marker-malignant patient-level CLDN4 %pos
  (gate: `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`).

GEO histology for GSE189357 (series matrix `!Sample_characteristics_ch1`, typo
"histolgical type"): TD1 IAC, TD2 IAC, TD3 MIA, TD4 MIA, TD5 AIS, TD6 MIA,
TD7 AIS, TD8 AIS, TD9 IAC. Three each of AIS / MIA / IAC (Zhu / Fan / Jiang,
PMID 36434043).
