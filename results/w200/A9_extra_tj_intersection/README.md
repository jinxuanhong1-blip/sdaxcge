# A9 extra-cohort supplement — TJ genes in TACSTD2-high vs low LUAD

**3 of 83 catalog genes pass in all three extra cohorts:
CLDN4, MICALL2, PARD6B.** 20 genes pass in ≥2/3.
OncoSG ∩ CPTAC is the same 3 genes
as the 3/3 set — CPTAC is the limiting cohort.

User A9 (intersection CLDN1/4/7, F11R, PARD3) is **taken as given**.
This table is **not** a check of that slide list. It asks a different
question: in three extra public LUAD RNA cohorts, which genes from a
published tight-junction catalog are up in TACSTD2-high tumors, and
which of those calls recur?

## Honest overlap

| Item | Count | Genes |
|------|------:|-------|
| TJ universe (pre-specified) | 83 | published GO:0005923 ∪ R-HSA-420029 ∪ JAM2 |
| Pass OncoSG LUAD | 21 / 82 measured | see per-cohort |
| Pass CPTAC LUAD RNA | 13 / 83 measured | see per-cohort |
| Pass GSE31210 LUAD | 26 / 77 measured | see per-cohort |
| Recur in **3/3** extra cohorts | 3 | CLDN4, MICALL2, PARD6B |
| Recur in **≥2/3** extra cohorts | 20 | CLDN4, MICALL2, PARD6B, AMOT, ANK3, CGN, CLDN1, CLDN12, CLDN3, CLDN6, CLDN7, CRB3, CXADR, CYTH3, EPCAM, F11R, MARVELD3, OCLN, PATJ, TJP3 |
| Pairwise OncoSG ∩ CPTAC | 3 | CLDN4, MICALL2, PARD6B |
| Pairwise OncoSG ∩ GSE31210 | 15 | CGN, CLDN1, CLDN12, CLDN3, CLDN4, CLDN6, CLDN7, CRB3, CXADR, EPCAM, MARVELD3, MICALL2, PARD6B, PATJ, TJP3 |
| Pairwise CPTAC ∩ GSE31210 | 8 | AMOT, ANK3, CLDN4, CYTH3, F11R, MICALL2, OCLN, PARD6B |

Empty cells are real negatives, not missing analyses.

Median-split sensitivity (does **not** define the table): 3/3 =
CLDN1, CLDN12, CLDN4, CRB3, MICALL2 (n=5); ≥2/3 n=20.
CLDN4 and MICALL2 are in both the tertile and median 3/3 sets. PARD6B is
tertile-only 3/3. The exact three-gene list is split-dependent; that is
reported rather than resolved by picking the larger set.

### Per-cohort passes (tertile, FDR across the TJ catalog)

- **OncoSG LUAD** (n=169, high=57, low=57): 82/83 TJ genes measured; **21 pass**. CGN, CLDN1, CLDN12, CLDN16, CLDN3, CLDN4, CLDN6, CLDN7, CLDN8, CLDN9, CRB3, CXADR, EPCAM, MARVELD3, MICALL2, PARD6A, PARD6B, PATJ, TJP3, VAPA, YBX3
- **CPTAC LUAD RNA** (n=110, high=37, low=37): 83/83 TJ genes measured; **13 pass**. AMOT, AMOTL2, ANK3, CLDN4, CYTH3, F11R, FRMD4B, MICALL2, OCLN, PARD6B, RAP2B, TJP1, TJP2
- **GSE31210 LUAD** (n=226, high=76, low=76): 77/83 TJ genes measured; **26 pass**. AMOT, ANK3, CGN, CLDN1, CLDN12, CLDN23, CLDN3, CLDN4, CLDN6, CLDN7, CRB3, CXADR, CYTH3, EPCAM, F11R, MARVELD2, MARVELD3, MICALL2, MPP7, OCLN, PARD3, PARD6B, PATJ, POF1B, STRN, TJP3

### A9 slide genes — annotation only, not the finding

These five genes were not used to define the universe or the overlap.
They are listed so a reader can see where they landed under the broader FDR.

- CLDN4: 3/3 extra cohorts (OncoSG pass, CPTAC pass, GSE31210 pass)
- CLDN1: 2/3 extra cohorts (OncoSG pass, CPTAC no-call, GSE31210 pass)
- CLDN7: 2/3 extra cohorts (OncoSG pass, CPTAC no-call, GSE31210 pass)
- F11R: 2/3 extra cohorts (OncoSG no-call, CPTAC pass, GSE31210 pass)
- PARD3: 1/3 extra cohorts (OncoSG no-call, CPTAC no-call, GSE31210 pass)

## What this is

- A supplement table of **recurrent TJ-catalog genes** in TACSTD2-high vs
  TACSTD2-low tumors in three extra public LUAD RNA cohorts.
- Call rule matches A9 in spirit (ρ > 0, tertile log2FC > 0, Welch BH-FDR < 0.05)
  with FDR across the full TJ catalog in that cohort, not across five slide genes.

## What this is not

- Not a replication test or revision of User A9.
- Not GSEA, not protein, not scRNA, not a TROP2-binding result.
- Not GSE72094 (Moffitt / US LUAD). The East-Asian GEO used here is GSE31210.
- OncoSG values are z-scores; the high-vs-low difference is in z units, not
  log2 RNA. Direction is still valid. Spearman is rank-invariant.

## Pre-specified methods

| Item | Choice |
|------|--------|
| Splitter | TACSTD2, within-cohort tertiles (high ≥ 2/3, low ≤ 1/3) |
| Universe | GO:0005923 ∪ Reactome R-HSA-420029 ∪ JAM2 (n=83) |
| Cohorts | OncoSG LUAD; CPTAC LUAD tumor RNA; GSE31210 primary tumors |
| Association | Spearman ρ vs continuous TACSTD2 |
| High-vs-low | Welch t-test on the cohort native scale |
| Call | ρ>0 AND log2FC>0 AND Welch BH-FDR<0.05 |
| FDR | BH across TJ genes **present in that cohort** |
| Recur | 3/3 and ≥2/3 extra cohorts (pairwise counts also reported) |
| Sensitivity | median split (labeled; does not define recurrence) |

## Sources

- OncoSG LUAD (Chen et al. 2020) via cBioPortal `luad_oncosg_2020`
- CPTAC LUAD RNA, LinkedOmics / S3 freeze v1.2 tumor RSEM log2
- GSE31210 (Okayama et al. *Cancer Res* 2012), GPL570, tumors only

## Reproduce

```bash
python3 scripts/w200/A9_extra_tj_intersection/download.py
python3 scripts/w200/A9_extra_tj_intersection/analyze.py
```
