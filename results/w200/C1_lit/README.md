# C1 public literature catalog

Wave: `w200`  
Claim slice: **C1** — TROP2 ADC internalization and CLDN4 trafficking  
Mode: **catalog only**. No new experiments. No invented numbers.

## Scope

Two literature arms, plus their documented intersection:

1. **TROP2 ADC internalization** — antibody/ADC binding, endocytosis, intracellular trafficking, lysosomal delivery, and payload release for TROP2-directed agents (RS7/hRS7, sacituzumab govitecan, datopotamab deruxtecan, RN927C/PF-06664178).
2. **CLDN4 trafficking** — internalization, endocytic route, recycling vs degradation, and phosphorylation-linked removal of claudin-4 from the plasma membrane / tight junction.
3. **TROP2–claudin intersection** — physical proximity or binding of TROP2/TACSTD2 to claudins, including CLDN4. This is catalogued separately because it is **not** the same as measured co-internalization of a TROP2 ADC with CLDN4.

## Extraction rules

- A number appears in this folder only if it was read in a retrieved abstract or full text.
- Review-cited numbers that were **not** re-read in the cited primary PDF are tagged `review_cited_unverified_in_primary`.
- Abstract-only access is tagged `abstract_only`.
- Full-text access is tagged `full_text`.
- Qualitative statements are allowed when the source states them; they are not converted into rates, half-lives, fold-changes, or IC50s.
- Missing table cells, supplementary figures not retrieved, and secondary-blog numbers are **not** filled in.
- Clinical efficacy (PFS/OS, ORR) is out of scope unless a paper uses it only to identify the ADC.

## Files

| File | Contents |
| --- | --- |
| `catalog.tsv` | One row per catalogued record (PMID, DOI, topic, access, quoted numbers) |
| `trop2_adc_internalization.md` | TROP2 antibody/ADC internalization notes |
| `cldn4_trafficking.md` | CLDN4 endocytosis / trafficking notes |
| `trop2_cldn_intersection.md` | TROP2–claudin binding / proximity notes |
| `gaps.md` | Explicitly not found; do not infer |

## Search snapshot

PubMed (2026-08-16, `esearch`):

- `TROP2 OR Trop-2 OR TACSTD2 AND (internalization OR endocytosis OR lysosome) AND (ADC OR antibody-drug OR hRS7 OR sacituzumab OR datopotamab)` → **45** hits (query also matches the journal abbreviation ADC; list was manually filtered).
- `claudin-4 OR CLDN4 AND (internalization OR endocytosis OR trafficking OR recycling)` → **50** hits.

This catalog is a curated subset of those hits plus core primary papers cited by them. It is not a systematic review and is not exhaustive.
