# RESULTS — STORY GAP: KL/STK11 cold + Tacstd2/Cldn4 high

## Purpose

PPT-ready **disease context** slide pack. Assemble already-locked public numbers so the story can say cold STK11/KL **and** Tacstd2/Cldn4-high **without fabricating human antigen elevation**.

## Verdict

| Limb | Public support | Where it lives |
|---|---|---|
| STK11/KL cold / ICI-resistant | **YES** | Skoulidis 2018 + MSK KRAS-restricted DCB (PR #49) |
| Tacstd2/Cldn4 high | **YES — mouse cell-line bulk lock** | GSE137244 Δ +3.238 / +5.570; MW p=0.00794 |
| Human STK11 → high TACSTD2/CLDN4 | **NO (reverse in TCGA bulk)** | PR #49 |
| Human STK11 scRNA antigen genotype | **NO powered claim** | GSE280232 (PR #113) |

## Deliverables

| File | Use |
|---|---|
| `PPT_EVIDENCE_TABLE.md` | Slide paste |
| `tables/evidence_table.tsv` | Full rows |
| `provenance.json` | PR / PMID citations |
| `figures/honesty_matrix.png` | Bulk vs scRNA honesty grid |
| `FINDING.md` | Narrative |

## Reproduce display only

No analysis rerun. Display:

```bash
python3 results/story_gap_stk11_cold_tacstd2/scripts/make_honesty_figure.py
```
