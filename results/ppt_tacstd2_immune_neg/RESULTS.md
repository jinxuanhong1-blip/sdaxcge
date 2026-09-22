# PPT corrected: TACSTD2-high ↔ immune NEGATIVE (public multi-cohort)

## Purpose

Build one PPT-ready evidence table for the **corrected opening claim**:

> **TROP2 / TACSTD2-high ↔ immune NEGATIVE** (fewer T/NK, lower CD8 / ImmuneScore).  
> **Not** drug resistance / ICB failure / Tacstd2↑ after therapy.

Required limbs (all public): **CosMx, concordant-4, TCGA, OncoSG, CPTAC, mouse public.**

No new downloads. No fabricated statistics. Display files:

| File | Role |
|---|---|
| `PPT_EVIDENCE_TABLE.md` | paste into slide |
| `evidence_table.tsv` | full rows + caveats |
| `provenance.json` | PR + path + exact fields |

## Overall calls (TACSTD2 only)

| Cohort | Supports immune-negative? |
|---|---|
| CosMx short-range CD8+NK / immune | **Yes** (10–20 µm; 8/8 & 5/5) |
| concordant-4 TACSTD2 vs T/NK | **Null** (ρ=−0.112); CLDN4 −0.531 is separate |
| TCGA CD8 ‖ ABSOLUTE | **Yes** (pooled −0.191; LUSC stronger) |
| TCGA ImmuneScore ‖ ABSOLUTE | **Partial** (LUSC −0.131; LUAD +0.072) |
| TCGA 8-cohort keratin-adj CD8 | **Weak** (ρ=−0.069; 6/8) |
| OncoSG CD8 / T / NK ‖ purity | **Yes** (−0.31 / −0.40 / −0.42) |
| CPTAC TROP2 protein | **Partial** (LUAD xCell yes; LSCC null) |
| Mouse public pool | **Weak / ns** (c=−0.277, p=0.137) |

## Explicitly out of frame

- **TISMO 49/64 Tacstd2↑ after ICB** — drug/ICB induction, not fewer T/NK (see excluded row X).
- Human ICI bulk response (GSE126044, GSE135222, etc.) — underpowered / mixed; not this slide.
- Visium same-spot correlation — not spatial exclusion.
- Private 8KL — not in this repo; not merged with public mouse.

## Relation to PR #735

PR #735 was a first funnel table. It omitted **concordant-4** and **mouse public**, and kept **TISMO** in the main table. This corrected pack fills the six required cohorts and removes drug-resistance framing.

## Reproduce display only

```bash
# No analysis rerun — numbers are copied from cited PRs.
python3 -c "print(open('results/ppt_tacstd2_immune_neg/evidence_table.tsv').readline())"
```
