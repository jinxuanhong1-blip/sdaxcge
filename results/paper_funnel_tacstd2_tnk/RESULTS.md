# Tacstd2-high → fewer T/NK: public paper-funnel evidence table

## Purpose

One PPT-ready corroboration table for the claim **Tacstd2 / TACSTD2-high associates with fewer T/NK**, using only already-computed public analyses:

| Required source | Used |
|---|---|
| CosMx He2022 | yes (PR #726) |
| TCGA TACSTD2 vs CD8, keratin-adjusted | yes (PR #593 primary; #705 sensitivity) |
| OncoSG | yes (PR #139) |
| CPTAC TROP2 protein | yes if present (PR #99 xCell; #721 CD8A protein) |
| TISMO | yes (PR #152 / #173 locked 49/64; #727 CD8 honesty) |

No new downloads. No fabricated statistics. Machine table: `evidence_table.tsv`. Slide paste: `PPT_EVIDENCE_TABLE.md`.

## Overall verdict

| Limb | Supports fewer T/NK? |
|---|---|
| CosMx short-range CD8+NK | **Yes** |
| OncoSG CD8 / IMSIG T+NK | **Yes** |
| TCGA keratin-adj CD8 | **Weak yes** (ρ≈−0.07; 6/8) |
| CPTAC TROP2 protein | **Partial** (LUAD xCell CD8); **No** for CD8A protein |
| TISMO | **No** for T/NK; **Yes** for Tacstd2↑ after ICB |

## Provenance

See `provenance.json`. All numeric cells cite a PR path that already exists on GitHub.

## Reproduce display only

```bash
# No analysis rerun. Render TSV → markdown if needed:
python3 -c "import csv; print(open('results/paper_funnel_tacstd2_tnk/evidence_table.tsv').read()[:200])"
```
