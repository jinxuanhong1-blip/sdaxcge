# Pair GSE189357 + GSE205335 — CLDN4-only REAL Slingshot / PAGA

ADDITIVE. **CLDN4 only.** The PR #459 pair that already differs (malignant CLDN4 %pos vs T/NK, n=31, ρ=−0.478) is taken as given and is **not** re-ranked.

- Epithelium only. No dual-high TACSTD2∩CLDN4.
- Root / Slingshot start cluster is **not** CLDN4-high.
- Programs along pseudotime: CLDN4, barrier/keratin (CLDN4 held out), IFN (Hallmark IFNα∩IFNγ).
- Inferential unit = patient (9 + 22). Done when `results/tables/lineage_table.tsv` exists.

```bash
pip install -r methods/pair_189357_205335_slingshot_cldn4/requirements.txt
python3 methods/pair_189357_205335_slingshot_cldn4/scripts/download.py
python3 methods/pair_189357_205335_slingshot_cldn4/scripts/extract_epithelium.py
python3 methods/pair_189357_205335_slingshot_cldn4/scripts/analyze.py
```
