# Disease-context story gap (public)

Assemble **KL/STK11 cold** + **Tacstd2/Cldn4 high** for one PPT slide without fabricating human antigen elevation.

| File | Role |
|---|---|
| `PPT_EVIDENCE_TABLE.md` | Slide paste |
| `FINDING.md` | Narrative |
| `RESULTS.md` | PR summary |
| `tables/evidence_table.tsv` | Machine table |
| `provenance.json` | Citations |
| `figures/honesty_matrix.png` | Bulk vs scRNA grid |

```bash
pip install -r requirements.txt
python3 scripts/make_honesty_figure.py
```
