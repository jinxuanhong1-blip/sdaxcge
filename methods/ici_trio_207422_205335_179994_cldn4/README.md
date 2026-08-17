# ICI trio CLDN4-only (GSE207422 + GSE205335 + GSE179994)

Additive CLDN4-only ICI-adjacent scRNA merge. No dual-high. Include 207422.

If GSE179994 has no public processed all-cell matrix <2 GB, the merge still
runs on GSE207422 + GSE205335 and says so.

```bash
python3 -m pip install -r methods/ici_trio_207422_205335_179994_cldn4/requirements.txt
python3 methods/ici_trio_207422_205335_179994_cldn4/download.py
python3 methods/ici_trio_207422_205335_179994_cldn4/analyze.py
```

See `FINDING.md` for honest n, the combo table, and the outgoing LR table.
