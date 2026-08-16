# OncoSG A1 rework — how to run

Self-contained recompute of TACSTD2 vs CD8 / GEP18 / immune after purity
on public cBioPortal `luad_oncosg_2020`.

```bash
pip install -r scripts/rework/OncoSG_A1/requirements.txt
python3 scripts/rework/OncoSG_A1/analyze.py \
  --cache-dir /tmp/oncosg_a1 \
  --out-dir results/rework/OncoSG_A1
```

The script downloads the revision-pinned cBioPortal Datahub files if they
are not already in `--cache-dir`. The 23 MB expression matrix is **not**
committed.
