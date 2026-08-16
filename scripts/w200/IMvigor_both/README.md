# w200 / IMvigor_both

Public IMvigor210 RNA: **TACSTD2 and CLDN4** (plus an exploratory both-high
AND gate) versus confirmed **ORR** and **OS**.

```bash
python3 -m pip install -r scripts/w200/IMvigor_both/requirements.txt
python3 scripts/w200/IMvigor_both/download.py
python3 scripts/w200/IMvigor_both/analyze.py
```

Outputs land in `results/w200/IMvigor_both/`. Requires R only for the
CountDataSet dump (`extract_cds.R`). Raw tarball is gitignored.
