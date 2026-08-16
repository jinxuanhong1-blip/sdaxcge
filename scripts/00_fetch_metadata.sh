#!/usr/bin/env bash
# Fetch recount3 SRA project indexes and per-study SRA metadata.
# This is how the local cache under $HOME/hunt_data was built. The analysis
# scripts read that cache; they do not need this script unless you rebuild it.
set -euo pipefail
ROOT="${HUNT_DATA:-$HOME/hunt_data}"
mkdir -p "$ROOT"
cd "$ROOT"
curl -sS --retry 3 -o human.recount_project.MD.gz \
  https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/sra/metadata/sra.recount_project.MD.gz
curl -sS --retry 3 -o mouse.recount_project.MD.gz \
  https://recount-opendata.s3.amazonaws.com/recount3/release/mouse/data_sources/sra/metadata/sra.recount_project.MD.gz
python3 - <<'PY'
import gzip, collections, os
from pathlib import Path
root = Path.cwd()
for org, base in [
    ("human", "https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/sra/metadata"),
    ("mouse", "https://recount-opendata.s3.amazonaws.com/recount3/release/mouse/data_sources/sra/metadata"),
]:
    c = collections.Counter()
    with gzip.open(root / f"{org}.recount_project.MD.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i = header.index("project")
        for line in fh:
            c[line.split("\t")[i]] += 1
    keep = [p for p, n in c.items() if n >= 6]
    (root / "meta" / org).mkdir(parents=True, exist_ok=True)
    with open(root / f"urls_{org}.txt", "w") as out:
        for p in sorted(keep):
            out.write(f'url = "{base}/{p[-2:]}/{p}/sra.sra.{p}.MD.gz"\n')
            out.write(f'output = "meta/{org}/{p}.MD.gz"\n')
    print(org, "studies>=6", len(keep))
PY
curl -sS --parallel --parallel-max 32 --retry 3 --retry-delay 2 --config urls_human.txt
curl -sS --parallel --parallel-max 32 --retry 3 --retry-delay 2 --config urls_mouse.txt
