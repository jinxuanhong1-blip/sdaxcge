#!/usr/bin/env python3
"""Parse GSE253013 series matrix: 89 lanes, 9 patients, no response labels."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default=Path("results/gse253013"))
    args = ap.parse_args()
    text = gzip.open(args.matrix, "rt").read().splitlines()
    rows: dict[str, list[str]] = {}
    for line in text:
        if line.startswith("!") and "\t" in line:
            k, *vals = line.split("\t")
            rows[k] = [v.strip('"') for v in vals]
    titles = rows["!Sample_title"]
    accs = rows["!Sample_geo_accession"]
    groups = []
    for k, v in rows.items():
        if k.startswith("!Sample_characteristics_ch1") and v and v[0].startswith("group:"):
            groups = [x.split(":", 1)[1].strip() for x in v]
            break
    if not groups:
        groups = ["Adjacent Non-Tumor" if "ANT" in t.upper() else "Tumor" for t in titles]

    out_rows = []
    for gsm, title, group in zip(accs, titles, groups):
        pid_m = re.search(r"(MRC0*\d+)", title.replace(" ", ""))
        patient = pid_m.group(1) if pid_m else "unknown"
        tissue = "ANT" if "ANT" in title.upper() or "Adjacent" in group else "Tumor"
        cd45 = "CD45+" in title
        out_rows.append(
            {
                "gsm": gsm,
                "title": title,
                "patient": patient,
                "tissue": tissue,
                "group": group,
                "cd45_sorted": cd45,
            }
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    tsv = args.outdir / "sample_metadata.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    response_keys = [
        k
        for k in rows
        if any(s in k.lower() for s in ("response", "mpr", "pcr", "treatment", "ici", "pd-1", "pdl"))
    ]
    audit = {
        "n_gsm": len(out_rows),
        "n_patients": len({r["patient"] for r in out_rows}),
        "patients": dict(Counter(r["patient"] for r in out_rows)),
        "tissue": dict(Counter(r["tissue"] for r in out_rows)),
        "cd45_sorted_lanes": sum(r["cd45_sorted"] for r in out_rows),
        "response_like_geo_keys": response_keys,
        "public_mpr_r_labels": False,
        "blood_only": False,
        "source": "GSE253013_series_matrix.txt.gz",
    }
    (args.outdir / "geo_label_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
