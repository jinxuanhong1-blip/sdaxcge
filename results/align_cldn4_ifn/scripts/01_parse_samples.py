#!/usr/bin/env python3
"""Extract per-sample metadata (title/source/characteristics/suppl) from GEO SOFT."""
import csv
import glob
import os
import re

META = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "meta"))
KEEP = ("title", "source_name_ch1", "characteristics_ch1", "description",
        "organism_ch1", "supplementary_file", "supplementary_file_1",
        "label_ch1", "label_ch2", "source_name_ch2", "characteristics_ch2")

rows = []
for path in sorted(glob.glob(os.path.join(META, "*.all.txt"))):
    gse = os.path.basename(path).split(".")[0]
    gsm = None
    cur = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("^SAMPLE"):
                if gsm:
                    rows.append((gse, gsm, cur))
                gsm = line.split("=")[1].strip()
                cur = {}
            elif line.startswith("!Sample_") and gsm:
                m = re.match(r"!Sample_(\w+) = (.*)", line.rstrip("\n"))
                if m and m.group(1) in KEEP:
                    k, v = m.group(1), m.group(2)
                    cur[k] = (cur.get(k, "") + " | " + v).strip(" |") if k in cur else v
            elif line.startswith("!sample_table_begin"):
                # skip the data block entirely (can be tens of MB)
                for l2 in fh:
                    if l2.startswith("!sample_table_end"):
                        break
    if gsm:
        rows.append((gse, gsm, cur))

out = os.path.join(META, "sample_design.tsv")
cols = ["gse", "gsm"] + list(KEEP)
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for gse, gsm, cur in rows:
        w.writerow([gse, gsm] + [cur.get(k, "") for k in KEEP])
print(f"wrote {out}: {len(rows)} samples")
for gse in sorted({r[0] for r in rows}):
    sub = [r for r in rows if r[0] == gse]
    print(f"\n=== {gse} (n={len(sub)}) {sub[0][2].get('organism_ch1','')}")
    for _, gsm, cur in sub:
        print(f"  {gsm}\t{cur.get('title','')}\t[{cur.get('source_name_ch1','')}]\t{cur.get('characteristics_ch1','')}")
