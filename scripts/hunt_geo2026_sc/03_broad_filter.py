#!/usr/bin/env python3
"""Filter the broad 2026 lung scRNA search for IO-related series not already listed."""
import csv
import json
import re

LISTED = {l.strip() for l in open("results/hunt_geo2026_sc/already_listed.txt") if l.strip()}
IO = re.compile(
    r"immunotherap|checkpoint|PD-1|PD1|PD-L1|PDL1|pembrolizumab|nivolumab|atezolizumab|"
    r"durvalumab|sintilimab|tislelizumab|camrelizumab|toripalimab|ipilimumab|cemiplimab|"
    r"ICI\b|ICB\b|anti-PD", re.I)
RESP = re.compile(r"respond|MPR|pCR|RECIST|non-respond|pathologic", re.I)

raw = json.load(open("results/hunt_geo2026_sc/geo_esearch_2026_broad.json"))
rows = []
for d in raw["docs"].values():
    acc = d.get("accession", "")
    if not acc.startswith("GSE"):
        continue
    title = d.get("title", "")
    summary = d.get("summary", "")
    text = f"{title} {summary}"
    io = bool(IO.search(text))
    rows.append({
        "accession": acc,
        "already_listed": acc in LISTED,
        "io_kw": io,
        "resp_kw": bool(RESP.search(text)),
        "n_samples": int(d.get("n_samples", 0) or 0),
        "pdat": d.get("pdat", ""),
        "gdstype": d.get("gdstype", ""),
        "suppfile": d.get("suppfile", ""),
        "title": title,
        "summary": summary[:500],
    })
rows.sort(key=lambda r: (not r["io_kw"], r["already_listed"], -r["n_samples"], r["accession"]))
with open("results/hunt_geo2026_sc/broad_2026_lung_scrna.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
new_io = [r for r in rows if r["io_kw"] and not r["already_listed"]]
print(f"broad 2026 GSE: {len(rows)}  IO: {sum(r['io_kw'] for r in rows)}  NEW-IO: {len(new_io)}")
print("\n=== NEW 2026 lung-scRNA with IO keywords ===")
for r in new_io:
    print(f"{r['accession']:12s} n={r['n_samples']:4d} resp={int(r['resp_kw'])} {r['pdat']}  {r['title'][:100]}")
