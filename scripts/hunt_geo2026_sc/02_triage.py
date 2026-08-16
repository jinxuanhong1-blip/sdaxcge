#!/usr/bin/env python3
"""Triage live GEO hits: drop already-listed series, keep GSE expression series,
score for 2026 + response labels + malignant/T-NK + usable supplementary files.
"""
import csv
import json
import re

RAW = "results/hunt_geo2026_sc/geo_esearch_raw.json"
LISTED = "results/hunt_geo2026_sc/already_listed.txt"
OUT = "results/hunt_geo2026_sc/candidates_new.tsv"

RESP = re.compile(
    r"respond|responder|non-respond|MPR|pCR|major pathologic|pathologic(al)? response|"
    r"RECIST|progressive disease|partial response|complete response|"
    r"resistan|sensitive|benefit|DCB|durable", re.I)
NEO = re.compile(r"neoadjuvant|pre[- ]?treatment|post[- ]?treatment|baseline", re.I)
MAL = re.compile(r"malignant|tumor cell|cancer cell|epithelial|tumor tissue", re.I)
TNK = re.compile(r"T cell|NK|lymphocyte|CD8|immune", re.I)
SUPP = re.compile(r"H5AD|H5|RDS|MTX|TAR|CSV|TXT|TSV|LOOM", re.I)
SCRNA = re.compile(r"single.?cell|scRNA|snRNA|10[xX]|CITE", re.I)


def main():
    listed = {l.strip() for l in open(LISTED) if l.strip()}
    raw = json.load(open(RAW))
    rows = []
    for uid, d in raw["docs"].items():
        acc = d.get("accession", "")
        if not acc.startswith("GSE"):
            continue
        gdstype = d.get("gdstype", "")
        title = d.get("title", "")
        summary = d.get("summary", "")
        text = f"{title} {summary}"
        pdat = d.get("pdat", "")
        year = int(pdat[:4]) if pdat[:4].isdigit() else 0
        n = int(d.get("n_samples", 0) or 0)
        supp = d.get("suppfile", "")
        already = acc in listed
        score = 0
        has_resp = bool(RESP.search(text))
        has_neo = bool(NEO.search(text))
        has_sc = bool(SCRNA.search(text) or "single cell" in gdstype.lower())
        score += 4 * (year >= 2026)
        score += 3 * has_resp + 2 * has_neo + 2 * has_sc
        score += 1 * bool(MAL.search(text)) + 1 * bool(TNK.search(text))
        score += 2 * bool(SUPP.search(supp))
        score += min(n, 40) / 20.0
        if already:
            score -= 20
        rows.append({
            "accession": acc,
            "already_listed": already,
            "year": year,
            "pdat": pdat,
            "score": round(score, 2),
            "has_response_kw": has_resp,
            "has_neoadjuvant_kw": has_neo,
            "has_scrna_kw": has_sc,
            "n_samples": n,
            "gdstype": gdstype,
            "suppfile": supp,
            "title": title,
            "summary": summary[:700],
            "ftplink": d.get("ftplink", ""),
        })
    rows.sort(key=lambda r: (-r["score"], r["accession"]))
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    new = [r for r in rows if not r["already_listed"]]
    y2026 = [r for r in new if r["year"] >= 2026]
    print(f"total GSE: {len(rows)}  already_listed: {sum(r['already_listed'] for r in rows)}  NEW: {len(new)}  NEW-2026: {len(y2026)}")
    print("\n=== NEW (not already listed) ===")
    for r in new:
        print(f"{r['accession']:12s} {r['pdat']} score={r['score']:5} resp={int(r['has_response_kw'])} "
              f"n={r['n_samples']:4d}  {r['title'][:95]}")


if __name__ == "__main__":
    main()
