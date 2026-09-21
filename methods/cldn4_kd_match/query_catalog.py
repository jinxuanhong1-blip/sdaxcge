#!/usr/bin/env python3
"""Query the DepMap no-captcha download catalog for Expression Public and CRISPR.

Portal page downloads sit behind a browser check. The catalog at
https://depmap.org/portal/api/no-captcha/download/files lists each release
file and whether a bulk URL is populated. Expression Public is the
protein-coding log2(TPM+1) matrix. It is baseline RNA of the models, not
RNA measured after a CLDN4 knockout.
"""
from __future__ import annotations

import csv
import urllib.request
from pathlib import Path

CATALOG_URL = "https://depmap.org/portal/api/no-captcha/download/files"
KEEP = {
    "CRISPRGeneEffect.csv",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "Model.csv",
}
PROTEOMICS_RELEASES = {
    "Proteomics",
    "Harmonized Public Proteomics 24Q4",
    "Harmonized Public Proteomics 26Q1",
}


def main() -> int:
    out = Path("methods/cldn4_kd_match/tables/expression_public_catalog.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(CATALOG_URL, headers={"User-Agent": "cldn4-kd-match/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines()))
    kept = []
    for row in rows:
        release = row["release"]
        name = row["filename"]
        take = (release.startswith("DepMap Public") and name in KEEP) or (
            release in PROTEOMICS_RELEASES
        )
        if not take:
            continue
        kept.append(
            {
                "release": release,
                "release_date": row["release_date"],
                "filename": name,
                "url_present": "yes" if row["url"].strip() else "no",
            }
        )
    cols = ["release", "release_date", "filename", "url_present"]
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept)
    print(f"wrote {len(kept)} rows to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
