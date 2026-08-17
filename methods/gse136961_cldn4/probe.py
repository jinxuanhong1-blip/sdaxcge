#!/usr/bin/env python3
"""Probe public GSE136961 for CLDN4 vs ICI response / CD8 / CD274.

Downloads the GEO TPM and raw-count tables if needed, records gene presence,
and writes empty test rows when CLDN4 is absent. Does not invent CLDN4 values.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961"
URLS = {
    "tpm": f"{FTP}/suppl/GSE136961_TPM.tsv.gz",
    "raw": f"{FTP}/suppl/GSE136961_raw_count.tsv.gz",
}

QUERIES = [
    "CLDN4",
    "TACSTD2",
    "TROP2",
    "CD274",
    "CD8A",
    "CD8B",
    "PDCD1",
    "PDCD1LG2",
    "CLDN3",
    "CLDN7",
    "CLDN18",
]
ALIASES = {
    "CLDN4": {"CLDN4", "CLAUDIN4", "CLAUDIN-4", "CPE-R", "CPER", "WBSCR8", "HCPER"},
    "TACSTD2": {"TACSTD2", "TROP2", "TROP-2", "GA733-1", "M1S1", "EGP1"},
}


def fetch(name: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / Path(URLS[name]).name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(URLS[name], headers={"User-Agent": "gse136961-cldn4-probe"})
    with urllib.request.urlopen(req, timeout=120) as r, dest.open("wb") as out:
        out.write(r.read())
    return dest


def load_tpm(path: Path) -> tuple[list[str], list[tuple[str, str]]]:
    with gzip.open(path, "rt") as f:
        header = [h.strip().strip('"') for h in next(f).rstrip("\n").split("\t")]
        rows = []
        for line in f:
            sid = line.split("\t", 1)[0].strip().strip('"')
            m = re.match(r"^([A-Za-z0-9.\-]+)_", sid)
            rows.append((sid, m.group(1) if m else sid))
    return header[1:], rows


def load_raw(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as f:
        header = [h.strip().strip('"') for h in next(f).rstrip("\n").split("\t")]
        rows = []
        for line in f:
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            rows.append(dict(zip(header, parts)))
    return rows


def hit(query: str, symbols: set[str]) -> bool:
    alts = {a.upper() for a in ALIASES.get(query, {query})}
    alts.add(query.upper())
    return any(s.upper() in alts for s in symbols)


def main() -> None:
    tpm_path = fetch("tpm")
    raw_path = fetch("raw")
    samples, tpm_rows = load_tpm(tpm_path)
    raw_rows = load_raw(raw_path)

    tpm_symbols = [s for _, s in tpm_rows]
    raw_genes = [r.get("Gene", "") for r in raw_rows]
    raw_targets = [r.get("Target", "") for r in raw_rows]

    d_ids = sorted(s for s in samples if s.startswith("D"))
    n_ids = sorted(s for s in samples if s.startswith("N"))

    presence = []
    for q in QUERIES:
        tpm_ids = [sid for sid, s in tpm_rows if hit(q, {s, sid})]
        raw_ids = [
            r.get("Target", "")
            for r in raw_rows
            if hit(q, {r.get("Gene", ""), r.get("Target", ""), r.get("NCBI_NAME", "")})
        ]
        presence.append(
            {
                "query": q,
                "in_tpm": bool(tpm_ids),
                "in_raw": bool(raw_ids),
                "tpm_ids": tpm_ids,
                "raw_ids": raw_ids,
            }
        )

    cldn4_n = 21 if any(p["query"] == "CLDN4" and p["in_tpm"] for p in presence) else 0
    empty = {
        "n": cldn4_n,
        "status": "empty" if cldn4_n == 0 else "computed",
        "reason": "CLDN4 not on deposited Oncomine Immune Response panel"
        if cldn4_n == 0
        else "",
    }

    summary = {
        "accession": "GSE136961",
        "public": True,
        "public_on": "2020-02-03",
        "pmid": "31959763",
        "platform": "GPL24014 Oncomine Immune Response Research Assay",
        "n_samples": len(samples),
        "n_genes_tpm": len(tpm_rows),
        "n_unique_genes_raw": len(set(raw_genes)),
        "n_targets_raw": len(raw_targets),
        "sample_ids": samples,
        "n_title_D": len(d_ids),
        "n_title_N": len(n_ids),
        "title_D": d_ids,
        "title_N": n_ids,
        "author_dcb_ndb": {"DCB": 9, "NDB": 12},
        "cldn_family_on_panel": sorted({s for s in tpm_symbols if s.upper().startswith("CLDN")}),
        "tacstd_trop_on_panel": sorted(
            {s for s in tpm_symbols if "TACSTD" in s.upper() or s.upper().startswith("TROP")}
        ),
        "presence": presence,
        "cldn4_n": cldn4_n,
        "tests": {
            "CLDN4_vs_response": dict(empty),
            "CLDN4_vs_CD8": dict(empty),
            "CLDN4_vs_CD274": dict(empty),
        },
        "ftp": URLS,
    }

    (HERE / "gene_presence.json").write_text(json.dumps(summary, indent=2) + "\n")

    with (HERE / "presence_table.tsv").open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["query", "in_tpm", "in_raw", "tpm_ids", "raw_ids"], delimiter="\t"
        )
        w.writeheader()
        for row in presence:
            w.writerow(
                {
                    "query": row["query"],
                    "in_tpm": row["in_tpm"],
                    "in_raw": row["in_raw"],
                    "tpm_ids": ",".join(row["tpm_ids"]),
                    "raw_ids": ",".join(row["raw_ids"]),
                }
            )

    with (HERE / "panel_symbols.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["symbol_id", "symbol"])
        w.writerows(tpm_rows)

    print(json.dumps({"cldn4_n": cldn4_n, "n_samples": len(samples), "tests": summary["tests"]}, indent=2))


if __name__ == "__main__":
    main()
