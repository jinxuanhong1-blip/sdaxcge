#!/usr/bin/env python3
"""Triage leftover 2026-only GEO lung ICI series for TACSTD2/CLDN4.

Scope: PDAT year == 2026. Prior analyzed 2024/2025 series (GSE261345,
GSE261348, GSE233203) are out of scope. Writes:
  results/w200/GEO_2026/leftover_2026_catalog.tsv
  results/w200/GEO_2026/series_matrix/<GSE>.characteristics.txt (new ones)
  results/w200/GEO_2026/leftover_2026_supp.tsv
"""
import csv
import gzip
import io
import json
import re
import time
import urllib.request
from pathlib import Path

RES = Path(__file__).resolve().parents[2] / "results" / "w200" / "GEO_2026"
META = RES / "series_matrix"
META.mkdir(parents=True, exist_ok=True)

# Already closed in the 2024–2025 leftover writeup; still catalogued here
# because they have 2026 PDAT, but they are not new analyses.
PREVIOUSLY_CLOSED = {
    "GSE292421",  # TME phenotype, unmappable IDs
    "GSE329813",  # spatial, no per-ROI response in prior check
}

RESP_KEYS = (
    "respon", "recist", "pfs", "overall surviv", "os time", "os event",
    "dcb", "benefit", "pcr", " mpr", "pathologic", "progression",
    "sensitiv", "resistan", "outcome", "efficacy", "best response",
    "irrecist", "clinical benefit",
)


def fetch_bytes(url: str, timeout=90):
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-2026"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise RuntimeError(f"{url}: {last}")


def series_matrix(gse: str):
    stub = gse[:-3] + "nnn"
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
           f"{gse}_series_matrix.txt.gz")
    try:
        raw = fetch_bytes(url)
    except Exception as e:  # noqa: BLE001
        return None, str(e)
    return gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode(
        "utf-8", "replace"), url


def parse_chars(txt: str):
    keys = []
    char_lines = []
    for ln in txt.splitlines():
        if ln.startswith(("!Sample_characteristics", "!Sample_title",
                          "!Sample_source_name", "!Sample_geo_accession",
                          "!Sample_description", "!Series_summary",
                          "!Series_overall_design", "!Series_title",
                          "!Series_type", "!Series_platform_id")):
            char_lines.append(ln)
        if ln.startswith("!Sample_characteristics"):
            for p in ln.split("\t")[1:]:
                p = p.strip().strip('"')
                if ":" in p:
                    keys.append(p.split(":", 1)[0].strip().lower())
                    break
    block = "\n".join(char_lines).lower()
    hits = sorted({k.strip() for k in RESP_KEYS if k in block})
    return sorted(set(keys)), hits, char_lines


def ftp_dir(gse: str):
    stub = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/"
    try:
        html = fetch_bytes(url, timeout=60).decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        return url, [], str(e)
    files = []
    for m in re.finditer(
            r'<a href="([^"/][^"]*)">[^<]+</a>\s*'
            r'([0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]+)\s+([0-9.]+[KMGT]?|-)',
            html):
        files.append((m.group(1), m.group(3)))
    if not files:
        for m in re.finditer(r'<a href="([^"/?][^"]*)">', html):
            if m.group(1) not in ("Parent Directory",):
                files.append((m.group(1), "?"))
    return url, files, ""


def classify(row, keys, hits, files, summary):
    acc = row["accession"]
    text = (row["title"] + " " + row["summary"] + " " + " ".join(keys)).lower()
    file_names = " ".join(n for n, _ in files).lower()
    reasons = []
    tissue = "unknown"
    if any(x in text for x in (
            "pbmc", "peripheral blood", "plasma", "serum", "circulating",
            "blood gene", "cell-free", "cfDNA", "whole blood")):
        tissue = "blood/fluid"
        reasons.append("blood/fluid assay")
    if any(x in text for x in (
            "cell line", "in vitro", "murine", "mouse", "mice", "xenograft",
            "knock-down", "knockout", "shrna", "sirna", "crispr")):
        tissue = "mechanistic/cell-line"
        reasons.append("cell-line/mechanistic")
    if "methylation" in row["gdstype"].lower() or "non-coding" in row["gdstype"].lower():
        reasons.append(f"assay={row['gdstype']}")
    if "atac" in text or "chip-seq" in text or "cut&run" in text:
        reasons.append("chromatin assay, not gene expression")
    if hits:
        reasons.append("response-like metadata: " + ",".join(hits))
    else:
        reasons.append("no response-like series-matrix field")
    if any(x in file_names for x in ("mtx", "h5", "csv", "xlsx", "txt", "tsv")):
        reasons.append("has processed supp")
    n = row.get("n_samples") or ""
    try:
        n_i = int(float(n))
    except ValueError:
        n_i = 0
    if n_i and n_i < 4:
        reasons.append(f"n={n_i} too small for group test")
    return tissue, "; ".join(reasons)


def main():
    cands = list(csv.DictReader((RES / "geo_search_candidates.tsv").open(),
                                delimiter="\t"))
    y2026 = [r for r in cands if str(r["pdat"]).startswith("2026")]
    print(f"2026 leftover candidates: {len(y2026)}")

    catalog = []
    supp_rows = []
    for i, r in enumerate(y2026, 1):
        gse = r["accession"]
        print(f"[{i}/{len(y2026)}] {gse}")
        txt, mtx_err = series_matrix(gse)
        keys, hits, char_lines = [], [], []
        status = "OK"
        if txt is None:
            status = "NO_SINGLE_MATRIX"
            keys, hits = [], []
        else:
            keys, hits, char_lines = parse_chars(txt)
            (META / f"{gse}.characteristics.txt").write_text("\n".join(char_lines))
        url, files, ferr = ftp_dir(gse)
        for name, size in files:
            supp_rows.append((gse, name, size, url + name))
        if ferr:
            supp_rows.append((gse, f"ERROR {ferr}", "", url))
        tissue, why = classify(r, keys, hits, files, r["summary"])
        catalog.append({
            "accession": gse,
            "pdat": r["pdat"],
            "n_samples": r["n_samples"],
            "gdstype": r["gdstype"],
            "title": r["title"],
            "matrix_status": status,
            "char_keys": ";".join(keys),
            "resp_hits": ";".join(hits),
            "n_supp_files": str(len(files)),
            "supp_names": ";".join(n for n, _ in files[:12]),
            "tissue_guess": tissue,
            "why": why,
            "previously_closed": "yes" if gse in PREVIOUSLY_CLOSED else "no",
        })
        time.sleep(0.2)

    cols = list(catalog[0].keys())
    with (RES / "leftover_2026_catalog.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in catalog:
            fh.write("\t".join(str(row[c]).replace("\t", " ") for c in cols) + "\n")
    with (RES / "leftover_2026_supp.tsv").open("w") as fh:
        fh.write("gse\tfile\tsize\turl\n")
        for row in supp_rows:
            fh.write("\t".join(row) + "\n")
    print("wrote", RES / "leftover_2026_catalog.tsv")
    print("response-like 2026 series:")
    for row in catalog:
        if row["resp_hits"] or row["accession"] in (
                "GSE329813", "GSE292299", "GSE317309", "GSE337519",
                "GSE311200", "GSE299684", "GSE305086"):
            print(f"  {row['accession']} hits={row['resp_hits']} "
                  f"tissue={row['tissue_guess']} why={row['why'][:140]}")


if __name__ == "__main__":
    main()
