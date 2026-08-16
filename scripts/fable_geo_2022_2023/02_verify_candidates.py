#!/usr/bin/env python3
"""
Verify the GEO search candidates.

For every candidate GSE from stage 01 we:
  1. Classify relevance from title + summary text:
        is_human, is_lung, is_ici, is_expression, mentions_outcome
  2. List the series supplementary files (names + sizes) from the NCBI HTTPS
     mirror of the GEO FTP tree, and compute the total processed-download size.
  3. Flag whether an *open, processed* supplementary download of < 2 GB exists.

Output: results/fable_geo_2022_2023/geo_verified.csv
        results/fable_geo_2022_2023/geo_suppl_files.csv (one row per file)
"""
import csv
import re
import sys
import time
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "results" / "fable_geo_2022_2023"
CAND = OUTDIR / "geo_search_candidates.csv"

FTP_HTTPS = "https://ftp.ncbi.nlm.nih.gov/geo/series"

LUNG_RE = re.compile(
    r"\b(lung|nsclc|non[- ]?small[- ]?cell|luad|lusc|sclc|small[- ]cell lung|"
    r"pulmonary|adenocarcinoma of the lung|lung squamous|lung adenocarcinoma)\b",
    re.I)
ICI_RE = re.compile(
    r"\b(immune checkpoint|checkpoint inhibitor|checkpoint blockade|"
    r"immunotherap|anti[- ]?pd[- ]?l?1|pd[- ]?1|pd[- ]?l1|ctla[- ]?4|"
    r"nivolumab|pembrolizumab|atezolizumab|durvalumab|cemiplimab|avelumab|"
    r"ipilimumab|tislelizumab|sintilimab|camrelizumab|toripalimab|ici\b)\b",
    re.I)
EXPR_RE = re.compile(r"(expression profiling|rna-?seq|microarray|transcriptom)", re.I)
OUTCOME_RE = re.compile(
    r"\b(response|responder|non[- ]?responder|resistan|sensitiv|"
    r"progression[- ]free|overall survival|\bpfs\b|\bos\b|\borr\b|"
    r"durable clinical benefit|clinical benefit|outcome|efficacy|"
    r"recist|prognos)\b", re.I)


class LinkSizeParser(HTMLParser):
    """Parse an Apache/nginx autoindex page for filenames and byte sizes."""

    def __init__(self):
        super().__init__()
        self.files = []  # (name, size_bytes_or_None)
        self._cur_href = None

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for k, v in attrs:
                if k == "href":
                    self._cur_href = v


def list_suppl(accession):
    """Return list of (filename, size_bytes) for a series suppl/ dir."""
    m = re.match(r"GSE(\d+)", accession)
    if not m:
        return []
    num = m.group(1)
    prefix = "GSE" + (num[:-3] + "nnn" if len(num) > 3 else "nnn")
    url = f"{FTP_HTTPS}/{prefix}/{accession}/suppl/"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            html = r.read().decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        return [("__ERROR__", str(e))]
    files = []
    # NCBI autoindex rows look like:
    # <a href="GSExxx_file.txt.gz">GSExxx_file.txt.gz</a>   2022-01-01 12:00  1.2M
    for mrow in re.finditer(
            r'href="([^"?/][^"]*)"[^<]*</a>\s*([0-9-]{10}\s+[0-9:]{5})?\s*'
            r'([0-9.]+[KMGT]?|-)?', html):
        name = mrow.group(1)
        size_txt = mrow.group(3)
        if name in ("../",) or name.endswith("/"):
            continue
        files.append((name, parse_size(size_txt)))
    # de-dup preserving order
    seen = set()
    out = []
    for n, s in files:
        if n not in seen:
            seen.add(n)
            out.append((n, s))
    return out


def parse_size(txt):
    if not txt or txt == "-":
        return None
    mult = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    m = re.match(r"([0-9.]+)([KMGT]?)", txt)
    if not m:
        return None
    val = float(m.group(1))
    return int(val * mult.get(m.group(2), 1))


PROC_EXT_RE = re.compile(
    r"(counts?|fpkm|tpm|rpkm|expression|matrix|normalized|processed|"
    r"\.csv|\.tsv|\.txt|\.xlsx?|\.rds|series_matrix)", re.I)
RAW_ONLY_RE = re.compile(r"(\.bam|\.fastq|\.fq|\.cel(\.gz)?$|\.idat)", re.I)


def main():
    df = pd.read_csv(CAND)
    rows = []
    file_rows = []
    for _, c in df.iterrows():
        acc = c["accession"]
        text = f"{c['title']} {c['summary']} {c['gdstype']}"
        taxon = str(c["taxon"])
        is_human = "Homo sapiens" in taxon
        is_lung = bool(LUNG_RE.search(text))
        is_ici = bool(ICI_RE.search(text))
        is_expr = bool(EXPR_RE.search(f"{text} {c['gdstype']}"))
        has_outcome = bool(OUTCOME_RE.search(text))

        files = list_suppl(acc)
        total = 0
        has_proc = False
        proc_size = 0
        file_names = []
        err = ""
        for name, size in files:
            if name == "__ERROR__":
                err = size
                continue
            file_names.append(name)
            if isinstance(size, int):
                total += size
            is_proc = bool(PROC_EXT_RE.search(name)) and not RAW_ONLY_RE.search(name)
            if is_proc:
                has_proc = True
                if isinstance(size, int):
                    proc_size += size
            file_rows.append({
                "accession": acc, "file": name,
                "size_bytes": size if isinstance(size, int) else "",
                "is_processed_candidate": is_proc,
            })

        rows.append({
            "accession": acc,
            "title": c["title"],
            "taxon": taxon,
            "n_samples": c["n_samples"],
            "pdat": c["pdat"],
            "gdstype": c["gdstype"],
            "is_human": is_human,
            "is_lung": is_lung,
            "is_ici": is_ici,
            "is_expression": is_expr,
            "mentions_outcome": has_outcome,
            "relevant": is_human and is_lung and is_ici and is_expr,
            "n_suppl_files": len(file_names),
            "suppl_total_bytes": total,
            "suppl_total_mb": round(total / 1024**2, 1) if total else 0,
            "has_processed_suppl": has_proc,
            "processed_bytes": proc_size,
            "under_2gb": (total < 2 * 1024**3) if total else True,
            "suppl_files": ";".join(file_names),
            "suppl_error": err,
            "summary": c["summary"],
        })
        print(f"{acc}: relevant={rows[-1]['relevant']} outcome={has_outcome} "
              f"files={len(file_names)} {rows[-1]['suppl_total_mb']}MB",
              file=sys.stderr)
        time.sleep(0.2)

    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "geo_verified.csv", index=False)
    pd.DataFrame(file_rows).to_csv(OUTDIR / "geo_suppl_files.csv", index=False)

    rel = out[out["relevant"]]
    print(f"\nRelevant (human+lung+ICI+expression): {len(rel)}/{len(out)}",
          file=sys.stderr)
    print(f"  ...also mentioning outcome: {int(rel['mentions_outcome'].sum())}",
          file=sys.stderr)
    print(f"  ...with processed suppl <2GB: "
          f"{int((rel['has_processed_suppl'] & rel['under_2gb']).sum())}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
