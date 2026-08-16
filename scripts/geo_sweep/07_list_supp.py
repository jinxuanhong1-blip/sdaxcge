#!/usr/bin/env python3
"""
Step 7a: list series-level supplementary files for every catalog series that
lacks an embedded expression table, so we can find processed gene-level tables
(counts / TPM / FPKM / normalized) to extract TACSTD2 / CLDN4 from.

Writes notes/geo_sweep/supp_listing.json.
"""
import json
import re
import time
import urllib.request
from pathlib import Path

NOTES = Path("notes/geo_sweep")
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

ROW = re.compile(
    r'href="(?P<name>[^"?/][^"]*)".*?(?P<date>\d{4}-\d{2}-\d{2} \d{2}:\d{2})'
    r'\s+(?P<size>[\d.]+[KMG]?)',
)
UNIT = {"K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "": 1}


def parse_size(s):
    m = re.match(r"([\d.]+)([KMG]?)", s.strip())
    return int(float(m.group(1)) * UNIT[m.group(2)]) if m else None


def suppl_url(acc):
    stub = acc[:-3] + "nnn"
    return f"{FTP}/{stub}/{acc}/suppl/"


def http(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-sweep/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            code = getattr(e, "code", None)
            if code == 404:
                return None
            last = e
            time.sleep(2 * (i + 1))
    return None


# filenames that look like processed gene-level expression tables
PROC = re.compile(
    r"count|cpm|tpm|fpkm|rpkm|normaliz|expression|expr|_mat|matrix|deseq|"
    r"genes?[._]|abundance|quant|rsem|salmon|featurecount|processed|log2|voom",
    re.IGNORECASE,
)
BAD = re.compile(r"\.(h5|h5ad|mtx|loom|idat|fcs|bam|bai|bigwig|bw|pdf|png|"
                 r"tar|fastq|fq)(\.gz)?$", re.IGNORECASE)


def main():
    catalog_lines = (NOTES / "catalog.tsv").read_text().splitlines()
    header = catalog_lines[0].split("\t")
    idx = {c: i for i, c in enumerate(header)}
    rows = [ln.split("\t") for ln in catalog_lines[1:]]

    out = []
    for r in rows:
        acc = r[idx["accession"]]
        if r[idx["has_expression_matrix"]] == "Y":
            continue
        url = suppl_url(acc)
        html = http(url)
        files = []
        if html:
            for m in ROW.finditer(html):
                name = m.group("name")
                if name in ("Parent Directory",) or name.endswith("/"):
                    continue
                files.append({
                    "name": name,
                    "size_bytes": parse_size(m.group("size")),
                    "url": url + name,
                    "candidate": bool(PROC.search(name) and not BAD.search(name)),
                })
        cand = [f for f in files if f["candidate"]]
        out.append({
            "accession": acc,
            "is_ici": r[idx["is_ici"]],
            "is_kd_ko": r[idx["is_kd_ko"]],
            "suppl_url": url,
            "n_files": len(files),
            "files": files,
            "candidates": cand,
        })
        print(f"{acc:12s} ici={r[idx['is_ici']]} kd={r[idx['is_kd_ko']]} "
              f"files={len(files)} candidates={len(cand)}")
        time.sleep(0.15)

    (NOTES / "supp_listing.json").write_text(json.dumps(out, indent=2))
    n_with_cand = sum(1 for o in out if o["candidates"])
    print(f"\n{n_with_cand}/{len(out)} non-embedded series have a candidate "
          f"processed expression file.")


if __name__ == "__main__":
    main()
