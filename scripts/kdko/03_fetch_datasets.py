#!/usr/bin/env python3
"""Download open processed matrices + full SOFT metadata for the verified KD/KO series.

Only open-access processed files are fetched. No SRA/FASTQ, no raw MS, nothing >2 GB.
GSE15212_RAW.tar (463 MB of Affymetrix .CEL) is deliberately skipped in favour of the
7.5 MB GEO series matrix (MAS5/RMA-processed by the submitters).

Writes into data/kdko/<GSE>/ and records checksums in results/kdko/downloads.tsv
"""
from __future__ import annotations

import csv
import hashlib
import os
import sys
import time
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "data", "kdko")
OUT = os.path.join(ROOT, "results", "kdko")
FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
MAX_BYTES = 2 * 1024 ** 3

# (GSE, [relative FTP paths]) -- every path was confirmed to return HTTP 200 by
# scripts/kdko/02 + a HEAD sweep before being written here.
TARGETS: list[tuple[str, list[str]]] = [
    ("GSE50927", [  # mouse lung Cldn4-/- vs WT, ventilator-induced lung injury
        "series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
        "series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtGenes.csv.gz",
        "series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz",
        "series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz",
    ]),
    ("GSE334497", [  # mouse 4T1 Trop2 KO vs WT tumours
        "series/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz",
    ]),
    ("GSE289287", [  # human T-47D Trop-2 KO (and DSG2 KO) cells + xenografts
        "series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz",
        "series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-DSG2KO_tumors_vs_WT.tsv.gz",
        "series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-DSG2KO_cells_vs_DSG2WT.tsv.gz",
    ]),
    ("GSE207704", [  # human T47D / MCF7 CLDN4-/- CRISPR KO
        "series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz",
    ]),
    ("GSE22493", [  # human SKOV-3 CLDN4 siRNA silencing (2-colour array)
        "series/GSE22nnn/GSE22493/matrix/GSE22493_series_matrix.txt.gz",
    ]),
    ("GSE15212", [  # human SW480 TACSTD2 siRNA knockdown screen (Affymetrix)
        "series/GSE15nnn/GSE15212/matrix/GSE15212_series_matrix.txt.gz",
    ]),
]

SOFT = {g: f"series/{g[:-3]}nnn/{g}/soft/{g}_family.soft.gz" for g, _ in TARGETS}
SOFT["GSE50927"] = "series/GSE50nnn/GSE50927/soft/GSE50927_family.soft.gz"
SOFT["GSE334497"] = "series/GSE334nnn/GSE334497/soft/GSE334497_family.soft.gz"
SOFT["GSE289287"] = "series/GSE289nnn/GSE289287/soft/GSE289287_family.soft.gz"
SOFT["GSE207704"] = "series/GSE207nnn/GSE207704/soft/GSE207704_family.soft.gz"
SOFT["GSE22493"] = "series/GSE22nnn/GSE22493/soft/GSE22493_family.soft.gz"
SOFT["GSE15212"] = "series/GSE15nnn/GSE15212/soft/GSE15212_family.soft.gz"


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD",
                                 headers={"User-Agent": "cldn4-trop2-kdko/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as fh:
            v = fh.headers.get("Content-Length")
            return int(v) if v else None
    except Exception:  # noqa: BLE001
        return None


def download(rel: str, dest: str) -> dict:
    url = f"{FTP}/{rel}"
    size = head_size(url)
    if size is not None and size > MAX_BYTES:
        print(f"SKIP (>2GB, {size}) {rel}")
        return {"path": rel, "status": "skipped_too_large", "bytes": size, "sha256": ""}
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"cached {os.path.basename(dest)}")
    else:
        for attempt in range(5):
            try:
                req = urllib.request.Request(
                    url, headers={"User-Agent": "cldn4-trop2-kdko/1.0"})
                with urllib.request.urlopen(req, timeout=300) as fh, \
                        open(dest + ".part", "wb") as out:
                    while chunk := fh.read(1 << 20):
                        out.write(chunk)
                os.replace(dest + ".part", dest)
                break
            except Exception as exc:  # noqa: BLE001
                sys.stderr.write(f"retry {attempt} {rel}: {exc}\n")
                if attempt == 4:
                    return {"path": rel, "status": f"failed: {exc}",
                            "bytes": 0, "sha256": ""}
                time.sleep(4 * 2 ** attempt)
    h = hashlib.sha256()
    with open(dest, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    nbytes = os.path.getsize(dest)
    print(f"  ok {os.path.basename(dest):58s} {nbytes/1e6:7.2f} MB")
    return {"path": rel, "status": "ok", "bytes": nbytes, "sha256": h.hexdigest()}


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    log = []
    for gse, files in TARGETS:
        d = os.path.join(DATA, gse)
        os.makedirs(d, exist_ok=True)
        print(f"=== {gse}")
        for rel in files + [SOFT[gse]]:
            rec = download(rel, os.path.join(d, os.path.basename(rel)))
            rec["gse"] = gse
            log.append(rec)
    with open(os.path.join(OUT, "downloads.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gse", "path", "status", "bytes", "sha256"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(log)
    bad = [r for r in log if r["status"] != "ok"]
    print(f"\n{len(log)-len(bad)}/{len(log)} downloaded ok")
    for r in bad:
        print("  !!", r["path"], r["status"])


if __name__ == "__main__":
    main()
