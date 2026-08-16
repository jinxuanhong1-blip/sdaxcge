#!/usr/bin/env python3
"""Download public GEO files used by the TROP2-ADC analog slice.

Public only. Does not download private SKB264 data. Raw files stay under
results/c_public_trop2_adc_analogs/raw/ (git-ignored).
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "results" / "c_public_trop2_adc_analogs" / "raw"

FILES = [
    (
        "GSE312098_gene_fpkm.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312098/suppl/GSE312098_gene_fpkm.txt.gz",
    ),
    (
        "GSE311016_gene_fpkm.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311016/suppl/GSE311016_gene_fpkm.txt.gz",
    ),
    (
        "GSE304294_gene_fpkm.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304294/suppl/GSE304294_gene_fpkm.txt.gz",
    ),
    (
        "GSE334497_normalized_counts.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz",
    ),
    (
        "GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz",
    ),
    (
        "GSE235812_RAW.tar",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235812/suppl/GSE235812_RAW.tar",
        "GSE235812/GSE235812_RAW.tar",
    ),
]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"exists {dest.name} ({dest.stat().st_size} bytes)")
        return
    print(f"GET {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"  -> {dest} {dest.stat().st_size}")


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    prov = []
    for item in FILES:
        name, url = item[0], item[1]
        rel = item[2] if len(item) > 2 else name
        dest = RAW / rel
        fetch(url, dest)
        prov.append(
            {
                "file": rel,
                "url": url,
                "bytes": dest.stat().st_size,
                "md5": md5(dest),
            }
        )
    out = ROOT / "results" / "c_public_trop2_adc_analogs" / "provenance.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"downloaded": prov}, indent=2) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
