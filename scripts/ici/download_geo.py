#!/usr/bin/env python3
"""Download only the small, open GEO processed files used in this analysis."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import urllib.request


MAX_BYTES = 2_000_000_000
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"

FILES = {
    "GSE126044": [
        f"{BASE}/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    ],
    "GSE135222": [
        f"{BASE}/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    ],
    "GSE166449": [
        f"{BASE}/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
    ],
    "GSE207422": [
        f"{BASE}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        f"{BASE}/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
    ],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def remote_size(url: str) -> int:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        return int(response.headers["Content-Length"])


def retrieve(url: str, destination: Path) -> None:
    size = remote_size(url)
    if size > MAX_BYTES:
        raise RuntimeError(f"Refusing file over 2 GB ({size} bytes): {url}")
    if destination.exists() and destination.stat().st_size == size:
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    urllib.request.urlretrieve(url, temporary)
    if temporary.stat().st_size != size:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Incomplete download: {url}")
    temporary.replace(destination)


def family_soft_url(accession: str) -> str:
    stem = accession[:-3] + "nnn"
    return f"{BASE}/{stem}/{accession}/soft/{accession}_family.soft.gz"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/ici/data"),
        help="Destination directory (default: results/ici/data)",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    manifest = []
    for accession, urls in FILES.items():
        for url in [family_soft_url(accession), *urls]:
            destination = args.output / url.rsplit("/", 1)[-1]
            retrieve(url, destination)
            manifest.append(
                {
                    "accession": accession,
                    "file": destination.name,
                    "bytes": destination.stat().st_size,
                    "sha256": sha256(destination),
                    "url": url,
                }
            )

    columns = ["accession", "file", "bytes", "sha256", "url"]
    with (args.output / "manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for record in manifest:
            handle.write("\t".join(str(record[column]) for column in columns) + "\n")


if __name__ == "__main__":
    main()
