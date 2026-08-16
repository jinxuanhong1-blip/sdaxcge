#!/usr/bin/env python3
"""Download the public GSE166449 files used by this analysis."""

from __future__ import annotations

import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
FILES = {
    "GSE166449_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/matrix/"
        "GSE166449_series_matrix.txt.gz"
    ),
    "GSE166449_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/soft/"
        "GSE166449_family.soft.gz"
    ),
    "GSE166449_Raw_gene_TPM_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/"
        "GSE166449_Raw_gene_TPM_matrix.txt.gz"
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    retrieved_at = datetime.now(timezone.utc).isoformat()
    rows = []
    for name, url in FILES.items():
        dest = DATA_DIR / name
        if not dest.exists():
            urllib.request.urlretrieve(url, dest)
        rows.append(
            {
                "file": f"data/{name}",
                "url": url,
                "bytes": dest.stat().st_size,
                "sha256": sha256(dest),
                "retrieved_or_verified_utc": retrieved_at,
            }
        )
    (HERE / "file_manifest.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
