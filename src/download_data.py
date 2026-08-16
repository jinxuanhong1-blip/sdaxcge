#!/usr/bin/env python3
"""Download the three public lung RNA-seq datasets used in this hunt.

All datasets are per-sample gene expression matrices for human lung:

  1. GTEx v8   - normal lung, gene TPM (GENCODE v26)
  2. TCGA-LUAD - lung adenocarcinoma, STAR TPM  (UCSC Xena / GDC hub, GENCODE v36)
  3. TCGA-LUSC - lung squamous cell carcinoma, STAR TPM (UCSC Xena / GDC hub)

Nothing here is committed to git (see .gitignore); everything is re-downloadable
from the stable public URLs below. A provenance record (URL, size, sha256,
timestamp) is written to data/provenance.json so results are auditable.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import urllib.request
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(HERE), "data")

DATASETS = {
    "gtex_lung": {
        "url": (
            "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
            "tpms-by-tissue/gene_tpm_2017-06-05_v8_lung.gct.gz"
        ),
        "filename": "gtex_v8_lung_gene_tpm.gct.gz",
        "description": "GTEx v8 normal lung, gene-level TPM per sample (GENCODE v26).",
    },
    "tcga_luad": {
        "url": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz",
        "filename": "tcga_luad_star_tpm.tsv.gz",
        "description": "TCGA-LUAD lung adenocarcinoma, STAR TPM (log2(tpm+1)); UCSC Xena GDC hub.",
    },
    "tcga_lusc": {
        "url": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.star_tpm.tsv.gz",
        "filename": "tcga_lusc_star_tpm.tsv.gz",
        "description": "TCGA-LUSC lung squamous cell carcinoma, STAR TPM (log2(tpm+1)); UCSC Xena GDC hub.",
    },
}


def sha256_of(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: str, retries: int = 4) -> None:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cldn-hunt/1.0"})
            with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as out:
                total = 0
                while True:
                    block = resp.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
                    total += len(block)
            print(f"    downloaded {total/1e6:.1f} MB")
            return
        except Exception as exc:  # noqa: BLE001 - report and retry with backoff
            wait = 4 * (2 ** attempt)
            print(f"    attempt {attempt+1} failed: {exc}; retrying in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"failed to download {url} after {retries} attempts")


def main() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    provenance = {"downloaded_at_utc": datetime.now(timezone.utc).isoformat(), "datasets": {}}
    for key, meta in DATASETS.items():
        dest = os.path.join(DATA_DIR, meta["filename"])
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"[{key}] already present: {dest}")
        else:
            print(f"[{key}] downloading {meta['url']}")
            download(meta["url"], dest)
        size = os.path.getsize(dest)
        print(f"[{key}] hashing {dest} ({size/1e6:.1f} MB) ...")
        provenance["datasets"][key] = {
            "url": meta["url"],
            "filename": meta["filename"],
            "description": meta["description"],
            "size_bytes": size,
            "sha256": sha256_of(dest),
        }
    prov_path = os.path.join(DATA_DIR, "provenance.json")
    with open(prov_path, "w") as fh:
        json.dump(provenance, fh, indent=2)
    print(f"wrote provenance -> {prov_path}")


if __name__ == "__main__":
    main()
