#!/usr/bin/env python3
"""A10: download public TCGA lung RNA and extract GRHL family + TROP2/claudins.

Claim A10 (as stated by the user):
    "GRHL1 is a public regulator/correlate of TACSTD2 (TROP2) and CLDN4 in lung."

This script pulls the two public UCSC Xena TCGA lung cohorts and writes a small
tidy expression table with only the genes we need. It does NOT filter or tune
anything toward a pre-specified answer.

Source (UCSC Xena, TCGA hub; RSEM gene-level, log2(norm_count+1)):
  TCGA LUAD (lung adenocarcinoma)
    https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz
    (redirects to)
    https://tcga-xena-hub.s3.dualstack.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz
  TCGA LUSC (lung squamous cell carcinoma)
    https://tcga-xena-hub.s3.dualstack.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap/HiSeqV2.gz

The matrix is gene x sample. TCGA sample-type is encoded in the barcode
(the 4th field, e.g. TCGA-05-4384-01 -> "01" primary tumor, "11" normal).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

RELEASE = "UCSC Xena TCGA hub, HiSeqV2 (RSEM gene-level, log2(norm_count+1))"
UA = (
    "sdaxcge-w200-A10-GRHL1/1.0 "
    "(reproducible public TCGA extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"
)

COHORTS = {
    "LUAD": "https://tcga-xena-hub.s3.dualstack.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
    "LUSC": "https://tcga-xena-hub.s3.dualstack.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap/HiSeqV2.gz",
}

# Primary target genes plus honest context (GRHL family + related epithelial TFs,
# and the claudins TROP2 is often discussed alongside).
TARGET_GENES = [
    "GRHL1",   # the claimed regulator
    "GRHL2",   # closest paralog; well-known claudin/epithelial regulator (context)
    "GRHL3",   # paralog (context)
    "TACSTD2", # TROP2
    "CLDN4",
    "CLDN3",   # context claudin
    "CLDN7",   # context claudin
    "ELF3",    # epithelial TF context
    "KLF5",    # epithelial TF context
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    print(f"GET {url} -> {dest}", flush=True)
    with urllib.request.urlopen(req, timeout=600) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def sample_type_from_barcode(barcode: str) -> str:
    """TCGA sample-type from the 4th barcode field.

    01 primary tumor, 02 recurrent, 03/04/... other; 1x = normal/control.
    We only need a tumor/normal/other split for honest cohort slicing.
    """
    parts = barcode.split("-")
    if len(parts) < 4:
        return "unknown"
    code = parts[3][:2]
    if code == "01":
        return "primary_tumor"
    if code == "02":
        return "recurrent_tumor"
    if code.startswith("1"):
        return "normal"
    return "other"


def extract_genes(gz_path: Path, cohort: str) -> tuple[list[str], dict[str, list[str]]]:
    """Return (sample_ids, {gene: [values...]}) for TARGET_GENES present in file."""
    wanted = set(TARGET_GENES)
    found: dict[str, list[str]] = {}
    samples: list[str] = []
    with gzip.open(gz_path, "rt", newline="") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in f:
            if not line:
                continue
            i = line.find("\t")
            gene = line[:i]
            if gene in wanted:
                found[gene] = line.rstrip("\n").split("\t")[1:]
                if len(found) == len(wanted):
                    # keep reading is unnecessary; break early
                    break
    missing = wanted - set(found)
    if missing:
        print(f"  [warn] {cohort}: genes not found in matrix: {sorted(missing)}", flush=True)
    return samples, found


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/xena_tcga_lung")
    p.add_argument("--out-dir", default="results/w200/A10_GRHL1")
    args = p.parse_args()

    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    tidy_rows: list[dict] = []
    files_meta: dict[str, dict] = {}

    for cohort, url in COHORTS.items():
        gz = cache / f"TCGA.{cohort}.HiSeqV2.gz"
        if not gz.exists() or gz.stat().st_size < 1_000_000:
            download(url, gz)
        samples, found = extract_genes(gz, cohort)
        files_meta[cohort] = {
            "url": url,
            "sha256": sha256_file(gz),
            "bytes": gz.stat().st_size,
            "n_samples": len(samples),
            "genes_found": sorted(found.keys()),
            "genes_missing": sorted(set(TARGET_GENES) - set(found)),
        }
        for j, s in enumerate(samples):
            row = {
                "sample": s,
                "cohort": cohort,
                "sample_type": sample_type_from_barcode(s),
            }
            for g in TARGET_GENES:
                if g in found:
                    v = found[g][j]
                    row[g] = float(v) if v not in ("", "NA", "NaN") else ""
                else:
                    row[g] = ""
            tidy_rows.append(row)

    # Write tidy CSV (small: ~1100 samples x ~12 cols).
    import csv

    cols = ["sample", "cohort", "sample_type"] + TARGET_GENES
    tidy_path = out / "expression_tcga_lung.csv"
    with tidy_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in tidy_rows:
            w.writerow(r)

    manifest = {
        "task": "A10_GRHL1",
        "release": RELEASE,
        "citation": (
            "Goldman MJ et al. (2020) Visualizing and interpreting cancer genomics "
            "data via the Xena platform. Nat Biotechnol. TCGA LUAD/LUSC HiSeqV2."
        ),
        "expression_scale": "log2(norm_count+1) RSEM gene-level",
        "target_genes": TARGET_GENES,
        "files": files_meta,
        "tidy_output": str(tidy_path),
        "n_rows": len(tidy_rows),
        "note": "Raw matrices cached locally only; repo stores the small gene extract.",
    }
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
