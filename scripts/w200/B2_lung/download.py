#!/usr/bin/env python3
"""Download DepMap 24Q4 public files and extract TACSTD2 + CLDN4 expression.

Primary source (Figshare+, not the bot-gated portal):
  DepMap, Broad (2024). DepMap 24Q4 Public.
  https://doi.org/10.25452/figshare.plus.27993248.v1
  Figshare article 27993248

Files:
  Model.csv
    https://ndownloader.figshare.com/files/51065297
  OmicsExpressionProteinCodingGenesTPMLogp1.csv  (log2(TPM+1))
    https://ndownloader.figshare.com/files/51065489

The full expression matrix is ~507 MB. This script streams it and writes only
the two requested genes plus ModelID so the repo never stores the matrix.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path

FIGSHARE_ARTICLE = "https://doi.org/10.25452/figshare.plus.27993248.v1"
MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
RELEASE = "DepMap Public 24Q4"
UA = "sdaxcge-w200-B2-lung/1.0 (reproducible public DepMap extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"

GENES = {
    "TACSTD2": ("TACSTD2", "4070"),
    "CLDN4": ("CLDN4", "1364"),
}


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


def match_gene_column(header: list[str], symbol: str, entrez: str) -> str:
    exact = f"{symbol} ({entrez})"
    if exact in header:
        return exact
    candidates = [
        h
        for h in header
        if h == symbol
        or h.startswith(f"{symbol} (")
        or h.startswith(f"{symbol}(")
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise SystemExit(f"Gene column not found for {symbol} ({entrez}). Header sample: {header[:8]}")
    raise SystemExit(f"Ambiguous columns for {symbol}: {candidates}")


def extract_two_genes(expr_path: Path, out_path: Path) -> dict:
    print(f"Extracting TACSTD2 + CLDN4 from {expr_path}", flush=True)
    with expr_path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        id_col = header[0]
        tac = match_gene_column(header, *GENES["TACSTD2"])
        cld = match_gene_column(header, *GENES["CLDN4"])
        i_tac = header.index(tac)
        i_cld = header.index(cld)
        n = 0
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["ModelID", "TACSTD2", "CLDN4", "TACSTD2_column", "CLDN4_column"])
            for row in reader:
                if not row:
                    continue
                w.writerow([row[0], row[i_tac], row[i_cld], tac, cld])
                n += 1
    return {
        "n_models_in_expression": n,
        "id_column": id_col,
        "TACSTD2_column": tac,
        "CLDN4_column": cld,
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/depmap_24q4")
    p.add_argument("--out-dir", default="results/w200/B2_lung")
    args = p.parse_args()

    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    model_path = cache / "Model.csv"
    expr_path = cache / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    two_gene_path = out / "depmap24q4_tacstd2_cldn4_all_models.csv"

    if not model_path.exists() or model_path.stat().st_size < 1000:
        download(MODEL_URL, model_path)
    if not expr_path.exists() or expr_path.stat().st_size < 1_000_000:
        download(EXPR_URL, expr_path)

    extract_meta = extract_two_genes(expr_path, two_gene_path)

    # Keep a working copy of Model.csv next to results (small).
    model_copy = out / "Model.csv"
    if not model_copy.exists() or model_copy.stat().st_size != model_path.stat().st_size:
        model_copy.write_bytes(model_path.read_bytes())

    manifest = {
        "release": RELEASE,
        "citation": "DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. https://doi.org/10.25452/figshare.plus.27993248.v1",
        "figshare_article": FIGSHARE_ARTICLE,
        "files": {
            "Model.csv": {
                "url": MODEL_URL,
                "figshare_file_id": 51065297,
                "sha256": sha256_file(model_path),
                "bytes": model_path.stat().st_size,
            },
            "OmicsExpressionProteinCodingGenesTPMLogp1.csv": {
                "url": EXPR_URL,
                "figshare_file_id": 51065489,
                "sha256": sha256_file(expr_path),
                "bytes": expr_path.stat().st_size,
                "note": "full matrix cached locally only; repo stores the two-gene extract",
            },
        },
        "extract": extract_meta,
        "expression_scale": "log2(TPM+1)",
        "expression_file": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    }
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
