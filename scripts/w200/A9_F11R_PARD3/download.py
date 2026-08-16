#!/usr/bin/env python3
"""A9: download the public lung matrices used for the F11R / PARD3 vs TACSTD2 test.

All sources are public and open access. Nothing here is controlled-access.

  TCGA-LUAD / TCGA-LUSC expression
    UCSC Xena, TCGA.{LUAD,LUSC}.sampleMap/HiSeqV2, gene x sample,
    unit log2(norm_count+1).
      https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz
      https://tcga.xenahubs.net/download/TCGA.LUSC.sampleMap/HiSeqV2.gz

  Tumour purity
    Aran, Sirota & Butte, Nat Commun 2015 (ncomms9971), Supplementary Data 1.
    Per-sample ESTIMATE / ABSOLUTE / LUMP / IHC / CPE purity.

The matrices are ~25 MB each and are NOT committed to git; rerun this script to
regenerate them. Downloads are idempotent and checksummed into
results/w200/A9_F11R_PARD3/data/SOURCES.json.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

DEPMAP_MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
DEPMAP_EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
DEPMAP_DOI = "https://doi.org/10.25452/figshare.plus.27993248.v1"
DEPMAP_GENES = {
    "TACSTD2": ("TACSTD2", "4070"),
    "F11R": ("F11R", "50848"),
    "PARD3": ("PARD3", "56288"),
    "CLDN1": ("CLDN1", "9076"),
    "CLDN4": ("CLDN4", "1364"),
    "CLDN7": ("CLDN7", "1366"),
    "EPCAM": ("EPCAM", "4072"),
    "KRT8": ("KRT8", "3856"),
    "PTPRC": ("PTPRC", "5788"),
    "CD8A": ("CD8A", "925"),
}

UA = (
    "sdaxcge-w200-A9-F11R-PARD3/1.0 (reproducible public lung extract; "
    "+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

FILES = {
    "TCGA.LUAD.HiSeqV2.gz": {
        "url": "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
        "source": "UCSC Xena TCGA.LUAD.sampleMap/HiSeqV2",
        "unit": "log2(norm_count+1)",
    },
    "TCGA.LUSC.HiSeqV2.gz": {
        "url": "https://tcga.xenahubs.net/download/TCGA.LUSC.sampleMap/HiSeqV2.gz",
        "source": "UCSC Xena TCGA.LUSC.sampleMap/HiSeqV2",
        "unit": "log2(norm_count+1)",
    },
    "Aran_CPE_purity.xlsx": {
        "url": (
            "https://media.springernature.com/original/springer-static/esm/"
            "art%3A10.1038%2Fncomms9971/MediaObjects/"
            "41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"
        ),
        "source": "Aran et al. Nat Commun 2015 ncomms9971 Supplementary Data 1",
        "unit": "purity fraction 0-1",
    },
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
    print(f"GET {url}\n -> {dest}", flush=True)
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
        if h == symbol or h.startswith(f"{symbol} (") or h.startswith(f"{symbol}(")
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise SystemExit(f"Gene column not found for {symbol} ({entrez})")
    raise SystemExit(f"Ambiguous columns for {symbol}: {candidates}")


def extract_depmap_panel(expr_path: Path, out_path: Path) -> dict:
    with expr_path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = {sym: match_gene_column(header, *pair) for sym, pair in DEPMAP_GENES.items()}
        idx = {sym: header.index(col) for sym, col in cols.items()}
        n = 0
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["ModelID", *DEPMAP_GENES.keys()])
            for row in reader:
                if not row:
                    continue
                w.writerow([row[0], *[row[idx[sym]] for sym in DEPMAP_GENES]])
                n += 1
    return {"n_models": n, "columns": cols}


def fetch_depmap(out: Path, cache: Path) -> dict:
    cache.mkdir(parents=True, exist_ok=True)
    model_path = cache / "Model.csv"
    expr_path = cache / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    if not model_path.exists() or model_path.stat().st_size < 1000:
        download(DEPMAP_MODEL_URL, model_path)
    if not expr_path.exists() or expr_path.stat().st_size < 1_000_000:
        download(DEPMAP_EXPR_URL, expr_path)
    extract = extract_depmap_panel(expr_path, out / "depmap24q4_panel_all_models.csv")
    (out / "Model.csv").write_bytes(model_path.read_bytes())
    return {
        "doi": DEPMAP_DOI,
        "model_sha256": sha256_file(model_path),
        "expr_bytes": expr_path.stat().st_size,
        "extract": extract,
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="results/w200/A9_F11R_PARD3/data")
    p.add_argument("--force", action="store_true", help="re-download even if present")
    p.add_argument("--with-depmap", action="store_true", help="also extract DepMap 24Q4 lung panel")
    p.add_argument("--depmap-cache", default="/tmp/depmap_24q4")
    args = p.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    manifest = {
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "files": {},
    }
    for name, meta in FILES.items():
        dest = out / name
        if dest.exists() and not args.force:
            print(f"exists, skipping: {dest}")
        else:
            download(meta["url"], dest)
        manifest["files"][name] = {
            **meta,
            "bytes": dest.stat().st_size,
            "sha256": sha256_file(dest),
        }

    if args.with_depmap:
        manifest["depmap"] = fetch_depmap(out, Path(args.depmap_cache))

    (out / "SOURCES.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
