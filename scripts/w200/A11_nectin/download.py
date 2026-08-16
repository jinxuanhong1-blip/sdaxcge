#!/usr/bin/env python3
"""Download public lung matrices and extract TACSTD2 + nectin genes.

Primary (tumors):
  UCSC Xena GDC hub STAR TPM (log2(TPM+1), GENCODE v36)
    TCGA-LUAD.star_tpm.tsv.gz
    TCGA-LUSC.star_tpm.tsv.gz
    gencode.v36.annotation.gtf.gene.probemap
  GDC open ABSOLUTE purity table
    https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5

Sensitivity (cell lines, optional):
  DepMap Public 24Q4 Figshare+ expression + Model.csv
  Streamed; only requested gene columns are written.

Raw matrices stay in a cache directory and are not committed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

GDC_HUB = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PURITY_URL = "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"
MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
UA = "sdaxcge-w200-A11-nectin/1.0 (public TCGA/DepMap extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"

# Current HGNC symbols (GENCODE v36 / DepMap 24Q4).
# Entrez used only to disambiguate DepMap "SYMBOL (entrez)" headers.
GENES = {
    "TACSTD2": "4070",
    "NECTIN1": "5818",
    "NECTIN2": "5819",
    "NECTIN3": "25945",
    "NECTIN4": "81607",
    "PVR": "5817",
    "CLDN4": "1364",
}

TCGA_FILES = {
    "TCGA-LUAD.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUAD.star_tpm.tsv.gz",
    "TCGA-LUSC.star_tpm.tsv.gz": f"{GDC_HUB}/TCGA-LUSC.star_tpm.tsv.gz",
    "gencode.v36.probemap": f"{GDC_HUB}/gencode.v36.annotation.gtf.gene.probemap",
    "tcga_absolute_purity.txt": PURITY_URL,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, timeout: int = 600) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)", flush=True)
        return
    print(f"[get ] {url} -> {dest}", flush=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    print(f"[done] {dest.name}  size={dest.stat().st_size:,}  sha256={sha256_file(dest)[:16]}...", flush=True)


def extract_tcga_genes(tpm_path: Path, probemap_path: Path, out_path: Path) -> dict:
    """Write a compact samples x genes table (log2(TPM+1)) for requested symbols."""
    pm = {}
    with open(probemap_path, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        id_i = header.index("id")
        gene_i = header.index("gene")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            pm[parts[id_i]] = parts[gene_i]

    wanted = set(GENES)
    # first pass: find Ensembl rows that map to wanted symbols
    keep_idx: dict[str, int] = {}
    means: dict[str, float] = {}
    with gzip.open(tpm_path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        for i, row in enumerate(reader):
            ens = row[0]
            sym = pm.get(ens)
            if sym not in wanted:
                continue
            vals = [float(x) if x not in ("", "NA") else float("nan") for x in row[1:]]
            mean = sum(v for v in vals if v == v) / max(1, sum(v == v for v in vals))
            if sym not in keep_idx or mean > means[sym]:
                keep_idx[sym] = i
                means[sym] = mean

    # second pass: collect the winning rows
    found: dict[str, list[float]] = {}
    with gzip.open(tpm_path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for i, row in enumerate(reader):
            ens = row[0]
            sym = pm.get(ens)
            if sym in keep_idx and keep_idx[sym] == i:
                found[sym] = [float(x) if x not in ("", "NA") else float("nan") for x in row[1:]]

    missing = sorted(wanted - set(found))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["sample"] + sorted(found))
        for j, sample in enumerate(samples):
            w.writerow([sample] + [found[g][j] for g in sorted(found)])
    return {
        "n_samples": len(samples),
        "genes_found": sorted(found),
        "genes_missing": missing,
        "winning_row_mean": {k: float(v) for k, v in means.items()},
    }


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
        raise SystemExit(f"DepMap column not found for {symbol} ({entrez}). Header sample: {header[:8]}")
    raise SystemExit(f"Ambiguous DepMap columns for {symbol}: {candidates}")


def extract_depmap(expr_path: Path, out_path: Path) -> dict:
    with expr_path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = {sym: match_gene_column(header, sym, entrez) for sym, entrez in GENES.items()}
        idx = {sym: header.index(col) for sym, col in cols.items()}
        n = 0
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["ModelID"] + list(GENES))
            for row in reader:
                if not row:
                    continue
                w.writerow([row[0]] + [row[idx[s]] for s in GENES])
                n += 1
    return {"n_models": n, "columns": cols}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/a11_nectin_cache")
    p.add_argument("--out-dir", default="results/w200/A11_nectin")
    p.add_argument("--skip-depmap", action="store_true")
    args = p.parse_args()

    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    manifest: dict = {
        "task": "A11_nectin",
        "genes": GENES,
        "tcga": {},
        "depmap": None,
    }

    for name, url in TCGA_FILES.items():
        dest = cache / name
        download(url, dest)
        manifest["tcga"][name] = {
            "url": url,
            "bytes": dest.stat().st_size,
            "sha256": sha256_file(dest),
        }

    pm = cache / "gencode.v36.probemap"
    for cohort in ("LUAD", "LUSC"):
        src = cache / f"TCGA-{cohort}.star_tpm.tsv.gz"
        dest = out / f"tcga_{cohort.lower()}_genes.csv"
        meta = extract_tcga_genes(src, pm, dest)
        manifest["tcga"][f"extract_{cohort}"] = meta
        print(f"[extract] {cohort}: {meta}", flush=True)

    # copy purity next to results (small)
    purity_src = cache / "tcga_absolute_purity.txt"
    purity_dst = out / "tcga_absolute_purity.txt"
    if not purity_dst.exists() or purity_dst.stat().st_size != purity_src.stat().st_size:
        purity_dst.write_bytes(purity_src.read_bytes())

    if not args.skip_depmap:
        try:
            model_path = cache / "Model.csv"
            expr_path = cache / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
            download(MODEL_URL, model_path)
            download(EXPR_URL, expr_path, timeout=1800)
            extract_meta = extract_depmap(expr_path, out / "depmap24q4_nectin_all_models.csv")
            model_copy = out / "Model.csv"
            if not model_copy.exists() or model_copy.stat().st_size != model_path.stat().st_size:
                model_copy.write_bytes(model_path.read_bytes())
            manifest["depmap"] = {
                "release": "DepMap Public 24Q4",
                "citation": "DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. https://doi.org/10.25452/figshare.plus.27993248.v1",
                "extract": extract_meta,
                "Model.csv_sha256": sha256_file(model_path),
            }
        except Exception as exc:
            manifest["depmap"] = {"error": str(exc), "note": "DepMap extract failed; TCGA analysis still proceeds."}
            print(f"[warn] DepMap skipped: {exc}", flush=True)

    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in manifest.items() if k != "tcga"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
