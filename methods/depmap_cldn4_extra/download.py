#!/usr/bin/env python3
"""Download DepMap Public 24Q4 files needed for the *extra* CLDN4 cut.

This slice is additive to methods/depmap_cldn4_ifn (CLDN4 RNA vs IFN/MHC/CD274
on all-lung n=214). Here we add CRISPR Chronos for CLDN4 and keep only the
RNA panel required to rebuild Hallmark IFN-γ / MHC-I / CD274.

Source (Figshare+, not the bot-gated portal):
  DepMap, Broad (2024). DepMap 24Q4 Public.
  https://doi.org/10.25452/figshare.plus.27993248.v1
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
CRISPR_URL = "https://ndownloader.figshare.com/files/51064667"
GMT_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
    "h.all.v2024.1.Hs.symbols.gmt"
)
RELEASE = "DepMap Public 24Q4"
DOI = "10.25452/figshare.plus.27993248.v1"
UA = (
    "sdaxcge-depmap-cldn4-extra/1.0 "
    "(public DepMap extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"
)

CORE_GENES = ["CLDN4", "TACSTD2", "CD274", "HLA-A", "HLA-B", "HLA-C", "B2M"]
HALLMARK_IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
CRISPR_GENES = ["CLDN4"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def parse_gmt_ifng(path: Path) -> list[str]:
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if parts and parts[0] == HALLMARK_IFNG and len(parts) >= 3:
            return [g for g in parts[2:] if g]
    raise SystemExit(f"{HALLMARK_IFNG} not found in {path}")


def symbol_of(col: str) -> str:
    return col.split(" (")[0].strip()


def extract_columns(src: Path, wanted: set[str], out_path: Path, id_name: str) -> dict:
    with src.open("r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
    keep_idx = [0]
    found: dict[str, str] = {}
    for i, col in enumerate(header[1:], start=1):
        sym = symbol_of(col)
        if sym in wanted and sym not in found:
            keep_idx.append(i)
            found[sym] = col
    missing = sorted(wanted - set(found))
    print(
        f"{src.name}: kept {len(found)}/{len(wanted)} columns; missing={missing[:12]}",
        flush=True,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    with src.open("r", newline="") as fh, out_path.open("w", newline="") as out:
        reader = csv.reader(fh)
        writer = csv.writer(out)
        header = next(reader)
        writer.writerow([id_name] + [symbol_of(header[i]) for i in keep_idx[1:]])
        for row in reader:
            writer.writerow([row[i] for i in keep_idx])
            n_rows += 1
    return {
        "n_models": n_rows,
        "n_genes_kept": len(found),
        "n_genes_wanted": len(wanted),
        "missing_genes": missing,
        "columns": found,
    }


def ensure(url: str, dest: Path, min_bytes: int) -> None:
    if dest.exists() and dest.stat().st_size >= min_bytes:
        print(f"OK exists {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    download(url, dest)
    if dest.stat().st_size < min_bytes:
        raise SystemExit(f"{dest} too small: {dest.stat().st_size}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/depmap_cldn4_extra")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    model_path = outdir / "Model.csv"
    gmt_path = outdir / "h.all.v2024.1.Hs.symbols.gmt"
    expr_full = outdir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    expr_slim = outdir / "expression_ifn_panel.csv"
    crispr_full = outdir / "CRISPRGeneEffect.csv"
    crispr_slim = outdir / "crispr_cldn4.csv"

    ensure(MODEL_URL, model_path, 100_000)
    ensure(GMT_URL, gmt_path, 5_000)

    ifng_genes = parse_gmt_ifng(gmt_path)
    wanted_rna = set(CORE_GENES) | set(ifng_genes)
    (outdir / "hallmark_ifng_genes.txt").write_text("\n".join(ifng_genes) + "\n")

    if not expr_slim.exists() or expr_slim.stat().st_size < 100_000:
        ensure(EXPR_URL, expr_full, 400_000_000)
        expr_meta = extract_columns(expr_full, wanted_rna, expr_slim, "ModelID")
        try:
            expr_full.unlink()
            print(f"Removed full matrix {expr_full.name}", flush=True)
        except OSError:
            pass
    else:
        print(f"OK exists {expr_slim} ({expr_slim.stat().st_size} bytes)", flush=True)
        with expr_slim.open() as fh:
            slim_header = next(csv.reader(fh))
        found_syms = slim_header[1:]
        expr_meta = {
            "n_models": sum(1 for _ in expr_slim.open()) - 1,
            "n_genes_kept": len(found_syms),
            "n_genes_wanted": len(wanted_rna),
            "missing_genes": sorted(wanted_rna - set(found_syms)),
            "columns": {s: s for s in found_syms},
        }

    if not crispr_slim.exists() or crispr_slim.stat().st_size < 1_000:
        ensure(CRISPR_URL, crispr_full, 300_000_000)
        crispr_meta = extract_columns(crispr_full, set(CRISPR_GENES), crispr_slim, "ModelID")
        try:
            crispr_full.unlink()
            print(f"Removed full matrix {crispr_full.name}", flush=True)
        except OSError:
            pass
    else:
        print(f"OK exists {crispr_slim} ({crispr_slim.stat().st_size} bytes)", flush=True)
        with crispr_slim.open() as fh:
            slim_header = next(csv.reader(fh))
        found_syms = slim_header[1:]
        crispr_meta = {
            "n_models": sum(1 for _ in crispr_slim.open()) - 1,
            "n_genes_kept": len(found_syms),
            "n_genes_wanted": len(CRISPR_GENES),
            "missing_genes": sorted(set(CRISPR_GENES) - set(found_syms)),
            "columns": {s: s for s in found_syms},
        }

    manifest = {
        "release": RELEASE,
        "doi": DOI,
        "model_url": MODEL_URL,
        "expression_url": EXPR_URL,
        "crispr_url": CRISPR_URL,
        "gmt_url": GMT_URL,
        "model_sha256": sha256_file(model_path),
        "expression_slim_sha256": sha256_file(expr_slim),
        "crispr_slim_sha256": sha256_file(crispr_slim),
        "gmt_sha256": sha256_file(gmt_path),
        "n_hallmark_ifng_genes_gmt": len(ifng_genes),
        "expression_extract": expr_meta,
        "crispr_extract": crispr_meta,
        "crispr_metric": "CRISPRGeneEffect.csv Chronos; more negative = stronger dependency",
    }
    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("release", "n_hallmark_ifng_genes_gmt")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
