#!/usr/bin/env python3
"""Download DepMap Public 24Q4 RNA + Model metadata and extract genes needed
for CLDN4 vs Hallmark IFN-γ / MHC-I / CD274 / TACSTD2.

Source (Figshare+, not the bot-gated portal):
  DepMap, Broad (2024). DepMap 24Q4 Public.
  https://doi.org/10.25452/figshare.plus.27993248.v1

Hallmark GMT: MSigDB 2024.1 Hs.
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
GMT_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
    "h.all.v2024.1.Hs.symbols.gmt"
)
RELEASE = "DepMap Public 24Q4"
DOI = "10.25452/figshare.plus.27993248.v1"
UA = "sdaxcge-depmap-cldn4-ifn/1.0 (public DepMap extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"

# Always extract these plus Hallmark IFN-γ members.
CORE_GENES = ["CLDN4", "TACSTD2", "CD274", "HLA-A", "HLA-B", "HLA-C", "B2M"]
HALLMARK_IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"


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


def extract_genes(expr_path: Path, wanted: set[str], out_path: Path) -> dict:
    with expr_path.open("r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
    id_col = header[0]
    keep_idx = [0]
    found: dict[str, str] = {}
    for i, col in enumerate(header[1:], start=1):
        sym = symbol_of(col)
        if sym in wanted and sym not in found:
            keep_idx.append(i)
            found[sym] = col
    missing = sorted(wanted - set(found))
    print(f"Expression columns kept: {len(found)} / {len(wanted)}; missing={missing[:12]}", flush=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    with expr_path.open("r", newline="") as fh, out_path.open("w", newline="") as out:
        reader = csv.reader(fh)
        writer = csv.writer(out)
        header = next(reader)
        writer.writerow(["ModelID"] + [symbol_of(header[i]) for i in keep_idx[1:]])
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/depmap_cldn4_ifn")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    model_path = outdir / "Model.csv"
    gmt_path = outdir / "h.all.v2024.1.Hs.symbols.gmt"
    expr_full = outdir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    expr_slim = outdir / "expression_ifn_panel.csv"

    if not model_path.exists() or model_path.stat().st_size < 100_000:
        download(MODEL_URL, model_path)
    else:
        print(f"OK exists {model_path} ({model_path.stat().st_size} bytes)", flush=True)

    if not gmt_path.exists() or gmt_path.stat().st_size < 5_000:
        download(GMT_URL, gmt_path)
    else:
        print(f"OK exists {gmt_path} ({gmt_path.stat().st_size} bytes)", flush=True)

    ifng_genes = parse_gmt_ifng(gmt_path)
    wanted = set(CORE_GENES) | set(ifng_genes)
    (outdir / "hallmark_ifng_genes.txt").write_text("\n".join(ifng_genes) + "\n")

    if not expr_slim.exists() or expr_slim.stat().st_size < 100_000:
        if not expr_full.exists() or expr_full.stat().st_size < 400_000_000:
            download(EXPR_URL, expr_full)
        else:
            print(f"OK exists {expr_full} ({expr_full.stat().st_size} bytes)", flush=True)
        extract_meta = extract_genes(expr_full, wanted, expr_slim)
        # Drop the 507 MB matrix after the slim extract so the workspace stays small.
        try:
            expr_full.unlink()
            print(f"Removed full matrix {expr_full.name}", flush=True)
        except OSError:
            pass
    else:
        print(f"OK exists {expr_slim} ({expr_slim.stat().st_size} bytes)", flush=True)
        # Reconstruct found/missing from slim header + wanted list.
        with expr_slim.open() as fh:
            slim_header = next(csv.reader(fh))
        found_syms = [c for c in slim_header[1:]]
        extract_meta = {
            "n_models": sum(1 for _ in expr_slim.open()) - 1,
            "n_genes_kept": len(found_syms),
            "n_genes_wanted": len(wanted),
            "missing_genes": sorted(wanted - set(found_syms)),
            "columns": {s: s for s in found_syms},
        }

    manifest = {
        "release": RELEASE,
        "doi": DOI,
        "model_url": MODEL_URL,
        "expression_url": EXPR_URL,
        "gmt_url": GMT_URL,
        "model_sha256": sha256_file(model_path),
        "expression_slim_sha256": sha256_file(expr_slim),
        "gmt_sha256": sha256_file(gmt_path),
        "n_hallmark_ifng_genes_gmt": len(ifng_genes),
        "extract": extract_meta,
    }
    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("release", "n_hallmark_ifng_genes_gmt")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
