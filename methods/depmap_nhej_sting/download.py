#!/usr/bin/env python3
"""Download DepMap Public 24Q4 columns for the lung NHEJ / cGAS–STING cut.

Source (Figshare, not the bot-gated portal):
  DepMap, Broad (2024). DepMap 24Q4 Public.
  https://doi.org/10.6084/m9.figshare.27993248

RNA is CCLE/DepMap log2(TPM+1) from
OmicsExpressionProteinCodingGenesTPMLogp1.csv.
Dependency is CRISPRGeneEffect.csv Chronos (more negative = stronger dependency).
Full matrices are streamed, reduced to the gene panel, then deleted.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.request
from pathlib import Path

RELEASE = "DepMap Public 24Q4"
DOI = "10.6084/m9.figshare.27993248"
MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
EXPR_URL = "https://ndownloader.figshare.com/files/51065489"
CRISPR_URL = "https://ndownloader.figshare.com/files/51064667"
SCREEN_URL = "https://ndownloader.figshare.com/files/51065804"
SCREEN_MAP_URL = "https://ndownloader.figshare.com/files/51065159"
UA = (
    "sdaxcge-depmap-nhej-sting/1.0 "
    "(public DepMap extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"
)

# HGNC symbols in 24Q4, with legacy aliases if a column uses the old name.
GENE_ALIASES = {
    "CLDN4": ["CLDN4"],
    "PRKDC": ["PRKDC"],
    "XRCC4": ["XRCC4"],
    "LIG4": ["LIG4"],
    "TP53BP1": ["TP53BP1"],
    "STING1": ["STING1", "TMEM173"],
    "CGAS": ["CGAS", "MB21D1"],
    "STAT1": ["STAT1"],
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, attempts: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, attempts + 1):
        tmp = dest.with_suffix(dest.suffix + ".part")
        try:
            print(f"GET {url} -> {dest} (attempt {attempt})", flush=True)
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            tmp.replace(dest)
            print(f"  wrote {dest.stat().st_size} bytes", flush=True)
            return
        except Exception as exc:
            if tmp.exists():
                tmp.unlink()
            if attempt == attempts:
                raise
            wait = 4 * 2 ** (attempt - 1)
            print(f"  failed ({exc}); retry in {wait}s", flush=True)
            time.sleep(wait)


def symbol_of(col: str) -> str:
    return col.split(" (")[0].strip()


def extract_columns(src: Path, wanted: set[str], out_path: Path) -> dict:
    alias_to_canon = {}
    for canon, aliases in GENE_ALIASES.items():
        for alias in aliases:
            alias_to_canon[alias] = canon

    with src.open("r", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
    found: dict[str, tuple[int, str]] = {}
    for i, col in enumerate(header[1:], start=1):
        canon = alias_to_canon.get(symbol_of(col))
        if canon in wanted and canon not in found:
            found[canon] = (i, col)
    missing = sorted(wanted - set(found))
    order = [g for g in GENE_ALIASES if g in found]
    print(f"{src.name}: kept {len(found)}/{len(wanted)}; missing={missing}", flush=True)

    n_rows = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with src.open("r", newline="") as fh, out_path.open("w", newline="") as out:
        reader = csv.reader(fh)
        writer = csv.writer(out)
        next(reader)
        writer.writerow(["ModelID"] + order)
        idxs = [found[g][0] for g in order]
        for row in reader:
            writer.writerow([row[0]] + [row[i] if i < len(row) else "" for i in idxs])
            n_rows += 1
    return {
        "n_models": n_rows,
        "n_genes_kept": len(found),
        "n_genes_wanted": len(wanted),
        "missing_genes": missing,
        "source_columns": {g: found[g][1] for g in order},
        "symbols_written": order,
    }


def stream_screen_prkdc(url: str, dest: Path) -> dict:
    """PRKDC is absent from integrated CRISPRGeneEffect in 24Q4.

    It remains in ScreenGeneEffect for Humagne-CD screens only. Keep that
    one column so the gap is documented rather than silently dropped.
    """
    import io

    print(f"STREAM {url} -> {dest} (PRKDC column only)", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    n_rows = 0
    n_finite = 0
    with urllib.request.urlopen(req, timeout=600) as resp:
        text = io.TextIOWrapper(resp, encoding="utf-8", newline="")
        reader = csv.reader(text)
        header = next(reader)
        idx = next(i for i, col in enumerate(header) if symbol_of(col) == "PRKDC")
        source_col = header[idx]
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["ScreenID", "PRKDC"])
            for row in reader:
                val = row[idx] if idx < len(row) else ""
                writer.writerow([row[0], val])
                n_rows += 1
                if val not in ("", "NA", "NaN", "nan"):
                    n_finite += 1
    return {
        "source": "ScreenGeneEffect.csv",
        "source_column": source_col,
        "n_screens": n_rows,
        "n_finite": n_finite,
        "note": "Integrated CRISPRGeneEffect.csv has no PRKDC/XRCC7 column in 24Q4. Finite values are Humagne-CD screens.",
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
    ap.add_argument("--outdir", default="data/depmap_nhej_sting")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    wanted = set(GENE_ALIASES)
    model_path = outdir / "Model.csv"
    expr_full = outdir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    expr_slim = outdir / "expression_panel.csv"
    crispr_full = outdir / "CRISPRGeneEffect.csv"
    crispr_slim = outdir / "crispr_panel.csv"

    ensure(MODEL_URL, model_path, 100_000)

    if expr_slim.exists() and expr_slim.stat().st_size > 10_000:
        print(f"OK exists {expr_slim}", flush=True)
        with expr_slim.open() as fh:
            header = next(csv.reader(fh))
        expr_meta = {
            "n_models": sum(1 for _ in expr_slim.open()) - 1,
            "n_genes_kept": len(header) - 1,
            "n_genes_wanted": len(wanted),
            "missing_genes": sorted(wanted - set(header[1:])),
            "source_columns": {g: g for g in header[1:]},
            "symbols_written": header[1:],
            "reused_slim": True,
        }
    else:
        ensure(EXPR_URL, expr_full, 100_000_000)
        expr_meta = extract_columns(expr_full, wanted, expr_slim)
        expr_full.unlink(missing_ok=True)
        print("Removed full expression matrix", flush=True)

    if crispr_slim.exists() and crispr_slim.stat().st_size > 10_000:
        print(f"OK exists {crispr_slim}", flush=True)
        with crispr_slim.open() as fh:
            header = next(csv.reader(fh))
        crispr_meta = {
            "n_models": sum(1 for _ in crispr_slim.open()) - 1,
            "n_genes_kept": len(header) - 1,
            "n_genes_wanted": len(wanted),
            "missing_genes": sorted(wanted - set(header[1:])),
            "source_columns": {g: g for g in header[1:]},
            "symbols_written": header[1:],
            "reused_slim": True,
        }
    else:
        ensure(CRISPR_URL, crispr_full, 100_000_000)
        crispr_meta = extract_columns(crispr_full, wanted, crispr_slim)
        crispr_full.unlink(missing_ok=True)
        print("Removed full CRISPR matrix", flush=True)

    screen_map = outdir / "CRISPRScreenMap.csv"
    screen_prkdc = outdir / "screen_prkdc.csv"
    ensure(SCREEN_MAP_URL, screen_map, 1_000)
    if screen_prkdc.exists() and screen_prkdc.stat().st_size > 1_000:
        print(f"OK exists {screen_prkdc}", flush=True)
        screen_meta = {"reused_slim": True, "source": "ScreenGeneEffect.csv"}
    else:
        screen_meta = stream_screen_prkdc(SCREEN_URL, screen_prkdc)

    manifest = {
        "release": RELEASE,
        "doi": DOI,
        "model_url": MODEL_URL,
        "expression_url": EXPR_URL,
        "crispr_url": CRISPR_URL,
        "screen_gene_effect_url": SCREEN_URL,
        "screen_map_url": SCREEN_MAP_URL,
        "expression_matrix": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "expression_units": "log2(TPM+1)",
        "crispr_matrix": "CRISPRGeneEffect.csv",
        "crispr_metric": "Chronos gene effect; more negative = stronger dependency",
        "gene_aliases": GENE_ALIASES,
        "model_sha256": sha256_file(model_path),
        "expression_slim_sha256": sha256_file(expr_slim),
        "crispr_slim_sha256": sha256_file(crispr_slim),
        "expression_extract": expr_meta,
        "crispr_extract": crispr_meta,
        "screen_map_sha256": sha256_file(screen_map),
        "screen_prkdc_sha256": sha256_file(screen_prkdc),
        "screen_prkdc_extract": screen_meta,
    }
    (outdir / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "release": RELEASE,
                "expression_missing": expr_meta["missing_genes"],
                "crispr_missing": crispr_meta["missing_genes"],
                "expression_n": expr_meta["n_models"],
                "crispr_n": crispr_meta["n_models"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
