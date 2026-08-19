#!/usr/bin/env python3
"""Download public KP/KL anti-PD1 processed RNA-seq matrices from GEO."""

from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "gemm"
UA = "tismo-icb-recompute/1.0"

# Focused public GEMM / GEMM-derived ICB series (processed matrices only).
FILES = [
    {
        "acc": "GSE114601",
        "model": "KP GEMM lung nodules",
        "note": "Adeegbe et al. KrasTrp53; Vehicle vs anti-PD-1 (n=2/arm)",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114601/suppl/GSE114601_counts.normalized.csv.gz",
        "name": "GSE114601_counts.normalized.csv.gz",
    },
    {
        "acc": "GSE169194",
        "model": "KPM (Kras;Trp53;Msh2) total viable",
        "note": "A2V ± anti-PD-1; no PD-1 monotherapy",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169194/suppl/GSE169194_annotated_log2_norm_count.txt.gz",
        "name": "GSE169194_annotated_log2_norm_count.txt.gz",
    },
    {
        "acc": "GSE157880",
        "model": "HKP1 (Kras/p53) orthotopic lung",
        "note": "IgG vs PD-1 ± RT; libraries are pooled mice",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157880/suppl/GSE157880_Bulk048.txt.gz",
        "name": "GSE157880_Bulk048.txt.gz",
    },
    {
        "acc": "GSE182228",
        "model": "Lkb1-deficient LUAD s.c. (KL)",
        "note": "anti-PD-1 ± palbociclib vs vehicle; RAW tar of per-sample files",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182228/suppl/GSE182228_RAW.tar",
        "name": "GSE182228_RAW.tar",
    },
]


def fetch(url: str) -> bytes:
    req = Request(url, headers={"User-Agent": UA})
    with urlopen(req, timeout=300) as resp:
        return resp.read()


def extract_genes(path: Path, genes: list[str]) -> Path | None:
    """Write a small gene-subset table so we do not commit huge matrices."""
    import pandas as pd

    name = path.name
    if name.endswith(".tar"):
        out_dir = path.parent / path.stem
        out_dir.mkdir(exist_ok=True)
        with tarfile.open(path) as tar:
            tar.extractall(out_dir)
        return out_dir

    if name.endswith(".csv.gz"):
        df = pd.read_csv(path, index_col=0)
    elif name.endswith(".txt.gz"):
        df = pd.read_csv(path, sep="\t", index_col=0)
    else:
        return None

    # gene symbols may be index or a column
    idx = df.index.astype(str)
    keep = df[idx.isin(genes)].copy()
    if keep.empty:
        for col in df.columns[:3]:
            if df[col].astype(str).isin(genes).any():
                keep = df[df[col].astype(str).isin(genes)].copy()
                break
        if keep.empty:
            # try case-insensitive on index
            low = {g.lower(): g for g in genes}
            mask = [str(x).split(".")[0].lower() in low for x in idx]
            keep = df.loc[mask].copy()
    slim = path.with_name(path.name.replace(".gz", "") + ".focus.tsv")
    # if still compressed name
    slim = Path(str(slim).replace(".csv.focus.tsv", ".focus.tsv").replace(".txt.focus.tsv", ".focus.tsv"))
    if name.endswith(".csv.gz"):
        slim = path.with_name(name.replace(".csv.gz", ".focus.tsv"))
    elif name.endswith(".txt.gz"):
        slim = path.with_name(name.replace(".txt.gz", ".focus.tsv"))
    keep.to_csv(slim, sep="\t")
    print(f"  focus table {slim} shape={keep.shape}")
    return slim


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    genes = ["Tacstd2", "Cldn4", "Cldn3", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
    manifest = []
    for spec in FILES:
        dest = DATA / spec["name"]
        print(f"[{spec['acc']}] {spec['name']}")
        if not dest.exists() or dest.stat().st_size < 100:
            blob = fetch(spec["url"])
            dest.write_bytes(blob)
            print(f"  wrote {dest} ({len(blob):,} bytes)")
        else:
            print(f"  cached {dest} ({dest.stat().st_size:,} bytes)")
        slim = extract_genes(dest, genes)
        rec = {**spec, "bytes": dest.stat().st_size, "focus": str(slim) if slim else None}
        manifest.append(rec)
    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
