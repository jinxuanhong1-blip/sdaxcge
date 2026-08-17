#!/usr/bin/env python3
"""Download processed GSE131907 + GSE205335 files (no EGA raw)."""

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]

FILES = [
    (
        "GSE131907",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    ),
    (
        "GSE131907",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz",
        "GSE131907_series_matrix.txt.gz",
    ),
    (
        "GSE131907",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    ),
    (
        "GSE205335",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "GSE205335_Lung_IO_CellIdentity.txt.gz",
    ),
    (
        "GSE205335",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz",
        "GSE205335_family.soft.gz",
    ),
    (
        "GSE205335",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "GSE205335_Lung_IO_UMI_matrix.rds.gz",
    ),
]


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"OK {dest.name} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse131907", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--gse205335", type=Path, default=Path("/tmp/gse205335"))
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    roots = {"GSE131907": args.gse131907, "GSE205335": args.gse205335}
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for cohort, url, name in FILES:
        dest = roots[cohort] / name
        fetch(url, dest)
        rows.append({"cohort": cohort, "file": name, "bytes": dest.stat().st_size, "path": str(dest)})
    (args.outdir / "file_manifest.tsv").write_text(
        "cohort\tfile\tbytes\tpath\n"
        + "".join(f"{r['cohort']}\t{r['file']}\t{r['bytes']}\t{r['path']}\n" for r in rows)
    )
    (args.outdir / "download.json").write_text(json.dumps(rows, indent=2))
    print(json.dumps({"n_files": len(rows), "ok": True}, indent=2), flush=True)


if __name__ == "__main__":
    main()
