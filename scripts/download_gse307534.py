#!/usr/bin/env python3
"""Download GSE307534 Visium processed spots (matrix + coordinates only)."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

KEEP_SUFFIXES = (
    "/filtered_feature_bc_matrix/features.tsv.gz",
    "/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "/filtered_feature_bc_matrix/matrix.mtx.gz",
    "/spatial/tissue_positions.csv",
    "/spatial/tissue_positions_list.csv",
    "/spatial/scalefactors_json.json",
)


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def section_ready(out_dir: Path) -> bool:
    mtx = out_dir / "filtered_feature_bc_matrix" / "matrix.mtx.gz"
    feat = out_dir / "filtered_feature_bc_matrix" / "features.tsv.gz"
    bc = out_dir / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"
    pos = out_dir / "spatial" / "tissue_positions.csv"
    pos_list = out_dir / "spatial" / "tissue_positions_list.csv"
    return mtx.is_file() and feat.is_file() and bc.is_file() and (pos.is_file() or pos_list.is_file())


def keep_member(name: str) -> bool:
    base = name.split("/")[-1]
    if base.startswith("._"):
        return False
    return name.endswith(KEEP_SUFFIXES)


def extract_needed(tar_path: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile() or not keep_member(member.name):
                continue
            parts = Path(member.name).parts
            if "filtered_feature_bc_matrix" in parts:
                dest = out_dir / "filtered_feature_bc_matrix" / Path(member.name).name
            elif "spatial" in parts:
                dest = out_dir / "spatial" / Path(member.name).name
            else:
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            src = tf.extractfile(member)
            if src is None:
                continue
            dest.write_bytes(src.read())


def download_one(row: dict[str, str], data_root: Path, tmp_root: Path) -> str:
    gsm = row["gsm"]
    label = row["label"]
    out_dir = data_root / f"{gsm}_{label}"
    if section_ready(out_dir):
        return f"SKIP {gsm}_{label}"
    tmp_root.mkdir(parents=True, exist_ok=True)
    tar_path = tmp_root / f"{gsm}_{label}.tar.gz"
    cmd = [
        "curl",
        "-fL",
        "--retry",
        "5",
        "--retry-delay",
        "4",
        "--retry-all-errors",
        "-o",
        str(tar_path),
        row["url"],
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        if tar_path.exists():
            tar_path.unlink()
        raise RuntimeError(f"{gsm}_{label} download failed: {proc.stderr[-400:]}")
    try:
        extract_needed(tar_path, out_dir)
    finally:
        if tar_path.exists():
            tar_path.unlink()
    if not section_ready(out_dir):
        raise RuntimeError(f"{gsm}_{label} missing matrix or positions after extract")
    return f"OK {gsm}_{label}"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, default=Path("scripts/sample_manifest.tsv"))
    p.add_argument("--data-root", type=Path, default=Path("data/gse307534"))
    p.add_argument("--tmp-root", type=Path, default=Path("/tmp/gse307534_tars"))
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--class-filter", default="", help="comma-separated: invasive,precursor,normal")
    args = p.parse_args()

    rows = load_manifest(args.manifest)
    if args.class_filter:
        keep = {x.strip() for x in args.class_filter.split(",") if x.strip()}
        rows = [r for r in rows if r["class"] in keep]
    args.data_root.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(download_one, row, args.data_root, args.tmp_root): row for row in rows}
        for fut in as_completed(futs):
            row = futs[fut]
            try:
                print(fut.result(), flush=True)
            except Exception as e:
                msg = f"ERR {row['gsm']}_{row['label']}: {e}"
                print(msg, flush=True)
                errors.append(msg)
    if errors:
        print(f"{len(errors)} failed", file=sys.stderr)
        return 1
    print(f"done n={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
