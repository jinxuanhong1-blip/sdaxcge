#!/usr/bin/env python3
"""Download public Visium LUAD/NSCLC matrices + coordinates only (no H&E).

Sources (open):
  - 10x CytAssist lung cancer demos (LUSC FFPE, 11 mm neuroendocrine)
  - GSE189487, GSE273378, GSE300676 (per-GSM processed files)
  - GSE307534 (per-GSM Space Ranger tars; extract matrix + positions)

No private KL. Images / BAM / FASTQ are not kept.
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tarfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

DATA = Path(os.environ.get("VISIUM_DATA", "/workspace/data"))
HERE = Path(__file__).resolve().parent
FTP_SAMP = "https://ftp.ncbi.nlm.nih.gov/geo/samples"
UA = ["-A", "Mozilla/5.0 visium-spatial-lag"]


def gsm_bucket(gsm: str) -> str:
    return gsm[:7] + "nnn"


def curl(url: str, dest: Path, retries: int = 5) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  exists {dest} ({dest.stat().st_size})", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    delay = 4
    last = None
    for attempt in range(1, retries + 1):
        cmd = [
            "curl",
            "-fL",
            "--retry",
            "5",
            "--retry-delay",
            "4",
            "--retry-all-errors",
            *UA,
            "-o",
            str(tmp),
            url,
        ]
        print(f"  GET {url} (try {attempt})", flush=True)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode == 0 and tmp.exists() and tmp.stat().st_size > 0:
            tmp.replace(dest)
            print(f"  wrote {dest} ({dest.stat().st_size})", flush=True)
            return
        last = proc.stderr[-400:] if proc.stderr else f"rc={proc.returncode}"
        if tmp.exists():
            tmp.unlink()
        time.sleep(delay)
        delay = min(delay * 2, 32)
    raise RuntimeError(f"failed {url}: {last}")


def extract_named(tar_path: Path, dest_dir: Path, keep_names: tuple[str, ...]) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            base = Path(m.name).name
            if base.startswith("._"):
                continue
            if base not in keep_names and not any(m.name.endswith("/" + n) for n in keep_names):
                continue
            dest = dest_dir / base
            if dest.exists() and dest.stat().st_size > 0:
                continue
            src = tf.extractfile(m)
            if src is None:
                continue
            dest.write_bytes(src.read())


GSE189 = [
    ("GSM5702473", "TD1", "IAC"),
    ("GSM5702474", "TD2", "IAC"),
    ("GSM5702475", "TD3", "MIA"),
    ("GSM5702476", "TD5", "AIS"),
    ("GSM5702477", "TD6", "MIA"),
    ("GSM5702478", "TD8", "AIS"),
]

GSE273 = [
    ("GSM8427428", "LM_SD_1216_1"),
    ("GSM8427429", "LM_SD_16"),
    ("GSM8427430", "LM_SD_11"),
    ("GSM8427431", "LM_SD_2"),
    ("GSM8427432", "LM_SD_3"),
    ("GSM8427433", "LM_SD_4"),
    ("GSM8427434", "LM_SD_5"),
    ("GSM8427435", "LM_SD_6"),
    ("GSM8427436", "LM_SD_7"),
    ("GSM8427437", "LM_SD_1216_8"),
    ("GSM8427438", "LM_SD_9"),
    ("GSM8427439", "LM_SD_10"),
    ("GSM8427440", "LM_SD_1216_12"),
    ("GSM8427441", "LM_SD_13"),
    ("GSM8427442", "LM_SD_1216_14"),
    ("GSM8427443", "LM_SD_15"),
]

GSE300 = [
    ("GSM9066288", "CRC_mPAP3_A"),
    ("GSM9066289", "CRC_mPAP3_B"),
    ("GSM9066290", "CRC_mPAP4_A"),
    ("GSM9066291", "CRC_mPAP4_B"),
    ("GSM9066292", "CRC_mPAP1_A"),
    ("GSM9066293", "CRC_mPAP1_B"),
    ("GSM9066294", "CRC_mPAP2_A"),
    ("GSM9066295", "CRC_mPAP2_B"),
]


def dl_10x() -> None:
    root = DATA / "10x_visium"
    specs = [
        (
            "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma",
            "2.0.0",
            "LUSC_FFPE",
        ),
        (
            "CytAssist_11mm_FFPE_Human_Lung_Cancer",
            "2.0.1",
            "NEC_11mm",
        ),
    ]
    for name, ver, sub in specs:
        h5 = root / f"{name}_filtered_feature_bc_matrix.h5"
        url_h5 = (
            f"https://cf.10xgenomics.com/samples/spatial-exp/{ver}/{name}/"
            f"{name}_filtered_feature_bc_matrix.h5"
        )
        curl(url_h5, h5)
        star = root / f"{name}_spatial.tar.gz"
        url_s = (
            f"https://cf.10xgenomics.com/samples/spatial-exp/{ver}/{name}/"
            f"{name}_spatial.tar.gz"
        )
        curl(url_s, star)
        sdir = root / sub / "spatial"
        extract_named(
            star,
            sdir,
            (
                "tissue_positions.csv",
                "tissue_positions_list.csv",
                "scalefactors_json.json",
            ),
        )


def dl_mtx_series() -> None:
    # GSE189487
    d = DATA / "geo" / "GSE189487"
    for gsm, sid, _stage in GSE189:
        for suf in (
            f"{gsm}_{sid}_barcodes.tsv.gz",
            f"{gsm}_{sid}_features.tsv.gz",
            f"{gsm}_{sid}_matrix.mtx.gz",
            f"{gsm}_{sid}_tissue_positions_list.csv.gz",
        ):
            url = f"{FTP_SAMP}/{gsm_bucket(gsm)}/{gsm}/suppl/{suf}"
            curl(url, d / suf)

    # GSE273378
    d = DATA / "geo" / "GSE273378"
    for gsm, sid in GSE273:
        for suf in (
            f"{gsm}_{sid}_barcodes.tsv.gz",
            f"{gsm}_{sid}_features.tsv.gz",
            f"{gsm}_{sid}_matrix.mtx.gz",
            f"{gsm}_{sid}_tissue_positions_list.csv.gz",
            f"{gsm}_{sid}_scalefactors_json.json.gz",
        ):
            url = f"{FTP_SAMP}/{gsm_bucket(gsm)}/{gsm}/suppl/{suf}"
            curl(url, d / suf)

    # GSE300676: filtered H5 + spatial tar (positions/scalefactors only kept)
    d = DATA / "geo" / "GSE300676"
    for gsm, sid in GSE300:
        h5n = f"{gsm}_{sid}_filtered_feature_bc_matrix.h5"
        curl(f"{FTP_SAMP}/{gsm_bucket(gsm)}/{gsm}/suppl/{h5n}", d / h5n)
        star = d / f"{gsm}_{sid}_spatial.tar.gz"
        curl(
            f"{FTP_SAMP}/{gsm_bucket(gsm)}/{gsm}/suppl/{gsm}_{sid}_spatial.tar.gz",
            star,
        )
        extract_named(
            star,
            d / f"{sid}_spatial",
            (
                "tissue_positions.csv",
                "tissue_positions_list.csv",
                "scalefactors_json.json",
            ),
        )


KEEP_307 = (
    "features.tsv.gz",
    "barcodes.tsv.gz",
    "matrix.mtx.gz",
    "tissue_positions.csv",
    "tissue_positions_list.csv",
    "scalefactors_json.json",
)


def section_ready_307(out_dir: Path) -> bool:
    mtx = out_dir / "filtered_feature_bc_matrix" / "matrix.mtx.gz"
    feat = out_dir / "filtered_feature_bc_matrix" / "features.tsv.gz"
    bc = out_dir / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"
    pos = out_dir / "spatial" / "tissue_positions.csv"
    pos_list = out_dir / "spatial" / "tissue_positions_list.csv"
    return mtx.is_file() and feat.is_file() and bc.is_file() and (pos.is_file() or pos_list.is_file())


def extract_307(tar_path: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:*") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            base = Path(member.name).name
            if base.startswith("._") or base not in KEEP_307:
                continue
            parts = Path(member.name).parts
            if "filtered_feature_bc_matrix" in parts:
                dest = out_dir / "filtered_feature_bc_matrix" / base
            elif "spatial" in parts or base.startswith("tissue_") or base == "scalefactors_json.json":
                dest = out_dir / "spatial" / base
            else:
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            src = tf.extractfile(member)
            if src is None:
                continue
            dest.write_bytes(src.read())


def dl_gse307534(class_filter: str, workers: int) -> None:
    man = HERE / "sample_manifest_gse307534.tsv"
    rows = list(csv.DictReader(man.open(), delimiter="\t"))
    if class_filter:
        keep = {x.strip() for x in class_filter.split(",") if x.strip()}
        rows = [r for r in rows if r["class"] in keep]
    root = DATA / "geo" / "GSE307534"
    tmp = Path("/tmp/gse307534_tars")
    tmp.mkdir(parents=True, exist_ok=True)

    def one(row: dict[str, str]) -> str:
        gsm, label = row["gsm"], row["label"]
        out_dir = root / f"{gsm}_{label}"
        if section_ready_307(out_dir):
            return f"SKIP {gsm}_{label}"
        tar_path = tmp / f"{gsm}_{label}.tar.gz"
        curl(row["url"], tar_path)
        try:
            extract_307(tar_path, out_dir)
        finally:
            if tar_path.exists():
                tar_path.unlink()
        if not section_ready_307(out_dir):
            raise RuntimeError(f"{gsm}_{label} missing matrix or positions")
        return f"OK {gsm}_{label}"

    errors = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(one, r): r for r in rows}
        for fut in as_completed(futs):
            row = futs[fut]
            try:
                print(fut.result(), flush=True)
            except Exception as e:
                errors.append(f"{row['gsm']}_{row['label']}: {e}")
                print(f"FAIL {row['gsm']}_{row['label']}: {e}", flush=True)
    if errors:
        print(f"{len(errors)} GSE307534 failures", flush=True)
        for e in errors:
            print(" ", e, flush=True)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--skip-10x", action="store_true")
    p.add_argument("--skip-geo", action="store_true")
    p.add_argument("--skip-307534", action="store_true")
    p.add_argument(
        "--gse307534-class",
        default="invasive,precursor",
        help="comma-separated class filter; empty = all",
    )
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    if not args.skip_10x:
        print("=== 10x demos ===", flush=True)
        dl_10x()
    if not args.skip_geo:
        print("=== GSE189487 / GSE273378 / GSE300676 ===", flush=True)
        dl_mtx_series()
    if not args.skip_307534:
        print("=== GSE307534 ===", flush=True)
        dl_gse307534(args.gse307534_class, args.workers)
    print("download done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
