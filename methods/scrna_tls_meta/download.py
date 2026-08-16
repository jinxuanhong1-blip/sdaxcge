#!/usr/bin/env python3
"""Download public processed lung scRNA files for the TLS / B meta.

Self-contained. Public GEO / GitHub only. No EGA FASTQ.
Filters are not tuned to a user-claimed number.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "scrna_tls_meta"

FILES = {
    "GSE131907": [
        (
            "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        ),
        (
            "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        ),
        (
            "GSE131907_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/GSE131907_series_matrix.txt.gz",
        ),
    ],
    "GSE207422": [
        (
            "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        ),
        (
            "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        ),
    ],
    "GSE241934": [
        (
            "GSE241934_IIT_Matrix.mtx.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Matrix.mtx.gz",
        ),
        (
            "GSE241934_IIT_Meta.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_Meta.txt.gz",
        ),
        (
            "GSE241934_IIT_barcodes.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_barcodes.tsv.gz",
        ),
        (
            "GSE241934_IIT_features.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_IIT_features.tsv.gz",
        ),
        (
            "GSE241934_Real_Matrix.mtx.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Matrix.mtx.gz",
        ),
        (
            "GSE241934_Real_Meta.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_Real_Meta.txt.gz",
        ),
        (
            "GSE241934_RWC_barcodes.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_barcodes.tsv.gz",
        ),
        (
            "GSE241934_RWC_features.tsv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RWC_features.tsv.gz",
        ),
    ],
    "GSE148071": [
        (
            "GSE148071_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/GSE148071_RAW.tar",
        ),
    ],
    "GSE253013": [
        (
            "GSE253013_all_luad_garnett_temp.rds.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/suppl/GSE253013_all_luad_garnett_temp.rds.gz",
        ),
        (
            "GSE253013_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/matrix/GSE253013_series_matrix.txt.gz",
        ),
    ],
    "GSE154826": [
        (
            "GSE154826_sample_annots.csv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/GSE154826_sample_annots.csv.gz",
        ),
    ],
}

LEADER = {
    "annots_list.csv": "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/annots_list.csv",
    "table_s1_sample_table.csv": "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/table_s1_sample_table.csv",
    "cell_metadata.csv": "https://raw.githubusercontent.com/effiken/Leader_et_al/master/input_tables/cell_metadata.csv",
}

GSE154826_BATCH_INDEX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/"


def sha256_file(path: Path, max_bytes: int | None = 32_000_000) -> str:
    h = hashlib.sha256()
    n = 0
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
            n += len(chunk)
            if max_bytes is not None and n >= max_bytes:
                return h.hexdigest() + f":partial_{n}"
    return h.hexdigest()


def wget(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    cmd = [
        "wget",
        "-c",
        "--tries=8",
        "--waitretry=8",
        "--timeout=60",
        "--progress=dot:giga",
        "-O",
        str(tmp),
        url,
    ]
    subprocess.check_call(cmd)
    tmp.rename(dest)


def list_154826_batches() -> list[tuple[str, str]]:
    html = urlopen(GSE154826_BATCH_INDEX, timeout=60).read().decode("utf-8", "replace")
    names = []
    for token in html.replace("<", " ").replace(">", " ").split():
        if token.startswith("GSE154826_amp_batch_ID_") and token.endswith(".tar.gz"):
            names.append(token)
    names = sorted(set(names))
    return [
        (n, f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/{n}")
        for n in names
    ]


def record(provenance: list[dict], cohort: str, name: str, url: str, dest: Path) -> None:
    st = dest.stat()
    provenance.append(
        {
            "cohort": cohort,
            "filename": name,
            "url": url,
            "bytes": st.st_size,
            "sha256_head": sha256_file(dest),
            "exists": True,
        }
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=DATA)
    ap.add_argument(
        "--cohorts",
        default="GSE131907,GSE207422,GSE241934,GSE148071,GSE154826,GSE253013",
    )
    args = ap.parse_args()
    want = {c.strip() for c in args.cohorts.split(",") if c.strip()}
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    prov: list[dict] = []

    for cohort, items in FILES.items():
        if cohort not in want:
            continue
        d = out / cohort
        d.mkdir(parents=True, exist_ok=True)
        for name, url in items:
            dest = d / name
            print(f"GET {cohort} {name}", flush=True)
            t0 = time.time()
            wget(url, dest)
            print(f"  {dest.stat().st_size} B in {time.time()-t0:.1f}s", flush=True)
            record(prov, cohort, name, url, dest)

    if "GSE154826" in want:
        lead = out / "leader"
        lead.mkdir(parents=True, exist_ok=True)
        for name, url in LEADER.items():
            dest = lead / name
            print(f"GET leader {name}", flush=True)
            wget(url, dest)
            record(prov, "GSE154826", name, url, dest)
        d = out / "GSE154826"
        d.mkdir(parents=True, exist_ok=True)
        for name, url in list_154826_batches():
            dest = d / name
            print(f"GET GSE154826 {name}", flush=True)
            t0 = time.time()
            wget(url, dest)
            print(f"  {dest.stat().st_size} B in {time.time()-t0:.1f}s", flush=True)
            record(prov, "GSE154826", name, url, dest)

    (out / "provenance_downloads.json").write_text(json.dumps(prov, indent=2) + "\n")
    print("wrote", out / "provenance_downloads.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
