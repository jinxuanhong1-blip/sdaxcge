#!/usr/bin/env python3
"""Download in-budget GSE131907 processed files; skip the 2.9 GB log2TPM text and EGA raw."""

from __future__ import annotations

import argparse
import json
import os
import urllib.request
from pathlib import Path

SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
SIZE_BUDGET = 2 * 1024**3

FILES = [
    ("suppl", "GSE131907_Lung_Cancer_Feature_Summary.xlsx", True),
    ("suppl", "GSE131907_Lung_Cancer_cell_annotation.txt.gz", True),
    ("matrix", "GSE131907_series_matrix.txt.gz", True),
    ("suppl", "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", True),
    ("suppl", "GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz", False),
    ("suppl", "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds.gz", False),
    ("suppl", "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz", False),
]


def head_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        return int(resp.headers["Content-Length"])


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parents[3] / "results" / "w200" / "A3_GSE131907",
    )
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for kind, name, wanted in FILES:
        url = (SUPP if kind == "suppl" else MATRIX) + name
        size = head_size(url)
        within = size < SIZE_BUDGET
        action = "download" if wanted and within else ("skip_over_2gb" if not within else "skip_unused")
        rows.append(
            {
                "file": name,
                "size_bytes": size,
                "size_gb": round(size / 1024**3, 3),
                "within_2gb_budget": within,
                "action": action,
                "url": url,
            }
        )
        if action == "download":
            print(f"GET {name} ({size} bytes)")
            fetch(url, args.workdir / name)
        else:
            print(f"SKIP {name} ({size} bytes; {action})")

    manifest_path = args.outdir / "file_manifest.tsv"
    with manifest_path.open("w") as fh:
        fh.write("file\tsize_bytes\tsize_gb\twithin_2gb_budget\taction\turl\n")
        for row in rows:
            fh.write(
                f"{row['file']}\t{row['size_bytes']}\t{row['size_gb']}\t"
                f"{row['within_2gb_budget']}\t{row['action']}\t{row['url']}\n"
            )

    umi = next(r for r in rows if r["file"].endswith("raw_UMI_matrix.txt.gz"))
    huge = next(r for r in rows if r["file"].endswith("log2TPM_matrix.txt.gz"))
    verdict = {
        "dataset": "GSE131907",
        "task": "W200-A3 analog: epithelial TACSTD2 vs T/NK (no MPR labels in this atlas)",
        "size_budget_gb": 2.0,
        "processed_umi_bytes": umi["size_bytes"],
        "processed_umi_gb": umi["size_gb"],
        "umi_within_budget": umi["within_2gb_budget"],
        "skipped_log2tpm_txt_bytes": huge["size_bytes"],
        "skipped_log2tpm_txt_gb": huge["size_gb"],
        "skipped_raw_ega": "EGAD00001005054",
        "decision": "ANALYZE" if umi["within_2gb_budget"] else "SKIP",
        "workdir": str(args.workdir),
    }
    with (args.outdir / "feasibility.json").open("w") as fh:
        json.dump(verdict, fh, indent=2)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
