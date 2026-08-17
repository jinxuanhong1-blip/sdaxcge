#!/usr/bin/env python3
"""Download in-budget GSE131907 processed files for LIANA/LR.

Uses the author cell annotation and the processed raw UMI matrix
(0.38 GB gzip). Skips the 2.86 GB log2TPM text and EGA raw FASTQ.
"""

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
MATRIX = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/matrix/"
SIZE_BUDGET = 2 * 1024**3

FILES = [
    ("suppl", "GSE131907_Lung_Cancer_cell_annotation.txt.gz", True),
    ("matrix", "GSE131907_series_matrix.txt.gz", True),
    ("suppl", "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", True),
    ("suppl", "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz", False),
]


def head_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as resp:
        return int(resp.headers["Content-Length"])


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"HAVE {dest.name} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET {url}", flush=True)
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"OK {dest.name} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
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
            fetch(url, args.workdir / name)
        else:
            print(f"SKIP {name} ({size} bytes; {action})", flush=True)

    manifest = args.outdir / "file_manifest.tsv"
    with manifest.open("w") as fh:
        fh.write("file\tsize_bytes\tsize_gb\twithin_2gb_budget\taction\turl\n")
        for row in rows:
            fh.write(
                f"{row['file']}\t{row['size_bytes']}\t{row['size_gb']}\t"
                f"{row['within_2gb_budget']}\t{row['action']}\t{row['url']}\n"
            )
    umi = next(r for r in rows if r["file"].endswith("raw_UMI_matrix.txt.gz"))
    verdict = {
        "dataset": "GSE131907",
        "size_budget_gb": 2.0,
        "processed_umi_gb": umi["size_gb"],
        "umi_within_budget": umi["within_2gb_budget"],
        "skipped_raw_ega": "EGAD00001005054",
        "decision": "ANALYZE" if umi["within_2gb_budget"] else "SKIP",
        "workdir": str(args.workdir),
    }
    (args.outdir / "feasibility.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2), flush=True)


if __name__ == "__main__":
    main()
