#!/usr/bin/env python3
"""Download in-budget GSE207422 processed files and paper supplements."""

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

SIZE_BUDGET = 2 * 1024**3
GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
ESM_BASE = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1186%2Fs13073-023-01164-9/MediaObjects/"
)

FILES = [
    {
        "file": "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "url": GEO_SUPPL + "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz",
        "role": "processed_umi",
        "wanted": True,
    },
    {
        "file": "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "url": GEO_SUPPL + "GSE207422_NSCLC_scRNAseq_metadata.xlsx",
        "role": "sample_metadata",
        "wanted": True,
    },
    {
        "file": "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
        "url": GEO_SUPPL + "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx",
        "role": "bulk_metadata_unused",
        "wanted": False,
    },
    {
        "file": "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        "url": GEO_SUPPL + "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz",
        "role": "bulk_unused",
        "wanted": False,
    },
    {
        "file": "13073_2023_1164_MOESM1_ESM.xlsx",
        "url": ESM_BASE + "13073_2023_1164_MOESM1_ESM.xlsx",
        "role": "paper_table_s1_clinical",
        "wanted": True,
    },
]


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length else None
    except Exception:
        return None


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", type=Path, default=Path("data/GSE207422"))
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parents[3] / "results" / "w200" / "A3_GSE207422_malig",
    )
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for spec in FILES:
        size = head_size(spec["url"])
        within = True if size is None else size < SIZE_BUDGET
        if spec["wanted"] and within:
            action = "download"
        elif not within:
            action = "skip_over_2gb"
        else:
            action = "skip_unused"
        dest = args.workdir / spec["file"]
        if dest.exists() and dest.stat().st_size > 0 and action == "download":
            action = "already_present"
            size = size or dest.stat().st_size
        rows.append(
            {
                "file": spec["file"],
                "role": spec["role"],
                "size_bytes": size,
                "size_gb": None if size is None else round(size / 1024**3, 3),
                "within_2gb_budget": within,
                "action": action,
                "url": spec["url"],
            }
        )
        if action == "download":
            print(f"GET {spec['file']} ({size} bytes)")
            fetch(spec["url"], dest)
        else:
            print(f"{action.upper()} {spec['file']} ({size} bytes)")

    manifest = args.outdir / "file_manifest.tsv"
    with manifest.open("w") as fh:
        fh.write("file\trole\tsize_bytes\tsize_gb\twithin_2gb_budget\taction\turl\n")
        for row in rows:
            fh.write(
                f"{row['file']}\t{row['role']}\t{row['size_bytes']}\t{row['size_gb']}\t"
                f"{row['within_2gb_budget']}\t{row['action']}\t{row['url']}\n"
            )

    umi = next(r for r in rows if r["role"] == "processed_umi")
    umi_bytes = umi["size_bytes"] or (args.workdir / umi["file"]).stat().st_size
    verdict = {
        "dataset": "GSE207422",
        "task": "W200-A3: malignant-only TACSTD2 NMPR vs MPR and vs T/NK",
        "size_budget_gb": 2.0,
        "processed_umi_bytes": umi_bytes,
        "processed_umi_gb": round(umi_bytes / 1024**3, 3),
        "umi_within_budget": umi_bytes < SIZE_BUDGET,
        "decision": "ANALYZE" if umi_bytes < SIZE_BUDGET else "CATALOG_ONLY",
        "author_cell_labels_public": False,
        "skipped_raw": "HRA001033 (GSA-Human controlled; not processed)",
        "workdir": str(args.workdir),
    }
    with (args.outdir / "feasibility.json").open("w") as fh:
        json.dump(verdict, fh, indent=2)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
