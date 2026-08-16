#!/usr/bin/env python3
"""Download TISMO metadata + official ICB gene CSVs used by this hunt.

Large Aliyun RDS matrices are optional (--rds). The committed analysis
uses the gene CSVs plus compact extracted tables under notes/.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tismo_client import (  # noqa: E402
    download_vivo_gene_csv,
    fetch_json,
    aliyun_download,
)

NOTES = ROOT / "notes" / "hunt_tismo_lung"
GENES = ["Tacstd2", "Cldn4", "Actb", "Cd8a", "Ifng", "Gzmb", "Cd274", "Epcam", "Krt8", "Cldn3", "Cldn7"]
META_ENDPOINTS = {
    "vivoMeta": "/tismo/metaData/vivoMeta",
    "vitroMeta": "/tismo/metaData/vitroMeta",
    "cellLineMeta": "/tismo/metaData/cellLineMeta",
    "getVivoTreatment": "/tismo/gene/getVivoTreatment",
    "getVivoCohort": "/tismo/gene/getVivoCohort",
    "getVitroTreatment": "/tismo/gene/getVitroTreatment",
    "getVitroCohort": "/tismo/gene/getVitroCohort",
}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rds", action="store_true", help="also download ~250 MB Aliyun RDS/annotation files")
    p.add_argument("--rds-dir", default="/tmp/tismo_data/aliyun")
    args = p.parse_args()

    raw = NOTES / "raw"
    raw.mkdir(parents=True, exist_ok=True)

    for name, path in META_ENDPOINTS.items():
        obj = fetch_json(path)
        (raw / f"{name}.json").write_text(json.dumps(obj, indent=2))
        print("meta", name, "status", obj.get("status"), "n", len(obj.get("data") or []))

    for gene in GENES:
        dest = raw / f"{gene}_vivo.csv"
        download_vivo_gene_csv(gene, dest)
        print("gene csv", gene, dest.stat().st_size)

    if args.rds:
        dest_dir = Path(args.rds_dir)
        for key in ["cellline", "vitro_ann", "vivo_ann", "immune", "vitro_expr", "vivo_expr"]:
            path = aliyun_download(key, dest_dir)
            print("aliyun", key, path, path.stat().st_size)


if __name__ == "__main__":
    main()
