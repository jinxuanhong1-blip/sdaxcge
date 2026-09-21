#!/usr/bin/env python3
"""Summarize the HTAN portal lung scRNA-seq file snapshot.

The snapshot ``notes/htan_cldn4_sting/lung_scrna_files.tsv`` was queried from
the public HTAN portal catalog (ClickHouse database ``htan_2026_932``, table
``files``) on 2026-09-21:

    SELECT atlas_name, level, downloadSource, FileFormat, assayName,
           Filename, synapseId, CellTotal,
           arrayStringConcat(PrimaryDiagnosis, ' | ') AS diagnosis,
           arrayStringConcat(organType, ' | ') AS organ
    FROM files
    WHERE has(organType, 'Lung') AND assayName = 'scRNA-seq'

Phase-2 database ``htan2_2026_926`` returned no lung RNA rows on that date.
Portal read credentials are not stored in this repo. Re-query the portal
explore catalog to refresh the TSV; do not commit secrets.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "notes" / "htan_cldn4_sting" / "lung_scrna_files.tsv"
OUT = ROOT / "notes" / "htan_cldn4_sting" / "inventory_summary.tsv"


def main():
    df = pd.read_csv(SRC, sep="\t")
    g = (
        df.groupby(["atlas_name", "level", "downloadSource", "FileFormat"], dropna=False)
        .size()
        .reset_index(name="n_files")
        .sort_values(["atlas_name", "level", "downloadSource", "FileFormat"])
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(OUT, sep="\t", index=False)
    print(g.to_string(index=False))
    meta = {
        "n_files": int(len(df)),
        "source_tsv": str(SRC.relative_to(ROOT)),
        "catalog": "HTAN portal files, database htan_2026_932, 2026-09-21",
        "phase2": "htan2_2026_926 lung RNA query returned 0 rows",
        "anonymous_synapse_download": "repo-prod file API returned HTTP 403 without a Synapse session",
    }
    (ROOT / "notes" / "htan_cldn4_sting" / "inventory_meta.json").write_text(
        json.dumps(meta, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
